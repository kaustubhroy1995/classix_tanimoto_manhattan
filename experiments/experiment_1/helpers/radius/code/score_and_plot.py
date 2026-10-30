"""Verify saved radius fits and plot ARI/AMI; never fit clustering estimators."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE / 'code'))
METHODS = {
    'classix_packed': ('CLASSIX_T', 'red', 'x', '--'),
    'dbscan': ('DBSCAN', 'blue', 'o', '-.'),
    'bitbirch': ('BitBIRCH', 'green', 'v', ':'),
    'bblean': ('BitBIRCH-Lean', '#9467bd', 's', '--'),
    'butina': ('Taylor-Butina', '#e68613', 'D', '-'),
}
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=HERE / 'runs/full')
    parser.add_argument('--output-dir', type=Path, default=HERE)
    args = parser.parse_args()
    run_dir, output = args.run_dir.resolve(), args.output_dir.resolve()
    assert not (run_dir / 'run.lock').exists(), 'Wait for the fitting process to finish.'
    import numpy as np
    import sklearn
    from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
    from state_checks import array_hash
    from prepare_dataset import digest
    plan = read_configuration(CONFIG)
    cases = sorted(plan['cases'], key=lambda case: case['x'])
    assert [case['x'] for case in cases] == [.3, .4, .5, .6, .7, .8]
    assert len({case['data_id'] for case in cases}) == 1
    generator = dict(num_clusters=10, pops=[50, 300], d=1024, n=1000, flip_prob=.3, seed=101)
    assert all(case['generator'] == generator for case in cases)
    directory = run_dir / 'data' / cases[0]['data_id']
    manifest = read(directory / 'manifest.json')
    assert manifest['shape'] == [10000, 1024] and manifest['parameters'] == generator
    for item in manifest['files'].values():
        array = np.load(directory / item['file'], mmap_mode='r', allow_pickle=False)
        assert digest(array) == item['sha256']
    truth = np.load(directory / 'truth.npy', allow_pickle=False)
    assert np.unique(truth, return_counts=True)[1].tolist() == [1000] * 10
    rows, evidence = [], []
    for case in cases:
        path = run_dir / 'results' / (case['case_id'] + '.json')
        raw = read(path)
        assert raw['completed'] and raw['pilot_rows'] is None and raw['case'] == case
        assert len(raw['runs']) == 5 and {run['variant'] for run in raw['runs']} == set(METHODS)
        for run in raw['runs']:
            assert run['status'] == 'ok' and not run['warmup'] and run['repeat'] == 0
            checked = run['validated']
            assert checked['logical_input_sha256'] == manifest['files']['int32']['sha256']
            assert run['ready']['rows'] == 10000 and run['ready']['features'] == 1024
            labels = np.load(run_dir / checked['output_relative_assignments_file'], allow_pickle=False)
            assert labels.shape == truth.shape and np.all(labels >= 0)
            mapping = {}
            canonical = np.array([mapping.setdefault(int(label), len(mapping)) for label in labels], dtype=np.int64)
            assert array_hash(canonical) == checked['canonical_partition_sha256']
            ari = float(adjusted_rand_score(truth, labels))
            ami = float(adjusted_mutual_info_score(truth, labels, average_method='arithmetic'))
            assert np.isclose(ari, checked['ari'], rtol=0, atol=1e-12)
            assert np.isfinite(ami)
            rows.append(dict(radius=case['x'], similarity=round(1-case['x'], 10), method=run['variant'],
                             ari=ari, ami=ami, clusters=int(len(np.unique(labels))), noise_points=0,
                             scored_rows=len(truth), seconds=run['fit_done']['clustering_seconds'],
                             canonical_partition_sha256=checked['canonical_partition_sha256'],
                             assignments_sha256=sha(run_dir / checked['output_relative_assignments_file'])))
        evidence.append(dict(path=str(path), sha256=sha(path)))
    assert len(rows) == 30
    output.mkdir(parents=True, exist_ok=True)
    result = dict(status='passed', generator=generator, radius_grid=[case['x'] for case in cases],
                  plan_sha256=sha(CONFIG), input_sha256=manifest['files']['int32']['sha256'],
                  sklearn_version=sklearn.__version__, ami_average_method='arithmetic',
                  all_rows_scored=True, run_records=evidence, records=rows)
    (output / 'scores.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                        'axes.titlesize': 10, 'axes.labelsize': 9,
                        'xtick.labelsize': 9, 'ytick.labelsize': 9,
                        'axes.grid': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.3), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.98, bottom=.18, top=.73, wspace=.15)
    handles = []
    for method, (label, color, marker, linestyle) in METHODS.items():
        selected = sorted((row for row in rows if row['method'] == method), key=lambda row: row['radius'])
        x = [row['radius'] for row in selected]
        for ax, metric in zip(axes, ('ari', 'ami')):
            # Open markers help show methods whose quality scores coincide.
            line, = ax.plot(x, [row[metric] for row in selected], label=label,
                            color=color, marker=marker, linestyle=linestyle,
                            linewidth=1.4, markersize=5,
                            markerfacecolor='none',
                            markeredgewidth=1.1)
            np.testing.assert_array_equal(line.get_ydata(), [row[metric] for row in selected])
        handles.append(Line2D([0], [0], label=label, color=color, marker=marker, linestyle=linestyle,
                              linewidth=1.4, markersize=5,
                              markerfacecolor='none'))
    for ax, metric in zip(axes, ('ARI', 'AMI')):
        ax.set_title(metric, pad=5)
        ax.set_xlabel('Distance threshold r')
        ax.set_xlim(.28, .82)
        ax.set_xticks([.3, .4, .5, .6, .7, .8])
        ax.set_ylim(-.04, 1.04)
        ax.set_yticks([0, .25, .5, .75, 1])
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.065, .995),
               ncol=3, fontsize=8, frameon=False, columnspacing=1.7, handlelength=2.4)
    assert not fig.texts
    plots = output / 'plots'
    plots.mkdir(exist_ok=True)
    artifacts = []
    for extension in ('pdf', 'png', 'svg'):
        path = plots / ('radius_ari_ami.' + extension)
        fig.savefig(path, dpi=300, facecolor='white')
        artifacts.append(dict(path=str(path), sha256=sha(path)))
    plt.close(fig)
    pairs = []
    for radius in result['radius_grid']:
        selected = {row['method']: row for row in rows if row['radius'] == radius}
        pairs.append(dict(radius=radius, identical_bitbirch_partitions=
                          selected['bitbirch']['canonical_partition_sha256'] == selected['bblean']['canonical_partition_sha256']))
    verification = dict(status='passed', full_method_settings=30, new_clustering_fits_in_this_script=0,
                        scored_rows_per_fit=10000, all_assignment_and_input_hashes_verified=True,
                        ari_independently_verified=True, ami_computed_from_saved_assignments=True,
                        all_methods_complete=True, score_sha256=sha(output / 'scores.json'),
                        bitbirch_comparison=pairs, artifacts=artifacts)
    (output / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n', encoding='utf-8')
    lines = ['| Radius | Similarity | Method | ARI | AMI | Clusters |', '|---|---|---|---:|---:|---:|']
    lines.extend(f"| {row['radius']:.1f} | {row['similarity']:.1f} | {METHODS[row['method']][0]} | {row['ari']:.4f} | {row['ami']:.4f} | {row['clusters']} |"
                 for row in rows)
    (output / 'scores.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'method_settings': len(rows), 'plot': str(plots / 'radius_ari_ami.pdf')}, indent=2))

if __name__ == '__main__':
    main()
