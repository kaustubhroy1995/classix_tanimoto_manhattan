"""Score saved Landrum assignments uniformly; never fit or tune a model.

Use --root STAGE --data-dir PREPARED --results-dir RESULTS. Optional verified
historical fits are supplied with --references JSON --reference-base WORKSPACE.
Fresh result files take precedence. Missing methods remain explicitly pending.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time

DATA_ID = 'landrum_pki_gt4_morgan1024'
METHODS = ('classix_packed', 'dbscan', 'butina', 'bitbirch', 'bblean')
NAMES = {'classix_packed': 'CLASSIX_T', 'dbscan': 'DBSCAN',
         'butina': 'Taylor-Butina', 'bitbirch': 'BitBIRCH', 'bblean': 'BitBIRCH-Lean'}
PARAMETER_KEYS = {'classix_packed': 'classix', 'dbscan': 'dbscan',
                  'butina': 'butina', 'bitbirch': 'bitbirch', 'bblean': 'bitbirch'}
SEED, MAX_QUERIES = 101, 2000


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def payload(array):
    return hashlib.sha256(memoryview(array).cast('B')).hexdigest()


def resolve(path, base):
    path = Path(path)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # Dataclass annotations require registration.
    spec.loader.exec_module(module)
    return module


def approved(method, parameters):
    if method == 'classix_packed':
        expected = {'radius': .3, 'minPts': 100, 'mergeScale': 1.33}
    elif method in ('bitbirch', 'bblean'):
        expected = {'threshold': 9 / 44, 'branching_factor': 50, 'merge_criterion': 'diameter'}
    elif method == 'dbscan':
        expected = {'eps': .3, 'min_samples': 50}
    elif method == 'butina':
        expected = {'cutoff': .3, 'reordering': False}
    else:
        return False
    return all((key in parameters and math.isclose(float(parameters[key]), value, rel_tol=0, abs_tol=1e-14))
               if isinstance(value, (int, float)) else parameters.get(key) == value
               for key, value in expected.items())


def from_run(raw_path, raw, run, origin, reference=None):
    case, method = raw['case'], run['variant']
    assert case['dataset'] == 'landrum' and case['data_id'] == DATA_ID
    assert case['rows'] == 90498 and case['features'] == 1024
    assert method in METHODS and run['repeat'] == 0 and not run['warmup']
    assert raw.get('pilot_rows') is None and run.get('ready', {}).get('pilot_rows') is None
    parameters = case[PARAMETER_KEYS[method]]
    assert approved(method, parameters), 'Unapproved fixed fast-method parameters'
    record = dict(dataset='landrum', method=method, status=run['status'], parameters=parameters,
                  origin=origin, raw_result=str(raw_path.resolve()), raw_result_sha256=sha(raw_path))
    if run['status'] != 'ok':
        record['reason'] = run.get('error', run.get('reason', run['status']))
        return record
    for key, value in run['ready']['parameters'].items():
        assert parameters[key] == value, (method, key)
    validated = run['validated']
    assignment = Path(validated['assignments_file'])
    if not assignment.is_file():
        relative = validated.get('output_relative_assignments_file')
        assert relative, 'Saved assignment cannot be located'
        assignment = raw_path.parent.parent / relative
    seconds = run['fit_done']['clustering_seconds']
    assert math.isfinite(seconds) and seconds > 0
    record.update(seconds=seconds, assignment=str(assignment.resolve()), validated=validated,
                  manifest_sha256=raw.get('data_manifest_sha256'), reference=reference,
                  source_protocol=raw.get('protocol', run.get('protocol', {})))
    return record


def collect(results_dir, references=None, reference_base=None):
    selected = {}
    if references is not None:
        base = reference_base or references.parent
        records = read(references)
        if isinstance(records, dict):
            records = records['records']
        for reference in records:
            method = reference['method']
            if (reference['dataset'] != 'landrum' or method not in METHODS
                    or not approved(method, reference['parameters'])):
                continue
            raw_path = resolve(reference['raw_result'], base)
            raw = read(raw_path)
            runs = [r for r in raw['runs'] if r['variant'] == method and r['repeat'] == 0 and not r['warmup']]
            assert len(runs) == 1, 'Ambiguous saved reference'
            record = from_run(raw_path, raw, runs[0], 'saved_existing_fit', reference)
            assert record['status'] == 'ok'
            assert math.isclose(record['seconds'], reference['seconds'], rel_tol=0, abs_tol=1e-12)
            assert method not in selected, 'Duplicate approved reference'
            selected[method] = record
    fresh = set()
    for path in sorted(results_dir.glob('*.json')):
        raw = read(path)
        if 'case' not in raw or 'runs' not in raw or raw['case'].get('dataset') != 'landrum':
            continue
        if raw.get('pilot_rows') is not None:
            continue
        for run in raw['runs']:
            method = run['variant']
            if method not in METHODS or run['warmup'] or run['repeat'] != 0:
                continue
            assert method not in fresh, 'More than one fresh fit per method; choose a single run directory'
            fresh.add(method)
            selected[method] = from_run(path, raw, run, 'fresh_fixed_fit')
    for method in METHODS:
        selected.setdefault(method, dict(dataset='landrum', method=method, status='pending',
                            reason='No measured full-dataset fit is available', origin='pending'))
    return selected


def canonical_hash(labels, np):
    mapping, canonical = {}, np.empty(len(labels), dtype=np.int64)
    for i, value in enumerate(labels):
        canonical[i] = mapping.setdefault(int(value), len(mapping))
    return payload(canonical)


def compact_labels(labels, np):
    keep = labels >= 0
    unique, inverse, counts = np.unique(labels[keep], return_inverse=True, return_counts=True)
    compact = np.full(len(labels), -1, dtype=np.int32)
    compact[keep] = inverse
    return compact, counts


def score_silhouette(reference, compact, native, samples_dir, identity, *, seed=SEED, maximum=MAX_QUERIES):
    """Uniform query sample; every retained reference row enters each distance mean."""
    import numpy as np
    keep = compact >= 0
    n_retained = int(keep.sum())
    counts = np.bincount(compact[keep]) if n_retained else np.empty(0, dtype=np.int64)
    k = len(counts)
    result = dict(status='undefined', value=None, reason=None, seed=seed, max_queries=maximum,
                  n_total=len(compact), n_retained=n_retained, k_retained=k,
                  coverage=n_retained / len(compact), singleton_clusters=int(np.count_nonzero(counts == 1)),
                  n_queries=0, query_sample_standard_error=None, full_query_census=False,
                  protocol='Seeded permutation of all row indices, retain non-noise indices, take first up to 2000; sort queries for evaluation. All retained rows are references; self excluded; singleton score zero.')
    if k < 2 or k >= n_retained:
        result['reason'] = 'fewer_than_two_retained_clusters' if k < 2 else 'all_retained_clusters_are_singletons'
        return result
    priority = np.random.default_rng(seed).permutation(len(compact))
    queries = np.sort(priority[keep[priority]][:min(maximum, n_retained)]).astype(np.int64)
    state = native.prepare_labels(compact, counts)
    silhouettes, within, nearest = [np.empty(len(queries), dtype=np.float64) for _ in range(3)]
    for start in range(0, len(queries), native.MAX_QUERIES):
        end = min(start + native.MAX_QUERIES, len(queries))
        indices = queries[start:end]
        distances = native.distance_block(reference, indices)
        s, a, b = native.reduce_block(distances, state, indices, check_distances=False)
        silhouettes[start:end], within[start:end], nearest[start:end] = s, a, b
    assert np.isfinite(silhouettes).all() and np.all(np.abs(silhouettes) <= 1 + 1e-12)
    m = len(queries)
    se = float(np.std(silhouettes, ddof=1) / np.sqrt(m) * np.sqrt(1 - m / n_retained)) if m < n_retained else 0.0
    samples_dir.mkdir(parents=True, exist_ok=True)
    path = samples_dir / f'{identity}.npz'
    np.savez_compressed(path, query_indices=queries, silhouette=silhouettes, within_mean=within, nearest_other_mean=nearest)
    result.update(status='ok', value=float(np.mean(silhouettes)), n_queries=m,
                  query_sample_standard_error=se, full_query_census=m == n_retained,
                  sample_file=str(path.resolve()), sample_file_sha256=sha(path),
                  query_indices_payload_sha256=payload(queries), sample_values_payload_sha256=payload(silhouettes))
    return result


def score_records(records, data_dir, native_code, chi_code, output_dir):
    import numpy as np
    import sklearn
    from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
    from threadpoolctl import threadpool_limits
    native = load_module('landrum_native_silhouette', native_code)
    chi = load_module('landrum_chi_metrics', chi_code)
    directory = data_dir / DATA_ID
    manifest_path = directory / 'manifest.json'
    manifest = read(manifest_path)
    assert manifest['data_id'] == DATA_ID and manifest['shape'] == [90498, 1024]
    assert manifest['parameters']['filter'] == {'column': 'pKi', 'operator': '>', 'threshold': 4.0}
    fp = manifest['parameters']['fingerprint']
    assert fp['radius'] == 2 and fp['fpSize'] == 1024
    arrays = {kind: np.load(directory / manifest['files'][kind]['file'], mmap_mode='r', allow_pickle=False)
              for kind in ('int32', 'packed_little', 'truth')}
    for kind, array in arrays.items():
        assert array.flags.c_contiguous and payload(array) == manifest['files'][kind]['sha256'], kind
    X, packed, truth = [arrays[k] for k in ('int32', 'packed_little', 'truth')]
    assert X.shape == (90498, 1024) and X.dtype == np.int32
    assert packed.shape == (len(X), 128) and packed.dtype == np.uint8 and truth.shape == (len(X),)
    input_hash, manifest_hash = payload(X), sha(manifest_path)
    output, cache, started = [], {}, time.perf_counter()
    with threadpool_limits(limits=1):
        packed_reference = native.prepare_packed(packed, n_features=1024)
        for method in METHODS:
            record = records[method]
            if record['status'] != 'ok':
                output.append(record)
                continue
            print(f'Scoring {method}: ARI, AMI, silhouette and CHI_T', flush=True)
            labels = np.load(record['assignment'], allow_pickle=False)
            assert labels.shape == (len(X),) and labels.dtype.kind in 'iu' and np.all(labels >= -1)
            labels_hash, file_hash = payload(np.ascontiguousarray(labels)), sha(record['assignment'])
            validated = record['validated']
            assert validated['logical_input_sha256'] == input_hash
            assert record['manifest_sha256'] == manifest_hash
            assert canonical_hash(labels, np) == validated['canonical_partition_sha256']
            reference = record['reference']
            if reference and reference.get('labels_sha256'):
                assert reference['labels_sha256'] == file_hash
            ari, ami = float(adjusted_rand_score(truth, labels)), float(adjusted_mutual_info_score(truth, labels, average_method='arithmetic'))
            for name, value in [('ari', ari), ('ami', ami)]:
                for prior in [validated, reference or {}]:
                    if name in prior:
                        assert math.isclose(prior[name], value, rel_tol=0, abs_tol=1e-12), (method, name)
            compact, counts = compact_labels(labels, np)
            identity = payload(compact)
            reused = identity in cache
            if not reused:
                quality = chi.evaluate(X, labels, min_sizes=(0,), pairwise=False, block_size=128)[0]
                silhouette = score_silhouette(packed_reference, compact, native, output_dir / 'samples', identity)
                assert quality['n_retained'] == silhouette['n_retained']
                assert quality['k_retained'] == silhouette['k_retained'] == len(counts)
                cache[identity] = (quality, silhouette, method)
            quality, silhouette, source_method = cache[identity]
            assert payload(np.ascontiguousarray(labels)) == labels_hash
            output.append(dict(dataset='landrum', method=method, status='ok', parameters=record['parameters'],
                               seconds=record['seconds'], ari=ari, ami=ami, quality=quality, silhouette_tanimoto=silhouette,
                               quality_reused_from_method=source_method if reused else None,
                               origin=record['origin'], raw_result=record['raw_result'], raw_result_sha256=record['raw_result_sha256'],
                               assignment=record['assignment'], assignment_file_sha256=file_hash, labels_payload_sha256=labels_hash,
                               canonical_partition_sha256=validated['canonical_partition_sha256'], logical_input_sha256=input_hash,
                               reference_scores_reproduced=bool(reference), source_protocol=record['source_protocol']))
        for kind, array in arrays.items():
            assert payload(array) == manifest['files'][kind]['sha256'], 'Scoring changed input'
    return dict(status='complete' if all(r['status'] == 'ok' for r in output) else 'partial', records=output,
                inputs=dict(data_id=DATA_ID, shape=list(X.shape), manifest=str(manifest_path.resolve()),
                            manifest_sha256=manifest_hash, logical_input_sha256=input_hash,
                            reference_categories=int(np.unique(truth).size), duplicates_retained=True,
                            packed_payload_sha256=manifest['files']['packed_little']['sha256'], truth_payload_sha256=manifest['files']['truth']['sha256']),
                protocol=dict(seed=SEED, max_silhouette_queries=MAX_QUERIES, filtering='All non-noise rows, including singleton clusters; no post-clustering size filter',
                              chi_definition=chi.DEFINITION_VERSION, chi_pairwise=False,
                              chi_code_sha256=sha(chi_code), native=native.native_metadata()),
                scoring_environment=dict(python=sys.version, platform=platform.platform(), numpy=np.__version__, sklearn=sklearn.__version__),
                scorer_source_sha256=sha(__file__), scoring_seconds=time.perf_counter()-started, fits_performed=0,
                timing_note='One measured fit per available row, no warmup. Saved timings remain historical. Loading, fingerprint conversion and scoring are excluded.',
                baseline_timing_note='DBSCAN uses an exact packed Tanimoto radius graph followed by scikit-learn DBSCAN with precomputed distances; graph construction and clustering are timed together. Taylor-Butina uses the exact streaming RDKit-distance adapter with reordering=False; all neighborhood work is timed.',
                external_score_note='ARI and arithmetic-normalized AMI use all 90,498 rows, with noise as one predicted category. The 35 supplied target categories are reference labels, not unique molecule classes.',
                quality_note='Silhouette and CHI_T use all non-noise clusters, including singletons, with no post-clustering size filter. Silhouette averages up to 2,000 seeded query rows against every retained reference row; singleton silhouettes are zero. CHI_T is evaluated on all retained rows using binary majority representatives. Undefined scores stay null.',
                uncertainty_note='The silhouette standard error describes query sampling only, with finite-population correction; it is not clustering-run variability.')


def table_rows(report):
    rows = []
    for row in report['records']:
        if row['status'] != 'ok':
            rows.append([NAMES[row['method']], row['status']] + ['—'] * 4)
            continue
        value = lambda score, digits: 'undefined' if score['value'] is None else f"{score['value']:.{digits}f}"
        rows.append([NAMES[row['method']], f"{row['seconds']:.3f}", f"{row['ari']:.4f}", f"{row['ami']:.4f}",
                     value(row['silhouette_tanimoto'], 4), value(row['quality']['chi_tanimoto'], 3)])
    return rows


def write_outputs(report, output_dir, plots=True):
    output_dir.mkdir(parents=True, exist_ok=True)
    headers, rows = ['Method', 'Runtime (s)', 'ARI', 'AMI', 'Silhouette', 'CHI_T'], table_rows(report)
    coverage = '; '.join(f"{NAMES[r['method']]} {100 * r['quality']['coverage']:.2f}%" for r in report['records'] if r['status'] == 'ok')
    params = 'CLASSIX_T: radius=0.3, minPts=100, mergeScale=1.33. BitBIRCH/Lean: similarity=9/44, branching_factor=50, diameter merging.'
    lines = ['# Greg Landrum molecular comparison', '', '| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |']
    lines += ['| ' + ' | '.join(row) + ' |' for row in rows]
    lines += ['', '90,498 molecules after pKi > 4; radius-2 Morgan fingerprints, 1,024 bits; 35 reference target categories.', '',
              params, '', 'Lower runtime and higher ARI, AMI, silhouette and CHI_T are preferable.', '',
              report['external_score_note'], '', report['quality_note'], '', 'Non-noise coverage: ' + coverage + '.', '',
              report['uncertainty_note'], '', report['timing_note'], '', report['baseline_timing_note'], '',
              'Timing origins: ' + '; '.join(f"{NAMES[r['method']]}: {r['origin']}" for r in report['records']) + '.']
    for row in report['records']:
        if row['method'] in ('dbscan', 'butina') and row.get('parameters'):
            lines += ['', f"{NAMES[row['method']]} parameters: `{json.dumps(row['parameters'], sort_keys=True)}`."]
        if row['status'] != 'ok':
            lines += ['', f"{NAMES[row['method']]}: {row['status']} — {row.get('reason', 'No successful assignment available')}."]
        else:
            for name, score in [('silhouette', row['silhouette_tanimoto']), ('CHI_T', row['quality']['chi_tanimoto'])]:
                if score['status'] != 'ok':
                    lines += ['', f"{NAMES[row['method']]} {name} is undefined: {score['reason']}."]
    lines += ['', 'Full precision, query-sample uncertainty, retained populations, parameters and evidence hashes are recorded in results_landrum.json.']
    artifacts = []
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
        fig, ax = plt.subplots(figsize=(12, 5.2))
        ax.axis('off')
        grid = ax.table(cellText=rows, colLabels=headers, cellLoc='center', bbox=[0, .31, 1, .61], colWidths=[.22, .16, .14, .14, .17, .17])
        grid.auto_set_font_size(False); grid.set_fontsize(11)
        for (r, c), cell in grid.get_celld().items():
            cell.set_edgecolor('#cccccc'); cell.set_linewidth(.4)
            cell.set_facecolor('#dddddd' if r == 0 else '#f2f4f6' if r % 2 else 'white')
            if r == 0: cell.set_text_props(weight='bold')
        ax.set_title('Greg Landrum: 90,498 molecules · Morgan radius 2, 1,024 bits', pad=12)
        ax.text(0, .255, 'ARI/AMI: all rows; 35 reference target categories. Internal scores: all non-noise clusters, including singletons.', fontsize=8, transform=ax.transAxes)
        ax.text(0, .21, 'Silhouette: up to 2,000 seeded queries against all retained rows. Coverage: ' + coverage + '.', fontsize=8, transform=ax.transAxes)
        ax.text(0, .165, 'CLASSIX r=.3, minPts=100, mergeScale=1.33; BitBIRCH/Lean similarity=9/44, branching factor=50, diameter merging.', fontsize=8, transform=ax.transAxes)
        ax.text(0, .12, 'DBSCAN: distance cutoff=.3, min_samples=50; Taylor-Butina: distance cutoff=.3, reordering=False.', fontsize=8, transform=ax.transAxes)
        ax.text(0, .075, 'DBSCAN: exact packed graph + scikit-learn; Butina: exact streaming RDKit adapter. All neighborhood work is timed.', fontsize=8, transform=ax.transAxes)
        origins = {r['origin'] for r in report['records'] if r['status'] == 'ok'}
        origin_note = 'One measured fit per row; saved and fresh timings are combined.' if len(origins) > 1 else 'One measured fit per row; timings are reused from verified saved fits.' if origins == {'saved_existing_fit'} else 'One measured fit per row; fresh reproduction timings.'
        pending_note = ' Pending rows have no measured scores.' if any(r['status'] == 'pending' for r in report['records']) else ''
        ax.text(0, .03, origin_note + pending_note, fontsize=8, transform=ax.transAxes)
        fig.tight_layout()
        plot_dir = output_dir / 'plots'; plot_dir.mkdir(parents=True, exist_ok=True)
        for extension in ('png', 'pdf', 'svg'):
            path = plot_dir / f'landrum_comparison_table.{extension}'
            fig.savefig(path, dpi=200, bbox_inches='tight'); artifacts.append(str(path.resolve()))
        plt.close(fig)
    (output_dir / 'table.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    report['artifacts'] = artifacts
    (output_dir / 'results_landrum.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--results-dir', type=Path, required=True)
    parser.add_argument('--references', type=Path)
    parser.add_argument('--reference-base', type=Path)
    parser.add_argument('--native-code', type=Path)
    parser.add_argument('--chi-code', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    native_code = args.native_code or args.root / 'code' / 'native_silhouette.py'
    chi_code = args.chi_code or args.root / 'code' / 'chi_tanimoto.py'
    output_dir = args.output_dir or args.root
    records = collect(args.results_dir, args.references, args.reference_base)
    report = score_records(records, args.data_dir, native_code, chi_code, output_dir)
    write_outputs(report, output_dir, plots=not args.no_plots)
    print(json.dumps({'status': report['status'], 'rows': len(report['records']), 'fits_performed': 0,
                      'output': str((output_dir / 'results_landrum.json').resolve())}), flush=True)


if __name__ == '__main__':
    main()
