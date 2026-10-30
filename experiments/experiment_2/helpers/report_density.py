"""Report completed Experiment 2 density baselines, without fitting a model.

Usage: python report_density.py --root RUN_ROOT [--no-plots]
Reads the frozen panel, standardized inputs, selected results and saved CLASSIX
grid. Requires the complete 20-dataset/6-variant result; never hides partial runs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import statistics

ALGORITHMS = ('classix', 'dbscan', 'optics')
METRICS = ('euclidean', 'manhattan')
NAMES = {'classix': 'CLASSIX', 'dbscan': 'DBSCAN', 'optics': 'OPTICS'}
DISPLAY = {'iris': 'Iris', 'banknote': 'Banknote', 'wine': 'Wine',
           'breast_cancer': 'Breast cancer', 'digits': 'Digits', 'seeds': 'Seeds',
           'glass': 'Glass', 'sonar': 'Sonar', 'vehicle': 'Vehicle',
           'segmentation': 'Segmentation', 'breast_tissue': 'Breast tissue',
           'leaf': 'Leaf', 'libras': 'Libras', 'vertebral': 'Vertebral',
           'wireless': 'Wireless', 'parkinsons': 'Parkinsons', 'coimbra': 'Coimbra',
           'raisin': 'Raisin', 'haberman': 'Haberman', 'musk': 'Musk'}
VARIANTS = [(algorithm, metric) for algorithm in ALGORITHMS for metric in METRICS]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, root, sources, *, jsonl=False):
    raw = path.read_bytes()
    sources[path.relative_to(root).as_posix()] = hashlib.sha256(raw).hexdigest()
    text = raw.decode('utf-8-sig')
    return [json.loads(line) for line in text.splitlines() if line.strip()] if jsonl else json.loads(text)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def key(row):
    return row['dataset'], row['algorithm'], row['metric']


def collect(root):
    import numpy as np
    sources = {}
    policy = read(root / 'protocol.json', root, sources)
    provenance = read(root / 'input_provenance.json', root, sources)
    selected = read(root / 'results/selected.json', root, sources)
    grid = read(root / 'results/classix_grid_prior.jsonl', root, sources, jsonl=True)
    datasets, panel = policy['datasets'], policy['panel']
    assert len(datasets) == len(set(datasets)) == 20
    assert len(panel) == len(set(panel)) == 10 and set(panel) <= set(datasets)
    expected = {(dataset, algorithm, metric) for dataset in datasets for algorithm, metric in VARIANTS}
    lookup = {key(row): row for row in selected}
    assert len(selected) == len(lookup) == 120 and set(lookup) == expected, 'Require all 120 selected configurations'
    metadata = {}
    for dataset in datasets:
        path = root / 'data' / f'{dataset}.npz'
        relative = path.relative_to(root).as_posix()
        digest = sha(path)
        assert digest == provenance['inputs'][relative], f'Changed input: {dataset}'
        sources[relative] = digest
        with np.load(path, allow_pickle=False) as data:
            X, y = data['X'], data['y']
            n, d = X.shape
            assert y.ndim == 1 and len(y) == n
            metadata[dataset] = dict(n=n, d=d, reference_classes=int(np.unique(y).size))
        for algorithm, metric in VARIANTS:
            row = lookup[dataset, algorithm, metric]
            assert row.get('status', 'ok') == 'ok'
            assert row['input_sha256'] == digest
            assert all(finite(row[name]) for name in ('ari', 'ami', 'fit_seconds_median'))
            assert -1 <= row['ari'] <= 1 and -1 <= row['ami'] <= 1
            times = row['timing_seconds']
            assert len(times) == 5 and all(finite(value) and value > 0 for value in times)
            assert math.isclose(row['fit_seconds_median'], statistics.median(times), rel_tol=1e-12, abs_tol=1e-15)
            assert 0 <= row['noise_points'] <= n
            assert math.isclose(row['retained_fraction'], 1 - row['noise_points'] / n, abs_tol=1e-12)
            label_path = root / row['labels_path']
            label_hash = sha(label_path)
            assert label_hash == row['labels_file_sha256']
            sources[label_path.relative_to(root).as_posix()] = label_hash
    iris = [row for row in grid if row['dataset'] == 'iris']
    assert iris and all(row.get('status', 'ok') == 'ok' for row in iris)
    return dict(policy=policy, metadata=metadata, selected=selected, lookup=lookup,
                iris_grid=iris, sources=sources)


def name(dataset):
    return DISPLAY.get(dataset, dataset.replace('_', ' ').capitalize())


def fmt(value, score=False):
    return f'{value:.4f}' if score else f'{value:.3g}'


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                       '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                      ['| ' + ' | '.join(str(cell).replace('|', '\\|') for cell in row) + ' |' for row in rows])


def compact_rows(report, datasets):
    rows = []
    for dataset in datasets:
        meta = report['metadata'][dataset]
        for metric in METRICS:
            cells = [f"{name(dataset)} ({meta['n']:,} x {meta['d']})", metric.capitalize()]
            for algorithm in ALGORITHMS:
                row = report['lookup'][dataset, algorithm, metric]
                cells += [fmt(row['ari'], True), fmt(row['ami'], True), fmt(1000 * row['fit_seconds_median'])]
            rows.append(cells)
    return rows


def export_csv(path, records):
    fields = list(dict.fromkeys(field for row in records for field in row))
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def records_for(report, datasets):
    rows = []
    for dataset in datasets:
        for algorithm, metric in VARIANTS:
            raw = report['lookup'][dataset, algorithm, metric]
            row = dict(dataset=dataset, **report['metadata'][dataset], algorithm=algorithm,
                       metric=metric, ari=raw['ari'], ami=raw['ami'],
                       fit_seconds_median=raw['fit_seconds_median'],
                       fit_ms_median=1000 * raw['fit_seconds_median'],
                       timing_repetitions=len(raw['timing_seconds']),
                       n_clusters=raw['n_clusters'], noise_points=raw['noise_points'],
                       retained_fraction=raw['retained_fraction'])
            for field in ('radius', 'min_pts', 'merge_scale', 'parameter_name', 'parameter',
                          'parameter_index', 'min_samples', 'selection_origin',
                          'partition_sha256', 'input_sha256', 'labels_file_sha256'):
                row[field] = raw.get(field, '')
            rows.append(row)
    return rows


def latex_table(report):
    lines = [r'\begin{table}[p]', r'\centering\color{red}',
             r'\begingroup\renewcommand{\baselinestretch}{1}\footnotesize',
             r'\setlength{\tabcolsep}{2.5pt}', r'\renewcommand{\arraystretch}{1.15}',
             r'\begin{tabular}{@{}llrrrrrrrrr@{}}', r'\toprule',
             r'Dataset & Metric & \multicolumn{3}{c}{CLASSIX} & \multicolumn{3}{c}{DBSCAN} & \multicolumn{3}{c}{OPTICS} \\',
             r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}',
             r'$(n,d)$ & & ARI & AMI & ms & ARI & AMI & ms & ARI & AMI & ms \\', r'\midrule']
    for index, dataset in enumerate(report['policy']['panel']):
        meta = report['metadata'][dataset]
        for metric in METRICS:
            first = name(dataset) if metric == 'euclidean' else f"$({meta['n']},{meta['d']})$"
            cells = [first, 'E' if metric == 'euclidean' else 'M']
            for algorithm in ALGORITHMS:
                row = report['lookup'][dataset, algorithm, metric]
                cells += [fmt(row['ari'], True), fmt(row['ami'], True), fmt(1000 * row['fit_seconds_median'])]
            lines.append(' & '.join(cells) + r' \\')
        if index != len(report['policy']['panel']) - 1:
            lines.append(r'\addlinespace[2pt]')
    lines += [r'\bottomrule', r'\end{tabular}', r'\endgroup', r'\DraftCaption',
              r'\caption{Original-feature comparison on the fixed illustrative ten-dataset panel. E and M denote Euclidean and Manhattan distance; $n$ and $d$ are sample and feature counts. Settings maximize in-sample ARI, and AMI uses the same assignments, including noise as one category. Times are medians of five fresh fits in milliseconds, excluding tuning and external preprocessing. The full twenty-dataset screen is in the Experiment~2 notebook.}',
              r'\label{tab:fulldimensiondata}\label{tab:originalfeatures}', r'\end{table}']
    return '\n'.join(lines) + '\n'


def outcome(a, b):
    # Comparisons use a tiny numerical-tie tolerance, not four-decimal rounding.
    return 'M' if b > a + 1e-12 else 'E' if a > b + 1e-12 else 'tie'


def summary(report):
    result = {'datasets': len(report['policy']['datasets']), 'selected_variants': len(report['selected']),
              'timed_fits': sum(len(row['timing_seconds']) for row in report['selected']),
              'panel': report['policy']['panel'], 'comparison_tolerance': 1e-12,
              'within_algorithm_metric_comparisons': {}, 'best_tested_variants': {},
              'median_time_ranges_ms': {}}
    for algorithm in ALGORITHMS:
        counts = {measure: Counter() for measure in ('ari', 'ami', 'joint')}
        for dataset in report['policy']['datasets']:
            e = report['lookup'][dataset, algorithm, 'euclidean']
            m = report['lookup'][dataset, algorithm, 'manhattan']
            decisions = [outcome(e[measure], m[measure]) for measure in ('ari', 'ami')]
            for measure, decision in zip(('ari', 'ami'), decisions):
                counts[measure][decision] += 1
            counts['joint'][decisions[0] if decisions[0] == decisions[1] and decisions[0] != 'tie' else 'mixed_or_tied'] += 1
        result['within_algorithm_metric_comparisons'][algorithm] = {measure: dict(count) for measure, count in counts.items()}
    for measure in ('ari', 'ami', 'fit_seconds_median'):
        counts = Counter()
        winners = {}
        for dataset in report['policy']['datasets']:
            values = {f'{algorithm}_{metric}': report['lookup'][dataset, algorithm, metric][measure] for algorithm, metric in VARIANTS}
            target = min(values.values()) if measure == 'fit_seconds_median' else max(values.values())
            tied = [variant for variant, value in values.items() if math.isclose(value, target, rel_tol=0, abs_tol=1e-12)]
            winners[dataset] = tied
            counts.update(tied)
        result['best_tested_variants'][measure] = {'co_winner_counts': dict(counts), 'by_dataset': winners,
                                                  'ties_credited_to_each_variant': True}
    for algorithm, metric in VARIANTS:
        values = [1000 * report['lookup'][dataset, algorithm, metric]['fit_seconds_median'] for dataset in report['policy']['datasets']]
        result['median_time_ranges_ms'][f'{algorithm}_{metric}'] = {'minimum': min(values), 'maximum': max(values)}
    return result


def write_tables(report, root):
    destination = root / 'tables'
    destination.mkdir(exist_ok=True)
    headers = ['Dataset (n x d)', 'Metric'] + [f'{NAMES[a]} {m}' for a in ALGORITHMS for m in ('ARI', 'AMI', 'ms')]
    tables = {}
    for stem, datasets in [('main_10', report['policy']['panel']), ('all_20', report['policy']['datasets'])]:
        content = table(headers, compact_rows(report, datasets))
        (destination / f'{stem}.md').write_text(content + '\n', encoding='utf-8')
        export_csv(destination / f'{stem}.csv', records_for(report, datasets))
        tables[stem] = content
    parameters = []
    for row in records_for(report, report['policy']['datasets']):
        setting = (f"radius={row['radius']:.8g}, minPts={row['min_pts']}, mergeScale={row['merge_scale']}"
                   if row['algorithm'] == 'classix' else f"{row['parameter_name']}={row['parameter']:.8g}, min_samples={row['min_samples']}")
        parameters.append([name(row['dataset']), NAMES[row['algorithm']], row['metric'], setting,
                           row['n_clusters'], row['noise_points'], f"{row['retained_fraction']:.4f}",
                           fmt(row['fit_ms_median'])])
    tables['parameters'] = table(['Dataset', 'Method', 'Metric', 'Selected setting', 'Clusters', 'Noise rows', 'Retained fraction', 'Median ms'], parameters)
    (destination / 'parameters.md').write_text(tables['parameters'] + '\n', encoding='utf-8')
    export_csv(destination / 'parameters.csv', records_for(report, report['policy']['datasets']))
    (destination / 'manuscript_table.tex').write_text(latex_table(report), encoding='utf-8')
    return tables


def render(report, root):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
    destination = root / 'plots'
    destination.mkdir(exist_ok=True)
    artifacts = []

    def save(fig, stem):
        for extension in ('pdf', 'png', 'svg'):
            path = destination / f'{stem}.{extension}'
            fig.savefig(path, dpi=200, bbox_inches='tight')
            artifacts.append(path.relative_to(root).as_posix())
        plt.close(fig)

    colors = {'classix': '#4477AA', 'dbscan': '#EE6677', 'optics': '#228833'}
    datasets = report['policy']['datasets']
    for measure, label, stem in [('ari', 'ARI (higher is better)', 'all20_ari'),
                                  ('ami', 'AMI (higher is better)', 'all20_ami'),
                                  ('fit_seconds_median', 'Median fit time (ms, log scale)', 'all20_runtime')]:
        fig, axes = plt.subplots(2, 1, figsize=(13, 8.2), sharey=True)
        for ax, group in zip(axes, [datasets[:10], datasets[10:]]):
            for variant_index, (algorithm, metric) in enumerate(VARIANTS):
                values = [report['lookup'][dataset, algorithm, metric][measure] for dataset in group]
                if measure == 'fit_seconds_median':
                    values = [1000 * value for value in values]
                x = np.arange(len(group)) + (variant_index - 2.5) * .13
                ax.bar(x, values, width=.12, label=f'{NAMES[algorithm]} {"E" if metric == "euclidean" else "M"}',
                       color=colors[algorithm], alpha=1 if metric == 'euclidean' else .55,
                       hatch='' if metric == 'euclidean' else '///', linewidth=.25, edgecolor=colors[algorithm])
            ax.set_xticks(range(len(group)), [name(dataset) for dataset in group], rotation=23, ha='right')
            ax.set_ylabel(label)
            ax.grid(axis='y', alpha=.22)
            ax.set_axisbelow(True)
            if measure == 'fit_seconds_median':
                ax.set_yscale('log')
            else:
                ax.axhline(0, color='#666666', linewidth=.7)
                ax.set_ylim(min(-.03, min(row[measure] for row in report['selected']) - .03), 1.05)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=6, frameon=False, bbox_to_anchor=(.5, .96))
        fig.suptitle('Experiment 2: complete twenty-dataset original-feature comparison', y=.995)
        footnote = ('E: Euclidean; M: Manhattan. Identical standardized original features; in-sample ARI-selected settings.\n'
                    'AMI uses those assignments, including noise; fit time is the median of five fresh fits, excluding tuning.')
        fig.text(.5, .015, footnote, ha='center', va='bottom', fontsize=9)
        fig.tight_layout(rect=(0, .075, 1, .915))
        save(fig, stem)

    fig, ax = plt.subplots(figsize=(13.5, 9.0))
    ax.axis('off')
    cells = compact_rows(report, report['policy']['panel'])
    for index, row in enumerate(cells):
        row[1] = 'E' if row[1] == 'Euclidean' else 'M'
        if index % 2:
            row[0] = ''
    headers = ['Dataset (n x d)', 'Metric'] + [f'{NAMES[a]}\n{m}' for a in ALGORITHMS for m in ('ARI', 'AMI', 'ms')]
    drawn = ax.table(cellText=cells, colLabels=headers, cellLoc='center',
                     colWidths=[.185, .05] + [.085] * 9, bbox=[0, .04, 1, .91])
    drawn.auto_set_font_size(False)
    drawn.set_fontsize(9)
    for (r, c), cell in drawn.get_celld().items():
        cell.set_linewidth(.35)
        cell.set_edgecolor('#cccccc')
        if r == 0:
            cell.set_facecolor('#dddddd')
            cell.set_text_props(weight='bold')
            cell.set_height(cell.get_height() * 1.25)
        else:
            cell.set_facecolor('#f3f5f7' if ((r - 1) // 2) % 2 == 0 else 'white')
            if c == 0:
                cell.set_text_props(ha='left')
    fig.suptitle('Experiment 2: fixed illustrative ten-dataset panel', y=.97)
    fig.text(.5, .025, 'E: Euclidean; M: Manhattan. Four-decimal scores at in-sample ARI-selected settings; time is median fit ms.\n'
             'AMI uses the same assignments; noise is one category. The complete twenty-dataset screen is retained.',
             ha='center', va='bottom', fontsize=9)
    fig.tight_layout(rect=(0, .07, 1, .96))
    save(fig, 'approved10_table')

    grid = report['iris_grid']
    minima = sorted({row['min_pts'] for row in grid})
    indices = sorted({row['radius_index'] for row in grid})
    arrays, radii, selections = {}, {}, {}
    assert len({row['merge_scale'] for row in grid}) == 1
    merge_scale = grid[0]['merge_scale']
    for metric in METRICS:
        records = [row for row in grid if row['method'] == metric]
        cells = {(row['min_pts'], row['radius_index']): row for row in records}
        assert len(cells) == len(records) == len(minima) * len(indices)
        arrays[metric] = np.array([[cells[m, j]['ari'] for j in indices] for m in minima])
        radii[metric] = [cells[minima[0], j]['radius'] for j in indices]
        chosen = report['lookup']['iris', 'classix', metric]
        matches = [row for row in records if row['min_pts'] == chosen['min_pts'] and
                   math.isclose(row['radius'], chosen['radius'], rel_tol=1e-12)]
        assert len(matches) == 1 and math.isclose(matches[0]['ari'], chosen['ari'], abs_tol=1e-12)
        assert math.isclose(chosen['ari'], float(arrays[metric].max()), abs_tol=1e-12)
        selections[metric] = matches[0]
    vmin = min(0., *(float(array.min()) for array in arrays.values()))
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.05), layout='constrained')
    for ax, metric in zip(axes, METRICS):
        image = ax.imshow(arrays[metric], origin='lower', aspect='auto', interpolation='nearest', cmap='viridis', vmin=vmin, vmax=1)
        ticks = np.unique(np.linspace(0, len(indices) - 1, 5).round().astype(int))
        ax.set_xticks(ticks, [f'{radii[metric][j]:.3f}' for j in ticks])
        ax.set_yticks(range(len(minima)), [str(value) for value in minima])
        ax.set_xlabel('Radius (geometric grid)')
        ax.set_ylabel('minPts')
        ax.set_title(f'{metric.capitalize()} CLASSIX')
        chosen = selections[metric]
        ax.scatter(indices.index(chosen['radius_index']), minima.index(chosen['min_pts']),
                   marker='*', s=85, facecolors='white', edgecolors='black', linewidths=.6)
    fig.colorbar(image, ax=axes, label='ARI', shrink=.92)
    save(fig, 'original_iris_radius_minpts')
    report['iris_heatmap'] = dict(saved_fits=len(grid), merge_scale=merge_scale,
                                  shared_color_limits=[vmin, 1], selected=selections,
                                  source='results/classix_grid_prior.jsonl', clustering_performed=False)
    return artifacts


def write_report(report, tables, root, figures):
    facts = summary(report)
    lines = ['# Experiment 2: original-feature CLASSIX, DBSCAN and OPTICS', '',
             'All six configurations use the same standardized original numeric features. DBSCAN and OPTICS are the scikit-learn implementations with Euclidean and Manhattan distance. UMAP data are not included. The ten-case panel is fixed from the earlier CLASSIX screen, before these baseline results; it is illustrative rather than a representative sample.', '',
             '## Fixed ten-dataset table', '', tables['main_10'], '',
             'E and M in the standalone table mean Euclidean and Manhattan. Scores use all rows, with noise as one category. Each setting maximizes in-sample ARI; AMI evaluates that identical assignment. Times are medians of five fresh fits, after one discarded warmup. Model-internal preprocessing and clustering are included; construction, input copies, external standardization, loading, tuning and scoring are excluded. The six variants are interleaved in a seeded order during timing.', '',
             '## Measured findings across the complete screen', '']
    for algorithm in ALGORITHMS:
        counts = facts['within_algorithm_metric_comparisons'][algorithm]['joint']
        lines.append(f"- {NAMES[algorithm]}: Manhattan has higher ARI and AMI on {counts.get('M', 0)} of 20 datasets; Euclidean on {counts.get('E', 0)}; {counts.get('mixed_or_tied', 0)} are mixed or tied.")
    lines += ['', 'These are descriptive counts of this screen, with a 1e-12 numerical tie tolerance. They do not establish significance, unseen-data performance or a metric-only cause. Comparisons use complete implementations with different normalization, sorting and execution kernels.', '',
              '## Complete twenty-dataset table', '', tables['all_20'], '',
              '## Selected settings and coverage', '', tables['parameters'], '',
              'CLASSIX retains the existing ARI-selected configurations and was retimed alongside the new baselines. Its original grid has 32 radii and six minPts values, with one prescribed boundary extension. DBSCAN uses 32 geometrically spaced, distance-calibrated eps values and min_samples {2,5,10,20,40,80}; OPTICS uses 32 geometric xi values from 0.01 to 0.3 and the same min_samples grid, with max_eps=infinity. Both baseline searches have one prescribed extension shared across metrics when a winner reaches a boundary. Full-precision ARI ties favor smaller min_samples and then the smaller parameter index.', '',
              'OPTICS search reuses each reachability graph over xi choices. Reported times come from fresh complete OPTICS.fit calls, including graph construction and extraction; search-extraction time is never presented as clustering runtime. Parameter grids compare different algorithms, not identical clustering objectives. Selected parameters are best tested settings, not global optima.', '',
              'The radius/minPts heatmap uses saved original-feature Iris CLASSIX fits with mergeScale 1.4 and a shared color scale; it performs no new clustering. Five-order CLASSIX checks remain in the original Experiment 2 report; no equivalent new order study is claimed for DBSCAN or OPTICS.', '',
              '## Files', '',
              'The tables directory contains full-precision CSV records, readable Markdown tables and a compact red LaTeX table. SUMMARY.json records all input hashes, selected records, dimensions, measured counts and per-dataset best-tested variants. Scores are not rounded before comparisons. The complete 20-dataset report retains all attempted cases.', '']
    for file in figures:
        if file.endswith('.png'):
            stem = file[:-4]
            lines.append(f'- [{Path(stem).name}]({file}) ([PDF]({stem}.pdf), [SVG]({stem}.svg))')
    (root / 'RESULTS.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    result = dict(status='complete', **facts, metadata=report['metadata'], selected=report['selected'],
                  sources=report['sources'], artifacts=figures, iris_heatmap=report.get('iris_heatmap'),
                  report_source_sha256=sha(__file__))
    (root / 'SUMMARY.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return facts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=DENSITY)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    report = collect(root)
    tables = write_tables(report, root)
    figures = [] if args.no_plots else render(report, root)
    for relative, digest in report['sources'].items():
        assert sha(root / relative) == digest, 'An input changed during report generation'
    facts = write_report(report, tables, root, figures)
    print(json.dumps(dict(status='complete', datasets=facts['datasets'], selected_variants=facts['selected_variants'],
                          timed_fits=facts['timed_fits'], plot_files=len(figures), root=str(root))))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
