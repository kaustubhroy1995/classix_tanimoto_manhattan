# Experiment settings

All clustering studies use seed 101. The notebooks retain the executed parameters and repetition schedules. CLASSIX_T is the packed implementation in `src/classix_t.py`, using `src/spmv/` for native operations.

## Synthetic data

The original generator uses `n` points per cluster and per-bit flip probability `flip_prob * prescribed_popcount / dimension`.

| Study | Total rows | Dimensions | CLASSIX radius | minPts |
| --- | --- | --- | --- | --- |
| Sample count | 10,000-90,000 | 1000 | 0.4 | 40 |
| Popcount separation | 10,000 | 1000 | 0.3 | 5 |
| Cluster count | 10,000 | 1000 | 0.4 | 5 |
| Dimension | 10,000 | 500-5000 | 0.4 | 5 |
| Earlier radius sweep in synthetic notebook | 10,000 | 1000 | 0.2-0.7 | 40 |
| Separate radius sensitivity notebook | 10,000 | 1024 | 0.3-0.8 | 1 |

CLASSIX uses `mergeScale=1/sqrt(radius)`. The separate radius study has ten clusters, popcounts between 50 and 300, `flip_prob=0.3`, DBSCAN `min_samples=1`, and DBSCAN epsilon/Taylor-Butina distance cutoff equal to radius. Both BitBIRCH implementations use similarity `1-radius`, branching factor 50 and diameter merging.

## Greg Landrum data

Use the supplied CSV with `pKi > 4`: 90,498 rows and 35 target categories. Morgan fingerprints have radius 2 and 1024 bits.

| Method | Parameters |
| --- | --- |
| CLASSIX_T | radius 0.3, minPts 100, mergeScale 1.33 |
| DBSCAN | epsilon 0.3, min_samples 50 |
| Taylor-Butina | distance cutoff 0.3, reordering disabled |
| BitBIRCH / BitBIRCH-Lean | similarity 9/44, branching factor 50, diameter merging |

No post-clustering size filter is applied to BitBIRCH assignments. ARI/AMI use all rows, with noise as one predicted category. Silhouette and CHI_T exclude noise. Silhouette uses at most 2000 seeded query points against every retained reference row.

Molecular DBSCAN times exact packed-bit radius-graph construction and scikit-learn clustering. Taylor-Butina uses the exact streaming adapter. Input loading and fingerprint construction are outside fit timing.

## Numeric data

The twenty real datasets retain their original features and share feature-wise standardization across methods. Radius/epsilon calibration is metric-specific. ARI selects settings; AMI evaluates those same assignments. Reported runtimes are medians of five fresh fits after a warmup.

The complete comparison includes Euclidean and Manhattan variants of CLASSIX, DBSCAN and OPTICS. `scripts/render_numeric_tables.py` formats the paper's four-method selection from completed results without fitting models.

## Execution

The notebooks create isolated environments with recorded dependencies and run timed fits in serial worker processes. Dataset generation, imports, compilation and scoring are outside clustering timing. Generated records retain per-fit settings and timing boundaries. Wall-clock times vary between machines and runs.

The probability notebook retains its formulas, grids and stored outputs. New Monte Carlo runs use seed 101; the original simulation seed was unrecorded.
