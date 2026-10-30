# Reproducing the paper experiments

The four notebooks below contain the experiment code, frozen designs and pinned
environments. They generate their working sources, inputs, assignments, records
and plots when run. Supplied plots are presentation artifacts, not evidence that
a notebook has been freshly executed. Reproduction concerns the recorded
experiments; it does not silently retune them to match later changes in prose.

## Paper-to-code map

| Paper output | Notebook / source | Generated output |
| --- | --- | --- |
| Figure 1, Manhattan search geometry | TikZ in the manuscript | No Python experiment |
| Figures 2–3, Tanimoto search-bound probabilities | Historical `classix_t_prob.ipynb` | `all_plot_calculation.pdf`, `plots/all_plots.pdf` |
| Figures 4–6, synthetic sample-count, popcount and dimension timings | [01_synthetic_benchmarks.ipynb](experiment_1/01_synthetic_benchmarks.ipynb) | `sample_count_runtime`, `popcount_separation_runtime`, `dimension_runtime` plots |
| Figure 7, radius ARI/AMI curves | [03_radius_sensitivity.ipynb](experiment_1/03_radius_sensitivity.ipynb) | `radius_ari_ami` plot |
| Table 1, synthetic settings | Frozen plans in the synthetic and radius notebooks | Per-case design records |
| Table 2, five-method Landrum comparison | [02_molecular_benchmarks.ipynb](experiment_1/02_molecular_benchmarks.ipynb) | Landrum comparison table and scored records |
| Tables 3–4, twenty numeric datasets | [01_euclidean_manhattan_comparison.ipynb](experiment_2/01_euclidean_manhattan_comparison.ipynb) | Selected-setting quality/runtime tables |

Experiment 2 additionally evaluates Euclidean DBSCAN/OPTICS and retains diagnostic
outputs. The paper's four-method view uses Manhattan CLASSIX, Euclidean CLASSIX,
Manhattan DBSCAN and Manhattan OPTICS. The report-only helper
[`render_numeric_tables.py`](../scripts/render_numeric_tables.py) extracts that
view from a completed notebook run; it performs no clustering.

The historical probability notebook is separate from the four self-contained
reruns. Its unused clustering/native imports were removed, output paths work
from the repository root or `experiments/`, and future Monte Carlo runs use seed
101. Its generator, probability formulas, numeric grids and stored outputs are
preserved. The original simulation seed was not recorded, so the historical
Monte Carlo noise is not promised to recur exactly. Only setup, definitions, a
tiny simulation and temporary plotting were checked; the full probability
calculation and simulation grids were not rerun for this release. Install
`numpy`, `matplotlib` and `tqdm` in the Jupyter environment for this notebook.

## Running

Use Windows x64, Python 3.12, Jupyter, and Visual Studio C++ Build Tools with the
x64 toolset and Windows SDK. Run a notebook from its containing directory and
read its setup/configuration cells before continuing. The Experiment 1 notebooks
create isolated environments, install pinned packages and compile their embedded
native kernels. Initial setup and the numeric dataset downloads require internet
access. Experiment 2 creates its own pinned scientific environment.

Start with the documented `pilot` or `check` mode where supplied. Molecular and
radius notebooks default to a full run; a full DBSCAN/Taylor–Butina comparison can
take a long time. Experiment 2 performs thousands of tuning fits. Keep the
20 GiB worker limits and serial execution policies specified in the notebooks.
Do not infer a full-study result from a successful pilot.

Each notebook documents its working directory and resume behaviour. Choose a new
work directory or `RUN_NAME` for an independent run, and do not mix pilot and full
records. Generated environments, data, assignments and result logs are excluded
from the repository. They are recreated by the notebooks.

For the numeric manuscript tables, after the complete density comparison and its
verification have finished, run from the repository root:

```powershell
python scripts/render_numeric_tables.py --run-dir experiments/experiment_2/_generated/reproduction_density_101/density_baselines
```

This writes `tables/manuscript/quality_table.tex`, `runtime_table.tex` and
`four_method_results.csv` under that generated run. Change the path if you chose
a different run name. It reads fresh selected results and checks their agreement
with the generated full-precision report and five-fit timing medians.

## Data, settings and scoring

Synthetic clustering datasets use seed 101 and the embedded original generator.
Its `n` parameter is points **per cluster**: ten clusters with `n=1000` give
10,000 total points. Seed bits are Bernoulli draws; the per-bit flip probability
is `flip_prob * prescribed_popcount / dimension`, not a fixed 30% bit-flip rate.
The notebooks preserve the generator's shuffle stream.

The sample-count study uses 10,000–90,000 total rows, dimension 1000,
CLASSIX radius 0.4, `minPts=40`, `mergeScale=1/sqrt(radius)`, DBSCAN `eps=0.4`,
`min_samples=5`, and Taylor–Butina distance cutoff 0.35. Popcount and dimension
studies use 10,000 total rows; their CLASSIX `minPts=5`, with radius 0.3 and 0.4,
respectively. The separate Figure 7 study uses 10,000 rows, dimension **1024**,
six radii 0.3–0.8, CLASSIX `minPts=1`, DBSCAN `min_samples=1`, and Taylor–Butina
**distance cutoff equal to radius**. Both BitBIRCH implementations use similarity
`1-radius`, branching factor 50 and diameter merging. Complete per-case settings
and repetition schedules are embedded in the notebooks.

Landrum uses the supplied `data/BL_Sets_data_rdkitblog.csv`, verified by hash,
filtered by `pKi > 4`: 90,498 rows and 35 reference target categories. Morgan
fingerprints have radius 2 and 1024 bits. CLASSIX uses `(radius, minPts,
mergeScale)=(0.3,100,1.33)`, DBSCAN `(eps,min_samples)=(0.3,50)`, Taylor–Butina
distance cutoff 0.3 without reordering, and both BitBIRCH implementations
similarity `9/44`. BitBIRCH assignments and internal scores have **no**
post-clustering minimum-size filter. All-row ARI/AMI include noise as one predicted
category; silhouette and CHI use non-noise rows. Silhouette uses at most 2000
seeded query points against all retained reference rows.

The twenty numeric datasets are independent real datasets, not randomly
generated data. All methods receive shared feature-wise `StandardScaler`
preprocessing, with identifier and label columns excluded as specified by each
loader. Radius/epsilon calibration is metric-specific. ARI selects settings and
AMI scores the same assignments; reported timing medians use five fresh fits
after a warmup. Seeds control stochastic steps, not the identity of real data.

Timing boundaries differ between studies and are recorded in their protocols.
The molecular DBSCAN includes exact packed-bit radius-graph construction before
scikit-learn DBSCAN, and Taylor–Butina uses the documented exact streaming
adapter. These should not be described as timings of an unmodified dense-input
baseline. Hardware and implementation differences affect timings; fresh runs
need not reproduce the original wall-clock values exactly.

## Manuscript reconciliation notes

The code preserves the experiments that produced the recorded plots and tables.
The accompanying manuscript currently contains differing descriptions:

- Its synthetic defaults state CLASSIX `minPts=1` and DBSCAN `eps=0.5`; the
  sample-count experiment used 40 and 0.4.
- Its Figure 7 description states dimension 1000 and Taylor–Butina cutoff
  `1-radius`; the producing study used 1024 and distance cutoff `radius`.
- Its Landrum description mentions a BitBIRCH size cutoff of 100; no such filter
  produced the reported table.
- The blanket perfect-recovery statements for the popcount/dimension studies
  are not supported for every method by the saved assignments and scores.

These notes identify description differences; they do not represent new
experiments or replacements for the reported measurements. Use the notebook
protocols when reproducing the recorded results.
