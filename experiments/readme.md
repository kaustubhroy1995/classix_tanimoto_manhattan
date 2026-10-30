# Experiments

Run the cells in order in these notebooks:

- [Synthetic scaling](experiment_1/01_synthetic_benchmarks.ipynb)
- [Greg Landrum molecular data](experiment_1/02_molecular_benchmarks.ipynb)
- [Radius sensitivity](experiment_1/03_radius_sensitivity.ipynb)
- [Euclidean and Manhattan comparison](experiment_2/01_euclidean_manhattan_comparison.ipynb)
- [Search-bound probabilities](classix_t_prob.ipynb)

Setup: [repository README](../README.md). Parameters and scoring: [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

The notebooks call ordinary Python helpers, the canonical `src/classix_t.py` and `src/classix_m.py`, and the compiled Euclidean CLASSIX package. Generated data, environments, labels and run records are ignored by Git. The supplied plots retain the recorded paper results; reruns write fresh results.
