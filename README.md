# CLASSIX in the Manhattan and Tanimoto distances

Code for *Fast and explainable clustering in the Manhattan and Tanimoto distance*.

## Setup

Windows x64, Python 3.12, Visual Studio C++ Build Tools and a Windows SDK:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_native.ps1
.\.venv\Scripts\python.exe -m jupyter lab
```

Open the binary experiment notebooks using this environment. Each notebook contains its data preparation, clustering calls, scoring and plotting code. Run its cells from top to bottom. It writes results under `_generated/` and figures under `plots/`.

| Notebook | Experiment |
| --- | --- |
| [Synthetic benchmarks](experiments/experiment_1/01_synthetic_benchmarks.ipynb) | Sample count, popcount, cluster count, dimension and radius scaling |
| [Greg Landrum](experiments/experiment_1/02_molecular_benchmarks.ipynb) | Five-method chemical fingerprint comparison |
| [Radius sensitivity](experiments/experiment_1/03_radius_sensitivity.ipynb) | Synthetic ARI and AMI curves |
| [Euclidean/Manhattan](experiments/experiment_2/01_euclidean_manhattan_comparison.ipynb) | Twenty original-feature datasets |
| [Search-bound probabilities](experiments/classix_t_prob.ipynb) | Probability calculations and simulations |

The numeric comparison uses the versions from its original runs. Open it in a separate environment:

```powershell
py -3.12 -m venv .venv-numeric
.\.venv-numeric\Scripts\python.exe -m pip install -r experiments/experiment_2/requirements.txt
.\.venv-numeric\Scripts\python.exe -m jupyter lab
```

Full benchmark runs, particularly DBSCAN and Taylor-Butina, can be lengthy. The notebooks time the same fitting operations as the recorded experiments, now directly in the notebook process. Fresh timings will vary. Supplied figures retain the reported measurements.

## Algorithms

```python
import sys
sys.path.insert(0, "src")
import numpy as np
from classix_t import CLASSIX_T

binary_rows = np.array([[1, 1, 0, 0], [1, 1, 0, 0],
                        [0, 0, 1, 1], [0, 0, 1, 1]], dtype=np.uint8)
packed = np.packbits(binary_rows, axis=1, bitorder="little")
model = CLASSIX_T(radius=0.3, minPts=1, mergeScale=1.4)
model.fit(packed)
labels = model.labels
```

`src/classix_t.py` is the packed Tanimoto implementation. Its C operations and compact adjacency are in `src/spmv/`. Input must be contiguous packed `uint8` rows with zero padding bits. `src/classix_m.py` supplies Manhattan CLASSIX. Euclidean CLASSIX comes from the compiled `classixclustering` package.

The comparison notebooks import the pinned third-party BitBIRCH library from `vendor/bitbirch/`. It retains its upstream license. Project code uses the [MIT license](LICENSE).
