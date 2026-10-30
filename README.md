# CLASSIX in the Manhattan and Tanimoto distances

Code for *Fast and explainable clustering in the Manhattan and Tanimoto distance*.

## Setup

Use Python 3.12 on Windows x64, with Visual Studio C++ Build Tools and a Windows SDK.
From the repository root in PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build_native.ps1
.\.venv\Scripts\python.exe .\scripts\smoke_classix.py
.\.venv\Scripts\python.exe -m jupyter lab
```

The build compiles the packed C operations in `src/spmv/`. Native binaries are generated locally. Linux/macOS builds have not been validated.

## Use

```python
import sys
sys.path.insert(0, "src")

import numpy as np
from classix_t import CLASSIX_T
from classix_m import CLASSIX_M

binary_rows = np.array([
    [1, 1, 0, 0],
    [1, 1, 0, 0],
    [0, 0, 1, 1],
    [0, 0, 1, 1],
], dtype=np.uint8)
packed_rows = np.packbits(binary_rows, axis=1, bitorder="little")

model = CLASSIX_T(radius=0.3, minPts=1, mergeScale=1.4)
model.fit(packed_rows)
print(model.labels)
```

`CLASSIX_T` uses packed binary fingerprints and native popcount/intersection operations. Pass nonempty fingerprints in a two-dimensional `uint8` array with little-endian bit packing and zero padding bits. The connection matrix uses packed bits, with sparse markers for special `minPts` edges. `CLASSIX_M` accepts floating-point feature matrices.

## Experiments

Open a notebook and run its cells in order. The setup cells install each study's recorded dependencies in an isolated environment. Full runs include lengthy DBSCAN and Taylor-Butina fits.

| Notebook | Results |
| --- | --- |
| [Synthetic benchmarks](experiments/experiment_1/01_synthetic_benchmarks.ipynb) | Sample count, popcount, cluster count and dimension scaling |
| [Greg Landrum benchmark](experiments/experiment_1/02_molecular_benchmarks.ipynb) | Five-method runtime, ARI, AMI, silhouette and CHI_T table |
| [Radius sensitivity](experiments/experiment_1/03_radius_sensitivity.ipynb) | Synthetic ARI and AMI curves |
| [Euclidean/Manhattan comparison](experiments/experiment_2/01_euclidean_manhattan_comparison.ipynb) | Twenty datasets; CLASSIX, DBSCAN and OPTICS |
| [Search-bound probabilities](experiments/classix_t_prob.ipynb) | Probability calculations and simulations |

Shared Python helpers are ordinary source files beside the notebooks. Tanimoto and Manhattan fits import the implementations in `src/`; Euclidean CLASSIX uses the compiled `classixclustering` package. Notebooks generate their datasets and results under ignored work directories; the Landrum source CSV is supplied in `data/`. See [experiment settings](experiments/REPRODUCIBILITY.md) for timing and scoring details.

## License

Project code is under the [MIT license](LICENSE). Vendored BitBIRCH retains its upstream license; installed comparison libraries retain their own licenses.
