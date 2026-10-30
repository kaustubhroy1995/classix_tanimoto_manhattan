# CLASSIX in the Manhattan and Tanimoto distances

Code supporting **Fast and explainable clustering in the Manhattan and Tanimoto distance**.

This revision includes packed binary Tanimoto clustering with native C popcount/intersection operations, compact group adjacency, the measured Manhattan implementation, and full reproduction notebooks.

## Packed Tanimoto: Windows setup

The native implementation is tested on **Windows x64 with Python 3.12**. Install Visual Studio C++ Build Tools (or Visual Studio with the C++ workload) and a Windows SDK. Linux/macOS support for this packed kernel has not been established; the earlier implementation is retained separately.

From the repository root in PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-packed.txt
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build_native.ps1
.\.venv\Scripts\python.exe .\scripts\smoke_packed.py
```

The build script locates an x64 MSVC toolchain, builds both native libraries and installs them beside the Python wrapper. Generated DLLs and build records are ignored by Git. Build outputs remain under `_build/native/`.

### Python API

```python
import sys
sys.path.insert(0, "src")

import numpy as np
from packed_tanimoto import CLASSIX_T_PACKED

binary_rows = np.array([
    [1, 1, 0, 0],
    [1, 1, 0, 0],
    [0, 0, 1, 1],
    [0, 0, 1, 1],
], dtype=np.uint8)
packed_rows = np.packbits(binary_rows, axis=1, bitorder="little")

model = CLASSIX_T_PACKED(radius=0.3, minPts=1, mergeScale=1.4)
model.fit(packed_rows)
print(model.labels)
```

Input is a contiguous two-dimensional `uint8` array of little-endian packed binary rows. Unused padding bits must be zero. `np.packbits` provides the required padding. The public class exposes labels as `model.labels`.

The adjacency implementation preserves logical values 0, 1 and 2 using one packed connection bitmatrix and sparse markers for `minPts` reassignment edges. Storage remains dependent on the number of aggregation groups; packing reduces the storage per connection.

The measured Manhattan implementation is available from `src/classix_m.py` as `CLASSIX_M`. Its only compatibility change is the use of `np.vstack` in place of `np.row_stack`.

## Reproduce the experiments

Install a notebook frontend in the Python 3.12 environment:

```powershell
.\.venv\Scripts\python.exe -m pip install jupyterlab matplotlib
.\.venv\Scripts\python.exe -m jupyter lab
```

Select that environment as the notebook kernel. Run notebooks from their own directories, or follow their repository-root setup instructions. The four curated notebooks contain the helper Python/C sources, fixed experimental designs and pinned environments needed for complete reruns:

| Notebook | Outputs |
| --- | --- |
| [`experiments/experiment_1/01_synthetic_benchmarks.ipynb`](experiments/experiment_1/01_synthetic_benchmarks.ipynb) | Synthetic scaling studies and Figures 4-6 |
| [`experiments/experiment_1/02_molecular_benchmarks.ipynb`](experiments/experiment_1/02_molecular_benchmarks.ipynb) | Five-method Greg Landrum comparison and Table 2 |
| [`experiments/experiment_1/03_radius_sensitivity.ipynb`](experiments/experiment_1/03_radius_sensitivity.ipynb) | Radius ARI/AMI study and Figure 7 |
| [`experiments/experiment_2/01_euclidean_manhattan_comparison.ipynb`](experiments/experiment_2/01_euclidean_manhattan_comparison.ipynb) | Twenty original-feature datasets, CLASSIX/DBSCAN/OPTICS comparisons and Tables 3-4 |

**Start with the notebook's `RUN_MODE="pilot"` or `"check"` where available.** Full DBSCAN and Taylor-Butina experiments are expensive. The full synthetic notebook generates substantial data and has no fit-time cutoff; numerical work is not performed merely by opening the notebook. The notebooks explain their full-run defaults and work directories.

The binary notebooks require Windows x64, Python 3.12, C++ build tools and package downloads. They create isolated, pinned run environments rather than using the root environment for measured fits. The numeric notebook also creates its own pinned environment. Generated code, datasets, labels, native binaries and records remain in ignored work directories; they are not bundled as historical run data.

The source chemical dataset is already in `data/BL_Sets_data_rdkitblog.csv`. Its LF line endings are preserved by `.gitattributes` because the molecular preparation verifies the source bytes.

See [`experiments/REPRODUCIBILITY.md`](experiments/REPRODUCIBILITY.md) for exact executed settings, timing/scoring definitions and distinctions between the manuscript descriptions and the frozen experiment protocols. Reproducers retain the executed settings; they do not silently retune them to match different prose.

### Generate the paper's numeric tables

After completing the numeric notebook, use its generated `density_baselines` result directory with:

```powershell
.\.venv\Scripts\python.exe .\scripts\render_numeric_tables.py --help
```

The helper formats the four variants reported in the paper from the complete six-variant study, without rerunning clustering or editing the manuscript.

## Earlier source and notebooks

The original `src/classix_t.py` is retained as the earlier **unpacked** CSR-based implementation. It loads `src/spmv/lib/spsubmatxvec.so` and has its own Unix-oriented build instructions in `src/spmv/`. It is separate from the packed implementation above.

The six notebooks directly under `experiments/` are historical studies. `classix_t_prob.ipynb` produces the probability-analysis figures; its setup no longer requires the unused clustering/native imports and seeded Monte Carlo reruns use 101. Historical stored simulation outputs did not record their seed, so exact historical sampling noise is not promised. The remaining original notebooks are not the reproduction entry points for the revised benchmark tables.

`requirements.txt` retains dependencies for the earlier studies. Use the isolated environments and instructions in the curated notebooks for the revised benchmarks.

## License and third-party methods

The project source is distributed under the root [MIT license](LICENSE). The packed source directory retains its license notice. BitBIRCH source embedded in the reproduction notebooks retains its supplied upstream notices; BBLean and the other comparison libraries are installed separately under their own licenses. Their notices are not replaced by this project's MIT license.
