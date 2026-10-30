# Packed Tanimoto CLASSIX

This directory contains the packed implementation measured in Experiment 1.
The algorithm and C sources retain their benchmarked numerical behavior.
`src/classix_t.py` is the earlier, unpacked implementation.

The implementation performs intersections and popcounts on packed fingerprints
in C. Its group graph stores one connection bit per entry and sparse markers
for minPts edges; together these preserve the logical adjacency values 0, 1, 2.
The connection buffer occupies `G * ceil(G / 8)` bytes for `G` aggregation
groups, plus the separate sparse-marker objects.

From the repository root on Windows x64:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_native.ps1
python scripts/smoke_packed.py
```

Building requires Visual Studio C++ Build Tools and a Windows SDK. The script
uses an x64 developer environment or discovers one with `vswhere.exe`. It
compiles both C sources, records the compiler, flags and hashes under
`_build/native/`, and installs the two generated DLLs beside `packed_native.py`.
Generated binaries are local build products, not repository inputs. Both DLLs
are needed because the existing wrapper loads both kernels at import time.
The current native implementation targets Windows x64; other platforms have
not been validated.

The Python runtime dependencies are NumPy, SciPy and tqdm. Use the notebook's
own pinned environment when reproducing a reported benchmark.

```python
import sys
import numpy as np

sys.path.insert(0, "src")
from packed_tanimoto import CLASSIX_T_PACKED

binary_rows = np.array([[1, 0, 1], [1, 1, 1], [0, 1, 0]], dtype=np.uint8)
packed_rows = np.packbits(binary_rows, axis=1, bitorder="little")
model = CLASSIX_T_PACKED(radius=0.3, minPts=1, mergeScale=1.4)
model.fit(packed_rows)
labels = model.labels
```

Input is a two-dimensional `uint8` array of packed rows. Padding bits in the
last byte must be zero; `np.packbits` supplies this padding. Use fingerprints
with at least one active bit: the benchmarked algorithm does not define a
special empty-fingerprint case. The implemented cutoff expressions require
`0 < radius < 1` and `0 < radius * mergeScale < 1`.

The smoke check uses 60 seed-101 observations, checks native intersections and
popcounts against NumPy, exercises the original CSR kernel and adjacency
markers, and compares the clustering partition with known separated groups.
It is a correctness check, not a performance benchmark.
