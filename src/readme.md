# Source

- `classix_t.py`: `CLASSIX_T`, using packed binary rows and native Tanimoto operations.
- `spmv/`: packed C operations, ctypes interface and compact merging adjacency.
- `classix_m.py`: `CLASSIX_M`, for real-valued Manhattan clustering.

Build the Tanimoto library with `scripts/build_native.ps1` from the repository root. Import algorithms with `src` on the Python path.
