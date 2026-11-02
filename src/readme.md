# Source

- `classix_t.py`: `CLASSIX_T`, Tanimoto operations, supports binary vector inputs.
- `spmv/`: C operations, ctypes interface for sparse submatrix - dense vector operations and adjacency computation.
- `classix_m.py`: `CLASSIX_M`, for Manhattan clustering.

Build the Tanimoto library with `scripts/build_native.ps1` from the repository root. Import algorithms with `src` on the Python path.
