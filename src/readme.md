# Source implementations

- `packed_tanimoto/`: the packed binary Tanimoto implementation used by the revised benchmarks, its compact adjacency representation and C kernels. Build it using `scripts/build_native.ps1` from the repository; see `packed_tanimoto/README.md` and the root README for its API.
- `classix_m.py`: the measured Manhattan implementation, with the NumPy stacking compatibility change.
- `classix_t.py` and `spmv/`: the earlier unpacked Tanimoto implementation and Unix-oriented CSR kernel. These are retained separately for historical compatibility.
- `setup.py`: the original spmv packaging file; it does not install the new packed implementation.

The curated reproduction notebooks freeze their benchmark source snapshots and pinned dependencies. They generate and compile their isolated run helpers rather than depending on locally modified algorithm code.
