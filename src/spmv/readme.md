# Native Tanimoto operations

`spsubmatxvec.c` computes row popcounts and packed binary submatrix-vector intersections. Inputs are contiguous `uint8` arrays produced by `np.packbits(rows, axis=1, bitorder="little")`; padding bits must be zero. Counts use `int32` outputs.

From the repository root on Windows x64:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_native.ps1
python scripts/smoke_classix.py
```

The build requires MSVC and a Windows SDK, and installs `spsubmatxvec.dll` here. `main.py` supplies the ctypes interface. `packed_adjacency.py` stores connections as bits with sparse markers for minPts edges.
