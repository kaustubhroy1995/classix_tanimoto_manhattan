"""Serial native packed Tanimoto/Jaccard distances and silhouette reduction.

Use prepare_packed once per dataset, prepare_labels once per partition, then
distance_block and reduce_block on at most 128 query indices per block. Every
retained reference row participates. Do not mutate prepared arrays afterwards.
No query sampling or clustering occurs in this module.
"""
from __future__ import annotations

import os
import ctypes as ct
from dataclasses import dataclass
import hashlib
from pathlib import Path
import platform

import numpy as np

ROOT = Path(os.environ["CLASSIX_EXPERIMENT_WORK"])
DLL_PATH = ROOT / "native" / ("silhouette_kernel.dll" if platform.system() == "Windows" else "libsilhouette_kernel.so")
MAX_QUERIES = 128
_library = None


def _native():
    global _library
    if _library is None:
        library = ct.CDLL(str(DLL_PATH))
        u8 = ct.POINTER(ct.c_uint8)
        i32 = ct.POINTER(ct.c_int32)
        i64 = ct.POINTER(ct.c_int64)
        f64 = ct.POINTER(ct.c_double)
        library.packed_row_popcounts.argtypes = [u8, ct.c_int64, ct.c_int64, i64]
        library.packed_row_popcounts.restype = ct.c_int
        library.tanimoto_distance_block.argtypes = [u8, i64, ct.c_int64, ct.c_int64, i64, ct.c_int64, f64]
        library.tanimoto_distance_block.restype = ct.c_int
        library.reduce_silhouette_block.argtypes = [f64, ct.c_int64, ct.c_int64, i32, i64, ct.c_int64, i64, f64, f64, f64]
        library.reduce_silhouette_block.restype = ct.c_int
        library.silhouette_has_hardware_popcount.argtypes = []
        library.silhouette_has_hardware_popcount.restype = ct.c_int
        _library = library
    return _library


def _pointer(array, ctype):
    return array.ctypes.data_as(ct.POINTER(ctype))


def _check(code, operation):
    if code:
        raise RuntimeError(f"Native {operation} failed with status {code}")


def _indices(indices, n):
    value = np.asarray(indices)
    if value.ndim != 1 or value.dtype.kind not in "iu":
        raise TypeError("Query indices must be a one-dimensional integer array")
    if len(value) > MAX_QUERIES:
        raise ValueError(f"At most {MAX_QUERIES} query rows are allowed per distance/reduction block")
    if len(value) and (np.any(value < 0) or np.any(value >= n)):
        raise IndexError("Query index outside the reference rows")
    return np.ascontiguousarray(value, dtype=np.int64)


def validate_packed(packed, n_features=None):
    """Validate packing once; optional true feature count checks zero tail bits."""
    value = np.asarray(packed)
    if value.dtype != np.uint8 or value.ndim != 2 or not value.flags.c_contiguous:
        raise TypeError("packed must be a C-contiguous two-dimensional uint8 array")
    if not all(value.shape):
        raise ValueError("packed must contain at least one row and one byte per row")
    if n_features is not None:
        if not isinstance(n_features, (int, np.integer)) or n_features < 1:
            raise ValueError("n_features must be a positive integer")
        if value.shape[1] != (int(n_features) + 7) // 8:
            raise ValueError("Packed byte width does not match n_features")
        used = int(n_features) % 8
        if used and np.any(value[:, -1] & np.uint8(255 ^ ((1 << used) - 1))):
            raise ValueError("Nonzero tail padding: expected little-endian packbits and zero unused bits")
    return value


@dataclass(frozen=True)
class PackedReference:
    packed: np.ndarray
    popcounts: np.ndarray
    n_features: int | None

    def distance_block(self, query_indices):
        return distance_block(self, query_indices)


@dataclass(frozen=True)
class LabelState:
    labels: np.ndarray
    counts: np.ndarray


def prepare_packed(packed, n_features=None):
    """Validate once and cache exact row popcounts; no full dataset copy."""
    value = validate_packed(packed, n_features)
    counts = np.empty(value.shape[0], dtype=np.int64)
    _check(_native().packed_row_popcounts(_pointer(value, ct.c_uint8), *value.shape,
                                         _pointer(counts, ct.c_int64)), "row popcounts")
    view = value.view()
    view.flags.writeable = False
    counts.flags.writeable = False
    return PackedReference(view, counts, int(n_features) if n_features is not None else None)


def prepare_labels(labels, counts=None):
    """Check compact 0..K-1 cluster IDs and -1 noise once per partition."""
    value = np.asarray(labels)
    if value.ndim != 1 or value.dtype.kind not in "iu" or not len(value):
        raise TypeError("labels must be a nonempty one-dimensional integer array")
    if np.any(value < -1) or np.any(value > np.iinfo(np.int32).max):
        raise ValueError("Labels must be -1 or nonnegative int32 cluster IDs")
    value = np.ascontiguousarray(value, dtype=np.int32)
    retained = value[value >= 0]
    if not len(retained):
        raise ValueError("No retained cluster members")
    unique, actual = np.unique(retained, return_counts=True)
    if not np.array_equal(unique, np.arange(len(unique), dtype=np.int32)):
        raise ValueError("Non-noise labels must be compact 0..K-1")
    actual = np.ascontiguousarray(actual, dtype=np.int64)
    if len(actual) < 2:
        raise ValueError("Silhouette requires at least two retained clusters")
    if counts is not None:
        supplied = np.asarray(counts)
        if supplied.ndim != 1 or supplied.dtype.kind not in "iu" or not np.array_equal(supplied, actual):
            raise ValueError("Cluster counts do not match retained reference labels")
    view = value.view()
    view.flags.writeable = False
    actual.flags.writeable = False
    return LabelState(view, actual)


def distance_block(packed, query_indices, *, n_features=None):
    """Return exact float64 Jaccard distances MÃ—N, including zero/zero = 0.

    Pass PackedReference to avoid revalidating packing and recomputing popcounts.
    Query indices may be out of order or repeated and refer to reference rows.
    """
    reference = packed if isinstance(packed, PackedReference) else prepare_packed(packed, n_features)
    if n_features is not None and reference.n_features != int(n_features):
        raise ValueError("Prepared feature count differs; prepare with the correct n_features first")
    n, width = reference.packed.shape
    queries = _indices(query_indices, n)
    output = np.empty((len(queries), n), dtype=np.float64)
    if len(queries):
        _check(_native().tanimoto_distance_block(_pointer(reference.packed, ct.c_uint8),
            _pointer(reference.popcounts, ct.c_int64), n, width, _pointer(queries, ct.c_int64),
            len(queries), _pointer(output, ct.c_double)), "distance block")
    return output


def reduce_block(distances, state, query_indices, *, check_distances=True):
    """Return (silhouette, a, b) against all retained reference rows.

    state must come from prepare_labels. check_distances=False skips a full D
    scan and is intended only for untouched distance_block output (or row subsets
    thereof) whose validity is already established. Shape/dtype/query/noise checks
    always run. Input arrays must remain unchanged after their preparation.
    """
    if not isinstance(state, LabelState):
        raise TypeError("state must be produced by prepare_labels")
    value = np.asarray(distances)
    queries = _indices(query_indices, len(state.labels))
    if value.dtype != np.float64 or value.shape != (len(queries), len(state.labels)) or not value.flags.c_contiguous:
        raise TypeError("distances must be C-contiguous float64 with shape (len(queries), len(labels))")
    if np.any(state.labels[queries] < 0):
        raise ValueError("Noise query rows must be removed before silhouette reduction")
    if check_distances and (not np.isfinite(value).all() or np.any(value < 0) or np.any(value > 1)):
        raise ValueError("Expected finite Tanimoto distances in [0,1]")
    silhouette, within, nearest = [np.empty(len(queries), dtype=np.float64) for _ in range(3)]
    if len(queries):
        _check(_native().reduce_silhouette_block(_pointer(value, ct.c_double), len(queries), len(state.labels),
            _pointer(state.labels, ct.c_int32), _pointer(state.counts, ct.c_int64), len(state.counts),
            _pointer(queries, ct.c_int64), _pointer(silhouette, ct.c_double), _pointer(within, ct.c_double),
            _pointer(nearest, ct.c_double)), "silhouette reduction")
    return silhouette, within, nearest


def silhouette_reduce(distances, labels, counts, query_indices, *, check_distances=True):
    """Convenience checked API; prepare_labels/reduce_block amortize validation."""
    state = labels if isinstance(labels, LabelState) else prepare_labels(labels, counts)
    if isinstance(labels, LabelState) and counts is not None and not np.array_equal(counts, labels.counts):
        raise ValueError("Supplied counts differ from prepared counts")
    return reduce_block(distances, state, query_indices, check_distances=check_distances)


def native_metadata():
    source = Path(__file__).with_name("silhouette_kernel.c")
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    return {"dll": str(DLL_PATH), "dll_sha256": sha(DLL_PATH), "c_source_sha256": sha(source),
            "wrapper_sha256": sha(Path(__file__)), "hardware_popcount": bool(_native().silhouette_has_hardware_popcount()),
            "threading": "serial C loops; no OpenMP or BLAS", "max_queries_per_block": MAX_QUERIES,
            "zero_union_distance": 0.0, "singleton_silhouette": 0.0,
            "reference_policy": "all non-noise reference members; self excluded by source row index"}
