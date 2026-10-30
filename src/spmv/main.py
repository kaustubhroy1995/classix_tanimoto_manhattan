"""Native packed-bit popcounts and binary submatrix-vector products."""
from ctypes import CDLL, c_int
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent
_packed = CDLL(str(_ROOT / "spsubmatxvec.dll"))
_int32_1d = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
_uint8_1d = np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1, flags="C_CONTIGUOUS")
_uint8_2d = np.ctypeslib.ndpointer(dtype=np.uint8, ndim=2, flags="C_CONTIGUOUS")
_packed.hardware_popcount_available.argtypes = []
_packed.hardware_popcount_available.restype = c_int
_packed.row_popcounts.argtypes = [_uint8_2d, c_int, c_int, _int32_1d]
_packed.row_popcounts.restype = None
_packed.packed_submatxvec.argtypes = [
    _uint8_2d, c_int, c_int, c_int, _uint8_1d, _int32_1d
]
_packed.packed_submatxvec.restype = None
HARDWARE_POPCOUNT = bool(_packed.hardware_popcount_available())
_INT_MAX = np.iinfo(np.int32).max


def _array(value, dtype, ndim, name):
    if not isinstance(value, np.ndarray):
        raise TypeError(f"{name} must be a numpy array")
    if value.dtype != dtype or value.ndim != ndim or not value.flags.c_contiguous:
        raise TypeError(f"{name} must be a C-contiguous {ndim}D {np.dtype(dtype)} array")
    return value


def _bounds(start, end, nrows):
    if not isinstance(start, (int, np.integer)) or not isinstance(end, (int, np.integer)):
        raise TypeError("row bounds must be integers")
    if not 0 <= start <= end <= nrows <= _INT_MAX:
        raise ValueError("invalid row bounds or too many rows for native int32 API")


def _packed_shape(data):
    _array(data, np.uint8, 2, "data")
    nrows, nbytes = data.shape
    if nrows > _INT_MAX or nbytes > _INT_MAX // 8:
        raise ValueError("packed shape exceeds native int32 API")
    return nrows, nbytes


def row_popcounts(data):
    nrows, nbytes = _packed_shape(data)
    output = np.empty(nrows, dtype=np.int32)
    _packed.row_popcounts(data, nrows, nbytes, output)
    return output


def packed_submatxvec(data, start, end, vector):
    nrows, nbytes = _packed_shape(data)
    _bounds(start, end, nrows)
    _array(vector, np.uint8, 1, "vector")
    if vector.size != nbytes:
        raise ValueError("packed vector length differs from packed matrix row width")
    output = np.empty(end - start, dtype=np.int32)
    _packed.packed_submatxvec(data, nbytes, start, end, vector, output)
    return output
