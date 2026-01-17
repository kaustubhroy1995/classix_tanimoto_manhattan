#!/usr/bin/env python
import platform
from ctypes import CDLL, c_int
import numpy as np
import pkg_resources
from scipy import sparse

sysname = platform.system()

lib_name = "spsubmatxvec.so"
lib_path_submat = pkg_resources.resource_filename('spmv', 'lib/{}'.format(lib_name))
testlib_submat = CDLL(lib_path_submat)
c_int_array = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
testlib_submat.multiply.argtypes = [c_int_array, c_int_array, c_int_array, c_int, c_int,
                                    c_int_array, c_int_array]

def spsubmatxvec(a_data, a_indptr, a_indices, start_row, end_row, V):
    O = np.zeros(end_row-start_row).astype(np.int32)
    status = testlib_submat.multiply(a_data, a_indptr, a_indices, start_row, end_row, V, O)
    return O