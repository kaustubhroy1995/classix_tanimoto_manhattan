"""Packed binary operations used by CLASSIX_T."""

from .main import HARDWARE_POPCOUNT, packed_submatxvec, row_popcounts

spsubmatxvec = packed_submatxvec

__all__ = ["HARDWARE_POPCOUNT", "packed_submatxvec", "spsubmatxvec", "row_popcounts"]
