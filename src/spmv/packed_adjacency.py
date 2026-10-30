"""CLASSIX adjacency with packed connections and sparse minPts edge markers.

Each logical entry is 0, 1, or 2. One bit records whether an edge exists; sparse
row-to-column sets mark edges whose logical value is 2. The connection buffer is
exactly G * ceil(G / 8) bytes for a G by G matrix. Sparse marker overhead is
separate and exposed as a Python-object size estimate. No dense graph is
allocated by construction, reads of individual edges, or writes. Decoding a row
returns an independent int8 array, not a writable view.

Supported writes match CLASSIX's existing operations:
    graph[i, columns] = 1
    graph[rows, j] = 1
    graph[i, j] = 2
The same forms accept 0 to clear edges. Repeated indices and different bits in
one byte are handled correctly using unbuffered NumPy ufunc.at updates.

Call iter_row_blocks(rows=...) for bounded logical hashing or validation.
Implicit conversion of the whole graph to a dense NumPy array is prohibited.
"""
import operator
import sys

import numpy as np


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, not a boolean")
    try:
        return operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error


class PackedAdjacency:
    """Logical int8 adjacency: one bit/entry plus sparse markers for value 2."""

    __slots__ = ("_n", "_row_bytes", "_connections", "_minpts")
    dtype = np.dtype(np.int8)
    ndim = 2

    def __init__(self, groups):
        groups = _integer(groups, "groups")
        if groups < 0:
            raise ValueError("groups must be nonnegative")
        self._n = groups
        self._row_bytes = (groups + 7) // 8
        self._connections = np.zeros((groups, self._row_bytes), dtype=np.uint8)
        self._minpts = {}

    @property
    def shape(self):
        return (self._n, self._n)

    @property
    def size(self):
        return self._n * self._n

    @property
    def nbytes(self):
        """Exact connection-buffer bytes; sparse marker objects are separate."""
        return self._connections.nbytes

    @property
    def connection_bits(self):
        """Read-only packed rows marking every connection, including minPts edges."""
        view = self._connections.view()
        view.flags.writeable = False
        return view

    @property
    def marker_count(self):
        """Number of directed edges marked as minPts edges."""
        return sum(len(columns) for columns in self._minpts.values())

    @property
    def marker_storage_bytes(self):
        """Estimated sparse-marker Python object bytes, for reporting only.

        Counts the dictionary, sets, row/column integers once by object identity.
        It excludes allocator overhead and may include shared/interned integers;
        therefore this is an estimate, not measured process memory. Cost is O(M)
        in stored directed markers, so avoid querying inside the clustering loop.
        """
        seen = set()
        def size_once(value):
            identity = id(value)
            if identity in seen:
                return 0
            seen.add(identity)
            return sys.getsizeof(value)
        total = size_once(self._minpts)
        for row, columns in self._minpts.items():
            total += size_once(row) + size_once(columns)
            total += sum(size_once(column) for column in columns)
        return total

    @property
    def logical_nbytes(self):
        """Bytes the equivalent dense int8 matrix would require."""
        return self.size

    def __len__(self):
        return self._n

    def __repr__(self):
        return f"PackedAdjacency(shape={self.shape}, nbytes={self.nbytes}, marker_count={self.marker_count})"

    def __array__(self, dtype=None, copy=None):
        raise TypeError("Whole-graph conversion is disabled; use iter_row_blocks(rows=...) or graph[i]")

    def _index(self, index):
        index = _integer(index, "index")
        if index < 0:
            index += self._n
        if not 0 <= index < self._n:
            raise IndexError("adjacency index out of bounds")
        return index

    def _selection(self, selection):
        if isinstance(selection, (int, np.integer)) and not isinstance(selection, (bool, np.bool_)):
            return self._index(selection)
        values = np.asarray(selection)
        if values.ndim != 1:
            raise TypeError("indices must be a scalar integer or a one-dimensional integer array")
        if values.size == 0:
            return np.empty(0, dtype=np.intp)
        if values.dtype.kind not in "iu":
            raise TypeError("indices must have an integer dtype")
        # Check unsigned values before casting to np.intp to prevent wraparound.
        if np.any(values >= self._n) or (values.dtype.kind == "i" and np.any(values < -self._n)):
            raise IndexError("adjacency index out of bounds")
        values = values.astype(np.intp, copy=False)
        if np.any(values < 0):
            values = np.where(values < 0, values + self._n, values)
        return values

    def _pair(self, key):
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError("edge access requires exactly two indices")
        row, column = (self._selection(index) for index in key)
        if isinstance(row, np.ndarray) and isinstance(column, np.ndarray):
            raise TypeError("at most one index may be a vector")
        return row, column

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            row = self._index(key)
            decoded = np.unpackbits(self._connections[row], bitorder="little", count=self._n).view(np.int8)
            markers = self._minpts.get(row)
            if markers:
                decoded[np.fromiter(markers, dtype=np.intp, count=len(markers))] = 2
            return decoded
        row, column = self._pair(key)
        byte_index = column >> 3
        offset = column & 7
        result = (self._connections[row, byte_index] >> offset) & 1
        if np.ndim(result) == 0:
            return np.int8(2 if column in self._minpts.get(row, ()) else result)
        result = np.asarray(result, dtype=np.int8)
        if not self._minpts:
            return result
        if isinstance(row, int):
            markers = self._minpts.get(row)
            if markers:
                mask = np.fromiter((int(c) in markers for c in column), dtype=bool, count=len(column))
                result[mask] = 2
        else:
            mask = np.fromiter((column in self._minpts.get(int(r), ()) for r in row), dtype=bool, count=len(row))
            result[mask] = 2
        return result

    def _write_marker(self, row, column, marked):
        if marked:
            markers = self._minpts.get(row)
            if markers is None:
                self._minpts[row] = {column}
            else:
                markers.add(column)
        else:
            markers = self._minpts.get(row)
            if markers is not None:
                markers.discard(column)
                if not markers:
                    del self._minpts[row]

    def __setitem__(self, key, value):
        value = _integer(value, "edge value")
        if not 0 <= value <= 2:
            raise ValueError("edge values must be 0, 1, or 2")
        row, column = self._pair(key)
        if isinstance(row, int) and isinstance(column, int):
            byte_index = column >> 3
            mask = np.uint8(1 << (column & 7))
            clear_mask = np.bitwise_not(mask)
            current = self._connections[row, byte_index]
            self._connections[row, byte_index] = current | mask if value else current & clear_mask
            self._write_marker(row, column, value == 2)
            return
        if isinstance(column, np.ndarray):
            byte_indices = column >> 3
            masks = np.left_shift(np.uint8(1), (column & 7).astype(np.uint8))
            target = (row, byte_indices)
        else:
            masks = np.uint8(1 << (column & 7))
            target = (row, column >> 3)
        if value:
            np.bitwise_or.at(self._connections, target, masks)
        else:
            np.bitwise_and.at(self._connections, target, np.bitwise_not(masks))
        if isinstance(row, int):
            if value == 2 and column.size:
                markers = self._minpts.get(row)
                if markers is None:
                    self._minpts[row] = set(map(int, column))
                else:
                    markers.update(map(int, column))
            elif row in self._minpts:
                markers = self._minpts[row]
                markers.difference_update(map(int, column))
                if not markers:
                    del self._minpts[row]
        elif value == 2 or self._minpts:
            # Normal vector writes during merging incur no Python per-edge loop
            # while the marker dictionary is empty. minPts writes are scalar in
            # CLASSIX; support vector marker writes for a complete tested API.
            for item in map(int, row):
                self._write_marker(item, column, value == 2)

    def iter_row_blocks(self, rows=64):
        """Yield dense int8 row chunks with shape at most (rows, G).

        One decoded int8 chunk is allocated. Callers should process chunks sequentially rather than
        collecting them, to keep memory bounded.
        """
        rows = _integer(rows, "rows")
        if rows <= 0:
            raise ValueError("rows must be positive")
        for start in range(0, self._n, rows):
            stop = min(start + rows, self._n)
            decoded = np.unpackbits(self._connections[start:stop], axis=1,
                                    bitorder="little", count=self._n).view(np.int8)
            for row in range(start, stop):
                markers = self._minpts.get(row)
                if markers:
                    decoded[row - start, np.fromiter(markers, dtype=np.intp, count=len(markers))] = 2
            yield decoded
