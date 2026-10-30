"""Adjusted Rand index from the contingency table, using Python integers."""
import numpy as np


def adjusted_rand_index(truth, labels):
    _, a = np.unique(truth, return_inverse=True)
    _, b = np.unique(labels, return_inverse=True)
    table = np.zeros((int(a.max())+1, int(b.max())+1), dtype=np.int64)
    np.add.at(table, (a, b), 1)
    def pairs(counts):
        return sum(int(x)*(int(x)-1)//2 for x in counts.flat)
    total = len(a)*(len(a)-1)//2
    cells, rows, cols = pairs(table), pairs(table.sum(axis=1)), pairs(table.sum(axis=0))
    denominator = total*(rows+cols)-2*rows*cols
    return 1.0 if denominator == 0 else 2*(total*cells-rows*cols)/denominator
