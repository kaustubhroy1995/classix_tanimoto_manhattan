"""Small correctness check for the published packed kernels and clustering API."""

import contextlib
import io
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from packed_tanimoto import CLASSIX_T_PACKED
from packed_tanimoto.packed_adjacency import PackedAdjacency
from packed_tanimoto.packed_native import (
    HARDWARE_POPCOUNT,
    original_spsubmatxvec,
    packed_submatxvec,
    row_popcounts,
    validate_original_csr,
)


def main():
    rng = np.random.default_rng(101)
    rows = np.zeros((60, 130), dtype=np.uint8)
    truth = np.repeat(np.arange(3), 20)
    for cluster in range(3):
        for row in range(cluster * 20, (cluster + 1) * 20):
            rows[row, cluster * 24:cluster * 24 + 12] = 1
            rows[row, rng.choice(np.arange(72, 130), size=2, replace=False)] = 1
    order = rng.permutation(len(rows))
    rows, truth = rows[order], truth[order]
    packed = np.packbits(rows, axis=1, bitorder="little")
    np.testing.assert_array_equal(row_popcounts(packed), rows.sum(axis=1))
    integer_rows = rows.astype(np.int32)
    for anchor in (0, 17, 59):
        expected = integer_rows[1:59] @ integer_rows[anchor]
        np.testing.assert_array_equal(
            packed_submatxvec(packed, 1, 59, packed[anchor]), expected
        )
    csr = sparse.csr_matrix(integer_rows)
    validate_original_csr(csr.data, csr.indptr, csr.indices, rows.shape[1])
    np.testing.assert_array_equal(
        original_spsubmatxvec(csr.data, csr.indptr, csr.indices, 1, 59, integer_rows[17]),
        integer_rows[1:59] @ integer_rows[17],
    )

    graph = PackedAdjacency(10)
    graph[0, [1, 2, 2, 9]] = 1
    graph[0, 9] = 2
    graph[9, 0] = 2
    assert graph[0, 1] == 1 and graph[0, 9] == graph[9, 0] == 2
    assert graph.marker_count == 2 and graph.nbytes == 20
    graph[0, 9] = 0
    assert graph[0, 9] == 0 and graph.marker_count == 1

    model = CLASSIX_T_PACKED(radius=0.3, minPts=1, mergeScale=1.4)
    # The research implementation prints timings/progress; keep the smoke output concise.
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        model.fit(packed)
    np.testing.assert_array_equal(
        model.labels[:, None] == model.labels[None, :],
        truth[:, None] == truth[None, :],
    )
    assert len(np.unique(model.labels)) == 3
    groups = len(model.splist)
    assert model.Adj.nbytes == groups * ((groups + 7) // 8)
    print(json.dumps({
        "status": "passed", "seed": 101, "rows": 60, "features": 130,
        "clusters": 3, "aggregation_groups": groups,
        "hardware_popcount": HARDWARE_POPCOUNT,
        "checks": ["row_popcounts", "packed_intersections", "original_csr_kernel",
                   "adjacency_markers", "known_cluster_partition"],
    }, indent=2))


if __name__ == "__main__":
    main()
