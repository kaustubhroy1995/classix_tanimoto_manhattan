"""Memory-bounded, exact static-order Butina using RDKit Tanimoto kernels.

This is an equivalent adapter, not RDKit's unchanged ClusterData function.
It preserves ClusterData(..., isDistData=True, reordering=False) cluster tuple
order for the same lower-triangle distances computed as 1 - RDKit similarity.
No threshold graph or quadratic distance array is stored. Auxiliary memory is
O(n), excluding the supplied fingerprints and returned cluster membership.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"

import math
from typing import Callable, Iterable

import numpy as np
from rdkit import DataStructs


def cluster_fingerprints(
    fingerprints: Iterable,
    dist_thresh: float,
    *,
    reordering: bool = False,
    progress: Callable[[dict], None] | None = None,
    progress_every: int = 1000,
) -> tuple[tuple[int, ...], ...]:
    """Return clusters with exactly RDKit's static neighbor-count tie ordering.

    ``dist_thresh`` is a Tanimoto distance in [0,1]. Reordering=True is deliberately
    unsupported. Progress callbacks receive small dictionaries after completed
    work units; an exception from the callback propagates (allowing cancellation).

    Work: n(n-1)/2 comparisons for fixed neighbor counts, followed by n comparisons
    per selected center whose original neighbor count exceeds one. Input row
    order is never changed. Neighbors are appended in ascending input index;
    candidate centers sort by descending (original neighbor count, input index).
    """
    if reordering is not False:
        raise ValueError('Only reordering=False is implemented and validated.')
    threshold = float(dist_thresh)
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('dist_thresh must be finite and in [0,1].')
    if not isinstance(progress_every, int) or progress_every < 1:
        raise ValueError('progress_every must be a positive integer.')
    fps = fingerprints if isinstance(fingerprints, list) else list(fingerprints)
    n = len(fps)

    def emit(stage, **values):
        if progress is not None:
            progress(dict(stage=stage, rows=n, **values))

    # Stock ClusterData sets diagonal distance to zero, including all-zero
    # fingerprints; count self explicitly rather than relying on self-similarity.
    counts = np.ones(n, dtype=np.int64)
    emit('count_neighbors', completed_rows=0, total_rows=max(0, n - 1))
    for i in range(1, n):
        similarities = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        distances = 1.0 - np.asarray(similarities, dtype=np.float64)
        near = distances <= threshold
        counts[:i] += near
        counts[i] += np.count_nonzero(near)
        if i % progress_every == 0 or i == n - 1:
            emit('count_neighbors', completed_rows=i, total_rows=n - 1,
                 evaluated_pairs=i * (i + 1) // 2)

    # The list scan is equivalent to stock sorted_indices.pop(0) when
    # reordering=False; no remaining-neighbor counts are recalculated.
    candidates = sorted(((int(count), index) for index, count in enumerate(counts)), reverse=True)
    seen = np.zeros(n, dtype=bool)
    clusters = []
    recomputed_centers = 0
    assigned = 0
    emit('assign_clusters', visited_candidates=0, assigned_points=0, clusters=0)
    for visited, (count, index) in enumerate(candidates, start=1):
        if not seen[index]:
            seen[index] = True
            if count <= 1:
                cluster = (index,)
            else:
                similarities = DataStructs.BulkTanimotoSimilarity(fps[index], fps)
                distances = 1.0 - np.asarray(similarities, dtype=np.float64)
                # flatnonzero has the same ascending order as stock np.where.
                neighbors = np.flatnonzero((distances <= threshold) & ~seen)
                seen[neighbors] = True
                cluster = (index, *neighbors.tolist())
                recomputed_centers += 1
            clusters.append(cluster)
            assigned += len(cluster)
        if visited % progress_every == 0 or visited == n:
            emit('assign_clusters', visited_candidates=visited, assigned_points=assigned,
                 clusters=len(clusters), recomputed_center_rows=recomputed_centers)
    assert assigned == n and bool(np.all(seen))
    emit('complete', clusters=len(clusters), assigned_points=assigned,
         lower_triangle_comparisons=n * (n - 1) // 2,
         recomputed_center_rows=recomputed_centers,
         recomputed_center_comparisons=n * recomputed_centers)
    return tuple(clusters)
