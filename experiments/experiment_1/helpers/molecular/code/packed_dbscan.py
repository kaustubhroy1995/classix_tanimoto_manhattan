"""Exact binary Jaccard radius graph followed by scikit-learn DBSCAN.

This is a distance-construction adapter, not a replacement clustering algorithm.
Every ordered pair is evaluated using the existing CLASSIX packed popcount
kernel. No approximate neighbors, data-dependent pruning, or tuning is used.
Graph construction and scikit-learn fitting both occur inside ``fit`` and must
both be included in reported clustering time. Input packing is separate.

The repository ``spmv`` module supplies the packed C kernel.
Input is a C-contiguous two-dimensional uint8 packed-bit array. Padding bits
must be zero. Bit order does not affect binary Jaccard distances.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"

from numbers import Integral, Real
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.cluster import DBSCAN

_KERNEL_DIRECTORY = REPO / "src"
if str(_KERNEL_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(_KERNEL_DIRECTORY))
from spmv import packed_submatxvec, row_popcounts


def _validate_input(packed, eps):
    if (not isinstance(packed, np.ndarray) or packed.ndim != 2
            or packed.dtype != np.uint8 or not packed.flags.c_contiguous):
        raise TypeError("packed must be a C-contiguous 2D uint8 ndarray")
    if not packed.shape[0] or not packed.shape[1]:
        raise ValueError("packed must contain at least one row and one byte")
    if not isinstance(eps, Real) or not np.isfinite(eps) or eps <= 0:
        raise ValueError("eps must be a finite positive distance")


def build_radius_graph(packed, eps, *, progress=None, progress_every=1000,
                       max_graph_bytes=None):
    """Return a complete symmetric CSR radius graph, including zero edges.

    Distances use ``mismatch / union``, exactly the binary SciPy Jaccard
    formula. Computing ``1 - intersection / union`` instead can change an
    edge at a floating-point threshold. Empty/empty fingerprints have distance
    zero. No sparse-zero elimination is performed: duplicate fingerprints and
    the explicit diagonal must remain neighbors.

    Memory is O(n + number_of_radius_edges); there is no dense n-by-n array.
    ``max_graph_bytes`` caps final CSR array storage, not total process memory.
    Construction temporarily retains row arrays and final CSR arrays together,
    and scikit-learn makes additional copies. A process-level memory limit is
    therefore still needed for a hard overall bound.
    """
    _validate_input(packed, eps)
    if not isinstance(progress_every, Integral) or progress_every < 1:
        raise ValueError("progress_every must be a positive integer")
    if max_graph_bytes is not None and (
            not isinstance(max_graph_bytes, Integral) or max_graph_bytes < 1):
        raise ValueError("max_graph_bytes must be a positive integer or None")

    n = len(packed)
    # int64 addition prevents overflow even at the native API's upper width.
    populations = row_popcounts(packed).astype(np.int64)
    row_indices, row_distances = [], []
    indptr = np.empty(n + 1, dtype=np.int32)
    indptr[0] = 0
    nnz = 0
    max_index = np.iinfo(np.int32).max
    for i in range(n):
        intersections = packed_submatxvec(packed, 0, n, packed[i])
        union = populations + populations[i] - intersections
        mismatch = union - intersections
        distances = np.zeros(n, dtype=np.float64)
        np.divide(mismatch, union, out=distances, where=union != 0)
        indices = np.flatnonzero(distances <= eps).astype(np.int32)
        values = distances[indices]
        # scikit-learn expects precomputed sparse rows sorted by distance.
        # Stable sorting preserves original row order for equal distances.
        order = np.argsort(values, kind="stable")
        indices, values = indices[order], values[order]
        nnz += len(indices)
        graph_bytes = int(nnz * 12 + (n + 1) * 4)
        if nnz > max_index:
            raise MemoryError("Radius graph exceeds int32 CSR edge capacity")
        if max_graph_bytes is not None and graph_bytes > max_graph_bytes:
            raise MemoryError(
                f"Radius graph needs at least {graph_bytes} bytes; "
                f"configured CSR cap is {max_graph_bytes} bytes")
        indptr[i + 1] = nnz
        row_indices.append(indices)
        row_distances.append(values)
        if progress is not None and ((i + 1) % progress_every == 0 or i + 1 == n):
            progress({"stage": "radius_graph", "rows_complete": i + 1,
                      "rows_total": n, "stored_edges": int(nnz),
                      "csr_bytes": graph_bytes})

    indices = np.concatenate(row_indices)
    values = np.concatenate(row_distances)
    del row_indices, row_distances
    graph = csr_matrix((values, indices, indptr), shape=(n, n), copy=False)
    stats = {"rows": n, "ordered_pair_evaluations": n * n,
             "stored_edges_including_diagonal": int(graph.nnz),
             "explicit_zero_edges": int(np.count_nonzero(graph.data == 0)),
             "csr_bytes": int(graph.data.nbytes + graph.indices.nbytes
                              + graph.indptr.nbytes),
             "distance_formula": "(pop_i + pop_j - 2*intersection) / union",
             "zero_union_distance": 0.0,
             "approximate_neighbors": False}
    return graph, stats


class PackedDBSCAN:
    """Exact binary Jaccard DBSCAN with packed C distance construction.

    ``metric='jaccard'`` describes the supplied packed features. The unchanged
    scikit-learn clustering stage receives their exact sparse radius graph with
    ``metric='precomputed'``. Only the single-threaded auto algorithm is used
    in this benchmark. ``labels_`` and ``core_sample_indices_`` are unmodified
    scikit-learn outputs. ``graph_stats_`` records construction diagnostics.
    """

    def __init__(self, eps=0.3, min_samples=50, metric="jaccard",
                 algorithm="auto", n_jobs=1, *, progress=None,
                 progress_every=1000, max_graph_bytes=None):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.algorithm = algorithm
        self.n_jobs = n_jobs
        self.progress = progress
        self.progress_every = progress_every
        self.max_graph_bytes = max_graph_bytes

    def fit(self, packed, y=None):
        _validate_input(packed, self.eps)
        if (not isinstance(self.min_samples, Integral)
                or isinstance(self.min_samples, bool) or self.min_samples < 1):
            raise ValueError("min_samples must be a positive integer")
        if self.metric != "jaccard" or self.algorithm != "auto" or self.n_jobs != 1:
            raise ValueError("This benchmark requires jaccard, auto, n_jobs=1")
        started = perf_counter()
        graph, self.graph_stats_ = build_radius_graph(
            packed, self.eps, progress=self.progress,
            progress_every=self.progress_every,
            max_graph_bytes=self.max_graph_bytes)
        graph_finished = perf_counter()
        self.model_ = DBSCAN(eps=self.eps, min_samples=self.min_samples,
                             metric="precomputed", algorithm="auto", n_jobs=1)
        self.model_.fit(graph)
        self.labels_ = self.model_.labels_
        self.core_sample_indices_ = self.model_.core_sample_indices_
        ended = perf_counter()
        self.graph_stats_.update({"graph_seconds": graph_finished - started,
                                  "sklearn_seconds": ended - graph_finished})
        if self.progress is not None:
            self.progress({"stage": "complete", **self.graph_stats_})
        return self

    def fit_predict(self, packed, y=None):
        return self.fit(packed, y).labels_
