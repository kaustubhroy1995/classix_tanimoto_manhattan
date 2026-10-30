"""Manhattan and compiled Euclidean CLASSIX adapters for Experiment 2.

This is a comparison of complete implementations, not a metric-only ablation:
normalization, ordering and native/Python execution paths also differ. The local
Manhattan uses the canonical repository source.
Progress rendering and standard output are suppressed without changing fits.
"""
from __future__ import annotations

import contextlib
import hashlib
import importlib
import importlib.metadata
import importlib.util
import inspect
import json
import os
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import platform
import sys
from time import perf_counter
import warnings

ROOT = RUN
SOURCE = MANHATTAN
_manhattan_module = None


def _method_name(method):
    aliases = {"manhattan": "manhattan", "classix_m": "manhattan", "euclidean": "euclidean",
               "classix_euclidean": "euclidean"}
    try:
        return aliases[str(method).lower()]
    except KeyError:
        raise ValueError("method must be 'manhattan' or 'euclidean'") from None


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _quiet_progress(iterable, *args, **kwargs):
    """Return the same iterable, removing only tqdm display bookkeeping."""
    return iterable


def _manhattan():
    global _manhattan_module
    if _manhattan_module is None:
        provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
        if _sha(SOURCE) != provenance["compatible_sha256"]:
            raise RuntimeError("Manhattan source differs from the run provenance")
        spec = importlib.util.spec_from_file_location("experiment_2_classix_m", SOURCE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.tqdm = _quiet_progress
        _manhattan_module = module
    return _manhattan_module.CLASSIX_M


def _euclidean():
    from classix import CLASSIX
    return CLASSIX


def preprocessing_scale(method, X):
    """Return the scalar divisor actually used by the selected fit path.

    Manhattan's zero median is returned as zero (the original has no fallback).
    Euclidean/PCA replaces a zero median centred L2 norm with one.
    """
    import numpy as np
    method = _method_name(method)
    data = np.asarray(X, dtype=np.float64)
    if data.ndim != 2 or not all(data.shape):
        raise ValueError("Expected a nonempty two-dimensional input")
    if method == "manhattan":
        return float(np.median(np.sum(data - data.min(axis=0), axis=1)))
    scale = float(np.median(np.linalg.norm(data - data.mean(axis=0), axis=1)))
    return 1.0 if scale == 0 else scale


def fit_model(method, X, radius, min_pts, merge_scale):
    """Fit once and return labels plus separately observed fit metadata.

    The external perf_counter interval encloses only model.fit(data). Estimator
    construction, a fresh C-contiguous float64 input copy, imports, thread-limit
    setup, scale reporting and label extraction are outside that interval.
    Internal centering/scaling/sorting/aggregation/merging/reallocation remain
    inside fit. The caller should use the same frozen input for both methods.
    """
    import numpy as np
    from threadpoolctl import threadpool_limits
    method = _method_name(method)
    data = np.array(X, dtype=np.float64, order="C", copy=True)
    if data.ndim != 2 or not all(data.shape) or not np.isfinite(data).all():
        raise ValueError("Expected a finite, nonempty two-dimensional input")
    cls = _manhattan() if method == "manhattan" else _euclidean()
    kwargs = dict(radius=float(radius), minPts=int(min_pts), mergeScale=float(merge_scale),
                  post_alloc=True, mergeTinyGroups=True, verbose=0)
    kwargs.update(sorting="popcount", group_merging="manhattan_distance") if method == "manhattan" else kwargs.update(
        sorting="pca", metric="euclidean", group_merging="distance")
    model = cls(**kwargs)
    with threadpool_limits(limits=1):
        scale = preprocessing_scale(method, data)
        with open(os.devnull, "w", encoding="utf-8") as sink, contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink), warnings.catch_warnings(record=True) as observed:
            warnings.simplefilter("always")
            started = perf_counter()
            model.fit(data)
            elapsed = perf_counter() - started
        labels = np.asarray(model.labels if method == "manhattan" else model.labels_, dtype=np.int64).copy()
    if labels.shape != (len(data),) or np.any(labels < -1):
        raise RuntimeError("Estimator returned invalid assignments")
    nonnoise = labels[labels >= 0]
    return {"labels": labels, "fit_seconds": float(elapsed),
            "n_groups": len(model.splist if method == "manhattan" else model.splist_),
            "n_clusters": int(np.unique(nonnoise).size), "noise_points": int(np.count_nonzero(labels == -1)),
            "preprocessing_scale": scale,
            "fitted_scale": float(model.dataScale_) if method == "euclidean" else scale,
            "warnings": [{"category": item.category.__name__, "message": str(item.message)} for item in observed]}


def native_details():
    """Import/construct only; report actual callable modules without fitting."""
    import numpy as np
    import scipy
    cls = _euclidean()
    model = cls(metric="euclidean", sorting="pca", group_merging="distance", verbose=0,
                mergeTinyGroups=True, post_alloc=True)
    kernels = {}
    for role, attribute in (("aggregation", "_aggregate"), ("distance_merging", "_distance_merge")):
        function = getattr(model, attribute)
        module = importlib.import_module(function.__module__)
        path = Path(module.__file__).resolve()
        kernels[role] = {"function": function.__name__, "module": function.__module__, "file": str(path),
                         "native": path.suffix.lower() in (".pyd", ".so"), "sha256": _sha(path)}
    _manhattan()
    return {"python": sys.version, "executable": sys.executable, "platform": platform.platform(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "classixclustering": importlib.metadata.version("classixclustering"),
            "euclidean_signature": str(inspect.signature(cls)), "euclidean_kernels": kernels,
            "all_euclidean_kernels_native": all(item["native"] for item in kernels.values()),
            "numpy_row_stack_available": hasattr(np, "row_stack"),
            "manhattan_source": json.loads(PROVENANCE.read_text(encoding="utf-8")),
            "methods_source_sha256": _sha(__file__),
            "timing_scope": "Only model.fit; own preprocessing included; construction, copy, label extraction excluded",
            "instrumentation": "threadpool_limits(1); stdout/stderr redirected to null; Manhattan tqdm replaced by identity iterable; warnings retained",
            "comparison_scope": "Complete implementations differ in metric, normalization, sorting and execution kernels"}
