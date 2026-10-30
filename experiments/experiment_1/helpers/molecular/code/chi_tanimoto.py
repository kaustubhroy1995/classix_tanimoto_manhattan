"""Tanimoto CHI only, preserving the audited binary-majority formula.

Extracted from experiment_1B/code/internal_quality_metrics.py. No DBI, Dunn,
pairwise cluster matrix, clustering or parameter selection is performed.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import numpy as np

DEFINITION_VERSION = 'tanimoto_quality_v1_1_all_clusters'


def score(value=None, reason=None):
    if reason is not None or value is None or not np.isfinite(value):
        return {'status': 'undefined', 'value': None, 'reason': reason or 'nonfinite_result'}
    return {'status': 'ok', 'value': float(value), 'reason': None}


def divide(numerator, denominator):
    a, b = np.broadcast_arrays(numerator, denominator)
    return np.divide(a, b, out=np.full(a.shape, np.nan, dtype=np.float64), where=b != 0)


def evaluate(X, labels, min_sizes=(0,), pairwise=False, block_size=128, checkpoint=None):
    """Return CHI and population metadata for each exclusive minimum size.

    Use min_sizes=(0,) to retain every non-noise cluster, including singletons.
    Representatives use binary majority, with half ties rounded to one. Label
    -1 is excluded as noise; all rows remain in the coverage denominator.
    Undefined zero-union/zero-within-dispersion cases preserve audited reasons.
    """
    if pairwise:
        raise ValueError('This module evaluates CHI only; pairwise must be False')
    X, labels = np.asarray(X), np.asarray(labels)
    if X.ndim != 2 or labels.shape != (len(X),) or labels.dtype.kind not in 'iu':
        raise ValueError('Expected a binary matrix and one integer label per row')
    if not len(X) or not X.shape[1] or np.any(labels < -1):
        raise ValueError('Empty input or invalid negative labels')
    if not isinstance(block_size, int) or block_size < 1:
        raise ValueError('block_size must be a positive integer')
    if not min_sizes or any(type(v) is not int or v < 0 for v in min_sizes):
        raise ValueError('Exclusive minimum sizes must be nonnegative integers')
    for start in range(0, len(X), 4096):
        block = X[start:start + 4096]
        if np.any((block != 0) & (block != 1)):
            raise ValueError('Fingerprints must be binary')
    used = np.flatnonzero(labels >= 0)
    unique, inverse, counts = np.unique(labels[used], return_inverse=True, return_counts=True)
    k_total, n_total = len(unique), len(X)
    n_noise = n_total - len(used)
    singleton_fraction = float(np.count_nonzero(counts == 1) / n_total)
    order = np.argsort(inverse, kind='stable')
    starts = np.r_[0, np.cumsum(counts)[:-1]] if k_total else np.array([], dtype=int)
    all_sums = (np.add.reduceat(X[used[order]], starts, axis=0, dtype=np.int64)
                if k_total else np.empty((0, X.shape[1]), dtype=np.int64))
    keep = counts > min(min_sizes)
    selected = np.flatnonzero(keep)
    remap = np.full(k_total, -1, dtype=int)
    remap[selected] = np.arange(len(selected))
    rows = used[keep[inverse]]
    groups = remap[inverse[keep[inverse]]]
    sizes, sums = counts[keep].astype(np.float64), all_sums[keep].astype(np.float64)
    del all_sums, order
    k = len(sizes)
    if n_total * n_total * X.shape[1] >= 2**53:
        raise ValueError('Input exceeds exact-integer range of float64 aggregate products')
    centres = (2 * sums >= sizes[:, None]).astype(np.float64)
    pops = centres.sum(axis=1)
    scatter_sq, bad_scatter = np.zeros(k), np.zeros(k, dtype=bool)
    for start in range(0, len(rows), 2048):
        if checkpoint:
            checkpoint()
        idx = groups[start:start + 2048]
        data = X[rows[start:start + 2048]]
        intersection = np.einsum('ij,ij->i', data, centres[idx], dtype=np.float64)
        distance = 1 - divide(intersection, data.sum(axis=1) + pops[idx] - intersection)
        valid = np.isfinite(distance)
        scatter_sq += np.bincount(idx[valid], weights=distance[valid]**2, minlength=k)
        bad_scatter[idx[~valid]] = True
    result = []
    for minimum in min_sizes:
        mask = sizes > minimum
        n, retained_k = int(sizes[mask].sum()), int(mask.sum())
        entry = {'min_size_exclusive': minimum, 'n_total': n_total, 'n_noise': n_noise,
                 'k_total': k_total, 'n_retained': n, 'k_retained': retained_k,
                 'coverage': n / n_total, 'singleton_fraction': singleton_fraction,
                 'singleton_cluster_fraction': float(np.count_nonzero(counts == 1) / k_total) if k_total else None,
                 'zero_centroids': int(np.count_nonzero(pops[mask] == 0)),
                 'definition_version': DEFINITION_VERSION}
        reason = 'fewer_than_two_retained_clusters' if retained_k < 2 else None
        entry['chi_tanimoto'] = score(reason=reason or 'not_evaluated')
        if not reason:
            global_centre = (2 * sums[mask].sum(axis=0) >= n).astype(np.float64)
            inter = centres[mask] @ global_centre
            between_dist = 1 - divide(inter, pops[mask] + global_centre.sum() - inter)
            wcss = float(scatter_sq[mask].sum())
            bcss = float(np.dot(sizes[mask], between_dist**2))
            if bad_scatter[mask].any() or not np.isfinite(between_dist).all():
                entry['chi_tanimoto'] = score(reason='zero_union_tanimoto')
            elif wcss == 0:
                entry['chi_tanimoto'] = score(reason='zero_within_dispersion')
            else:
                entry['chi_tanimoto'] = score(bcss * (n - retained_k) / (wcss * (retained_k - 1)))
            entry.update(within_squared_distance_sum=wcss,
                         between_squared_distance_sum=bcss if np.isfinite(bcss) else None)
        result.append(entry)
    return result
