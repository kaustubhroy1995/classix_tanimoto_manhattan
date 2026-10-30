"""Seeded serial adapter for the canonical clustering implementations."""
import numpy as np
from methods import fit_model as _fit, preprocessing_scale, native_details

def fit_model(method, X, radius, min_pts, merge_scale):
    state = np.random.get_state()
    np.random.seed(101)
    try:
        return _fit(method, X, radius, min_pts, merge_scale)
    finally:
        np.random.set_state(state)
