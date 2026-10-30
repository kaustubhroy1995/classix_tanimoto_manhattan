"""Validation after measurement, using bounded temporary buffers for hashes."""
import hashlib
import json
import numpy as np


def array_hash(array):
    array = np.ascontiguousarray(array)
    return hashlib.sha256(memoryview(array).cast("B")).hexdigest()


def logical_hash(data, n_features, packing=None):
    result = hashlib.sha256()
    for start in range(0, len(data), 256):
        block = data[start:start+256]
        if packing:
            block = np.unpackbits(block, axis=1, count=n_features, bitorder=packing)
        block = np.ascontiguousarray(block, dtype=np.int32)
        result.update(memoryview(block).cast("B"))
    return result.hexdigest()


def state_fingerprint(model, n_features, packed=False):
    fields = ("labels", "aggregation_labels", "splist", "group_sizes", "group_labels",
              "initial_merging_labels", "ind", "unsort_ind", "Adj")
    state = {}
    for name in fields:
        value = getattr(model, name)
        if name == "Adj" and hasattr(value, "iter_row_blocks"):
            digest = hashlib.sha256()
            for block in value.iter_row_blocks(rows=max(1, min(128, (8*1024*1024)//max(1, value.shape[1])))):
                digest.update(memoryview(np.ascontiguousarray(block)).cast("B"))
            state[name] = {"shape": list(value.shape), "dtype": "int8", "sha256": digest.hexdigest()}
            continue
        array = np.asarray(value)
        state[name] = {"shape": list(array.shape), "dtype": str(array.dtype), "sha256": array_hash(array)}
    state["logical_centers"] = logical_hash(model.spdata, n_features, "little" if packed else None)
    state["aggregation_dict"] = hashlib.sha256(json.dumps(model.aggregation_dict, sort_keys=True, default=lambda v: int(v)).encode()).hexdigest()
    return state
