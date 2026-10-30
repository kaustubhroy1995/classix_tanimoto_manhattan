"""Generate shared data afresh; validate cached arrays against this run's hashes."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import argparse
import gc
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from generate_data_snapshot import generate_data
HERE = Path(__file__).resolve().parent
ROOT = work_root(SOURCE_ROOT)

def digest(array):
    result = hashlib.sha256()
    row_bytes = array.dtype.itemsize * int(np.prod(array.shape[1:], dtype=np.int64))
    rows = max(1, (8 * 1024 * 1024) // max(1, row_bytes))
    for start in range(0, len(array), rows):
        result.update(memoryview(np.ascontiguousarray(array[start:start + rows])).cast("B"))
    return result.hexdigest()

def validate(directory, metadata, case):
    assert metadata["data_id"] == case["data_id"]
    assert metadata["parameters"] == case["generator"]
    assert metadata["shape"] == [case["rows"], case["features"]]
    assert metadata["generator_snapshot_sha256"] == hashlib.sha256((COMMON/"generate_data_snapshot.py").read_bytes()).hexdigest()
    n, d, k = case["rows"], case["features"], case["generator"]["num_clusters"]
    shapes = {"int32": (n,d), "int64": (n,d), "packed_little": (n,(d+7)//8),
              "packed_big": (n,(d+7)//8), "truth": (n,), "permutation": (n,), "centers": (k,d)}
    assert set(metadata["files"]) == set(shapes)
    for name, shape in shapes.items():
        rec = metadata["files"][name]
        assert Path(rec["file"]).name == rec["file"]
        array = np.load(directory/rec["file"], allow_pickle=False, mmap_mode="r")
        assert array.shape == shape and str(array.dtype) == rec["dtype"] and array.nbytes == rec["bytes"]
        assert digest(array) == rec["sha256"], f"Cached array changed: {name}"
    return len(shapes)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_id")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    plan = read_configuration(args.config)
    case = next(c for c in plan["cases"] if c["data_id"] == args.data_id)
    directory = args.data_dir.resolve()/args.data_id
    assert directory.is_relative_to((ROOT/"runs").resolve())
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory/"manifest.json"
    if manifest_path.exists():
        metadata = json.loads(manifest_path.read_text())
        print(json.dumps({"reusing": args.data_id, "array_hashes_verified": validate(directory, metadata, case)}), flush=True)
        return
    params = case["generator"]
    start = time.perf_counter()
    data, truth = generate_data(**params)
    centers = data[np.arange(params["num_clusters"])*params["n"]].copy()
    permutation = np.arange(len(data))
    np.random.shuffle(permutation)
    data, truth = np.ascontiguousarray(data[permutation]), truth[permutation]
    assert data.shape == (case["rows"], case["features"])
    assert np.all((data == 0) | (data == 1))
    generation_seconds = time.perf_counter()-start
    files = {}
    def save(name, array):
        np.save(directory/f"{name}.npy", array)
        files[name] = {"file": name+".npy", "bytes": array.nbytes,
                       "dtype": str(array.dtype), "sha256": digest(array)}
    save("int32", data)
    wide = data.astype(np.int64)
    save("int64", wide)
    del wide
    gc.collect()
    for order in ("little", "big"):
        save("packed_"+order, np.packbits(data, axis=1, bitorder=order))
    save("truth", truth)
    save("permutation", permutation)
    save("centers", centers)
    row_pops = data.sum(axis=1)
    metadata = {"data_id": args.data_id, "parameters": params, "shape": list(data.shape),
        "files": files, "zero_fraction": float(1-data.mean()), "zero_rows": int(np.sum(row_pops == 0)),
        "minimum_popcount": int(row_pops.min()), "maximum_popcount": int(row_pops.max()),
        "center_popcounts": centers.sum(axis=1).tolist(), "generation_seconds": generation_seconds,
        "source_notebook_sha256": plan["source_notebook_sha256"],
        "generator_snapshot_sha256": hashlib.sha256((COMMON/"generate_data_snapshot.py").read_bytes()).hexdigest(),
        "shuffle": "np.random.shuffle after generator, continuing its seeded NumPy RNG"}
    assert metadata["zero_rows"] == 0, "Tanimoto protocol requires nonempty vectors."
    validate(directory, metadata, case)
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(metadata, indent=2)+"\n")
    temporary.replace(manifest_path)
    print(json.dumps({"generated": args.data_id, "shape": metadata["shape"]}), flush=True)

if __name__ == "__main__":
    main()
