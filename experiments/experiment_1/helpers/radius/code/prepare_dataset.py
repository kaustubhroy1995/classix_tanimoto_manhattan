"""Generate one shared dataset from the notebook function, outside fit workers."""
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


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_recorded_data(directory, metadata, reference, case):
    """Check fresh or cached arrays against the frozen raw-array hashes."""
    require(metadata["data_id"] == reference["data_id"] == case["data_id"], "Dataset identity differs from recorded data")
    require(metadata["parameters"] == reference["parameters"] == case["generator"], "Generator parameters differ from recorded data")
    require(metadata["shape"] == reference["shape"] == [case["rows"], case["features"]], "Dataset shape differs from recorded data")
    require(metadata["generator_snapshot_sha256"] == reference["generator_snapshot_sha256"]
            == hashlib.sha256((COMMON/"generate_data_snapshot.py").read_bytes()).hexdigest(),
            "Generator source differs from the recorded snapshot")
    require(metadata["source_notebook_sha256"] == reference["source_notebook_sha256"], "Dataset source provenance differs")
    names = {"int32", "int64", "packed_little", "packed_big", "truth", "permutation", "centers"}
    require(set(metadata["files"]) == set(reference["files"]) == names, "Dataset file inventory differs")
    n, d, k = case["rows"], case["features"], case["generator"]["num_clusters"]
    shapes = {"int32": (n, d), "int64": (n, d), "packed_little": (n, (d+7)//8),
              "packed_big": (n, (d+7)//8), "truth": (n,), "permutation": (n,), "centers": (k, d)}
    for name in sorted(names):
        record, expected = metadata["files"][name], reference["files"][name]
        require(all(record[key] == expected[key] for key in ("file", "dtype", "bytes", "sha256")),
                f"Manifest differs from recorded array: {name}")
        require(Path(record["file"]).name == record["file"], "Array filename must be local")
        array = np.load(directory/record["file"], allow_pickle=False, mmap_mode="r")
        require(array.shape == shapes[name] and str(array.dtype) == expected["dtype"]
                and array.nbytes == expected["bytes"], f"Actual array shape/dtype/bytes differ: {name}")
        require(digest(array) == expected["sha256"], f"Actual raw-array SHA256 differs from recorded data: {name}")
        del array
    return len(names)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_id")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    plan = read_configuration(args.config)
    case = next(c for c in plan["cases"] if c["data_id"] == args.data_id)
    # New 10k datasets have no historical 25k reference hashes.
    directory = args.data_dir.resolve() / args.data_id
    if not directory.is_relative_to((ROOT/"runs").resolve()):
        raise ValueError("Generated data must stay under experiment_1A/runs")
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory/"manifest.json"
    if manifest_path.exists():
        prior = json.loads(manifest_path.read_text())
        verified = verify_recorded_data(directory, prior, prior, case)
        print(json.dumps({"reusing": args.data_id, "fresh_raw_array_hashes_verified": verified}), flush=True)
        return
    params = case["generator"]
    start = time.perf_counter()
    data, truth = generate_data(**params)
    centers = data[np.arange(params["num_clusters"])*params["n"]].copy()
    permutation = np.arange(len(data))
    np.random.shuffle(permutation)
    data, truth = np.ascontiguousarray(data[permutation]), truth[permutation]
    assert data.shape == (case["rows"], case["features"])
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
        packed = np.packbits(data, axis=1, bitorder=order)
        save("packed_"+order, packed)
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
    verified = verify_recorded_data(directory, metadata, metadata, case)
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(metadata, indent=2)+"\n")
    temporary.replace(manifest_path)
    print(json.dumps({"generated": args.data_id, "shape": metadata["shape"], "zeros": metadata["zero_fraction"],
                      "fresh_raw_array_hashes_verified": verified}), flush=True)


if __name__ == "__main__":
    main()
