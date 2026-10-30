"""Generate fresh Morgan inputs from hash-verified original text sources.

Full mode preserves original source order, duplicate molecules and the Landrum
pKi>4 filter. Pilot mode uses a declared prefix and cannot become a full run.
No previous prepared arrays, manifests, results or experiment folders are read.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

ROOT = work_root(SOURCE_ROOT)
SPECS = {
    "landrum": ("landrum_pki_gt4_morgan1024", "BL_Sets_data_rdkitblog.csv", 90498,
                "ddc8dbf5d92b1f8255ef210456c008dd5bf91866bc9f93870592ebcb5868ba1d"),
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def array_sha(array):
    h = hashlib.sha256()
    for start in range(0, len(array), 1024):
        h.update(memoryview(array[start:start+1024]).cast("B"))
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=SPECS, required=True)
    parser.add_argument("--pilot-rows", type=int)
    args = parser.parse_args()
    assert not list((ROOT / "runs").rglob("run.lock")), "Preparation must precede timed stages"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    import numpy as np
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    data_id, filename, expected_rows, expected_hash = SPECS[args.dataset]
    source = ROOT / "data" / "source" / filename
    assert sha(source) == expected_hash, "Original text input hash differs"
    assert args.pilot_rows is None or 0 < args.pilot_rows <= expected_rows
    started = time.perf_counter()
    smiles, row_ids, targets = [], [], []
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        assert {"", "smiles", "target_chembl_id", "label", "pKi"} <= set(reader.fieldnames)
        for index, row in enumerate(reader):
            pki = float(row["pKi"])
            assert np.isfinite(pki)
            # Preserve retained source order and duplicate rows.
            # Invalid retained molecules stop preparation rather than being dropped.
            if pki > 4:
                smiles.append(row["smiles"])
                targets.append(row["target_chembl_id"])
                row_ids.append(index)
        assert index + 1 == 91663
    target_ids = sorted(set(targets))
    assert len(target_ids) == 35 and len(smiles) == expected_rows
    mapping = {target: index for index, target in enumerate(target_ids)}
    truth = np.asarray([mapping[target] for target in targets], dtype=np.int32)
    count = args.pilot_rows or expected_rows
    directory = ROOT / "data" / "prepared" / data_id
    manifest_path = directory / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        assert manifest["shape"] == [count, 1024] and manifest["source"]["sha256"] == expected_hash
        assert manifest["pilot_rows"] == args.pilot_rows, "Use different work directories for pilot/full mode"
        for info in manifest["files"].values():
            array = np.load(directory / info["file"], mmap_mode="r", allow_pickle=False)
            assert list(array.shape) == info["shape"] and str(array.dtype) == info["dtype"]
            assert array.nbytes == info["bytes"] and array_sha(array) == info["sha256"]
        print(json.dumps({"dataset": args.dataset, "reused_verified_arrays": True, "rows": count}), flush=True)
        return
    assert not directory.exists() or not any(directory.iterdir()), "Incomplete preparation requires inspection"
    directory.mkdir(parents=True, exist_ok=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    data = np.empty((count, 1024), dtype=np.uint8)
    for index, text in enumerate(smiles[:count]):
        molecule = Chem.MolFromSmiles(text.strip())
        assert molecule is not None, f"Invalid retained SMILES at row {row_ids[index]}"
        data[index] = generator.GetFingerprintAsNumPy(molecule)
        if (index + 1) % 10000 == 0:
            print(json.dumps({"dataset": args.dataset, "fingerprints": index+1, "total": count}), flush=True)
    assert np.all((data == 0) | (data == 1))
    files = {}
    def save_array(name, array):
        path = directory / (name + ".npy")
        np.save(path, array, allow_pickle=False)
        files[name] = {"file": path.name, "shape": list(array.shape), "dtype": str(array.dtype),
                       "bytes": int(array.nbytes), "sha256": array_sha(array)}
    save_array("int32", data.astype(np.int32))
    save_array("int64", data.astype(np.int64))
    for order in ("little", "big"):
        packed = np.packbits(data, axis=1, bitorder=order)
        assert np.array_equal(np.unpackbits(packed, axis=1, bitorder=order), data)
        save_array("packed_"+order, packed)
    save_array("source_row_ids", np.asarray(row_ids[:count], dtype=np.int64))
    if truth is not None:
        save_array("truth", truth[:count])
    pops = data.sum(axis=1)
    assert not np.any(pops == 0), "The current study requires nonempty fingerprints"
    parameters = {"fingerprint": {"method": "rdFingerprintGenerator.GetMorganGenerator", "radius": 2,
                                  "fpSize": 1024, "other_options": "RDKit defaults"},
                  "order": "original source order; duplicates retained", "seed": 101,
                  "seed_note": "Morgan generation is deterministic; this seed is used for fit reproducibility"}
    if args.dataset == "landrum":
        parameters["filter"] = {"column": "pKi", "operator": ">", "threshold": 4.0}
    manifest = {"data_id": data_id, "shape": [count, 1024], "files": files, "parameters": parameters,
                "source": {"file": "data/source/"+filename, "sha256": expected_hash},
                "zero_rows": 0, "zero_fraction": float(1-data.mean()), "pilot_rows": args.pilot_rows,
                "target_label_mapping": mapping, "generator_info": generator.GetInfoString(),
                "prepared_utc": datetime.now(timezone.utc).isoformat(), "preparation_seconds": time.perf_counter()-started,
                "preparation_source_sha256": sha(Path(__file__)),
                "environment": {"python": sys.version, "numpy": np.__version__, "rdkit": rdkit.__version__, "platform": platform.platform()}}
    manifest_path.write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"dataset": args.dataset, "rows": count, "prepared_manifest_sha256": sha(manifest_path)}), flush=True)


if __name__ == "__main__":
    main()
