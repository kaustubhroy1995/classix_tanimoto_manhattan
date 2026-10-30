"""Prepare the fixed 20-dataset Experiment 2B roster without clustering.

Primary inputs are the original numeric features with a shared StandardScaler.
All twenty datasets are reconstructed from public source files or pinned sklearn data.
No earlier experimental artifacts are read.
All current dataset source bytes and transformations are recorded separately.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.resources
import io
import json
import os
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import platform
import shutil
import sys
import time
import zipfile

ROOT = RUN
DATA = ROOT / "data"
ROSTER = {
 "iris": (150,4), "banknote": (1372,4), "wine": (178,13), "breast_cancer": (569,30),
 "digits": (1797,64), "seeds": (210,7), "glass": (214,9), "sonar": (208,60),
 "vehicle": (846,18), "segmentation": (2310,19), "breast_tissue": (106,9),
 "leaf": (340,14), "libras": (360,90), "vertebral": (310,6), "wireless": (2000,7),
 "parkinsons": (195,22), "coimbra": (116,9), "raisin": (900,7), "haberman": (306,3), "musk": (476,166)}
EXPECTED_CLASSES = dict(zip(ROSTER, (3,2,3,2,10,3,6,2,4,7,6,30,15,3,4,2,2,2,2,2)))
ORIGINAL_EIGHT = tuple(ROSTER)[:8]
NEW_UCI_IDS = {"vehicle":149, "segmentation":50, "breast_tissue":192, "leaf":288,
 "libras":181, "vertebral":212, "wireless":422, "parkinsons":174,
 "coimbra":451, "raisin":850, "haberman":43, "musk":74}
PACKAGED = {"wine": ("wine_data.csv",109), "breast_cancer": ("breast_cancer.csv",17), "digits": ("digits.csv.gz",80)}
UCI = {
 "seeds":{"id":236,"url":"https://archive.ics.uci.edu/ml/machine-learning-databases/00236/seeds_dataset.txt","features":["area","perimeter","compactness","kernel_length","kernel_width","asymmetry_coefficient","kernel_groove_length"]},
 "glass":{"id":42,"url":"https://archive.ics.uci.edu/ml/machine-learning-databases/glass/glass.data","features":["RI","Na","Mg","Al","Si","K","Ca","Ba","Fe"]},
 "sonar":{"id":151,"url":"https://archive.ics.uci.edu/ml/machine-learning-databases/undocumented/connectionist-bench/sonar/sonar.all-data","features":[f"Attribute{i}" for i in range(1,61)]}}

# Expected source-derived array hashes, not clustering results.
BASELINE_ARRAY_HASHES = {'iris': {'X_raw': '58280bb559506628655ae68c493688079f736a0f3142f8b526ee3983cae22030', 'X': 'ac959912d5365b262b1e09d3ae4df1e12804ceef971852204a8c210f3c263212', 'y': '734737f11162991213afc3a96aa5e8823b2d9cdcfb8992df1e85a0bae44b2395', 'row_ids': 'cfe9ef49abc35a06b2e2eea71e8b6a3f9b7874159db5c2b63e734cfee3cec739'}, 'banknote': {'X_raw': '51d78de1717ccad9798c90cf61c561c9524baddd098421db30f5a4f1c8f71974', 'X': '4eb4422a6b79e1e5203ce95ef0a2291498dfb1e39994d39785e6d51a0cd7fb5a', 'y': 'a12511f8e5a48806a5f14263769e005e1511e23e9ede3088dfd97657c117155c', 'row_ids': 'a6fbccd5a811a05ea65c110c56dfbeed81b5091048ab0558d96755a2e063aa9e'}, 'wine': {'X_raw': '8edcf3903afd97c64d51e0212eb10b213f7943da650574d1c055b836c5c35d37', 'X': 'd7f713b79c8aa9e6e5a8452681aaf68ab2fcc7fd28ea9b6e18294d670697ada8', 'y': 'b6514e07a6d27854e0da000ca82d8d8cc3d0cb904c51600344ca34f1ee04884b', 'row_ids': '36f97ef14cd7769a635ee76380be5a8d1bb87771207cd1f1072ea5bf8ec39e5f'}, 'breast_cancer': {'X_raw': '6b202a2072f9a0385f405a8f8605b1b06f6f36ae6d23d9cd6cbbc0974a416bc7', 'X': '2d05022a09664922654ffcb11d398f2623746f69b3c0e1cde03bd2f7ad9072e3', 'y': '1763c0f8cd3f454b9eda72031113884593bb382411c0cd511e3ec6e7ac85cddb', 'row_ids': 'a7a0c058c22cc0b88ea47e9eb61e5cfbeeed1769206accbf6a59978d8111e91f'}, 'digits': {'X_raw': '20def7f70a702f0af9732fbba4375e147a7d54fe70d8c45569b8e7c1c7010c10', 'X': 'dddd7515d18f8dfd1b0dff6b52dd0e469a7ebf6c872b4b994bbcbc7452d80829', 'y': 'a3c91c262eddcf7ba8f0e37507c30284493c9b20412ffe4af30d536401f7ba21', 'row_ids': '67a20cc9a4089e17305c2563b9f49897c4fbe4cee8eaac614bccebe6aa3704b9'}, 'seeds': {'X_raw': 'b370a0d4eae52f0c3f9e082ff742b1ad3ee16ee3e0eb0dbcb76cee89814648ea', 'X': '733f671f01e4ceed73bc7f86065feff8f04f69bdcb31d9aa5eaa41d364bb80e4', 'y': '7b0d3c10c58a02be51a128c452c83e99b66b897ab8accf693992a457a9144b12', 'row_ids': '4698e7cf925510c88263784a9478f0d35d4f7cb64144ca887a8d66fdd260dbb7'}, 'glass': {'X_raw': 'ec3bb3cbe02e01600db223bfebc8573a269585bfbf07dd37d7e4c25ec77396f6', 'X': '145ed3a9b294d83e7943eae92ba86d13a60f90bb15c74e7b384679102af3bbaa', 'y': '6bcfea2c4f6eca88cd729d30b681da511fb953d7f450d0f4ea4d7da33306c5ff', 'row_ids': '9084ecbea3e0b5015a563f787755eec578aa3f4c2f6357b80de20637d4165637'}, 'sonar': {'X_raw': '0d6ab10c8d741bdda9f7525380780a205260255ac763c9e103fe0da3e41deec8', 'X': '22d6b48b2f2c0439b2521d8170cce0fdf1f21bad77d53ca57f740b469ddf8476', 'y': '4b4149bb5d933462f19a37b5386358f2a88358de662446ebd41a80c521216b74', 'row_ids': '714c9068559488fc60c7475070dfe380606647aec734e2a2af423170b89470ac'}, 'vehicle': {'X_raw': '5d549b80abb6a5fb231b5616172799b0077fbeab7ce692d489c5409fa15c3d39', 'X': '7c87afaa09d7b42ed7d1762cade5ffdcca65d581533107dbacf59ded6574cb11', 'y': 'bf17e0c8a44165630aaa717587fb6921ae341dd2677474fbb741feddaf869d38', 'row_ids': 'b7c06cb9ed106f60239ce9237cc4173c558b3a6126b628bac8de08ff85c0b6f5'}, 'segmentation': {'X_raw': '28985038149ca8d7cf628ec3e9d6b0be8f4731d9949506f6d00958cfefcb6f79', 'X': '1d1eebefeddccc9e0cda359a4b0694c005b9c141292b1cfa0092430f381208fd', 'y': '1d0c9899beaa354cba666b9ff4005dea2313668453ebf58728650d0bb578f5ce', 'row_ids': '2f53e9d0a369fcd8b53e14caf1ee7144e7dd8e97568893a208aeec078f01a77b'}, 'breast_tissue': {'X_raw': '528536d56346822cf1ab5c78ed09c6d206b7b5d841e488dc9158981ff7aeb4b7', 'X': '73c8480f02706b81f33e13dce4cc5134bf2bb53cac17752112b308d38333f507', 'y': 'cadbc4a65950e9b07f7250f14d814c7ae4ae63aa731bd874d70950cd3c3994b2', 'row_ids': '5e8ecb03daec6f350174da46332d25c5c10e581b1c59bba81bcd655595a9a2d6'}, 'leaf': {'X_raw': '0609fa9c2d3c02ea28647828aee8bd3528f4f271412792f31a3305a4139c1d55', 'X': '44c652309b607e0d5df814ee53c1bc06c5316330460c2f6c7818f845b9580aa0', 'y': 'f028a1779a3ffde991a70418aac4d323a5196f248a2bac4ef9e6e5adce64522b', 'row_ids': '894f2a35f85b5c32ee16b8bbf8aa7736c0f880be11ce0fe8f603aaa8fd21dccd'}, 'libras': {'X_raw': 'b3bf6652437cd4e73b934d713110534fb1467881592d92417e886ecb038a156f', 'X': '770776db23d8926691477128a806732e1675ac0251982537340c14933b70c58a', 'y': '708d64755fbdf194880cc70da0a3b475dab8e7afad3eb5e099c9ba34d19bafa0', 'row_ids': 'c5b5d89abc3aae5b9debf1091964dedb9e385c3348ba8746fecd1165ae058289'}, 'vertebral': {'X_raw': 'b2b5476e272380fec75990026cef94e9e966525f69d551327a8f0441a15db81a', 'X': '112bc4b7e369d0b3115225aca4155b1fcf2ba2ec422a58b909afa44f52be6ac4', 'y': '8517ad01791428d38f130203ab3cb1122968517052129d2f17f3c279a5aaefda', 'row_ids': '82e1526c66e17bf611e5b28ed09fb561d7f66ab6cd61065f02a4445816ec2c85'}, 'wireless': {'X_raw': 'bda2583cd6a2038ca221b71e5fd5ee78e0a6a2ce0e1f6cebe971310bb2ea71b8', 'X': '3c57a1ff89d37a446382e611513cc5dd689aa7b567c806618bc9992f3f32ff8d', 'y': '143cfdb0f701f726c0005f442b1ca9fc1570f4cee6d7d7e416c9746f8d6f014a', 'row_ids': '55f385cf2332d9056aaed6f496e7bebd2df52c6a9547ce2144b309432d4b0290'}, 'parkinsons': {'X_raw': '95465a3e4bd04587305b328ac7491fa81d747462a43540f5c674c4fb2e6f75e0', 'X': 'e8321e11ab2e1c888c18a67b039db45212c2103d141e71017c6e5296ed49a642', 'y': '7a81563ca99d9cd4e31a44451d5a5a1e4d3f2457f27f0672540efbb0c7ca2cf4', 'row_ids': 'ccac2e3783be8b3f2ca1239aff37da890641e0f0fd4f9256ab25d439d74a7954'}, 'coimbra': {'X_raw': '862ed9f2e288027019c73c5bfdf1e1a2b6609fbf8cf7002f2b862f3b76588b67', 'X': '9fc6a5a0043deecf3134844de595e36cbd625831818628d2d8db457577b76931', 'y': '28133af969a893818cded36a57650aaf495b83fafccb93b7c3a3124c9762b6e2', 'row_ids': 'b761f1f690ecb0af7ec17dc78b466b57b15a19107d9dac195322a4ffa4b27b82'}, 'raisin': {'X_raw': 'efb5e08d8b85a0d27dd0c7d1ada03c582bc015cfc0e5bcb180072801a58e4de9', 'X': '55c46a6a9e58986c86c647c30bbe3ffd16c482591ea2a24d8c039a7633306770', 'y': 'e397e17482579b50885ea8e6641d8f515a1fb78d124c31955d13e0d3bf5f8245', 'row_ids': '1227ffb5bb6f93f86cb920f1c01e4cea108f97044322266d2e8d889f936425d5'}, 'haberman': {'X_raw': '1e4872ba6ffdb17171b0ce59407c3e08de8d4f89ced01e4adb14632256ee5345', 'X': '9f21eeebcc4fd29d8f4bb34bba36f6ed3303fb033c42ae6a1358ce8d370d4c51', 'y': '98af84f7a082ef1254a4f6d05c4ab42f836e7ad29548829047ee0a44eb7c0539', 'row_ids': '4cb2c706a494a4ebec8bf21da38da43e4404bdc72f58fdaf97d590c732d0dee8'}, 'musk': {'X_raw': '95e7b38177932feaf39bb8e0d95f7e9b899532f0b1d12b65d90c540e1fbc7d6b', 'X': '0a6ef346647054b8aad3f9387801bb53487522c7f6e3db869b7415a26f9414cd', 'y': 'ba237fa68f8c7eed0840c15028acf2f73ff2607474e8692e4903729ee1cac47b', 'row_ids': '8fe69d200611eab87b450fa19281654c4197603602d9be801cecc4ccdd5e9876'}}

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def array_record(array):
    import numpy as np
    return {"shape": list(array.shape), "dtype": str(array.dtype), "bytes": int(array.nbytes),
            "sha256": hashlib.sha256(memoryview(np.ascontiguousarray(array)).cast("B")).hexdigest()}

def save_json(value, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)

def download(url, path):
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    metadata_path = path.with_suffix(path.suffix + ".source.json")
    if path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        assert metadata["requested_url"] == url and sha(path) == metadata["sha256"]
        return metadata
    path.parent.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=Retry(total=4, backoff_factor=1,
                  status_forcelist=(429, 500, 502, 503, 504), allowed_methods=("GET",))))
    with session.get(url, timeout=(30, 90)) as response:
        response.raise_for_status()
        temporary = path.with_suffix(path.suffix + ".download")
        temporary.write_bytes(response.content)
        metadata = {"requested_url": url, "resolved_url": response.url,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "sha256": sha(temporary), "bytes": temporary.stat().st_size}
        temporary.replace(path)
    save_json(metadata, metadata_path)
    return metadata

def load_packaged(name):
    import numpy as np
    from sklearn import datasets
    bundle = getattr(datasets, "load_" + name)()
    filename, uci_id = PACKAGED[name]
    raw_path = DATA / "raw" / name / filename
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    content = importlib.resources.files("sklearn.datasets.data").joinpath(filename).read_bytes()
    if raw_path.exists():
        assert raw_path.read_bytes() == content
    else:
        raw_path.write_bytes(content)
    description_path = raw_path.parent / "description.txt"
    if description_path.exists():
        assert description_path.read_text(encoding="utf-8") == bundle.DESCR
    else:
        description_path.write_text(bundle.DESCR, encoding="utf-8")
    X = np.ascontiguousarray(bundle.data, dtype=np.float64)
    labels = np.asarray(bundle.target, dtype=np.int64)
    return X, labels, np.arange(len(X), dtype=np.int64), {
        "loader": f"sklearn.datasets.load_{name}", "sklearn_version": importlib.metadata.version("scikit-learn"),
        "source_url": f"https://scikit-learn.org/1.6/modules/generated/sklearn.datasets.load_{name}.html",
        "original_uci_url": f"https://archive.ics.uci.edu/dataset/{uci_id}",
        "packaged_file": str(raw_path.relative_to(DATA)), "packaged_file_sha256": sha(raw_path),
        "description_file": str(description_path.relative_to(DATA)), "description_sha256": sha(description_path),
        "feature_names": list(map(str, bundle.feature_names)), "target_names": list(map(str, bundle.target_names)),
        "removed_identifier_columns": [],
        "subset_note": "All 1,797 8x8 images shipped by sklearn; the UCI optical digits test subset, not MNIST or the complete UCI training-plus-test corpus." if name == "digits" else "All rows provided by the pinned sklearn bundled dataset."}

def load_legacy_uci(name):
    import numpy as np
    settings = UCI[name]
    folder = DATA / "raw" / name
    metadata_path = folder / "uci_metadata.json"
    metadata_source = download(f"https://archive.ics.uci.edu/api/dataset?id={settings['id']}", metadata_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))["data"]
    assert metadata["uci_id"] == settings["id"]
    path = folder / settings["url"].rsplit("/", 1)[1]
    source = download(settings["url"], path)
    rows = [line.split(",") if name in ("glass", "sonar") else line.split()
            for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    expected_columns = ROSTER[name][1] + 1 + int(name == "glass")
    assert len(rows) == ROSTER[name][0]
    assert all(len(row) == expected_columns for row in rows), "Unexpected source columns; no silent parser repair."
    start = int(name == "glass")
    X = np.asarray([[float(item) for item in row[start:-1]] for row in rows], dtype=np.float64)
    labels = np.asarray([row[-1] for row in rows])
    ids = [row[0] for row in rows] if start else None
    source_metadata = {
        "uci_id": settings["id"], "source_url": metadata["repository_url"],
        "source_metadata": metadata_source, "data_source": source,
        "uci_citation": metadata.get("citation"), "uci_doi": metadata.get("dataset_doi"),
        "feature_names": settings["features"],
        "parser": "Comma-separated original UCI rows" if name in ("glass", "sonar") else "Whitespace split on each nonempty original UCI line; repeated tabs/spaces treated as one separator",
        "removed_identifier_columns": ["Id_number"] if name == "glass" else [],
        "source_identifiers": ids,
        "data_note": "Original numeric measurements retained; no feature engineering."}
    if name == "glass":
        assert list(map(int, ids)) == list(range(1, 215))
        source_metadata["class_note"] = "Six observed glass classes; source class 4 has no rows and is not fabricated."
    return X, labels, np.arange(len(X), dtype=np.int64), source_metadata

def prepare_fresh(name):
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    started = time.perf_counter()
    versions = {key: importlib.metadata.version(key) for key in ("numpy", "scipy", "scikit-learn", "requests")}
    assert versions["scikit-learn"] == "1.6.1"
    artifact, manifest_path = DATA / f"{name}.npz", DATA / f"{name}.manifest.json"
    source_hash = sha(__file__)
    if artifact.exists() or manifest_path.exists():
        assert artifact.is_file() and manifest_path.is_file(), "Inspect incomplete preparation before reuse."
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["preparation_source_sha256"] == source_hash and manifest["versions"] == versions
        assert sha(artifact) == manifest["npz"]["sha256"]
        with np.load(artifact, allow_pickle=False) as arrays:
            assert set(arrays.files) == set(manifest["arrays"])
            for key in arrays.files:
                assert array_record(arrays[key]) == manifest["arrays"][key]
        print(json.dumps({"dataset": name, "reused_verified": True}), flush=True)
        return manifest
    loader = load_parent if name in ("iris", "banknote") else (load_packaged if name in PACKAGED else load_uci)
    X_raw, source_labels, row_ids, provenance = loader(name)
    assert X_raw.shape == ROSTER[name] and source_labels.shape == (len(X_raw),)
    joint_valid = np.isfinite(X_raw).all(axis=1) & np.asarray([str(value).strip() not in ("", "?", "nan", "None") for value in source_labels])
    assert joint_valid.all(), "Unexpected missing values: this fixed complete-data roster requires inspection before joint filtering."
    classes, y = np.unique(source_labels, return_inverse=True)
    assert len(classes) == EXPECTED_CLASSES[name], "Unexpected class count."
    y = np.ascontiguousarray(y, dtype=np.int64)
    X_raw = np.ascontiguousarray(X_raw, dtype=np.float64)
    raw_hash = array_record(X_raw)["sha256"]
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    X = np.ascontiguousarray(scaler.fit_transform(X_raw), dtype=np.float64)
    assert array_record(X_raw)["sha256"] == raw_hash and np.isfinite(X).all()
    unique_values = [np.unique(X_raw[:, index]) for index in range(X_raw.shape[1])]
    constants = [index for index, values in enumerate(unique_values) if len(values) == 1]
    two_level = [index for index, values in enumerate(unique_values) if len(values) == 2]
    zero_one = [index for index, values in enumerate(unique_values) if set(values.tolist()).issubset({0.0, 1.0})]
    arrays = {"X_raw": X_raw, "X": X, "y": y, "row_ids": row_ids}
    DATA.mkdir(parents=True, exist_ok=True)
    temporary = artifact.with_suffix(".npz.tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(artifact)
    labels, counts = np.unique(y, return_counts=True)
    manifest = {
        "schema": 1, "dataset": name, "study": "experiment_2b_original_features_standardized_primary",
        "created_utc": datetime.now(timezone.utc).isoformat(), "preparation_source_sha256": source_hash,
        "source": provenance, "versions": versions, "python": sys.version, "platform": platform.platform(),
        "row_order": "Original source/package order, unchanged; no shuffled or selected subset.",
        "label_mapping": {str(label): index for index, label in enumerate(classes)},
        "label_counts": {str(int(label)): int(count) for label, count in zip(labels, counts)},
        "missing_values": {"policy": "Check features and targets jointly; reject unexpected missing rows rather than silently changing the fixed roster.", "rows_removed": 0},
        "preprocessing": {"primary_input": "X", "raw_input": "X_raw", "embedding": None, "binarization": False,
            "description": "Identical StandardScaler transformation of all original numeric features supplied to both methods. Fit on all rows without labels; no held-out generalization claim. This deliberately differs from the earlier unscaled Iris/Banknote UMAP experiment.",
            "scaler": "sklearn.preprocessing.StandardScaler", "parameters": scaler.get_params(),
            "mean": scaler.mean_.tolist(), "variance_ddof0": scaler.var_.tolist(), "scale": scaler.scale_.tolist(),
            "n_samples_seen": int(scaler.n_samples_seen_), "constant_features_retained": constants},
        "feature_characteristics": {"unique_value_counts": [len(values) for values in unique_values],
            "constant_columns": constants, "two_level_columns": two_level,
            "two_level_values": {str(index): unique_values[index].tolist() for index in two_level},
            "zero_one_columns_including_constants": zero_one, "all_features_zero_one": len(zero_one) == X_raw.shape[1],
            "note": "Original feature cardinalities are reported, not an assumption that every numeric column is continuous. Discrete, two-level and constant columns remain included."},
        "arrays": {key: array_record(value) for key, value in arrays.items()},
        "npz": {"file": artifact.name, "sha256": sha(artifact), "bytes": artifact.stat().st_size},
        "preparation_seconds": time.perf_counter() - started,
        "clustering_timing": "Source loading, validation and StandardScaler preprocessing are outside subsequent clustering timings."}
    save_json(manifest, manifest_path)
    print(json.dumps({"dataset": name, "shape": list(X.shape), "classes": len(classes),
                      "constant_columns": constants, "two_level_columns": two_level,
                      "manifest": str(manifest_path)}), flush=True)
    return manifest


def metadata(name, identifier):
    path = DATA / "raw" / name / "uci_metadata.json"
    record = download(f"https://archive.ics.uci.edu/api/dataset?id={identifier}", path)
    value = json.loads(path.read_text(encoding="utf-8"))["data"]
    assert value["uci_id"] == identifier
    return value, record


def load_parent(name):
    """Fresh public-source fallback for the original Iris/Banknote arrays."""
    import numpy as np
    import pandas as pd
    identifier = {"iris": 53, "banknote": 267}[name]
    info, metadata_record = metadata(name, identifier)
    path = DATA / "raw" / name / "data.csv"
    url = f"https://archive.ics.uci.edu/static/public/{identifier}/data.csv"
    assert info["data_url"] == url
    data_record = download(url, path)
    frame = pd.read_csv(path)
    features = [item["name"] for item in info["variables"] if item["role"] == "Feature"]
    assert len(info["target_col"]) == 1
    target = info["target_col"][0]
    assert not frame[features + [target]].isna().any().any()
    return np.asarray(frame[features], dtype=np.float64), frame[target].to_numpy(), np.arange(len(frame), dtype=np.int64), {
        "loader": "Public UCI CSV; same rows and features as the earlier verified preparation",
        "uci_id": identifier, "source_url": info["repository_url"],
        "metadata_source": metadata_record, "data_source": data_record,
        "feature_names": features, "target_name": target, "removed_identifier_columns": []}


def load_uci(name):
    import numpy as np
    import pandas as pd
    if name in UCI:
        return load_legacy_uci(name)
    identifier = NEW_UCI_IDS[name]
    info, metadata_record = metadata(name, identifier)
    folder = DATA / "raw" / name
    provenance = {"uci_id": identifier, "source_url": info["repository_url"],
                  "metadata_source": metadata_record, "metadata_instance_count": info["num_instances"],
                  "metadata_feature_count": info["num_features"], "removed_identifier_columns": []}
    legacy = "https://archive.ics.uci.edu/ml/machine-learning-databases/"
    if name == "vehicle":
        rows, sources, ranges = [], [], []
        for suffix in "abcdefghi":
            filename = f"xa{suffix}.dat"
            path = folder / filename
            sources.append(download(legacy + "statlog/vehicle/" + filename, path))
            block = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
            assert len(block) == 94 and all(len(row) == 19 for row in block)
            ranges.append({"file": filename, "start_inclusive": len(rows), "stop_exclusive": len(rows) + len(block)})
            rows.extend(block)
        X = np.asarray([row[:-1] for row in rows], dtype=np.float64)
        labels = np.asarray([row[-1] for row in rows])
        provenance.update(data_sources=sources, concatenation_order=ranges,
            feature_names=[item["name"] for item in info["variables"] if item["role"] == "Feature"],
            source_note="Original nine files in alphabetical order, 94 rows each. UCI metadata reports 946 but the complete files contain 846. The API CSV has a shifted row/spurious class 204 and was not used; no source row was repaired or dropped.")
    elif name == "segmentation":
        rows, sources, ranges = [], [], []
        headers = []
        for filename, expected in (("segmentation.data", 210), ("segmentation.test", 2100)):
            path = folder / filename
            sources.append(download(legacy + "image/" + filename, path))
            lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
                     if line.strip() and not line.startswith(";")]
            headers.append(lines[0].split(","))
            block = [line.split(",") for line in lines[1:]]
            assert len(block) == expected and all(len(row) == 20 for row in block)
            ranges.append({"file": filename, "start_inclusive": len(rows), "stop_exclusive": len(rows) + len(block)})
            rows.extend(block)
        assert headers[0] == headers[1] and len(headers[0]) == 19
        X = np.asarray([row[1:] for row in rows], dtype=np.float64)
        labels = np.asarray([row[0] for row in rows])
        provenance.update(data_sources=sources, concatenation_order=ranges, feature_names=headers[0],
            source_note="Use all 2,310 published observations, training file first then test file, without shuffling. The API CSV has only the 210 training rows. These are small regions from seven source images, not independent new scenes. The constant region-pixel-count feature remains included.")
    elif name == "breast_tissue":
        path = folder / "BreastTissue.xls"
        source = download(legacy + "00192/BreastTissue.xls", path)
        frame = pd.read_excel(path, sheet_name="Data", engine="xlrd")
        assert list(frame.columns[:2]) == ["Case #", "Class"] and frame.shape == (106, 11)
        features = list(frame.columns[2:])
        X, labels = frame[features].to_numpy(dtype=np.float64), frame["Class"].to_numpy()
        provenance.update(data_source=source, sheet="Data", feature_names=features,
            removed_identifier_columns=["Case #"], source_identifiers=frame["Case #"].astype(int).tolist(),
            parser="pandas.read_excel(engine='xlrd', sheet_name='Data')",
            xlrd_version=importlib.metadata.version("xlrd"),
            source_note="All original six tissue classes retained; no merging to the alternate four-class task.")
    elif name == "leaf":
        archive = folder / "leaf.zip"
        source = download(legacy + "00288/leaf.zip", archive)
        with zipfile.ZipFile(archive) as bundle:
            members = [item for item in bundle.infolist() if Path(item.filename).name.lower() == "leaf.csv"]
            assert len(members) == 1
            member = members[0]
            assert member.file_size < 1024**2
            content = bundle.read(member)
        path = folder / "leaf.csv"
        if path.exists():
            assert path.read_bytes() == content
        else:
            path.write_bytes(content)
        frame = pd.read_csv(path, header=None)
        assert frame.shape == (340, 16)
        labels = frame.iloc[:, 0].to_numpy()
        X = frame.iloc[:, 2:].to_numpy(dtype=np.float64)
        features = ["eccentricity", "aspect_ratio", "elongation", "solidity", "stochastic_convexity",
                    "isoperimetric_factor", "maximal_indentation_depth", "lobedness", "average_intensity",
                    "average_contrast", "smoothness", "third_moment", "uniformity", "entropy"]
        provenance.update(data_source=source, archive_member=member.filename, member_sha256=sha(path),
            member_bytes=len(content), feature_names=features, removed_identifier_columns=["Specimen Number"],
            source_identifiers=frame.iloc[:, 1].astype(int).tolist(),
            source_note="Use the published 340-row numeric leaf.csv only; no feature extraction from images. The accompanying image collection description mentions 40 species, while this numeric table contains 30 observed species. Archive images are not extracted or used.")
    elif name in ("libras", "wireless"):
        suffix = "libras/movement_libras.data" if name == "libras" else "00422/wifi_localization.txt"
        path = folder / suffix.rsplit("/", 1)[1]
        source = download(legacy + suffix, path)
        frame = pd.read_csv(path, header=None, sep="," if name == "libras" else r"\s+")
        X, labels = frame.iloc[:, :-1].to_numpy(dtype=np.float64), frame.iloc[:, -1].to_numpy()
        features = ([f"frame_{i // 2 + 1}_{'x' if i % 2 == 0 else 'y'}" for i in range(90)] if name == "libras"
                    else [f"wifi_signal_{i + 1}" for i in range(7)])
        provenance.update(data_source=source, feature_names=features,
            source_note="All published movement_libras.data rows; retain its provided 90 normalized trajectory coordinates and 15 gesture labels, not one of the smaller alternate files." if name == "libras" else "All 2,000 original received-signal-strength vectors and four room labels.")
    else:
        path = folder / "data.csv"
        url = f"https://archive.ics.uci.edu/static/public/{identifier}/data.csv"
        assert info["data_url"] == url
        source = download(url, path)
        frame = pd.read_csv(path)
        assert len(info["target_col"]) == 1
        target = info["target_col"][0]
        removed = ["name"] if name == "parkinsons" else (["molecule_name", "conformation_name"] if name == "musk" else [])
        features = [column for column in frame.columns if column not in [target] + removed]
        X, labels = frame[features].to_numpy(dtype=np.float64), frame[target].to_numpy()
        provenance.update(data_source=source, feature_names=features, target_name=target,
                          removed_identifier_columns=removed)
        if removed:
            provenance["source_identifiers"] = {key: frame[key].astype(str).tolist() for key in removed}
        if name == "parkinsons":
            provenance["source_note"] = "195 voice recordings. The description says 31 people; names have 32 distinct subject prefixes. Metadata says 197 recordings but the source has 195. All 22 voice measures retained. Duplicate metadata header names are disambiguated by pandas suffixes, never collapsed. Repeated recordings are not independent patients."
            provenance["subject_count"] = int(frame["name"].str.rsplit("_", n=1).str[0].nunique())
        if name == "musk":
            assert features == [f"f{i}" for i in range(1, 167)]
            provenance["source_note"] = "476 conformations of 92 molecules, with molecule-level musk labels. Exclude both string identifiers even though the UCI API marks them as features; keep all 166 numeric descriptors. Conformations are not independent molecules."
            provenance["molecule_count"] = int(frame["molecule_name"].nunique())
        if name == "vertebral":
            provenance["source_note"] = "The original three-class task: Hernia, Spondylolisthesis and Normal; no merging to a binary target."
    assert np.asarray(X).shape == ROSTER[name]
    return np.ascontiguousarray(X, dtype=np.float64), np.asarray(labels), np.arange(len(X), dtype=np.int64), provenance


def verify_arrays(name, manifest, artifact):
    import numpy as np
    assert sha(artifact) == manifest["npz"]["sha256"]
    with np.load(artifact, allow_pickle=False) as stored:
        assert set(stored.files) == {"X_raw", "X", "y", "row_ids"}
        arrays = {key: stored[key] for key in stored.files}
    for key, value in arrays.items():
        assert array_record(value) == manifest["arrays"][key]
    assert arrays["X"].shape == ROSTER[name] and np.isfinite(arrays["X"]).all()
    assert arrays["X_raw"].shape == ROSTER[name] and np.isfinite(arrays["X_raw"]).all()
    assert len(np.unique(arrays["y"])) == EXPECTED_CLASSES[name]
    assert np.array_equal(arrays["row_ids"], np.arange(len(arrays["X"])))
    if name in BASELINE_ARRAY_HASHES:
        assert {key: array_record(value)["sha256"] for key, value in arrays.items()} == BASELINE_ARRAY_HASHES[name]
    return arrays


def validate_neighbors(name, arrays, manifest):
    """Label-free checks of duplicate rows and the fifth non-self neighbour."""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    from threadpoolctl import threadpool_limits
    X = arrays["X"]
    result = {"dataset": name, "status": "ok", "preparation_source_sha256": sha(__file__),
              "X_sha256": manifest["arrays"]["X"]["sha256"], "rows": len(X), "features": X.shape[1],
              "finite": bool(np.isfinite(X).all()), "duplicate_rows_beyond_first": len(X) - len(np.unique(X, axis=0)),
              "definition": "For each input row, the sixth sorted neighbour including self is the fifth other-row neighbour; duplicates are retained. This check uses no labels.",
              "metrics": {}}
    with threadpool_limits(limits=1):
        for metric in ("euclidean", "manhattan"):
            model = NearestNeighbors(n_neighbors=6, metric=metric, algorithm="brute", n_jobs=1).fit(X)
            fifth = model.kneighbors(X, return_distance=True)[0][:, 5]
            median = float(np.median(fifth))
            assert np.isfinite(fifth).all() and median > 0, "Degenerate fifth-neighbour radius."
            result["metrics"][metric] = {"median_fifth_neighbor": median, "zero_fifth_neighbor_rows": int(np.count_nonzero(fifth == 0))}
    save_json(result, DATA / f"{name}.validation.json")
    return result




def prepare(name):
    artifact = DATA / f"{name}.npz"
    manifest = prepare_fresh(name)
    arrays = verify_arrays(name, manifest, artifact)
    validate_neighbors(name, arrays, manifest)
    return manifest


def main():
    global DATA
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", choices=tuple(ROSTER), nargs="+", default=tuple(ROSTER))
    parser.add_argument("--data-dir", type=Path, default=DATA)
    args = parser.parse_args()
    DATA = args.data_dir.resolve()
    assert len(args.datasets) == len(set(args.datasets))
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    for name in args.datasets:
        prepare(name)


if __name__ == "__main__":
    main()
