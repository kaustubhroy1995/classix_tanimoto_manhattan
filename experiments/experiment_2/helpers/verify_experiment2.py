"""Portable post-run numerical audit for the twenty-dataset original-feature screen.

No fits, downloads, or historical workspace dependencies. Strict mode validates
fresh-run identity and evidence. --numerical-only is for an explicitly labelled
packaging validation fixture; it does not certify a fresh benchmark execution.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import sys
import time
import traceback
import zipfile

ROOT = RUN
METHODS = ("manhattan", "euclidean")
ANCHORS = ("iris", "banknote", "wine", "breast_cancer", "digits", "seeds", "glass", "sonar")
ROSTER = ANCHORS + ("vehicle", "segmentation", "breast_tissue", "leaf", "libras", "vertebral", "wireless", "parkinsons", "coimbra", "raisin", "haberman", "musk")
SHAPES = dict(zip(ROSTER, ((150,4),(1372,4),(178,13),(569,30),(1797,64),(210,7),(214,9),(208,60),
                         (846,18),(2310,19),(106,9),(340,14),(360,90),(310,6),(2000,7),(195,22),(116,9),(900,7),(306,3),(476,166))))
CLASSES = dict(zip(ROSTER, (3,2,3,2,10,3,6,2,4,7,6,30,15,3,4,2,2,2,2,2)))


def require(condition, message):
    if not condition:
        raise ValueError(message)

def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))

def lines(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    for attempt in range(40):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt == 39:
                raise
            time.sleep(.05)

def process_alive(pid):
    if int(pid) == os.getpid():
        return True
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.GetExitCodeProcess.restype = wintypes.BOOL
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        handle = kernel.OpenProcess(0x1000, False, int(pid))
        if not handle:
            return ctypes.get_last_error() == 5
        try:
            code = wintypes.DWORD()
            require(kernel.GetExitCodeProcess(handle, ctypes.byref(code)), "Cannot establish producer process status")
            return code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(int(pid), 0)
        return True
    except ProcessLookupError:
        return False

def close(actual, expected, message):
    require(isinstance(actual, (float, int)) and not isinstance(actual, bool)
            and math.isfinite(actual) and math.isfinite(float(expected))
            and math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12), message)

def array_record(array, np):
    return dict(shape=list(array.shape), dtype=str(array.dtype), bytes=int(array.nbytes),
                sha256=hashlib.sha256(memoryview(np.ascontiguousarray(array)).cast("B")).hexdigest())

def normalization(method, X, np):
    if method == "manhattan":
        return float(np.median(np.sum(X-X.min(axis=0), axis=1)))
    return float(np.median(np.linalg.norm(X-X.mean(axis=0), axis=1))) or 1.0

def rank(rows):
    return min(rows, key=lambda row: (-row["ari"], row["min_pts"], row["radius_index"]))

def scientific(row):
    return {key: value for key, value in row.items() if key not in ("origin", "labels_path")}

def row_checks(row, n):
    require(row["status"] == "ok", "Non-success result cannot be counted as completed")
    for key in ("ari", "ami"):
        require(math.isfinite(row[key]) and -1.00000001 <= row[key] <= 1.00000001, "Invalid score: "+key)
    require(math.isfinite(row["fit_seconds"]) and row["fit_seconds"] > 0, "Invalid fit time")
    require(1 <= row["n_groups"] <= n and 0 <= row["n_clusters"] <= n and 0 <= row["noise_points"] <= n, "Invalid group/cluster/noise counts")
    require(isinstance(row["warnings"], list) and row["origin"] in ("imported", "fresh"), "Missing fit warning/origin metadata")

def labels_from(row, n, results, seen, np):
    path = (ROOT/row["labels_path"]).resolve()
    require(path.is_relative_to((results/"assignments").resolve()), "Assignment escapes selected results directory")
    require(path not in seen, "Assignment path reused")
    seen.add(path)
    with np.load(path, allow_pickle=False) as arrays:
        require(set(arrays.files) == {"labels"}, "Unexpected assignment fields")
        labels = arrays["labels"]
    require(labels.shape == (n,) and labels.dtype == np.int64 and np.all(labels >= -1), "Invalid saved assignments")
    require(np.unique(labels[labels >= 0]).size == row["n_clusters"] and int(np.sum(labels == -1)) == row["noise_points"], "Saved partition counts differ")
    return labels, dict(file=path.relative_to(ROOT).as_posix(), sha256=sha(path), array_sha256=array_record(labels, np)["sha256"])

def scores(row, truth, labels, metrics):
    for key, function in metrics.items():
        close(row[key], function(truth, labels), "Saved "+key+" differs")

def verify_plan(dataset, X, results, np):
    from scipy.spatial.distance import pdist, squareform
    design = read(results/"plans"/(dataset+".json"))
    require(design["dataset"] == dataset and design["kind"] == "real"
            and design["input_sha256"] == sha(ROOT/"data"/(dataset+".npz"))
            and design["radius_count"] == 32 and design["min_pts"] == [0,5,10,20,40,80]
            and design["merge_scale"] == 1.4 and set(design["methods"]) == set(METHODS), "Dataset plan differs")
    for method in METHODS:
        distances = pdist(X, metric="cityblock" if method == "manhattan" else "euclidean")
        # Diagonal zero occupies index0; sorted index5 is the fifth other row.
        near = float(np.median(np.sort(squareform(distances), axis=1)[:,5]))
        far = float(np.median(distances[distances > 0]))
        scale = normalization(method,X,np)
        require(0 < .35*near < .9*far and scale > 0, "Invalid unchanged calibration interval")
        geometry = design["methods"][method]
        for key, value in (("median_5nn",near),("median_pair",far),("physical_low",.35*near),("physical_high",.9*far),("internal_scale",scale)):
            close(geometry[key],value,"Label-free radius calibration differs: "+dataset+"/"+method+"/"+key)
        require(np.allclose(geometry["radii"],np.geomspace(.35*near,.9*far,32)/scale,rtol=1e-12,atol=1e-12), "Radius grid differs")
    return design

def verify_grid(dataset, X, rows, design, results, np, warnings):
    initial_keys = {(method,index,minimum) for method in METHODS for index in range(32) for minimum in (0,5,10,20,40,80)}
    initial = [row for row in rows if (row["method"],row["radius_index"],row["min_pts"]) in initial_keys]
    require(len(initial) == 384, "Initial grid incomplete/duplicated")
    winners = {method:rank([row for row in initial if row["method"] == method]) for method in METHODS}
    radius_edge = any(row["radius_index"] in (0,31) for row in winners.values())
    min_edge = any(row["min_pts"] == 80 for row in winners.values())
    indices = list(range(-4,36)) if radius_edge else list(range(32))
    minima = [0,5,10,20,40,80]+([120,160] if min_edge else [])
    extension = read(results/"plans"/(dataset+"_extension.json"))
    require(extension["radius_edge"] == radius_edge and extension["min_edge"] == min_edge
            and extension["indices"] == indices and extension["min_pts"] == minima
            and {method:scientific(row) for method,row in extension["initial_winners"].items()} == {method:scientific(row) for method,row in winners.items()}, "Predeclared symmetric extension differs")
    expected = {(method,index,minimum) for method in METHODS for index in indices for minimum in minima}
    actual = [(row["method"],row["radius_index"],row["min_pts"]) for row in rows]
    require(len(actual) == len(set(actual)) == len(expected) and set(actual) == expected, "Final equal-budget grid incomplete/duplicated")
    for row in rows:
        row_checks(row,len(X))
        require(row["kind"] == "real" and row["merge_scale"] == 1.4, "Grid parameters differ")
        geometry = design["methods"][row["method"]]
        radius = geometry["physical_low"]*(geometry["physical_high"]/geometry["physical_low"])**(row["radius_index"]/31)/geometry["internal_scale"]
        close(row["radius"],radius,"Recorded radius differs")
        for field in ("preprocessing_scale","fitted_scale"):
            close(row[field],geometry["internal_scale"],"Recorded normalization differs")
    for method in METHODS:
        winner = rank([row for row in rows if row["method"] == method])
        if winner["radius_index"] in (min(indices),max(indices)) or winner["min_pts"] == max(minima):
            warnings.append(dict(dataset=dataset,method=method,radius_index=winner["radius_index"],min_pts=winner["min_pts"],reason="Winner remains at outer tested boundary; best only within tested grid"))
    return dict(total_fits=len(rows),fits_per_method=len(rows)//2,radius_extension=radius_edge,min_pts_extension=min_edge,final_radius_count=len(indices),final_min_pts=minima)

def illustrative_panel(selected, permutations):
    chosen = {(row["dataset"],row["method"]):row for row in selected}
    orders = {(row["dataset"],row["method"],row["replicate"]):row for row in permutations}
    classifications, ranges = {}, {}
    for dataset in ROSTER:
        pairs = [(chosen[dataset,"manhattan"],chosen[dataset,"euclidean"])] + [(orders[dataset,"manhattan",i],orders[dataset,"euclidean",i]) for i in range(5)]
        deltas = [(a["ari"]-b["ari"],a["ami"]-b["ami"]) for a,b in pairs]
        classifications[dataset] = ("manhattan" if all(a>0 and b>0 for a,b in deltas) else "euclidean" if all(a<0 and b<0 for a,b in deltas) else "mixed")
        ranges[dataset] = dict(classification=classifications[dataset],original_delta_ari=deltas[0][0],original_delta_ami=deltas[0][1],
            permutation_delta_ari_range=[min(x[0] for x in deltas[1:]),max(x[0] for x in deltas[1:])],
            permutation_delta_ami_range=[min(x[1] for x in deltas[1:]),max(x[1] for x in deltas[1:])])
    panel = list(ANCHORS)
    for direction in ("manhattan","euclidean"):
        candidates = sorted(name for name in ROSTER if name not in ANCHORS and classifications[name] == direction)
        if candidates:
            panel.append(candidates[0])
    for group in ("mixed",None):
        for name in sorted(set(ROSTER)-set(panel)):
            if len(panel) == 10:
                break
            if group is None or classifications[name] == group:
                panel.append(name)
    require(len(panel) == 10 and len(set(panel)) == 10, "Panel fallback did not produce ten distinct datasets")
    return dict(dataset_ids=panel,classification_counts=dict(Counter(classifications.values())),datasets=ranges,
                interpretation="Outcome-stratified illustrative panel; all twenty datasets remain the screened population.")

def checked_download(directory, record):
    path = directory/record["requested_url"].rsplit("/",1)[1]
    # Metadata API URLs do not have a useful filename; callers supply that file.
    require(sha(path) == record["sha256"] and path.stat().st_size == record["bytes"], "Cached source bytes changed")
    require(read(path.with_suffix(path.suffix+".source.json")) == record, "Download provenance changed")
    return path

def mapped_labels(values, mapping, np):
    keys = list(mapping)
    require(sorted(mapping.values()) == list(range(len(keys))), "Invalid zero-based target mapping")
    try:
        numeric = [float(key) for key in keys]
    except ValueError:
        numeric = None
    require(keys == sorted(keys,key=float) if numeric is not None else keys == sorted(keys), "Target mapping order differs from sorted source classes")
    encoded = []
    for value in values:
        key = str(value)
        if key not in mapping and numeric is not None:
            candidates = [candidate for candidate in keys if float(candidate) == float(value)]
            require(len(candidates) == 1, "Ambiguous numeric target mapping")
            key = candidates[0]
        require(key in mapping, "Source target absent from recorded mapping")
        encoded.append(mapping[key])
    require(set(encoded) == set(mapping.values()), "Unused/fabricated target class")
    return np.asarray(encoded,dtype=np.int64)

def parse_source(dataset, manifest, data_root, np):
    source = manifest["source"]
    directory = data_root/"raw"/dataset
    features = source["feature_names"]
    require(len(features) == SHAPES[dataset][1] and len(set(features)) == len(features), "Feature roster is incomplete/duplicated")
    metadata_record = source.get("metadata_source",source.get("source_metadata"))
    if metadata_record:
        path = directory/"uci_metadata.json"
        require(sha(path) == metadata_record["sha256"]
                and read(path.with_suffix(path.suffix+".source.json")) == metadata_record, "UCI metadata evidence changed")
        metadata = read(path)["data"]
        if "uci_id" in source:
            require(metadata["uci_id"] == source["uci_id"], "UCI dataset identity differs")
    if dataset in ("wine","breast_cancer","digits"):
        path = data_root/source["packaged_file"]
        require(sha(path) == source["packaged_file_sha256"]
                and sha(data_root/source["description_file"]) == source["description_sha256"], "Packaged source changed")
        if dataset == "digits":
            with gzip.open(path,"rt") as stream:
                rows = np.loadtxt(stream,delimiter=",")
        else:
            rows = np.loadtxt(path,delimiter=",",skiprows=1)
        X, labels = rows[:,:-1],rows[:,-1].astype(np.int64)
    elif dataset in ("vehicle","segmentation"):
        expected_names = [f"xa{suffix}.dat" for suffix in "abcdefghi"] if dataset == "vehicle" else ["segmentation.data","segmentation.test"]
        blocks, ranges, headers = [],[],[]
        for record, filename in zip(source["data_sources"],expected_names):
            path = checked_download(directory,record)
            require(path.name == filename, "Concatenation file order differs")
            if dataset == "vehicle":
                block = [line.split() for line in path.read_text().splitlines() if line.strip()]
                require(len(block) == 94 and all(len(row) == 19 for row in block), "Vehicle source rows/columns differ")
            else:
                text_rows = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith(";")]
                headers.append(text_rows[0].split(","))
                block = list(csv.reader(text_rows[1:]))
                require(len(block) == (210 if filename.endswith(".data") else 2100)
                        and all(len(row) == 20 for row in block), "Segmentation source rows/columns differ")
            start = sum(len(item) for item in blocks)
            ranges.append(dict(file=filename,start_inclusive=start,stop_exclusive=start+len(block)))
            blocks.append(block)
        require(len(source["data_sources"]) == len(expected_names) and ranges == source["concatenation_order"], "Concatenation metadata differs")
        rows = [row for block in blocks for row in block]
        X = np.asarray([row[:-1] if dataset == "vehicle" else row[1:] for row in rows],dtype=np.float64)
        labels = [row[-1] if dataset == "vehicle" else row[0] for row in rows]
        if dataset == "segmentation":
            require(headers[0] == headers[1] == features, "Segmentation feature order differs")
        else:
            require(features == [item["name"] for item in metadata["variables"] if item["role"] == "Feature"], "Vehicle feature order differs")
    elif dataset == "breast_tissue":
        import xlrd
        path = checked_download(directory,source["data_source"])
        sheet = xlrd.open_workbook(path).sheet_by_name("Data")
        headers = sheet.row_values(0)
        rows = [sheet.row_values(i) for i in range(1,sheet.nrows)]
        require(headers[:2] == ["Case #","Class"] and headers[2:] == features
                and [int(row[0]) for row in rows] == source["source_identifiers"]
                and source["removed_identifier_columns"] == ["Case #"], "Breast Tissue identifiers/feature selection differs")
        X,labels = np.asarray([row[2:] for row in rows],dtype=np.float64),[row[1] for row in rows]
    elif dataset == "leaf":
        path = checked_download(directory,source["data_source"])
        with zipfile.ZipFile(path) as archive:
            content = archive.read(source["archive_member"])
        require(len(content) == source["member_bytes"] and hashlib.sha256(content).hexdigest() == source["member_sha256"]
                and content == (directory/"leaf.csv").read_bytes(), "Leaf archive/member differs")
        rows = list(csv.reader(content.decode("utf-8-sig").splitlines()))
        require(all(len(row) == 16 for row in rows) and [int(row[1]) for row in rows] == source["source_identifiers"]
                and source["removed_identifier_columns"] == ["Specimen Number"], "Leaf identifier exclusion differs")
        X,labels = np.asarray([row[2:] for row in rows],dtype=np.float64),[int(row[0]) for row in rows]
    elif dataset in ("seeds","glass","sonar","libras","wireless"):
        path = checked_download(directory,source["data_source"])
        rows = [line.split() if dataset in ("seeds","wireless") else next(csv.reader([line])) for line in path.read_text().splitlines() if line.strip()]
        start = int(dataset == "glass")
        require(all(len(row) == SHAPES[dataset][1]+1+start for row in rows), "Delimited source columns differ")
        X,labels = np.asarray([row[start:-1] for row in rows],dtype=np.float64),[row[-1] for row in rows]
        if dataset == "glass":
            require([row[0] for row in rows] == source["source_identifiers"]
                    and [int(row[0]) for row in rows] == list(range(1,215)), "Glass identifier removal differs")
    else:
        path = checked_download(directory,source["data_source"])
        with path.open(encoding="utf-8-sig",newline="") as stream:
            table = list(csv.reader(stream))
        headers,rows = table[0],table[1:]
        # Match pandas' recorded suffixes, retaining all duplicate-named columns.
        used,renamed = {},[]
        for name in headers:
            count = used.get(name,0)
            renamed.append(name if count == 0 else name+"."+str(count))
            used[name] = count+1
        target = source["target_name"]
        removed = source["removed_identifier_columns"]
        require(features == [name for name in renamed if name not in [target]+removed], "Numeric feature/target/ID column selection differs")
        require(target in renamed and not set(features)&set(removed), "Target/identifier leaked into features")
        X = np.asarray([[float(row[renamed.index(name)]) for name in features] for row in rows],dtype=np.float64)
        labels = [row[renamed.index(target)] for row in rows]
        if removed:
            ids = {name:[row[renamed.index(name)] for row in rows] for name in removed}
            require(ids == source["source_identifiers"], "Source identifiers/repeated-observation mapping differs")
            if dataset == "parkinsons":
                require(len({value.rsplit("_",1)[0] for value in ids["name"]}) == source["subject_count"], "Parkinsons recorded identifier-prefix count differs")
            if dataset == "musk":
                require(len(set(ids["molecule_name"])) == source["molecule_count"] == 92, "Musk molecule count differs")
    require(X.shape == SHAPES[dataset], "Independent source row/column shape differs")
    return X,mapped_labels(labels,manifest["label_mapping"],np)

def check_report(results,panel,selected,permutations):
    path = results.parent/"SUMMARY.json"
    if not path.exists():
        return dict(status="not_present",note="Scientific verification does not require rendering first.")
    summary = read(path)
    if summary.get("experiment") == "2":
        approved = list(ANCHORS)+["breast_tissue","haberman"]
        audit = summary["panel_selection_audit"]
        require(summary["roster"] == list(ROSTER) and summary["approved_panel"] == approved
                and summary["all_selected_scores_available"] is True
                and audit["complete"] is True
                and audit["rule_selected_ids_this_run"] == panel["dataset_ids"]
                and audit["matches_approved_panel"] == (panel["dataset_ids"] == approved),
                "Portable fixed-panel or selection-rule audit differs")
        chosen = {(row["dataset"],row["method"]):row for row in selected}
        labels = {"manhattan":"M consistent in 5/5","euclidean":"E consistent in 5/5","mixed":"mixed/tied/order-sensitive"}
        require(len(summary["records"]) == 20 and {row["dataset"] for row in summary["records"]} == set(ROSTER),
                "Portable report omits/duplicates a screened dataset")
        for row in summary["records"]:
            dataset = row["dataset"]
            require(row["shape"] == list(SHAPES[dataset])
                    and row["ordering_classification"] == labels[panel["datasets"][dataset]["classification"]],
                    "Portable report dataset/order classification differs")
            for method in METHODS:
                require(row["methods"][method]["selected"] == chosen[dataset,method],"Portable table's selected scores differ")
        require(summary["ordering_classification_counts"] == {labels[key]:value for key,value in panel["classification_counts"].items()},
                "Portable report screen counts differ")
        require(summary["report_source_sha256"] in {sha(file) for file in HELPERS.glob("*.py")},"Portable report producer source differs")
        for source in summary["sources"]:
            require(sha(ROOT/source["file"]) == source["sha256"],"Portable report input snapshot changed")
        return dict(status="passed",summary_sha256=sha(path),approved_panel=approved,
                    recomputed_rule_panel=panel["dataset_ids"],matches_approved_panel=audit["matches_approved_panel"])
    require(summary["all_twenty_complete"] is True and summary["panel_unfilled_slots"] == 0
            and [row["dataset"] for row in summary["panel_selection"]] == panel["dataset_ids"], "Report panel differs from independent frozen-rule selection")
    require(summary["report_source_sha256"] in {sha(path) for path in HELPERS.glob("*.py")}, "Report producer source differs")
    chosen = {(row["dataset"],row["method"]):row for row in selected}
    orders = {(row["dataset"],row["method"],row["replicate"]):row for row in permutations}
    class_names = {"manhattan":"M consistent in 5/5", "euclidean":"E consistent in 5/5", "mixed":"mixed/tied/order-sensitive"}
    for row in summary["full_screen"]:
        dataset = row["dataset"]
        require(row["complete"] is True, "Report marks a complete dataset incomplete")
        for method in METHODS:
            require(row["methods"][method]["selected"] == chosen[dataset,method], "Report selected result differs")
        for key,value in (("ari",panel["datasets"][dataset]["original_delta_ari"]),("ami",panel["datasets"][dataset]["original_delta_ami"])):
            close(row["original"]["delta"][key],value,"Report original quality delta differs")
        require(row["classification"] == class_names[panel["datasets"][dataset]["classification"]], "Report ordering classification differs")
        require(len(row["permutations"]) == 5 and {item["replicate"] for item in row["permutations"]} == set(range(5)), "Report permutation records incomplete")
        for item in row["permutations"]:
            for method in METHODS:
                require(item["records"][method] == orders[dataset,method,item["replicate"]], "Report permutation record differs")
            for metric in ("ari","ami"):
                close(item["delta"][metric],item["records"]["manhattan"][metric]-item["records"]["euclidean"][metric],"Report paired quality delta differs")
        for metric in ("ari","ami"):
            bounds = panel["datasets"][dataset]["permutation_delta_"+metric+"_range"]
            require(row["permutation_delta_ranges"][metric] == dict(n=5,min=bounds[0],max=bounds[1]), "Report permutation range differs")
    require(len(summary["full_screen"]) == 20 and {row["dataset"] for row in summary["full_screen"]} == set(ROSTER), "Report omits/duplicates screened datasets")
    require(summary["screen_classification_counts"] == {class_names[key]:value for key,value in panel["classification_counts"].items()}, "Report classification counts differ")
    for source in summary["sources"]:
        require(sha(ROOT/source["file"]) == source["sha256"], "Report input snapshot changed")
    return dict(status="passed",summary_sha256=sha(path),panel_dataset_ids=panel["dataset_ids"])


def local_source_hashes():
    files = [*sorted(HELPERS.glob("*.py")), MANHATTAN]
    return {path.relative_to(REPO).as_posix():sha(path) for path in files}


def identity_checks(results,numerical_only):
    policy = read(ROOT/"protocol.json")
    require(policy["datasets"] == list(ROSTER) and policy["anchors"] == list(ANCHORS)
            and policy["seed"] == 101 and policy["radius_count"] == 32
            and policy["min_pts"] == [0,5,10,20,40,80] and policy["merge_scale"] == 1.4
            and policy["permutations"] == policy["selected_timing_repetitions"] == 5, "Scientific protocol differs")
    hashes = local_source_hashes()
    if numerical_only:
        return dict(mode="numerical_only",runtime_identity_verified=False,fresh_execution_verified=False,
                    note="Packaging fixture only: runtime/source identity and fresh-execution claims are deliberately not certified.",
                    current_local_source_sha256=hashes)
    manifest = read(results/"run_manifest.json")
    require(manifest["protocol_sha256"] == sha(ROOT/"protocol.json")
            and manifest["run_source_sha256"] in hashes.values()
            and manifest["fit_adapter_sha256"] == sha(HELPERS/"fit_adapter.py")
            and manifest["roster"] == list(ROSTER) and manifest.get("imported_mode",False) is False,
            "Fresh runtime/source/protocol identity differs")
    sys.path.insert(0,str(HELPERS))
    from fit_adapter import native_details
    native = native_details()  # Constructor/import only; no clustering.
    require(native == manifest["native"] and native["all_euclidean_kernels_native"] is True,
            "Native method identity/fallback differs")
    require(native["methods_source_sha256"] == sha(HELPERS/"methods.py"), "Fit adapter source identity differs")
    compatible = MANHATTAN
    info = native["manhattan_source"]
    require(sha(compatible) == info["compatible_sha256"], "Manhattan source fingerprint differs")
    for role,kernel in native["euclidean_kernels"].items():
        require(kernel["native"] is True and sha(kernel["file"]) == kernel["sha256"], "Compiled kernel differs: "+role)
    expected = {"numpy":"2.2.6","scipy":"1.15.3","scikit-learn":"1.6.1","classixclustering":"1.5.1",
                "pandas":"2.2.3","xlrd":"2.0.1","threadpoolctl":"3.7.0"}
    for name,version in expected.items():
        require(importlib.metadata.version(name) == version,"Pinned numerical package differs: "+name)
    require(not (results/"imported_provenance.json").exists(), "Fresh reproduction contains historical import evidence")
    return dict(mode="strict_fresh",runtime_identity_verified=True,fresh_execution_verified=True,native=native,
                versions=expected,run_manifest_sha256=sha(results/"run_manifest.json"),current_local_source_sha256=hashes)


def data_checks(dataset,np,numerical_only):
    from sklearn.preprocessing import StandardScaler
    directory = ROOT/"data"
    manifest_path,path = directory/(dataset+".manifest.json"),directory/(dataset+".npz")
    manifest = read(manifest_path)
    require(manifest["dataset"] == dataset and sha(path) == manifest["npz"]["sha256"]
            and path.stat().st_size == manifest["npz"]["bytes"], "Prepared archive hash differs")
    require(not (directory/(dataset+".import.json")).exists(), "Portable data must be freshly prepared, not imported from another workspace")
    if not numerical_only:
        require(manifest["preparation_source_sha256"] in local_source_hashes().values(), "Preparation source identity differs")
    for package,version in manifest["versions"].items():
        require(importlib.metadata.version(package) == version,"Prepared package version differs")
    with np.load(path,allow_pickle=False) as archive:
        require(set(archive.files) == set(manifest["arrays"]) == {"X_raw","X","y","row_ids"},"Prepared array fields differ")
        arrays = {key:archive[key] for key in archive.files}
    for key,value in arrays.items():
        require(array_record(value,np) == manifest["arrays"][key],"Prepared logical array hash differs")
    X,raw,y,ids = (arrays[key] for key in ("X","X_raw","y","row_ids"))
    require(X.shape == raw.shape == SHAPES[dataset] and X.dtype == raw.dtype == np.float64
            and y.shape == ids.shape == (len(X),) and y.dtype == ids.dtype == np.int64
            and np.isfinite(X).all() and np.isfinite(raw).all()
            and np.array_equal(ids,np.arange(len(X),dtype=np.int64)),"Dataset dimensions/type/row order differ")
    expected_raw,expected_y = parse_source(dataset,manifest,directory,np)
    require(np.allclose(raw,expected_raw,rtol=1e-13,atol=1e-13) and np.array_equal(y,expected_y),
            "Independent raw-source feature/target alignment differs: "+dataset)
    prep = manifest["preprocessing"]
    require(prep["embedding"] is None and prep["binarization"] is False
            and prep["parameters"] == dict(copy=True,with_mean=True,with_std=True)
            and manifest["missing_values"]["rows_removed"] == 0,"Unexpected preprocessing/row filtering")
    scaler = StandardScaler().fit(raw)
    require(np.array_equal(scaler.transform(raw),X),"StandardScaler matrix differs")
    for key,value in (("mean",scaler.mean_),("variance_ddof0",scaler.var_),("scale",scaler.scale_)):
        require(np.array_equal(prep[key],value),"StandardScaler fitted parameter differs")
    require(prep["n_samples_seen"] == len(X),"Scaler sample count differs")
    labels,counts = np.unique(y,return_counts=True)
    require(len(labels) == CLASSES[dataset] and manifest["label_counts"] == {str(int(label)):int(count) for label,count in zip(labels,counts)},
            "Observed class counts differ")
    feature_info = manifest["feature_characteristics"]
    unique = [np.unique(raw[:,i]) for i in range(raw.shape[1])]
    constant = [i for i,values in enumerate(unique) if len(values)==1]
    two = [i for i,values in enumerate(unique) if len(values)==2]
    binary = [i for i,values in enumerate(unique) if set(values.tolist()).issubset({0.,1.})]
    require(feature_info["unique_value_counts"] == [len(values) for values in unique]
            and feature_info["constant_columns"] == prep["constant_features_retained"] == constant
            and feature_info["two_level_columns"] == two and feature_info["zero_one_columns_including_constants"] == binary,
            "Numeric-feature cardinality metadata differs")
    validation_path = directory/(dataset+".validation.json")
    if validation_path.exists():
        validation = read(validation_path)
        require(validation["status"] == "ok" and validation["X_sha256"] == manifest["arrays"]["X"]["sha256"]
                and validation["duplicate_rows_beyond_first"] == len(X)-len(np.unique(X,axis=0)),"Data readiness record differs")
    return X,y,dict(rows=len(X),features=X.shape[1],classes=len(labels),npz_sha256=sha(path),
        manifest_sha256=sha(manifest_path),raw_source_alignment=True,exact_standard_scaler=True,
        raw_decimal_comparison_tolerance="rtol=atol=1e-13; prepared payload and source hashes are exact",
        duplicate_rows=len(X)-len(np.unique(X,axis=0)),removed_identifiers=manifest["source"].get("removed_identifier_columns",[]))


def main(argv=None):
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--results-dir",type=Path)
    parser.add_argument("--numerical-only",action="store_true",help="Validate a packaging fixture; does NOT certify runtime identity or a fresh run.")
    args = parser.parse_args(argv)
    ROOT = args.root.resolve()
    results = args.results_dir.resolve() if args.results_dir else ROOT/"results"
    require(results.is_relative_to(ROOT) and results != ROOT,"Results must lie inside the generated working directory")
    completion = read(results/"progress.json")
    require(completion["stage"] == "complete","Full screen has not completed")
    if not args.numerical_only:
        require(not process_alive(completion["pid"]),"Producer still alive; do not overlap verification with fitting")
    for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
        os.environ[key] = "1"
    started = time.perf_counter()
    report = dict(schema=1,status="running",mode="numerical_only" if args.numerical_only else "strict_fresh",
        generated_utc=datetime.now(timezone.utc).isoformat(),results_directory=str(results),issues=[],warnings=[],datasets={},selected=[],permutations=[],
        limitations=["No clustering or download is performed by this verifier.",
            "Full grid and individual timing-repeat labels are not saved; retained selected/permutation partitions are rescored. Timing repeat equality relies on producer assertions.",
            "Intermediate permutation-order labels are not retained; the order hash, restored assignments, and unmapping source are checked.",
            "The twenty-dataset screen and illustrative ten-case panel are exploratory and outcome-stratified. Permutations test ordering, not independent sampling.",
            "Complete implementations differ in normalization, sorting, native execution and metric."])
    provenance = {}
    try:
        import numpy as np
        from sklearn.metrics import adjusted_rand_score,adjusted_mutual_info_score
        from threadpoolctl import threadpool_limits
        metrics = dict(ari=adjusted_rand_score,ami=adjusted_mutual_info_score)
        report["runtime"] = identity_checks(results,args.numerical_only)
        grid,selected,permutations = lines(results/"grid.jsonl"),read(results/"selected.json"),lines(results/"permutation.jsonl")
        require({row["dataset"] for row in grid} == set(ROSTER),"Grid roster differs")
        require(len(selected) == 40 and {(row["dataset"],row["method"]) for row in selected} == {(dataset,method) for dataset in ROSTER for method in METHODS},"Selected schedule incomplete/duplicated")
        require(len(permutations) == 200 and {(row["dataset"],row["method"],row["replicate"]) for row in permutations} == {(dataset,method,repeat) for dataset in ROSTER for method in METHODS for repeat in range(5)},"Permutation schedule incomplete/duplicated")
        if not args.numerical_only:
            require(all(row["origin"] == "fresh" for row in grid+selected+permutations),"Strict reproduction contains imported observations")
        chosen = {(row["dataset"],row["method"]):row for row in selected}
        seen = set()
        with threadpool_limits(limits=1):
            for dataset in ROSTER:
                X,y,evidence = data_checks(dataset,np,args.numerical_only)
                design = verify_plan(dataset,X,results,np)
                rows = [row for row in grid if row["dataset"] == dataset]
                evidence["grid"] = verify_grid(dataset,X,rows,design,results,np,report["warnings"])
                report["datasets"][dataset] = evidence
                references = {}
                for method in METHODS:
                    row = chosen[dataset,method]
                    winner = rank([item for item in rows if item["method"] == method])
                    require(all(row[key] == value for key,value in winner.items()),"Selected setting differs from ARI optimum/tie-break")
                    row_checks(row,len(X))
                    labels,assignment = labels_from(row,len(X),results,seen,np)
                    scores(row,y,labels,metrics)
                    times = row["timing_seconds"]
                    require(len(times) == 5 and all(math.isfinite(value) and value>0 for value in times),"Five valid measured timing values required")
                    close(row["fit_seconds_median"],np.median(times),"Timing median differs")
                    references[method] = labels
                    report["selected"].append(dict(dataset=dataset,method=method,origin=row["origin"],ari=row["ari"],ami=row["ami"],timing_repeats=5,**assignment))
                for row in [item for item in permutations if item["dataset"] == dataset]:
                    method,repeat = row["method"],row["replicate"]
                    order = np.random.default_rng(np.random.SeedSequence([101,ROSTER.index(dataset),repeat,867])).permutation(len(X))
                    digest = hashlib.sha256(order.tobytes()).hexdigest()
                    require(row["permutation_sha256"] == digest,"Paired seeded permutation differs")
                    require(all(row[key] == chosen[dataset,method][key] for key in ("radius","radius_index","min_pts","merge_scale")),"Permutation parameters were retuned")
                    row_checks(row,len(X))
                    for key in ("preprocessing_scale","fitted_scale"):
                        close(row[key],normalization(method,X[order],np),"Permutation normalization differs")
                    labels,assignment = labels_from(row,len(X),results,seen,np)
                    scores(row,y,labels,metrics)
                    close(row["stability_ari"],adjusted_rand_score(references[method],labels),"Permutation stability differs")
                    report["permutations"].append(dict(dataset=dataset,method=method,replicate=repeat,origin=row["origin"],ari=row["ari"],ami=row["ami"],stability_ari=row["stability_ari"],permutation_sha256=digest,**assignment))
        require(seen == {path.resolve() for path in (results/"assignments").glob("*.npz")},"Saved assignments missing/unreferenced")
        require(completion["datasets"] == 20 and completion["selected_configurations"] == 40
                and completion["grid_fits"] == len(grid) and completion["permutation_fits"] == 200
                and completion["grid_origins"] == dict(Counter(row["origin"] for row in grid))
                and completion["permutation_origins"] == dict(Counter(row["origin"] for row in permutations)),"Completion counters differ")
        report["panel"] = illustrative_panel(selected,permutations)
        report["report_check"] = check_report(results,report["panel"],selected,permutations)
        failures = lines(results/"failures.jsonl") if (results/"failures.jsonl").exists() else []
        report["saved_failures"] = failures
        report["recorded_fit_warnings"] = [dict(dataset=row["dataset"],method=row["method"],warning=warning) for row in grid+permutations for warning in row["warnings"]]
        report.update(status="passed",counts=dict(datasets=20,grid_fits=len(grid),grid_origins=dict(Counter(row["origin"] for row in grid)),
            selected_configurations=40,selected_measured_times=200,permutation_fits=200,permutation_origins=dict(Counter(row["origin"] for row in permutations)),
            saved_partitions_verified=len(seen),saved_failures=len(failures),recorded_fit_warnings=len(report["recorded_fit_warnings"])))
        provenance["source_sha256"] = local_source_hashes()
        provenance["source_sha256"]["protocol.json"] = sha(ROOT/"protocol.json")
        provenance["data_sha256"] = {path.relative_to(ROOT).as_posix():sha(path) for path in sorted((ROOT/"data").rglob("*")) if path.is_file()}
        names = ["grid.jsonl","selected.json","permutation.jsonl","progress.json"]
        names += [name for name in ("run_manifest.json","readiness.json","failures.jsonl") if (results/name).exists()]
        paths = [results/name for name in names]+sorted((results/"plans").glob("*.json"))+sorted(seen)
        provenance["result_sha256"] = {path.relative_to(ROOT).as_posix():sha(path) for path in paths}
    except Exception as error:
        report["status"] = "failed"
        report["issues"].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
    report["verification_seconds"] = time.perf_counter()-started
    save(results/"verification.json",report)
    provenance.update(schema=1,generated_utc=report["generated_utc"],verification_status=report["status"],mode=report["mode"],
        verification_file=(results/"verification.json").relative_to(ROOT).as_posix(),verification_sha256=sha(results/"verification.json"))
    save(results.parent/"provenance.json",provenance)
    print(json.dumps(dict(status=report["status"],mode=report["mode"],counts=report.get("counts"),warnings=len(report["warnings"]),issues=report["issues"])))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
