"""Independently audit completed Experiment 2C records; never fit any estimator.

Run only after timing finishes. Recomputes input calibration, selected scores,
partitions, graph extraction and grid selection without importing the producer.
"""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import sys
import traceback

import numpy as np
import scipy
from scipy.spatial.distance import pdist, squareform
import sklearn
from sklearn.cluster import cluster_optics_xi
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROSTER = ["iris", "banknote", "wine", "breast_cancer", "digits", "seeds", "glass", "sonar",
          "vehicle", "segmentation", "breast_tissue", "leaf", "libras", "vertebral", "wireless",
          "parkinsons", "coimbra", "raisin", "haberman", "musk"]
PANEL = ROSTER[:8] + ["breast_tissue", "haberman"]
METRICS = ["euclidean", "manhattan"]
ALGORITHMS = ["dbscan", "optics"]
MINIMA = [2, 5, 10, 20, 40, 80]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(a, b, message, *, atol=1e-12, rtol=1e-12):
    require(np.isfinite(a) and np.isfinite(b) and np.isclose(a, b, atol=atol, rtol=rtol), message)


def canonical(labels):
    """First-occurrence ranks for non-noise groups, with noise kept at -1."""
    labels = np.asarray(labels)
    mask = labels >= 0
    output = np.full(labels.shape, -1, dtype="<i8")
    if mask.any():
        _, first, inverse = np.unique(labels[mask], return_index=True, return_inverse=True)
        ranks = np.empty(len(first), dtype=np.int64)
        ranks[np.argsort(first)] = np.arange(len(first))
        output[mask] = ranks[inverse]
    return output


def partition_hash(labels):
    return hashlib.sha256(canonical(labels).tobytes()).hexdigest()


def key(row):
    return row["dataset"], row["algorithm"], row["metric"], row["min_samples"], row["parameter_index"]


def winner(rows):
    return sorted(rows, key=lambda row: (-row["ari"], row["min_samples"], row["parameter_index"]))[0]


def load_labels(path, n):
    with np.load(path, allow_pickle=False) as archive:
        labels = archive["labels"].copy()
    require(labels.shape == (n,) and labels.dtype.kind in "iu", f"Invalid labels in {path}")
    require(np.all(labels >= -1), f"Unexpected negative labels in {path}")
    return labels


def rescore(row, y, labels):
    close(row["ari"], adjusted_rand_score(y, labels), "Selected ARI differs")
    close(row["ami"], adjusted_mutual_info_score(y, labels, average_method="arithmetic"), "Selected AMI differs")
    require(row["noise_points"] == np.count_nonzero(labels == -1), "Noise count differs")
    require(row["n_clusters"] == np.unique(labels[labels >= 0]).size, "Non-noise cluster count differs")
    close(row["retained_fraction"], np.mean(labels >= 0), "Coverage differs")
    require(row["partition_sha256"] == partition_hash(labels), "Canonical partition hash differs")


def check_row(row, n):
    require(row["status"] == "ok", "Grid contains failed status requiring review")
    require(row["algorithm"] in ALGORITHMS and row["metric"] in METRICS, "Unknown grid method")
    require(isinstance(row["min_samples"], int) and 2 <= row["min_samples"] <= n, "Invalid min_samples")
    require(isinstance(row["parameter_index"], int), "Noninteger parameter index")
    require(np.isfinite(row["parameter"]) and row["parameter"] > 0, "Invalid grid parameter")
    if row["algorithm"] == "optics":
        require(row["parameter"] < 1 and row["parameter_name"] == "xi", "Invalid xi")
        require("search_extraction_seconds" in row and "search_fit_seconds" not in row, "OPTICS search timing mislabeled")
        seconds = row["search_extraction_seconds"]
    else:
        require(row["parameter_name"] == "eps" and "search_fit_seconds" in row, "DBSCAN search timing mislabeled")
        seconds = row["search_fit_seconds"]
    require(np.isfinite(seconds) and seconds > 0, "Invalid search time")
    for name in ("ari", "ami"):
        require(np.isfinite(row[name]) and -1.000000001 <= row[name] <= 1.000000001, f"Invalid grid {name}")
    noise, groups = row["noise_points"], row["n_clusters"]
    require(isinstance(noise, int) and isinstance(groups, int) and 0 <= noise <= n and 0 <= groups <= n-noise, "Impossible grid counts")
    require((groups == 0) == (noise == n), "Zero non-noise groups disagree with noise count")
    close(row["retained_fraction"], (n-noise)/n, "Grid coverage differs")
    require(len(row["partition_sha256"]) == 64, "Malformed grid partition hash")
    if noise == n:
        close(row["ari"], 0.0, "All-noise ARI is not zero")
        close(row["ami"], 0.0, "All-noise AMI is not zero")
        require(row["partition_sha256"] == partition_hash(np.full(n, -1, dtype=np.int64)), "All-noise partition hash differs")


def check_protocol(root, prior_root, policy, provenance, report):
    require(policy["datasets"] == ROSTER and policy["panel"] == PANEL, "Dataset roster/panel differs")
    require(policy["metrics"] == METRICS and policy["seed"] == 101, "Metric/seed policy differs")
    require(policy["min_samples"] == MINIMA and policy["candidate_count"] == 32, "Initial grid policy differs")
    op = policy["optics"]
    require(op["xi_low"] == .01 and op["xi_high"] == .3 and op["cluster_method"] == "xi", "OPTICS grid differs")
    require(op["max_eps"] == "infinity" and op["memory"] is None and op["min_cluster_size"] is None
            and op["predecessor_correction"] is True, "OPTICS extraction policy differs")
    require(op["n_jobs"] == policy["dbscan"]["n_jobs"] == 1, "Thread policy differs")
    require(op["algorithm"] == policy["dbscan"]["algorithm"] == "auto", "Neighbor search differs")
    manifest = read(root/"results/run_manifest.json")
    require(manifest["protocol_sha256"] == sha(root/"protocol.json"), "Protocol changed after start")
    require(manifest["runner_sha256"] == sha(HELPERS/"run_density.py"), "Runner changed after start")
    for name, value in (("numpy", np.__version__), ("scipy", scipy.__version__), ("sklearn", sklearn.__version__)):
        require(manifest[name] == value, f"Verifier environment differs for {name}")
    require(manifest["native"]["methods_source_sha256"] == sha(HELPERS/"methods.py"), "CLASSIX adapter source differs")
    for record in manifest["sklearn_source"].values():
        require(sha(record["path"]) == record["sha256"], "Recorded sklearn source changed")
    for record in manifest["native"]["euclidean_kernels"].values():
        require(record["native"] and sha(record["file"]) == record["sha256"], "Euclidean native kernel differs")
    manhattan = read(PROVENANCE)
    require(sha(MANHATTAN) == manhattan["compatible_sha256"], "Manhattan snapshot differs")
    require(sha(root/"results/classix_selected_prior.json") == provenance["prior_classix_selected_sha256"], "Imported CLASSIX selection changed")
    require(sha(prior_root/"results/selected.json") == provenance["prior_classix_selected_sha256"], "Original CLASSIX selection differs")
    require(sha(prior_root/"protocol.json") == provenance["prior_protocol_sha256"], "Prior protocol changed")
    if "prior_classix_assignments" in provenance:
        snapshots = provenance["prior_classix_assignments"]
        require(len(snapshots) == 40, "Prior assignment snapshot map incomplete")
        for relative, expected_hash in snapshots.items():
            snapshot = (root/relative).resolve()
            require(snapshot.is_relative_to((root/"results/classix_assignments_prior").resolve()), "Prior snapshot path escapes its directory")
            require(sha(snapshot) == expected_hash, "Prior assignment snapshot hash differs")
    readiness = read(root/"results/readiness.json")
    require(readiness["status"] == "passed" and len(readiness["records"]) == 12, "Readiness record incomplete")
    expected = {(d, a, m) for d in ("iris", "segmentation", "musk") for a in ALGORITHMS for m in METRICS}
    require({(r["dataset"], r["algorithm"], r["metric"]) for r in readiness["records"]} == expected, "Readiness schedule differs")
    require(all(r["unscored"] and r["seconds"] > 0 for r in readiness["records"]), "Readiness not marked unscored")
    report["source_audit"] = {"runner_sha256": manifest["runner_sha256"], "protocol_sha256": manifest["protocol_sha256"],
                              "native_euclidean_kernels_verified": True, "sklearn_source_hashes_verified": True,
                              "runner_read_only_review": "Fresh selected fits include graph construction; search extraction timings have separate fields. Producer checks repeated partitions. No estimator fitted by verifier."}


def check_data(root, prior_root, dataset, provenance):
    path = root/"data"/(dataset+".npz")
    require(sha(path) == provenance["inputs"]["data/"+dataset+".npz"] == sha(prior_root/"data"/(dataset+".npz")), "Frozen input bytes differ: "+dataset)
    with np.load(path, allow_pickle=False) as data:
        X, y, raw = data["X"].copy(), data["y"].copy(), data["X_raw"].copy()
    require(X.dtype == np.float64 and X.ndim == 2 and y.shape == (len(X),) and np.isfinite(X).all(), "Invalid frozen input: "+dataset)
    require(raw.shape == X.shape and np.array_equal(StandardScaler().fit_transform(raw), X), "StandardScaler input differs: "+dataset)
    require(len(np.unique(y)) > 1, "Unexpected one-class reference data")
    design = read(root/"results/plans"/(dataset+".json"))
    require(design["input_sha256"] == sha(path) and design["shape"] == list(X.shape), "Calibration input differs")
    for metric in METRICS:
        distances = pdist(X, metric="cityblock" if metric == "manhattan" else "euclidean")
        dense = squareform(distances)
        # Remove only self; retain duplicate observations as other neighbors.
        np.fill_diagonal(dense, np.inf)
        near = float(np.median(np.partition(dense, 4, axis=1)[:, 4]))
        far = float(np.median(distances[distances > 0]))
        expected = dict(median_5nn=near, median_pair=far, eps_low=.35*near, eps_high=.9*far)
        for name, value in expected.items():
            close(design["methods"][metric][name], value, "Calibration differs: "+dataset+"/"+metric+"/"+name)
    return X, y, design


def check_grid(root, dataset, X, rows, design, report):
    n = len(X)
    for row in rows:
        check_row(row, n)
        low, high = (.01, .3) if row["algorithm"] == "optics" else (design["methods"][row["metric"]]["eps_low"], design["methods"][row["metric"]]["eps_high"])
        close(row["parameter"], low*(high/low)**(row["parameter_index"]/31), "Wrong metric/parameter grid value")
    expected_keys = set()
    stats = {}
    for algorithm in ALGORITHMS:
        subset = [r for r in rows if r["algorithm"] == algorithm]
        initial = [r for r in subset if 0 <= r["parameter_index"] < 32 and r["min_samples"] in MINIMA]
        require(len(initial) == 2*32*6, "Initial grid incomplete")
        initial_winners = {m: winner([r for r in initial if r["metric"] == m]) for m in METRICS}
        edge = any(r["parameter_index"] in (0, 31) for r in initial_winners.values())
        min_edge = any(r["min_samples"] == 80 for r in initial_winners.values())
        indices = list(range(-4, 36)) if edge else list(range(32))
        minima = MINIMA + ([value for value in (120, 160) if value <= n] if min_edge else [])
        extension = read(root/"results/plans"/(dataset+"_"+algorithm+"_extension.json"))
        require(extension["parameter_edge"] == edge and extension["min_samples_edge"] == min_edge
                and extension["indices"] == indices and extension["min_samples"] == minima, "Extension policy differs")
        require(extension["initial_winners"] == initial_winners, "Recorded extension winners differ")
        expected_keys.update((dataset, algorithm, metric, minimum, index) for metric in METRICS for minimum in minima for index in indices)
        stats[algorithm] = dict(candidates=len(subset), parameter_extension=edge, min_samples_extension=min_edge,
                                indices=indices, min_samples=minima)
        for metric in METRICS:
            selected = winner([r for r in subset if r["metric"] == metric])
            if selected["parameter_index"] in (min(indices), max(indices)) or selected["min_samples"] == max(minima):
                report["warnings"].append(dict(dataset=dataset, algorithm=algorithm, metric=metric,
                                               reason="Selected setting remains at an outer tested boundary; best within the tested grid only."))
    require({key(r) for r in rows} == expected_keys and len(rows) == len(expected_keys), "Full grid has missing/extra/duplicate candidates")
    return stats


def check_selected(root, prior_root, dataset, X, y, selected, prior, grid, report):
    n = len(X)
    require(len(selected) == 6 and {(r["algorithm"], r["metric"]) for r in selected} == {(a, m) for a in ["classix"]+ALGORITHMS for m in METRICS}, "Selected schedule differs")
    for row in selected:
        require(row["input_sha256"] == sha(root/"data"/(dataset+".npz")), "Selected input hash differs")
        path = (root/row["labels_path"]).resolve()
        require(path.is_relative_to((root/"results/assignments").resolve()), "Selected labels path escapes assignments")
        require(sha(path) == row["labels_file_sha256"], "Saved selected label hash differs")
        labels = load_labels(path, n)
        rescore(row, y, labels)
        times = row["timing_seconds"]
        require(len(times) == 5 and all(np.isfinite(t) and t > 0 for t in times), "Five positive selected times required")
        close(row["fit_seconds_median"], float(np.median(times)), "Selected median timing differs", atol=1e-15)
        if row["algorithm"] == "classix":
            old = prior[dataset, row["metric"]]
            require(all(row[name] == old[name] for name in ("radius", "min_pts", "merge_scale")), "CLASSIX configuration retuned")
            require(row["selection_origin"] == "existing_Experiment2_grid", "CLASSIX selection origin differs")
            close(row["prior_fit_seconds_median"], old["fit_seconds_median"], "Prior timing provenance differs")
            old_relative = Path(old["labels_path"].replace("\\", "/"))
            original_path = (prior_root/old_relative).resolve()
            require(original_path.is_relative_to((prior_root/"results/assignments").resolve()), "Prior labels escape assignments")
            snapshot_path = root/"results/classix_assignments_prior"/old_relative.name
            old_path = snapshot_path if snapshot_path.exists() else original_path
            if snapshot_path.exists() and original_path.exists():
                require(sha(snapshot_path) == sha(original_path), "Copied prior CLASSIX assignment differs from original")
            old_labels = load_labels(old_path, n)
            require(np.array_equal(canonical(labels), canonical(old_labels)), "Fresh CLASSIX partition changed: "+dataset+"/"+row["metric"])
            close(old["ari"], row["ari"], "Original CLASSIX ARI changed")
            close(old["ami"], row["ami"], "Original CLASSIX AMI changed")
            report["prior_classix_partitions"].append(dict(dataset=dataset, metric=row["metric"],
                                                           original_labels_sha256=sha(old_path), partition_equal=True))
        else:
            best = winner([r for r in grid if r["algorithm"] == row["algorithm"] and r["metric"] == row["metric"]])
            require(key(row) == key(best) and all(row[name] == best[name] for name in ("parameter", "partition_sha256", "n_clusters", "noise_points")), "Selected density winner differs")
            close(row["ari"], best["ari"], "Selected density ARI differs")
            close(row["ami"], best["ami"], "Selected density AMI differs")
            require(row["selection_origin"] == "new_Experiment2_density_grid", "Density origin differs")
            if row["algorithm"] == "optics":
                graph_path = root/"results/graphs"/f"{dataset}_{row['metric']}_{row['min_samples']}.npz"
                with np.load(graph_path, allow_pickle=False) as data:
                    graph = {name: data[name].copy() for name in data.files}
                extracted, _ = cluster_optics_xi(reachability=graph["reachability"], predecessor=graph["predecessor"],
                                                 ordering=graph["ordering"], min_samples=row["min_samples"], xi=row["parameter"],
                                                 min_cluster_size=None, predecessor_correction=True)
                require(np.array_equal(canonical(extracted), canonical(labels)), "Saved graph extraction differs from fresh selected OPTICS fit")
        report["selected"].append(dict(dataset=dataset, algorithm=row["algorithm"], metric=row["metric"], ari=row["ari"],
                                       ami=row["ami"], n_clusters=row["n_clusters"], noise_points=row["noise_points"],
                                       median_seconds=row["fit_seconds_median"], partition_sha256=row["partition_sha256"],
                                       independent_scores_verified=True, labels_file_sha256=sha(path)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DENSITY)
    parser.add_argument("--prior-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    output = args.output or root/"results/verification.json"
    report = dict(status="running", verified_utc=datetime.now(timezone.utc).isoformat(), root=str(root),
                  verifier_sha256=sha(__file__), issues=[], warnings=[], datasets={}, selected=[], prior_classix_partitions=[],
                  limitations=["No estimator fit was performed by this verifier. Fresh fit timing boundaries and equality across the five repetitions are producer assertions, supported by code review; individual repetition labels are not retained.",
                               "Full grids were checked for completeness, bounds, selection and all-noise consistency. Individual nonselected grid assignments are not retained and were not independently rescored."])
    try:
        progress = read(root/"results/progress.json")
        require(progress["stage"] == "complete", "Timing run is not complete; do not verify during active timing.")
        frozen_paths = [root/"results/grid.jsonl", root/"results/selected.json", root/"protocol.json", HELPERS/"run_density.py"]
        before = {str(path): sha(path) for path in frozen_paths}
        policy, provenance = read(root/"protocol.json"), read(root/"input_provenance.json")
        prior_root = args.prior_root or Path(provenance["prior_root"])
        if not prior_root.is_absolute():
            prior_root = root/prior_root
        prior_root = prior_root.resolve()
        report["prior_root"] = str(prior_root)
        check_protocol(root, prior_root, policy, provenance, report)
        rows = [json.loads(line) for line in (root/"results/grid.jsonl").read_text().splitlines() if line.strip()]
        selected = read(root/"results/selected.json")
        prior_rows = read(root/"results/classix_selected_prior.json")
        prior = {(r["dataset"], r["method"]): r for r in prior_rows}
        require(len(prior) == len(prior_rows) == 40, "Prior CLASSIX schedule incomplete")
        require(len({key(row) for row in rows}) == len(rows), "Duplicate full-grid candidates")
        require({r["dataset"] for r in rows} == set(ROSTER), "Grid roster differs")
        require(len(selected) == 120 and len({(r["dataset"], r["algorithm"], r["metric"]) for r in selected}) == 120, "Selected total differs")
        require(selected == [r for d in ROSTER for r in read(root/"results/selected"/(d+".json"))], "Combined and per-dataset selected files differ")
        graph_total = 0
        with threadpool_limits(limits=1):
            for dataset in ROSTER:
                X, y, design = check_data(root, prior_root, dataset, provenance)
                current_grid = [r for r in rows if r["dataset"] == dataset]
                grid_report = check_grid(root, dataset, X, current_grid, design, report)
                # Check all cached graph files, not only the selected min_samples.
                for metric in METRICS:
                    for minimum in grid_report["optics"]["min_samples"]:
                        path = root/"results/graphs"/f"{dataset}_{metric}_{minimum}.npz"
                        meta = read(path.with_suffix(".json"))
                        require(meta["input_sha256"] == design["input_sha256"] and meta["graph_sha256"] == sha(path), "Graph cache hash differs")
                        require((meta["dataset"], meta["metric"], meta["min_samples"]) == (dataset, metric, minimum), "Graph identity differs")
                        require(np.isfinite(meta["search_graph_fit_seconds"]) and meta["search_graph_fit_seconds"] > 0, "Invalid graph fit time")
                        with np.load(path, allow_pickle=False) as graph:
                            require(set(graph.files) == {"ordering", "reachability", "predecessor", "core_distances"}, "Graph payload differs")
                            require(all(graph[name].shape == (len(X),) for name in graph.files), "Graph array shape differs")
                            require(np.array_equal(np.sort(graph["ordering"]), np.arange(len(X))), "Invalid OPTICS ordering")
                        graph_total += 1
                current_selected = [r for r in selected if r["dataset"] == dataset]
                check_selected(root, prior_root, dataset, X, y, current_selected, prior, current_grid, report)
                report["datasets"][dataset] = dict(shape=list(X.shape), input_sha256=design["input_sha256"],
                                                  grid=grid_report, selected_variants=len(current_selected), input_and_scaler_verified=True)
        require(before == {str(path): sha(path) for path in frozen_paths}, "Producer records changed during verification")
        require(progress["complete_datasets"] == 20 and progress["selected_variants"] == 120 and progress["grid_candidates"] == len(rows), "Final progress record differs")
        require(len(report["prior_classix_partitions"]) == 40, "Not all prior CLASSIX partitions checked")
        report["counts"] = dict(datasets=20, selected_variants=120, fresh_recorded_fit_times=600, classix_partition_checks=40,
                                density_selected_winners=80, selected_optics_graph_partition_checks=40, cached_graphs_verified=graph_total,
                                grid_candidates=len(rows), grid_by_algorithm=dict(Counter(r["algorithm"] for r in rows)),
                                all_noise_grid_candidates=sum(r["noise_points"] == report["datasets"][r["dataset"]]["shape"][0] for r in rows),
                                all_noise_selected=sum(r["noise_points"] == report["datasets"][r["dataset"]]["shape"][0] for r in selected))
        report["record_hashes"] = before
        report["status"] = "passed"
    except Exception as exc:
        report["status"] = "failed"
        report["issues"].append({"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "output": str(output), "counts": report.get("counts"), "issues": report["issues"]}, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
