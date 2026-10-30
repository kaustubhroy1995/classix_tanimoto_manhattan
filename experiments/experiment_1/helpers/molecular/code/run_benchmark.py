"""Fresh, serial Windows chemical benchmarks with preserved records.

Default --check validates the selected design, prepared array headers, dependencies
and native imports, without fitting. --run executes only explicit case schedules.
Writes stay beneath runs/<name>. PAUSE stops at the next fit boundary;
STOP_CURRENT preserves and interrupts the current attempt. Remove the control
file and repeat the same command to resume. No fit deadline is imposed by default.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import queue
import statistics
import subprocess
import sys
import threading
import time
import traceback

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = work_root(SOURCE_ROOT)
BASE = ROOT
RUN_DIR = None
PLAN_PATH = CONFIG
DATA_DIR = BASE / "data" / "prepared"
METHODS = ("dbscan", "classix_packed", "bitbirch", "bblean", "butina")
PACKAGES = ("numpy", "scipy", "scikit-learn", "rdkit", "bblean", "threadpoolctl", "tqdm")
TERMINAL = {"timeout", "memory_limit", "algorithm_error", "worker_error",
            "startup_timeout", "validation_timeout", "validation_error"}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(value, path):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def case_methods(case):
    return case.get("methods", list(case["method_repeat_ids"]))


def safe_identifier(value):
    return isinstance(value, str) and bool(value) and all(char.isalnum() or char in "_-" for char in value)


def load_design():
    plan = read_configuration(PLAN_PATH)
    require(plan.get("methods") and len(plan["methods"]) == len(set(plan["methods"]))
            and set(plan["methods"]) <= set(METHODS), "Plan methods must be recognized and unique")
    require(plan.get("seed", 101) == 101, "Fixed experiment seed must be101")
    cases = plan["cases"]
    require(cases and len(cases) == len({case["case_id"] for case in cases}), "Cases must be nonempty and unique")
    plan.setdefault("fit_timeout_seconds", None)
    plan.setdefault("job_memory_limit_bytes", 20 * 1024**3)
    require(plan["fit_timeout_seconds"] is None or plan["fit_timeout_seconds"] > 0, "Invalid fit deadline")
    require(isinstance(plan["job_memory_limit_bytes"], int) and 0 < plan["job_memory_limit_bytes"] <= 20 * 1024**3,
            "Job memory limit must be positive and no greater than 20 GiB")
    if "butina" in plan["methods"]:
        require(plan.get("butina_implementation") == "streaming_exact_static_rdkit", "Explicitly label the streaming exact static-order Butina implementation")
    parameter_keys = {"classix_packed": "classix", "classix_euclidean": "classix_euclidean",
                      "bitbirch": "bitbirch", "bblean": "bitbirch", "dbscan": "dbscan", "butina": "butina"}
    for case in cases:
        require(safe_identifier(case["case_id"]) and safe_identifier(case["data_id"]), "Unsafe case or dataset ID")
        require(all(isinstance(case[key], int) and case[key] > 0 for key in ("rows", "features")), "Invalid dimensions")
        require(type(case.get("seed", plan.get("seed", 101))) is int and case.get("seed", plan.get("seed", 101)) == 101, "Fixed-case seed must be101")
        methods = case_methods(case)
        require(methods and len(methods) == len(set(methods)) and set(methods) <= set(plan["methods"]), "Invalid case method subset")
        require(set(case["method_repeat_ids"]) == set(methods), "Schedule keys must match case methods")
        for method in methods:
            require(isinstance(case[parameter_keys[method]], dict), "Missing method parameters")
            ids = case["method_repeat_ids"][method]
            require(isinstance(ids, list) and ids and all(type(value) is int and value >= -1 for value in ids)
                    and ids == sorted(set(ids)) and any(value >= 0 for value in ids), "Invalid repetition schedule")
            measured = [value for value in ids if value >= 0]
            require(measured == list(range(len(measured))), "Measured repeats must start at zero without gaps")
    return plan


def validate_prepared(case, methods):
    """Read small manifests and NPY headers only; no data generation or array scan."""
    import numpy as np
    directory = DATA_DIR / case["data_id"]
    manifest_path = directory / "manifest.json"
    manifest = read(manifest_path)
    require(manifest["shape"] == [case["rows"], case["features"]], "Prepared dimensions differ from case")
    if "data_id" in manifest:
        require(manifest["data_id"] == case["data_id"], "Prepared dataset ID differs")
    kinds = {"classix_packed": "packed_little", "classix_euclidean": "float64" if "float64" in manifest["files"] else "int32",
             "bitbirch": "int64", "bblean": "packed_big", "dbscan": "packed_little", "butina": "packed_little"}
    required = {kinds[method] for method in methods} | {"int32"}
    if manifest["files"].get("truth"):
        required.add("truth")
    file_records = {}
    for kind in sorted(required):
        record = manifest["files"][kind]
        path = (directory / record["file"]).resolve()
        require(path.is_relative_to(directory.resolve()) and path.is_file(), "Missing or unsafe prepared array path")
        with path.open("rb") as stream:
            version = np.lib.format.read_magic(stream)
            require(version in ((1, 0), (2, 0)), "Expected a numeric NPY v1/v2 array")
            reader = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
            shape, fortran, dtype = reader(stream)
            offset = stream.tell()
        require(not dtype.hasobject and not fortran, "Arrays must be C-contiguous and not use pickle")
        expected_shape = (case["rows"],) if kind == "truth" else (case["rows"], (case["features"]+7)//8 if kind.startswith("packed_") else case["features"])
        require(shape == expected_shape, "Prepared array shape differs: " + kind)
        expected_dtype = {"int32": "int32", "int64": "int64", "float64": "float64", "packed_little": "uint8", "packed_big": "uint8"}.get(kind)
        require(expected_dtype is None or dtype == np.dtype(expected_dtype), "Prepared dtype differs: " + kind)
        payload = int(dtype.itemsize)
        for size in shape:
            payload *= size
        require(path.stat().st_size == offset + payload and record["bytes"] == payload, "Prepared array size differs")
        require(isinstance(record.get("sha256"), str) and len(record["sha256"]) == 64, "Missing logical array hash")
        file_records[kind] = {"shape": list(shape), "dtype": str(dtype), "payload_bytes": payload}
    return {"data_id": case["data_id"], "manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
            "truth_available": bool(manifest["files"].get("truth")), "arrays": file_records,
            "note": "Header and manifest validation; full logical input hash is checked after each fit."}


def select_cases(plan, requested):
    known = {case["case_id"] for case in plan["cases"]}
    require(not requested or set(requested) <= known, "Unknown case ID; --check lists available cases")
    return [case for case in plan["cases"] if not requested or case["case_id"] in requested]


def output_root(value):
    root = Path(value).resolve()
    allowed = (ROOT / "runs").resolve()
    require(root != allowed and root.is_relative_to(allowed),
            "--output-dir must be a named folder beneath runs; recorded evidence is read-only")
    return root


def environment_check(plan, cases, methods):
    """Import selected dependencies, validate input headers; never fit."""
    require(os.name == "nt" and platform.machine().lower() in ("amd64", "x86_64") and sys.maxsize > 2**32,
            "Use 64-bit x86 Windows with the recorded native binaries")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    versions = {name: importlib.metadata.version(name) for name in PACKAGES}
    sys.path.insert(0, str(REPO / "src"))
    sys.path.insert(0, str(VENDOR))
    modules = {"numpy", "scipy", "memory_win", "job_limit"}
    lookup = {"classix_packed": ("classix_t",), "dbscan": ("sklearn.cluster", "packed_dbscan"),
              "bitbirch": ("bitbirch.bitbirch",), "bblean": ("bblean",),
              "butina": ("rdkit.DataStructs", "rdkit.ML.Cluster.Butina", "streaming_butina")}
    for method in methods:
        modules.update(lookup[method])
    for name in sorted(modules):
        importlib.import_module(name)
    if "bblean" in methods:
        similarity = importlib.import_module("bblean.similarity")
        require(similarity._jt_sim_arr_vec_packed.__module__ == "bblean._cpp_similarity", "BBLean native packed-bit kernel required")
    schedule_counts = {method: {"warmups": sum(-1 in case["method_repeat_ids"].get(method, []) for case in cases),
        "measured": sum(sum(value >= 0 for value in case["method_repeat_ids"].get(method, [])) for case in cases)} for method in methods}
    prepared = [validate_prepared(case, [method for method in methods if method in case_methods(case)]) for case in cases]
    hardware = None
    if "classix_packed" in methods or "dbscan" in methods:
        hardware = bool(importlib.import_module("spmv").HARDWARE_POPCOUNT)
    return {"status": "passed", "prepared_data": prepared, "python": sys.version, "executable": sys.executable,
            "platform": platform.platform(), "machine": platform.machine(), "processor": platform.processor(), "packages": versions,
            "hardware_popcount": hardware, "selected_cases": len(cases), "selected_methods": methods,
            "case_ids": [case["case_id"] for case in cases], "selected_datasets": len({case["data_id"] for case in cases}),
            "repetitions": schedule_counts, "fit_timeout_seconds": plan["fit_timeout_seconds"],
            "job_memory_limit_bytes": plan["job_memory_limit_bytes"],
            "note": "Imports/header checks only; no fingerprint generation or clustering."}


def protocol(plan, environment, pilot_rows=None):
    names = ("run_benchmark.py", "worker.py", "memory_win.py", "job_limit.py", "quality.py",
             "state_checks.py", "packed_dbscan.py", "streaming_butina.py")
    files = [HERE / name if (HERE / name).exists() else COMMON / name for name in names]
    files += [path for path in (REPO / "src").rglob("*") if path.is_file() and path.suffix in (".py", ".dll", ".c")]
    external = [path for path in (VENDOR).rglob("*") if path.is_file() and path.suffix in (".py", ".pyd", ".dll")]
    # Record installed native and Python implementations without copying vendor trees.
    import bblean, sklearn, rdkit
    for module in (bblean, sklearn, rdkit):
        directory = Path(module.__file__).resolve().parent
        if module.__name__ == "bblean":
            external += [path for path in directory.rglob("*") if path.is_file() and path.suffix in (".py", ".pyd", ".dll")]
        elif module.__name__ == "sklearn":
            external += [path for path in (directory/"cluster").glob("*dbscan*") if path.is_file() and path.suffix in (".py", ".pyd")]
        else:
            external += [path for path in (directory/"DataStructs").glob("*") if path.is_file() and path.suffix in (".py", ".pyd")]
            external.append(directory/"ML/Cluster/Butina.py")
    return {"kind": "landrum_fixed_benchmark", "schema": 1, "plan_sha256": sha(PLAN_PATH),
            "source_sha256": {path.relative_to(REPO).as_posix(): sha(path) for path in sorted(set(files))},
            "external_source_sha256": {str(path): sha(path) for path in sorted(set(external))},
            "prepared_manifests": {item["data_id"]: item["manifest_sha256"] for item in environment["prepared_data"]},
            "packages": environment["packages"], "python": sys.version, "platform": platform.platform(),
            "fit_timeout_seconds": plan["fit_timeout_seconds"], "job_memory_limit_bytes": plan["job_memory_limit_bytes"],
            "threads": 1, "pilot_rows": pilot_rows,
            "repetitions": "one measured fit per method, no warmup; diagnostic prefix" if pilot_rows else "explicit per-case method_repeat_ids",
            "timing_scope": plan.get("timed_operations", "Estimator construction, fit and complete assignment extraction; input loading, fingerprint representation preparation and validation excluded."),
            "input_order": "Prepared dataset order, unchanged", "default_random_seed": plan.get("seed", 101),
            "butina_implementation": plan.get("butina_implementation"),
            "dbscan_implementation": "sklearn_precomputed_exact_packed_radius_graph",
            "dbscan_note": "All pairwise packed Jaccard radius-graph construction and sklearn DBSCAN fitting are included in fit timing; no approximate neighbors.",
            "butina_note": "Streaming exact fixed-order Butina neighbor-count/assignment implementation, not stock RDKit condensed-distance storage. No new Butina optimization is claimed for the other methods."}


def execute(case, variant, repeat, plan, output_directory, pilot_rows=None):
    labels_path = output_directory / f"{variant}_{repeat}.npy"
    command = [sys.executable, str(HERE/"worker.py"), "--case", case["case_id"], "--variant", variant, "--output", str(labels_path), "--data-dir", str(DATA_DIR), "--config", str(PLAN_PATH)]
    if pilot_rows:
        command += ["--rows", str(pilot_rows)]
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", bufsize=1, env=environment, creationflags=subprocess.CREATE_NO_WINDOW)
    messages = queue.Queue()
    def reader():
        for line in process.stdout:
            messages.put(line.rstrip())
        messages.put(None)
    threading.Thread(target=reader, daemon=True).start()
    result = {"case_id": case["case_id"], "variant": variant, "repeat": repeat, "warmup": repeat < 0,
              "status": "running", "logs": [], "limits": {"fit_seconds": plan["fit_timeout_seconds"], "job_memory_bytes": plan["job_memory_limit_bytes"]}}
    phase, phase_start = "startup", time.perf_counter()
    job, handle = None, None
    sampled = {"samples": 0, "target_interval_seconds": .02, "max_interval_seconds": 0,
               "peak_private_commit_bytes": 0, "peak_working_set_bytes": 0}
    last_sample_time = None
    def sample():
        nonlocal last_sample_time
        now, values = time.perf_counter(), snapshot(handle)
        if last_sample_time is not None:
            sampled["max_interval_seconds"] = max(sampled["max_interval_seconds"], now-last_sample_time)
        last_sample_time = now
        sampled["samples"] += 1
        for metric in ("private_commit", "working_set"):
            key = "peak_"+metric+"_bytes"
            sampled[key] = max(sampled[key], values[metric+"_bytes"])
    heartbeat = time.perf_counter()
    last_control_check = time.perf_counter() - 2
    try:
        while True:
            elapsed = time.perf_counter()-phase_start
            if job is not None and time.perf_counter() - last_control_check >= 1:
                last_control_check = time.perf_counter()
                if (RUN_DIR / "STOP_CURRENT").exists():
                    result["status"] = "interrupted"
                    result["interruption_phase"] = phase
                    result["elapsed_in_phase_seconds"] = elapsed
                    break
            if time.perf_counter() - heartbeat >= 15:
                save({"case": case["case_id"], "variant": variant, "repeat": repeat,
                      "phase": phase, "elapsed_in_phase_seconds": elapsed,
                      "sampled_memory": sampled, "updated_utc": datetime.now(timezone.utc).isoformat()}, RUN_DIR/"active_fit.json")
                heartbeat = time.perf_counter()
            deadline = plan["fit_timeout_seconds"] if phase == "fit" else 300
            if deadline is not None and elapsed > deadline:
                result["status"] = "timeout" if phase == "fit" else phase+"_timeout"
                result["timeout_phase"] = phase
                result["elapsed_in_phase_seconds"] = elapsed
                break
            if phase == "fit" and handle:
                try:
                    sample()
                except OSError:
                    pass  # Read queued final/error output if the worker just exited.
            try:
                line = messages.get(timeout=.02 if phase == "fit" else .1)
            except queue.Empty:
                if process.poll() is not None and not messages.qsize():
                    break
                continue
            if line is None:
                break
            if not line.startswith("BENCH_JSON "):
                result["logs"].append(line)
                continue
            event = json.loads(line[len("BENCH_JSON "):])
            kind = event.pop("event")
            result[kind] = event
            if kind == "ready":
                # Monitor/constrain the actual Python PID, not its venv launcher.
                actual_pid = event["pid"]
                job = assign_job(pid=actual_pid, max_bytes=plan["job_memory_limit_bytes"])
                handle = open_process(actual_pid)
                result["monitored_pid"] = actual_pid
                result["launcher_pid"] = process.pid
                sample()
                phase, phase_start = "fit", time.perf_counter()
                process.stdin.write("GO\n")
                process.stdin.flush()
            elif kind == "fit_done":
                sample()
                phase, phase_start = "validation", time.perf_counter()
                process.stdin.write("VALIDATE\n")
                process.stdin.flush()
            elif kind == "fit_error":
                is_memory = "MemoryError" in event["error_type"] or "Unable to allocate" in event["message"]
                result["status"] = "memory_limit" if is_memory else "algorithm_error"
            elif kind == "validation_error":
                result["status"] = "validation_error"
            elif kind == "validated":
                result["status"] = "ok"
        if result["status"] == "running":
            result["status"] = "worker_error"
    finally:
        if job is not None:
            close_job(job)
        if process.poll() is None:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
        result["exit_code"] = process.wait(timeout=10)
        if handle:
            close_process(handle)
    if sampled["samples"]:
        result["sampled_memory"] = sampled
    if "fit_done" in result:
        result["native_peaks_attained_in_fit_window"] = {metric:
            result["fit_done"]["after_fit"]["lifetime_peak_"+metric+"_bytes"] > result["ready"]["before_fit"]["lifetime_peak_"+metric+"_bytes"]
            for metric in ("private_commit", "working_set")}
    save({"case": case["case_id"], "variant": variant, "repeat": repeat,
          "finished": True, "status": result["status"],
          "updated_utc": datetime.now(timezone.utc).isoformat()}, RUN_DIR/"active_fit.json")
    return result


def summarize(runs, methods):
    output = []
    for method in methods:
        all_runs = [r for r in runs if r["variant"] == method]
        failures = [r for r in all_runs if r["status"] != "ok"]
        measured = [r for r in all_runs if r["status"] == "ok" and not r["warmup"]]
        row = {"method": method, "status": failures[-1]["status"] if failures else "ok" if measured else "pending", "measured_repeats": len(measured),
               "warning_runs": sum(bool(r["logs"]) for r in all_runs)}
        if measured:
            values = [r["fit_done"]["clustering_seconds"] for r in measured]
            row.update(median_seconds=statistics.median(values), min_seconds=min(values), max_seconds=max(values), all_seconds=values,
                ami_values=sorted({r["validated"].get("ami") for r in measured if r["validated"].get("ami") is not None}),
                ari_values=sorted({r["validated"]["ari"] for r in measured if r["validated"]["ari"] is not None}), cluster_counts=sorted({r["validated"]["clusters"] for r in measured}),
                aggregation_groups=sorted({r["validated"]["aggregation_groups"] for r in measured if r["validated"]["aggregation_groups"] is not None}),
                median_peak_private_bytes=statistics.median(r["fit_done"]["after_fit"]["lifetime_peak_private_commit_bytes"] for r in measured),
                median_peak_working_set_bytes=statistics.median(r["fit_done"]["after_fit"]["lifetime_peak_working_set_bytes"] for r in measured),
                private_peak_verified_in_fit_window=all(r["native_peaks_attained_in_fit_window"]["private_commit"] for r in measured),
                working_set_peak_verified_in_fit_window=all(r["native_peaks_attained_in_fit_window"]["working_set"] for r in measured))
            row["truth_available"] = all(r["validated"].get("truth_available", False) for r in measured)
            row["ari_unavailable_reasons"] = sorted({r["validated"]["ari_unavailable_reason"] for r in measured if r["validated"].get("ari_unavailable_reason")})
            row["noise_counts"] = sorted({r["validated"]["noise_points"] for r in measured})
        output.append(row)
    return output


def saved_path(checked, kind):
    relative = checked.get("output_relative_" + kind)
    return RUN_DIR / relative if relative else Path(checked[kind])


def method_finished(runs, ids):
    if not runs:
        return False
    sequence = [run["repeat"] for run in runs]
    require(sequence == ids[:len(sequence)], "Saved repeats are not a unique ordered prefix of the recorded schedule")
    require(all(run["warmup"] == (run["repeat"] == -1) for run in runs), "Warmup flag differs")
    failed = [index for index, run in enumerate(runs) if run["status"] != "ok"]
    if failed:
        require(failed == [len(runs) - 1] and runs[-1]["status"] in TERMINAL,
                "Invalid terminal failure or additional attempts after failure")
        return True
    return sequence == ids


def verify_run(run, current_protocol, partitions, state_reference):
    require(run["protocol"] == current_protocol, "Existing run uses a different reproduction protocol; select a new output folder")
    if run["status"] != "ok":
        return state_reference
    checked = run["validated"]
    require(saved_path(checked, "assignments_file").is_file(), "Saved assignments are missing")
    partition = checked["canonical_partition_sha256"]
    require(partitions.setdefault(run["variant"], partition) == partition, "Partition changed between repetitions")
    if checked["state_file"]:
        state = read(saved_path(checked, "state_file"))
        if state_reference is None:
            state_reference = state
        require(state == state_reference, "CLASSIX state changed between repetitions")
        checked["classix_state_repeats_consistent"] = True
    return state_reference


@contextmanager
def prevent_idle_sleep():
    """Match the historical runner's Windows power request outside fit timing."""
    import ctypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.SetThreadExecutionState.argtypes = [ctypes.c_uint32]
    kernel.SetThreadExecutionState.restype = ctypes.c_uint32
    if not kernel.SetThreadExecutionState(0x80000001):
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        yield
    finally:
        kernel.SetThreadExecutionState(0x80000000)


def controlled_pause(case_id=None):
    if any((RUN_DIR / name).exists() for name in ("PAUSE", "STOP_CURRENT")):
        save({"paused": True, "case": case_id, "utc": datetime.now(timezone.utc).isoformat()}, RUN_DIR / "live_progress.json")
        print("Paused. Remove PAUSE/STOP_CURRENT in this output folder and rerun the same command to resume.", flush=True)
        return True
    return False


def run_cases(plan, cases, methods, current_protocol, pilot_rows=None):
    results_dir = RUN_DIR / "results"
    results_dir.mkdir(exist_ok=True)
    for case_index, case in enumerate(cases):
        case_id = case["case_id"]
        available_methods = case_methods(case)
        selected_methods = [method for method in methods if method in available_methods]
        schedules = {method: [0] for method in available_methods} if pilot_rows else case["method_repeat_ids"]
        path = results_dir / (case_id + ".json")
        if path.exists():
            result = read(path)
            require(result["case"] == case and result["protocol"] == current_protocol,
                    "Saved case or protocol differs; use a new output folder")
        else:
            result = {"case": case, "runs": [], "pilot_rows": pilot_rows,
                      "measured_repeats_requested": {method: sum(value >= 0 for value in ids) for method, ids in schedules.items()}, "protocol": current_protocol,
                      "method_repeat_ids": schedules, "reproduction": True}
        prepared = validate_prepared(case, selected_methods)
        require(result.get("data_manifest_sha256", prepared["manifest_sha256"]) == prepared["manifest_sha256"], "Prepared manifest changed since the saved results")
        require(result["method_repeat_ids"] == schedules, "Saved repetition schedule differs")
        partitions, state_reference = {}, None
        require(all(run["variant"] in available_methods for run in result["runs"]), "Unexpected saved method")
        for run in result["runs"]:
            state_reference = verify_run(run, current_protocol, partitions, state_reference)
        finished = {method: method_finished([run for run in result["runs"] if run["variant"] == method], schedules[method])
                    for method in available_methods}
        if all(finished[method] for method in selected_methods):
            print("REUSE " + case_id, flush=True)
            continue
        if controlled_pause(case_id):
            return False
        result["data_manifest"] = prepared["manifest"]
        result["data_manifest_sha256"] = prepared["manifest_sha256"]
        result["output_relative_data_manifest"] = Path(os.path.relpath(prepared["manifest"], RUN_DIR)).as_posix()
        save(result, path)
        existing = {(run["variant"], run["repeat"]) for run in result["runs"]}
        for repeat in sorted({value for method in selected_methods for value in schedules[method]}):
            offset = max(repeat, 0) % len(selected_methods)
            for variant in selected_methods[offset:] + selected_methods[:offset]:
                if finished[variant] or repeat not in schedules[variant] or (variant, repeat) in existing:
                    continue
                if controlled_pause(case_id):
                    return False
                stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                assignments = results_dir / (case_id + "_assignments") / f"{variant}_{repeat}_{stamp}"
                save({"case_index": case_index + 1, "case_total": len(cases), "case": case_id,
                      "variant": variant, "repeat": repeat, "started_utc": datetime.now(timezone.utc).isoformat()},
                     RUN_DIR / "live_progress.json")
                run = execute(case, variant, repeat, plan, assignments, pilot_rows)
                run["protocol"] = current_protocol
                if "validated" in run:
                    checked = run["validated"]
                    for kind in ("assignments_file", "state_file"):
                        if checked.get(kind):
                            checked["output_relative_" + kind] = Path(checked[kind]).relative_to(RUN_DIR).as_posix()
                if run["status"] == "interrupted":
                    attempts = RUN_DIR / "interrupted_attempts"
                    attempts.mkdir(exist_ok=True)
                    save(run, attempts / f"{case_id}_{variant}_{repeat}_{stamp}.json")
                    save({"paused": True, "interrupted_current_fit": True, "case": case_id,
                          "variant": variant, "repeat": repeat}, RUN_DIR / "live_progress.json")
                    print("Interrupted attempt preserved; rerun after removing the control file to retry.", flush=True)
                    return False
                try:
                    state_reference = verify_run(run, current_protocol, partitions, state_reference)
                except Exception as error:
                    run["status"] = "validation_error"
                    run["validation_error"] = {"error_type": type(error).__name__, "message": str(error)}
                result["runs"].append(run)
                existing.add((variant, repeat))
                finished[variant] = method_finished([value for value in result["runs"] if value["variant"] == variant], schedules[variant])
                result["summary"] = summarize(result["runs"], available_methods)
                result["completed"] = all(finished.values())
                save(result, path)
                print(json.dumps({"case": case_id, "variant": variant, "repeat": repeat,
                                  "status": run["status"], "seconds": run.get("fit_done", {}).get("clustering_seconds")}), flush=True)
                if run["status"] == "validation_error":
                    raise RuntimeError("Validation failed; preserved evidence requires review before continuing")
        result["summary"] = summarize(result["runs"], available_methods)
        result["completed"] = all(finished.values())
        save(result, path)
    save({"complete": True, "selected_cases": len(cases), "selected_methods": methods,
          "note": "Completion concerns this selection; per-case completed flags cover the methods explicitly assigned to that case."}, RUN_DIR / "live_progress.json")
    return True


def main():
    global RUN_DIR, PLAN_PATH, DATA_DIR, open_process, close_process, snapshot, assign_job, close_job
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Check headers, dependencies and design without fitting; default")
    mode.add_argument("--run", action="store_true", help="Run/resume the explicitly scheduled selected cases")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=BASE / "data" / "prepared")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "runs" / "benchmark")
    parser.add_argument("--cases", nargs="+", help="Case IDs listed by --check")
    parser.add_argument("--methods", nargs="+", choices=METHODS)
    parser.add_argument("--rows", "--pilot-rows", dest="pilot_rows", type=int,
                        help="Diagnostic prefix only: one measured fit, no warmup; requires --run and --cases")
    args = parser.parse_args()
    PLAN_PATH, DATA_DIR = args.config.resolve(), args.data_dir.resolve()
    require(DATA_DIR == (BASE / "data" / "prepared").resolve(), "Only the prepared Landrum inputs from this generated work tree are supported")
    plan = load_design()
    require(not args.methods or set(args.methods) <= set(plan["methods"]), "Requested method is absent from this plan")
    methods = [method for method in plan["methods"] if not args.methods or method in args.methods]
    require(methods, "No requested methods are present in this plan")
    cases = [case for case in select_cases(plan, args.cases) if set(case_methods(case)) & set(methods)]
    require(cases, "No selected case uses the requested methods")
    if args.pilot_rows is not None:
        require(args.run and args.cases and 0 < args.pilot_rows <= min(case["rows"] for case in cases),
                "--rows requires --run, explicit --cases, and a positive prefix within each dataset")
    RUN_DIR = output_root(args.output_dir)
    environment = environment_check(plan, cases, methods)
    if args.pilot_rows:
        environment["pilot_rows"] = args.pilot_rows
        environment["repetitions"] = {method: {"warmups": 0, "measured": sum(method in case_methods(case) for case in cases)} for method in methods}
        environment["pilot_note"] = "Functional smoke check only. Full arrays load before slicing, so pilot memory does not describe a standalone small dataset."
    print(json.dumps(environment, indent=2), flush=True)
    if not args.run:
        return 0
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    lock = RUN_DIR / "run.lock"
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(descriptor, json.dumps({"pid": os.getpid(), "utc": datetime.now(timezone.utc).isoformat()}).encode())
    os.close(descriptor)
    try:
        current_protocol = protocol(plan, environment, args.pilot_rows)
        protocol_file = RUN_DIR / "protocol.json"
        if protocol_file.exists():
            require(read(protocol_file) == current_protocol, "Output folder uses a different protocol; choose a new run name")
        else:
            save(current_protocol, protocol_file)
        save(environment, RUN_DIR / "environment.json")
        from memory_win import open_process, close_process, snapshot
        from job_limit import assign_job, close_job
        with prevent_idle_sleep():
            environment["prevents_idle_sleep_during_run"] = True
            save(environment, RUN_DIR / "environment.json")
            run_cases(plan, cases, methods, current_protocol, args.pilot_rows)
        return 0
    except BaseException as error:
        save({"status": "controller_error", "error_type": type(error).__name__, "message": str(error),
              "traceback": traceback.format_exc(), "utc": datetime.now(timezone.utc).isoformat()}, RUN_DIR / "controller_error.json")
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
