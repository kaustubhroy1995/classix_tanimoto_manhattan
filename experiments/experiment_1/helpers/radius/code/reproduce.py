"""Fresh Windows reproductions; recorded evidence is never used as new results.

Default: inspect the plan, inputs to generation, dependencies and native imports.
Use --prepare to generate data, or --run to execute the selected cases. Every
write stays beneath experiment_1A/runs/<name>. PAUSE stops at the next fit
boundary; STOP_CURRENT interrupts the current attempt and leaves it retriable.
Remove the control file and repeat the same command to resume.
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
RUN_DIR = None
METHODS = ("dbscan", "classix_packed", "bitbirch", "bblean", "butina")
PACKAGES = ("numpy", "scipy", "scikit-learn", "joblib", "cloudpickle", "narwhals",
            "rdkit", "bblean", "threadpoolctl", "tqdm")
TERMINAL = {"timeout", "memory_limit", "algorithm_error", "worker_error",
            "startup_timeout", "validation_timeout", "validation_error"}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(value, path):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for attempt in range(40):
        try:
            temporary.replace(path)
            break
        except PermissionError:
            if attempt == 39:
                raise
            time.sleep(.05)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_design():
    plan = read_configuration(CONFIG)
    policy = read_repetitions(CONFIG)
    require(tuple(plan["methods"]) == METHODS, "Expected the five recorded methods")
    cases = plan["cases"]
    require(len(cases) == len({case["case_id"] for case in cases}) == 6,
            "Expected the six-setting radius plan")
    require(len({case["data_id"] for case in cases}) == 1, "Expected one shared dataset")
    require("sample_count_100000" not in {case["case_id"] for case in cases},
            "The cancelled 100k sample case must remain excluded")
    require(plan["fit_timeout_seconds"] is None and plan["job_memory_limit_bytes"] == 20 * 1024**3,
            "The recorded resource policy is no fit cutoff and a 20-GiB job limit")
    for case in cases:
        generator = case["generator"]
        require(generator["seed"] == 101 and case["rows"] == generator["n"] * generator["num_clusters"],
                "Dataset seed or total size differs")
        if case["study"] == "radius":
            require(case["rows"] == 10000 and .3 <= case["x"] <= .8, "Radius dataset differs")
        for method in METHODS:
            ids = policy["schedules"][case["case_id"]][method]
            require(ids == sorted(set(ids)) and ids and any(value >= 0 for value in ids)
                    and all(value in (-1, 0, 1, 2) for value in ids), "Invalid recorded repetition schedule")
    return plan, policy


def select_cases(plan, requested):
    known = {case["case_id"] for case in plan["cases"]}
    require(not requested or set(requested) <= known, "Unknown case ID; --check lists available cases")
    return [case for case in plan["cases"] if not requested or case["case_id"] in requested]


def output_root(value):
    root = Path(value).resolve()
    allowed = (ROOT / "runs").resolve()
    require(root != allowed and root.is_relative_to(allowed),
            "--output-dir must be a named folder beneath experiment_1A/runs; recorded evidence is read-only")
    return root


def environment_check(plan, policy, cases, methods):
    """Import dependencies and native bindings, without constructing estimators."""
    require(os.name == "nt", "The recorded native kernels and memory protocol require 64-bit Windows")
    require(platform.machine().lower() in ("amd64", "x86_64") and sys.maxsize > 2**32,
            "Use a 64-bit x86 Windows Python with the included native binaries")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    versions = {name: importlib.metadata.version(name) for name in PACKAGES}
    sys.path.insert(0, str(REPO / "src"))
    sys.path.insert(0, str(VENDOR))
    for name in ("numpy", "scipy", "sklearn.cluster", "rdkit.DataStructs", "rdkit.ML.Cluster.Butina",
                 "bblean", "classix_t", "bitbirch.bitbirch", "memory_win", "job_limit"):
        importlib.import_module(name)
    similarity = importlib.import_module("bblean.similarity")
    require(similarity._jt_sim_arr_vec_packed.__module__ == "bblean._cpp_similarity",
            "BBLean must use its native packed-bit kernel")
    native = importlib.import_module("spmv")
    schedule_counts = {method: {
        "warmups": sum(-1 in policy["schedules"][case["case_id"]][method] for case in cases),
        "measured": sum(sum(value >= 0 for value in policy["schedules"][case["case_id"]][method]) for case in cases),
    } for method in methods}
    return {"status": "passed", "python": sys.version, "executable": sys.executable,
            "platform": platform.platform(), "machine": platform.machine(), "packages": versions,
            "hardware_popcount": bool(native.HARDWARE_POPCOUNT), "selected_cases": len(cases),
            "selected_methods": methods, "case_ids": [case["case_id"] for case in cases],
            "selected_datasets": len({case["data_id"] for case in cases}), "repetitions": schedule_counts,
            "fit_timeout_seconds": plan["fit_timeout_seconds"], "job_memory_limit_bytes": plan["job_memory_limit_bytes"],
            "note": "Imports only; no dataset generation or clustering fits were run."}


def protocol(plan, environment, pilot_rows=None):
    files = [HERE / name if (HERE / name).exists() else COMMON / name for name in ("reproduce.py", "worker.py", "prepare_dataset.py",
             "generate_data_snapshot.py", "quality.py", "state_checks.py", "memory_win.py", "job_limit.py")]
    files += [path for path in (REPO / "src").rglob("*") if path.is_file() and path.suffix in (".py", ".dll", ".c")]
    files += [path for path in (VENDOR).rglob("*.py")]
    return {"kind": "standalone_fresh_reproduction", "schema": 1,
            "plan_sha256": sha(CONFIG),
            "repetition_amendment_sha256": sha(CONFIG),
            "source_sha256": {path.relative_to(REPO).as_posix(): sha(path) for path in sorted(set(files))},
            "packages": environment["packages"], "python": sys.version,
            "platform": platform.platform(), "fit_timeout_seconds": plan["fit_timeout_seconds"],
            "job_memory_limit_bytes": plan["job_memory_limit_bytes"], "threads": 1, "pilot_rows": pilot_rows,
            "repetitions": "one measured fit per method, no warmup; diagnostic pilot" if pilot_rows else "recorded per-case schedules",
            "timing_scope": plan["timed_operations"]}


def execute(case, variant, repeat, plan, output_directory, pilot_rows=None):
    labels_path = output_directory / f"{variant}_{repeat}.npy"
    command = [sys.executable, str(HERE/"worker.py"), "--case", case["case_id"], "--variant", variant, "--output", str(labels_path), "--data-dir", str(RUN_DIR/"data"), "--config", str(CONFIG)]
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
        row = {"method": method, "status": failures[-1]["status"] if failures else "ok", "measured_repeats": len(measured),
               "warning_runs": sum(bool(r["logs"]) for r in all_runs)}
        if measured:
            values = [r["fit_done"]["clustering_seconds"] for r in measured]
            row.update(median_seconds=statistics.median(values), min_seconds=min(values), max_seconds=max(values), all_seconds=values,
                ari_values=sorted({r["validated"]["ari"] for r in measured}), cluster_counts=sorted({r["validated"]["clusters"] for r in measured}),
                aggregation_groups=sorted({r["validated"]["aggregation_groups"] for r in measured if r["validated"]["aggregation_groups"] is not None}),
                median_peak_private_bytes=statistics.median(r["fit_done"]["after_fit"]["lifetime_peak_private_commit_bytes"] for r in measured),
                median_peak_working_set_bytes=statistics.median(r["fit_done"]["after_fit"]["lifetime_peak_working_set_bytes"] for r in measured),
                private_peak_verified_in_fit_window=all(r["native_peaks_attained_in_fit_window"]["private_commit"] for r in measured),
                working_set_peak_verified_in_fit_window=all(r["native_peaks_attained_in_fit_window"]["working_set"] for r in measured))
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


def prepare(case):
    command = [sys.executable, str(HERE / "prepare_dataset.py"), case["data_id"],
               "--data-dir", str(RUN_DIR / "data"), "--config", str(CONFIG)]
    subprocess.run(command, check=True, cwd=ROOT, creationflags=subprocess.CREATE_NO_WINDOW,
                   env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))


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


def run_cases(plan, policy, cases, methods, current_protocol, pilot_rows=None):
    results_dir = RUN_DIR / "results"
    results_dir.mkdir(exist_ok=True)
    for case_index, case in enumerate(cases):
        case_id = case["case_id"]
        schedules = {method: [0] for method in METHODS} if pilot_rows else policy["schedules"][case_id]
        path = results_dir / (case_id + ".json")
        if path.exists():
            result = read(path)
            require(result["case"] == case and result["protocol"] == current_protocol,
                    "Saved case or protocol differs; use a new output folder")
        else:
            result = {"case": case, "runs": [], "pilot_rows": pilot_rows,
                      "measured_repeats_requested": 1 if pilot_rows else plan["measured_repeats"], "protocol": current_protocol,
                      "method_repeat_ids": schedules, "reproduction": True}
        partitions, state_reference = {}, None
        require(all(run["variant"] in METHODS for run in result["runs"]), "Unexpected saved method")
        for run in result["runs"]:
            state_reference = verify_run(run, current_protocol, partitions, state_reference)
        finished = {method: method_finished([run for run in result["runs"] if run["variant"] == method], schedules[method])
                    for method in METHODS}
        if all(finished[method] for method in methods):
            print("REUSE " + case_id, flush=True)
            continue
        if controlled_pause(case_id):
            return False
        prepare(case)
        result["data_manifest"] = str(RUN_DIR / "data" / case["data_id"] / "manifest.json")
        result["output_relative_data_manifest"] = (Path("data") / case["data_id"] / "manifest.json").as_posix()
        save(result, path)
        existing = {(run["variant"], run["repeat"]) for run in result["runs"]}
        for repeat in sorted({value for method in methods for value in schedules[method]}):
            offset = max(repeat, 0) % len(methods)
            for variant in methods[offset:] + methods[:offset]:
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
                result["summary"] = summarize(result["runs"], METHODS)
                result["completed"] = all(finished.values())
                save(result, path)
                print(json.dumps({"case": case_id, "variant": variant, "repeat": repeat,
                                  "status": run["status"], "seconds": run.get("fit_done", {}).get("clustering_seconds")}), flush=True)
                if run["status"] == "validation_error":
                    raise RuntimeError("Validation failed; preserved evidence requires review before continuing")
        result["summary"] = summarize(result["runs"], METHODS)
        result["completed"] = all(finished.values())
        save(result, path)
    save({"complete": True, "selected_cases": len(cases), "selected_methods": methods,
          "note": "Completion concerns this selection; per-case completed flags cover all five methods."}, RUN_DIR / "live_progress.json")
    return True


def main():
    global RUN_DIR, open_process, close_process, snapshot, assign_job, close_job
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Check dependencies and design without generating data or fitting; default")
    mode.add_argument("--prepare", action="store_true", help="Generate shared data for selected cases only")
    mode.add_argument("--run", action="store_true", help="Run/resume selected cases using the recorded repetition schedule")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "runs" / "reproduction")
    parser.add_argument("--cases", nargs="+", help="Case IDs listed by --check")
    parser.add_argument("--methods", nargs="+", choices=METHODS)
    parser.add_argument("--pilot-rows", type=int, help="Diagnostic prefix only: one measured fit per method, no warmup; requires --run and --cases")
    args = parser.parse_args()
    plan, policy = load_design()
    cases = select_cases(plan, args.cases)
    if args.pilot_rows is not None:
        require(args.run and args.cases and 0 < args.pilot_rows <= min(case["rows"] for case in cases),
                "--pilot-rows requires --run, explicit --cases, and a positive prefix within each dataset")
    methods = [method for method in METHODS if not args.methods or method in args.methods]
    RUN_DIR = output_root(args.output_dir)
    environment = environment_check(plan, policy, cases, methods)
    if args.pilot_rows:
        environment["pilot_rows"] = args.pilot_rows
        environment["repetitions"] = {method: {"warmups": 0, "measured": len(cases)} for method in methods}
        environment["pilot_note"] = "Functional smoke check only. Full arrays load before slicing, so pilot memory does not describe a standalone small dataset."
    print(json.dumps(environment, indent=2), flush=True)
    if not (args.prepare or args.run):
        return 0
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    lock = RUN_DIR / "run.lock"
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(descriptor, json.dumps({"pid": os.getpid(), "utc": datetime.now(timezone.utc).isoformat()}).encode())
    os.close(descriptor)
    try:
        current_protocol = protocol(plan, environment, args.pilot_rows)
        protocol_file = RUN_DIR / ("preparation_protocol.json" if args.prepare else "protocol.json")
        if protocol_file.exists():
            require(read(protocol_file) == current_protocol, "Output folder uses a different protocol; choose a new run name")
        else:
            save(current_protocol, protocol_file)
        save(environment, RUN_DIR / "environment.json")
        if args.prepare:
            done = set()
            for case in cases:
                if controlled_pause(case["case_id"]):
                    return 0
                if case["data_id"] not in done:
                    prepare(case)
                    done.add(case["data_id"])
            return 0
        from memory_win import open_process, close_process, snapshot
        from job_limit import assign_job, close_job
        with prevent_idle_sleep():
            environment["prevents_idle_sleep_during_run"] = True
            save(environment, RUN_DIR / "environment.json")
            run_cases(plan, policy, cases, methods, current_protocol, args.pilot_rows)
        return 0
    except BaseException as error:
        save({"status": "controller_error", "error_type": type(error).__name__, "message": str(error),
              "traceback": traceback.format_exc(), "utc": datetime.now(timezone.utc).isoformat()}, RUN_DIR / "controller_error.json")
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
