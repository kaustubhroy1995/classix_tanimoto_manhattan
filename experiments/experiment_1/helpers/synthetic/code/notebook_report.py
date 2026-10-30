from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import json
import math
from pathlib import Path
import re
import statistics
import textwrap

MIB = 1024 ** 2


METHODS = {
    "classix_packed": ("CLASSIX_T (packed)", "red", "x", "--"),
    "dbscan": ("DBSCAN", "blue", "o", "-."),
    "bitbirch": ("BitBIRCH", "green", "v", ":"),
    "bblean": ("BitBIRCH-Lean (reference)", "#9467bd", "s", "--"),
    "butina": ("Taylor-Butina", "#e68613", "D", "-"),
}


STUDIES = (
    ("sample_count", "Scaling with number of data points", "Number of points", "log"),
    ("popcount_separation", "Scaling with separability of data", "Max pop of seed vector", "log"),
    ("cluster_count", "Scaling with number of clusters", "Number of clusters", "linear"),
    ("dimension", "Scaling with dimensionality", "Number of dimensions", "linear"),
    ("radius", "Sensitivity to Tanimoto radius", "Tanimoto radius r", "log"),
)


WARNING = re.compile(r"\b[A-Za-z]*Warning\b")


STATUS_MARKS = {
    "memory_limit": ("M", "memory limit"),
    "timeout": ("T", "fit timeout"),
    "pending": ("P", "pending"),
    "running": ("P", "pending"),
    "incomplete_repeats": ("I", "incomplete repeats"),
    "startup_timeout": ("ST", "startup timeout"),
    "validation_timeout": ("VT", "validation timeout"),
    "algorithm_error": ("A", "algorithm error"),
    "validation_error": ("V", "validation error"),
    "worker_error": ("W", "worker error"),
}


SHORT_NAMES = {"classix_packed": "CLASSIX_T", "dbscan": "DBSCAN", "bitbirch": "BitBIRCH",
               "bblean": "Lean", "butina": "Butina"}


def status_mark(status):
    return STATUS_MARKS.get(status, ("E", "other error"))


def method_record(case, raw, method, expected_repeats, expected_ids=None):
    """Retain censored/partial evidence but require the full protocol for curves."""
    expected_ids = list(expected_ids) if expected_ids is not None else [-1] + list(range(expected_repeats))
    if expected_ids != sorted(set(expected_ids)) or not expected_ids:
        raise ValueError("Invalid expected repetition IDs")
    expected_measured_ids = [value for value in expected_ids if value >= 0]
    expected_warmup_ids = [value for value in expected_ids if value < 0]
    if not expected_measured_ids:
        raise ValueError("A successful schedule must contain a measured fit")
    expected_repeats = len(expected_measured_ids)
    runs = [run for run in raw.get("runs", []) if run["variant"] == method]
    if any(run.get("reused_from") for run in runs):
        raise ValueError("The final seed-101 benchmark cannot reuse older runs")
    if any(run["repeat"] not in expected_ids or bool(run["warmup"]) != (run["repeat"] < 0) for run in runs):
        raise ValueError(f"Run falls outside the approved schedule: {case['case_id']} / {method}")
    measured = [run for run in runs if run["status"] == "ok" and not run["warmup"]]
    warmups = [run for run in runs if run["status"] == "ok" and run["warmup"]]
    failures = [run for run in runs if run["status"] != "ok"]
    repeats = sorted(run["repeat"] for run in measured)
    full_repeats = repeats == expected_measured_ids and sorted(run["repeat"] for run in warmups) == expected_warmup_ids
    status = failures[-1]["status"] if failures else (
        "ok" if full_repeats else "incomplete_repeats" if raw.get("completed") else "pending"
    )
    eligible = status == "ok" and full_repeats
    warnings = [line for run in runs for line in run.get("logs", []) if WARNING.search(line)]
    record = {
        "case_id": case["case_id"], "study": case["study"], "x": case["x"],
        "method": method, "status": status, "eligible_for_curves": eligible,
        "expected_repeats": expected_repeats, "measured_repeats": len(measured),
        "scheduled_repeat_ids": expected_ids, "expected_measured_repeat_ids": expected_measured_ids,
        "expected_warmup_repeat_ids": expected_warmup_ids,
        "measured_repeat_ids": repeats, "expected_discarded_warmups": len(expected_warmup_ids),
        "discarded_successful_warmups": len(warmups),
        "timing_statistic": None, "runtime_variability_estimated": False,
        "warning_runs": sum(any(WARNING.search(line) for line in run.get("logs", [])) for run in runs),
        "warning_messages": sorted(set(warnings)),
        "observed_measured_seconds": [], "observed_ari": [],
        "median_seconds": None, "min_seconds": None, "max_seconds": None,
        "median_ari": None, "ari_values": [], "cluster_counts": [], "noise_counts": [],
        "aggregation_groups": [], "adjacency_connection_buffer_bytes": [],
        "minpts_marker_counts": [], "minpts_marker_object_estimates_bytes": [],
        "median_os_peak_private_mib": None, "median_os_peak_working_set_mib": None,
        "median_prefit_private_mib": None, "median_extra_private_mib": None,
        "private_peak_verified_in_fit_window": None,
        "working_set_peak_verified_in_fit_window": None,
        "partition_stable_across_successful_repeats": None,
        "terminal_evidence": [],
    }
    partitions, logical_inputs = set(), set()
    for run in measured:
        timing = run["fit_done"]["clustering_seconds"]
        quality = run["validated"]["ari"]
        if not math.isfinite(timing) or timing <= 0 or not math.isfinite(quality):
            raise ValueError(f"Invalid successful result for {case['case_id']} / {method}")
        if run["validated"].get("input_unchanged") is False:
            raise ValueError("Successful run changed its input")
        record["observed_measured_seconds"].append(timing)
        record["observed_ari"].append(quality)
        partitions.add(run["validated"]["canonical_partition_sha256"])
        logical_inputs.add(run["validated"]["logical_input_sha256"])
    if len(partitions) > 1 or len(logical_inputs) > 1:
        raise ValueError(f"Repeated fits have changed input/partition: {case['case_id']} / {method}")
    if measured:
        record["partition_stable_across_successful_repeats"] = True if len(measured) > 1 else None
        record["logical_input_sha256"] = next(iter(logical_inputs))
    for run in failures:
        evidence = {key: run[key] for key in (
            "repeat", "warmup", "status", "exit_code", "timeout_phase", "elapsed_in_phase_seconds",
            "limits", "fit_error", "validation_error", "sampled_memory"
        ) if key in run}
        record["terminal_evidence"].append(evidence)
    if not eligible:
        return record
    times = record["observed_measured_seconds"]
    record.update(median_seconds=statistics.median(times), min_seconds=min(times), max_seconds=max(times),
                  median_ari=statistics.median(record["observed_ari"]), ari_values=sorted(set(record["observed_ari"])),
                  timing_statistic="single fit" if len(times) == 1 else "median of measured fits",
                  runtime_variability_estimated=len(times) > 1)
    for destination, source in (
        ("cluster_counts", "clusters"), ("noise_counts", "noise_points"),
        ("aggregation_groups", "aggregation_groups"), ("adjacency_connection_buffer_bytes", "adjacency_bytes"),
        ("minpts_marker_counts", "minpts_marker_count"),
        ("minpts_marker_object_estimates_bytes", "minpts_marker_storage_bytes"),
    ):
        record[destination] = sorted({run["validated"][source] for run in measured
                                      if run["validated"].get(source) is not None})
    for metric, name in (("private_commit", "private"), ("working_set", "working_set")):
        key = "lifetime_peak_" + metric + "_bytes"
        peaks = [run["fit_done"]["after_fit"][key] for run in measured]
        record[f"median_os_peak_{name}_mib"] = statistics.median(peaks) / MIB
        record[f"{name}_peak_verified_in_fit_window"] = all(
            run["fit_done"]["after_fit"][key] > run["ready"]["before_fit"][key] for run in measured
        )
    baseline = [run["ready"]["before_fit"]["private_commit_bytes"] for run in measured]
    record["median_prefit_private_mib"] = statistics.median(baseline) / MIB
    if record["private_peak_verified_in_fit_window"]:
        record["median_extra_private_mib"] = statistics.median(
            run["fit_done"]["after_fit"]["lifetime_peak_private_commit_bytes"]
            - run["ready"]["before_fit"]["private_commit_bytes"] for run in measured
        ) / MIB
    return record


def measurement_note_lines(summary, study, metric):
    """State the actual repetition policy without implying single-run variation."""
    prefix = "PARTIAL RESULTS; " if not summary["complete_case_results"] else ""
    if summary.get("repetition_policy"):
        groups = {}
        for method in METHODS:
            counts = tuple(sorted({row["expected_repeats"] for row in summary["records"]
                                   if row["study"] == study and row["method"] == method}))
            groups.setdefault(counts, []).append(SHORT_NAMES[method])
        lines = []
        for counts, names in groups.items():
            statistic = " / ".join("1 fit" if count == 1 else f"median of {count} fits" for count in counts)
            lines.append("/".join(names) + ": " + statistic + ".")
        lines[0] = prefix + "Seed 101; " + lines[0]
    else:
        counts = sorted({protocol.get("measured_repeats", summary["plan"]["measured_repeats"])
                         for protocol in summary["recorded_protocols"]})
        repetitions = "/".join(str(value) for value in counts) or str(summary["plan"]["measured_repeats"])
        lines = [prefix + f"Seed 101; {repetitions} measured repetitions; medians."]
    if study == "radius":
        cases = [case for case in summary["plan"]["cases"] if case["study"] == "radius"]
        sizes = sorted({case["rows"] for case in cases})
        dimensions = sorted({case["features"] for case in cases})
        if len(sizes) != 1 or len(dimensions) != 1:
            raise ValueError("Radius figure requires one fixed dataset size and dimension")
        lines.append(f"Radius study: {sizes[0]:,} points; {dimensions[0]:,} dimensions; same dataset at every radius.")
    if study == "sample_count" and "sample_count_100000" in summary.get("excluded_case_ids", []):
        lines.append("The cancelled 100,000-point setting is excluded from this study.")
    if metric == "memory":
        lines.append("Process peaks include input and imports.")
    return lines


def figures(summary, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory
    import numpy as np

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
        "axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    destination.mkdir(parents=True, exist_ok=True)
    lookup = {(row["case_id"], row["method"]): row for row in summary["records"]}
    artifacts = []
    for study, title, xlabel, time_scale in STUDIES:
        cases = sorted((case for case in summary["plan"]["cases"] if case["study"] == study), key=lambda case: case["x"])
        xs = np.array([case["x"] for case in cases], dtype=float)
        rows_by_method = {method: [lookup[(case["case_id"], method)] for case in cases] for method in METHODS}
        unavailable_rows = [(method, rows) for method, rows in rows_by_method.items()
                            if any(not row["eligible_for_curves"] for row in rows)]
        used_statuses = sorted({status_mark(row["status"]) for _, rows in unavailable_rows
                               for row in rows if not row["eligible_for_curves"]})
        key_lines = textwrap.wrap("Status: " + "; ".join(f"{mark} = {label}" for mark, label in used_statuses),
                                  width=104) if used_statuses else []
        span = float(xs[-1] - xs[0])
        padding = .05 * span if span else max(1, abs(float(xs[0])) * .05)
        for metric, key, ylabel in (
            ("runtime", "median_seconds", "Time taken (in seconds)"),
            ("ari", "median_ari", "Adjusted Rand index"),
            ("memory", "median_os_peak_private_mib", "Peak process private memory (MiB)"),
        ):
            clean_sample_runtime = study == "sample_count" and metric == "runtime"
            clean_popcount_runtime = study == "popcount_separation" and metric == "runtime"
            clean_dimension_runtime = study == "dimension" and metric == "runtime"
            clean_runtime = clean_sample_runtime or clean_popcount_runtime or clean_dimension_runtime
            note_lines = measurement_note_lines(summary, study, metric)
            note_top = .018 + .023 * (len(note_lines) - 1)
            key_top = note_top + .030 + .023 * max(0, len(key_lines) - 1)
            status_bottom = key_top + .034
            main_bottom = (.18 + .028 * len(unavailable_rows) + .023 * max(0, len(key_lines) - 1)
                           + .023 * (len(note_lines) - 1))
            fig, ax = plt.subplots(figsize=(6.4, 4.0) if clean_runtime else (6.4, 4.8))
            fig.subplots_adjust(left=.14, right=.97, top=.9, bottom=.15 if clean_runtime else main_bottom)
            handles, finite_values = [], []
            for method, (label, color, marker, linestyle) in METHODS.items():
                if clean_runtime:
                    label = label.replace(" (packed)", "").replace(" (reference)", "")
                if method == "dbscan" and study == "popcount_separation":
                    linestyle = "-." if clean_popcount_runtime else "--"
                rows = rows_by_method[method]
                ys = np.array([row[key] if row["eligible_for_curves"] else np.nan for row in rows], dtype=float)
                finite_values.extend(ys[np.isfinite(ys)].tolist())
                ax.plot(xs, ys, color=color, linestyle=linestyle, marker=marker, linewidth=1.5, markersize=6)
                handles.append(Line2D([0], [0], label=label, color=color, linestyle=linestyle, marker=marker,
                                      linewidth=1.5, markersize=6))
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            # Keep failed endpoints visible; numerical curves cannot set limits.
            plot_padding = 4000 if clean_sample_runtime else (20 if clean_popcount_runtime else padding)
            ax.set_xlim(float(xs[0]) - plot_padding, float(xs[-1]) + plot_padding)
            ax.set_xticks(xs)
            ax.tick_params(axis="x", labelsize=10 if clean_popcount_runtime else (8 if study == "sample_count" else 9))
            if metric == "runtime":
                ax.set_yscale(time_scale)
                if not finite_values:
                    ax.set_ylim((.01, 10) if time_scale == "log" else (0, 1))
                elif time_scale == "linear":
                    ax.set_ylim(bottom=0)
            elif metric == "ari":
                ax.set_ylim(min(-.05, min(finite_values, default=0) - .05), 1.055)
            else:
                ax.set_ylim(bottom=0)
            location = "upper right" if study == "popcount_separation" else "upper left"
            legend_location = "upper left" if clean_runtime and not clean_sample_runtime else ("best" if finite_values else location)
            ax.legend(handles=handles, loc=legend_location, fontsize=8, frameon=True)
            if clean_popcount_runtime and finite_values:
                # Reserve enough log-axis headroom for the fixed upper-left legend.
                fig.canvas.draw()
                lower = ax.get_ylim()[0]
                box = ax.get_legend().get_window_extent(fig.canvas.get_renderer())
                fraction = (box.y0 - 6 * fig.dpi / 100 - ax.bbox.y0) / ax.bbox.height
                needed = lower * (max(finite_values) / lower) ** (1 / fraction)
                unit = 10 ** math.floor(math.log10(needed))
                ax.set_ylim(top=math.ceil(needed / unit) * unit)
            if not finite_values:
                ax.text(.5, .43, "No completed measurements", transform=ax.transAxes,
                        ha="center", va="center", fontsize=10, color="#777777")
                if metric == "memory":
                    ax.set_ylim(0, 1)
            if summary["is_pilot"]:
                ax.text(.5, .5, f"PILOT: {summary['pilot_rows']} rows per setting\nFull-data x labels show intended settings",
                        transform=ax.transAxes, ha="center", va="center", fontsize=11, color="#777777", alpha=.8)
            if not clean_runtime:
                for index, line in enumerate(note_lines):
                    fig.text(.5, note_top - .023 * index, line, ha="center", fontsize=7)
                if unavailable_rows:
                    # These are figure text rows, not a second plot or y observations.
                    # x follows the main data axis, while y is in figure coordinates.
                    status_transform = blended_transform_factory(ax.transData, fig.transFigure)
                    for index, (method, rows) in enumerate(unavailable_rows):
                        y = status_bottom + .028 * (len(unavailable_rows) - 1 - index)
                        color = METHODS[method][1]
                        fig.text(.133, y, SHORT_NAMES[method], ha="right", va="center", fontsize=7, color=color)
                        for x, row in zip(xs, rows):
                            if not row["eligible_for_curves"]:
                                mark, _ = status_mark(row["status"])
                                fig.text(x, y, mark, transform=status_transform, ha="center", va="center",
                                         fontsize=7.5, color=color, fontweight="bold")
                    for index, line in enumerate(key_lines):
                        fig.text(.5, key_top - .023 * index, line, ha="center", fontsize=7)
                else:
                    fig.text(.5, note_top + .037, "Fixed nominal settings; different methods use different clustering criteria.", ha="center", fontsize=7)
            paths = []
            for extension in ("png", "pdf", "svg"):
                path = destination / f"{study}_{metric}.{extension}"
                fig.savefig(path, dpi=300, facecolor="white")
                paths.append(str(path.resolve()))
            plt.close(fig)
            artifacts.append({"study": study, "metric": metric, "paths": paths,
                              "yscale": time_scale if metric == "runtime" else "linear",
                              "measurement_notes": note_lines,
                              "x_limits": [float(xs[0]) - plot_padding, float(xs[-1]) + plot_padding],
                              "status_rows": [{"method": method, "points": [
                                  {"x": row["x"], "status": row["status"], "mark": status_mark(row["status"])[0]}
                                  for row in rows if not row["eligible_for_curves"]]}
                                  for method, rows in unavailable_rows]})
    return artifacts


def collect_fresh(root, run_directory):
    plan = read_configuration(CONFIG)
    policy = read_repetitions(CONFIG)
    records, missing, protocols = [], [], []
    pilot_values = set()
    raws = {}
    for case in plan["cases"]:
        path = run_directory/"results"/(case["case_id"]+".json")
        raw = json.loads(path.read_text()) if path.is_file() else {}
        raws[case["case_id"]] = raw
        if raw:
            assert raw["case"] == case
            pilot_values.add(raw.get("pilot_rows"))
            protocols.append(raw["protocol"])
        if not raw.get("completed"):
            missing.append(case["case_id"])
    assert len(pilot_values) <= 1
    pilot_rows = next(iter(pilot_values)) if pilot_values else None
    for case in plan["cases"]:
        raw = raws[case["case_id"]]
        for method in METHODS:
            ids = [0] if pilot_rows else policy["schedules"][case["case_id"]][method]
            records.append(method_record(case, raw, method, len([i for i in ids if i >= 0]), ids))
    return {"plan": plan, "records": records, "complete_case_results": not missing,
            "missing_or_unfinished_cases": missing, "repetition_policy": policy,
            "excluded_case_ids": ["sample_count_100000"], "recorded_protocols": protocols,
            "is_pilot": bool(pilot_rows), "pilot_rows": pilot_rows,
            "data_origin": "Fresh generation and measurements from this notebook; no historical result values loaded."}


def verify_fresh(root, run_directory):
    import numpy as np
    from sklearn.metrics import adjusted_rand_score
    import hashlib
    plan = read_configuration(CONFIG)
    policy = read_repetitions(CONFIG)
    from reproduce import method_finished
    successful, failures, missing = [], [], []
    for case in plan["cases"]:
        path = run_directory/"results"/(case["case_id"]+".json")
        if not path.is_file():
            missing.append(case["case_id"])
            continue
        raw = json.loads(path.read_text())
        assert raw["case"] == case
        manifest = json.loads((run_directory/"data"/case["data_id"]/"manifest.json").read_text())
        truth = np.load(run_directory/"data"/case["data_id"]/manifest["files"]["truth"]["file"])
        for method in METHODS:
            runs = [r for r in raw["runs"] if r["variant"] == method]
            ids = [0] if raw.get("pilot_rows") else policy["schedules"][case["case_id"]][method]
            if not method_finished(runs, ids):
                missing.append(case["case_id"]+"/"+method)
            for run in runs:
                assert run["protocol"] == raw["protocol"]
                if run["status"] != "ok":
                    failures.append({"case_id": case["case_id"], "method": method, "status": run["status"]})
                    continue
                checked = run["validated"]
                relative = checked["output_relative_assignments_file"]
                label_path = (run_directory/relative).resolve()
                assert label_path.is_relative_to(run_directory.resolve())
                labels = np.load(label_path, allow_pickle=False)
                expected_n = raw.get("pilot_rows") or case["rows"]
                assert labels.shape == (expected_n,)
                score = adjusted_rand_score(truth[:len(labels)], labels)
                assert math.isclose(score, checked["ari"], rel_tol=1e-12, abs_tol=1e-12)
                _, counts = np.unique(labels[labels >= 0], return_counts=True)
                assert len(counts) == checked["clusters"]
                assert int(np.count_nonzero(labels < 0)) == checked["noise_points"]
                mapping = {}
                canonical = np.asarray([mapping.setdefault(int(v), len(mapping)) for v in labels], dtype=np.int64)
                assert hashlib.sha256(memoryview(canonical).cast("B")).hexdigest() == checked["canonical_partition_sha256"]
                if not raw.get("pilot_rows"):
                    assert checked["input_unchanged"] is True
                    assert checked["logical_input_sha256"] == manifest["files"]["int32"]["sha256"]
                successful.append({"case_id": case["case_id"], "method": method, "repeat": run["repeat"], "ari": score})
    return {"successful_assignments_checked": len(successful), "successful": successful,
            "terminal_failures": failures, "missing_or_incomplete": missing,
            "ari_reference": "sklearn.metrics.adjusted_rand_score; freshly generated truth"}
