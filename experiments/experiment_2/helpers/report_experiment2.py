"""Portable Experiment 2 report for twenty original-feature real datasets.

Usage: python report_experiment2.py --root RUN_ROOT [--no-plots]
Inputs are RUN_ROOT/data/*.manifest.json and RUN_ROOT/results/*.json[l].
Outputs are RUN_ROOT/RESULTS.md, SUMMARY.json, tables/ and plots/.
This module imports no project modules and never fits or prepares data.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE

ROSTER = ("iris", "banknote", "wine", "breast_cancer", "digits", "seeds", "glass", "sonar",
          "vehicle", "segmentation", "breast_tissue", "leaf", "libras", "vertebral", "wireless",
          "parkinsons", "coimbra", "raisin", "haberman", "musk")
ANCHORS = ROSTER[:8]
PANEL10 = (*ANCHORS, "breast_tissue", "haberman")
METHODS = ("euclidean", "manhattan")
METHOD_NAMES = {"euclidean": "Euclidean CLASSIX", "manhattan": "Manhattan CLASSIX"}
COLORS = {"euclidean": "#2166ac", "manhattan": "#b35806"}
GOOD = {"ok", "success", "successful", "complete", "completed"}
QUALITY_HEADERS = ("Dataset", "Points × features", "Euclidean ARI", "Manhattan ARI", "Euclidean AMI", "Manhattan AMI")


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def successful(row):
    return bool(row) and str(row.get("status", "ok")).lower() in GOOD


def score(row, name):
    value = row.get(name) if successful(row) else None
    return value if finite(value) else None


def display_name(name):
    return str(name).replace("_", " ").capitalize()


def number(value, digits=4):
    return f"{value:.{digits}f}" if finite(value) else "—"


def shape_text(shape):
    return " × ".join(f"{value:,}" if isinstance(value, int) else str(value) for value in shape) if shape else "unknown"


def clean_json(value):
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(item) for item in value]
    return None if isinstance(value, float) and not math.isfinite(value) else value


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path, root, provenance, *, jsonl=False, default=None):
    if not path.exists():
        return [] if default is None else default
    content = path.read_bytes()
    provenance.append({"file": path.relative_to(root).as_posix(), "sha256": hashlib.sha256(content).hexdigest()})
    text = content.decode("utf-8-sig")
    return [json.loads(line) for line in text.splitlines() if line.strip()] if jsonl else json.loads(text)


def latest(rows, fields):
    found = {}
    for row in rows:
        found[tuple(row.get(field) for field in fields)] = row
    return list(found.values())


def table(headers, rows):
    result = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    result.extend("| " + " | ".join(str(item).replace("|", "\\|").replace("\n", " ") for item in row) + " |" for row in rows)
    return "\n".join(result)


def outcome(euclidean, manhattan):
    e = [score(euclidean, name) for name in ("ari", "ami")]
    m = [score(manhattan, name) for name in ("ari", "ami")]
    if not all(finite(value) for value in e + m):
        return "unavailable"
    if all(a > b for a, b in zip(m, e)):
        return "M"
    if all(a < b for a, b in zip(m, e)):
        return "E"
    return "mixed/tied"


def current_classification(dataset, selected, permutations):
    original = outcome(selected.get((dataset, "euclidean"), {}), selected.get((dataset, "manhattan"), {}))
    repeated = [outcome(permutations.get((dataset, "euclidean", repeat), {}),
                        permutations.get((dataset, "manhattan", repeat), {})) for repeat in range(5)]
    if original == "unavailable" or "unavailable" in repeated:
        return "incomplete"
    if original == "M" and repeated == ["M"] * 5:
        return "M consistent in 5/5"
    if original == "E" and repeated == ["E"] * 5:
        return "E consistent in 5/5"
    return "mixed/tied/order-sensitive"


def apply_original_selection_rule(classifications):
    """Audit the old rule for this run; never replace the approved fixed panel."""
    chosen = list(ANCHORS)
    new = sorted(set(ROSTER) - set(ANCHORS))
    for category in ("M consistent in 5/5", "E consistent in 5/5"):
        eligible = [name for name in new if classifications[name] == category]
        if eligible:
            chosen.append(eligible[0])
    for mixed_only in (True, False):
        for name in new:
            if len(chosen) == 10:
                break
            if name in chosen or classifications[name] == "incomplete":
                continue
            if mixed_only and classifications[name] != "mixed/tied/order-sensitive":
                continue
            chosen.append(name)
    return chosen


def source_url(manifest):
    source = manifest.get("source", {})
    for item in (source, source.get("original_source", {})):
        for key in ("original_uci_url", "source_url", "repository_url"):
            value = item.get(key)
            if isinstance(value, str) and value.startswith(("https://", "http://")):
                return value
        if isinstance(item.get("uci_id"), int):
            return f"https://archive.ics.uci.edu/dataset/{item['uci_id']}"
    return None


def boundaries(chosen, grid):
    if not successful(chosen):
        return []
    radii = sorted({row["radius"] for row in grid if finite(row.get("radius"))})
    minima = sorted({row["min_pts"] for row in grid if finite(row.get("min_pts"))})
    flags = []
    if len(radii) > 1:
        if chosen.get("radius") == radii[0]:
            flags.append("lowest tested radius")
        if chosen.get("radius") == radii[-1]:
            flags.append("highest tested radius")
    if len(minima) > 1 and chosen.get("min_pts") == minima[-1]:
        flags.append("highest tested minPts")
    return flags


def collect(root):
    root = Path(root).resolve()
    provenance = []
    results = root / "results"
    selected_raw = read(results / "selected.json", root, provenance)
    grid_raw = read(results / "grid.jsonl", root, provenance, jsonl=True)
    permutations_raw = read(results / "permutation.jsonl", root, provenance, jsonl=True)
    failure_log = read(results / "failures.jsonl", root, provenance, jsonl=True)
    if any(not isinstance(rows, list) for rows in (selected_raw, grid_raw, permutations_raw, failure_log)):
        raise ValueError("Expected selected.json as a list and JSONL files as record streams")
    selected_rows = latest(selected_raw, ("dataset", "method"))
    grid = latest(grid_raw, ("dataset", "method", "radius_index", "min_pts"))
    permutations_rows = latest(permutations_raw, ("dataset", "method", "replicate"))
    selected = {(row.get("dataset"), row.get("method")): row for row in selected_rows}
    permutations = {(row.get("dataset"), row.get("method"), row.get("replicate")): row for row in permutations_rows}
    for path in sorted((results / "plans").glob("*.json")) + sorted((results / "grid_specs").glob("*.json")):
        read(path, root, provenance, default={})
    for name in ("run_manifest.json", "imported_provenance.json", "grid_specs.json"):
        read(results / name, root, provenance, default={})
    protocol = read(root / "protocol.json", root, provenance, default={})
    if protocol.get("datasets") is not None and tuple(protocol["datasets"]) != ROSTER:
        raise ValueError("The supplied dataset roster differs from the approved twenty-dataset screen")
    records, classifications = [], {}
    for dataset in ROSTER:
        manifest = read(root / "data" / f"{dataset}.manifest.json", root, provenance, default={})
        shape = manifest.get("arrays", {}).get("X", manifest.get("arrays", {}).get("X_raw", {})).get("shape", manifest.get("shape"))
        methods = {}
        for method in METHODS:
            chosen = selected.get((dataset, method), {})
            fits = [row for row in grid if row.get("dataset") == dataset and row.get("method") == method]
            methods[method] = {"selected": chosen, "successful_grid_settings": sum(successful(row) for row in fits),
                               "grid_settings": len(fits), "boundary_flags": boundaries(chosen, fits)}
        classification = current_classification(dataset, selected, permutations)
        classifications[dataset] = classification
        records.append({"dataset": dataset, "shape": shape, "source_url": source_url(manifest),
                        "methods": methods, "ordering_classification": classification})
    failures = [{"stage": stage, "record": row} for stage, saved in
                (("grid", grid_raw), ("selected", selected_raw), ("permutation", permutations_raw))
                for row in saved if not successful(row)]
    failures.extend({"stage": "failure_log", "record": row} for row in failure_log)
    historical_rule = apply_original_selection_rule(classifications)
    audit_complete = all(value != "incomplete" for value in classifications.values())
    counts = {name: {"saved_records": len(saved), "current_records": len(active),
                     "successful_records": sum(successful(row) for row in active),
                     "origins": dict(Counter(row.get("origin", "not recorded") for row in active))}
              for name, saved, active in (("grid", grid_raw, grid), ("selected", selected_raw, selected_rows),
                                          ("permutation", permutations_raw, permutations_rows))}
    return {"schema_version": 1, "experiment": "2", "roster": list(ROSTER), "approved_panel": list(PANEL10),
            "records": records, "counts": counts, "all_selected_scores_available": all(
                finite(score(selected.get((name, method), {}), metric_name))
                for name in ROSTER for method in METHODS for metric_name in ("ari", "ami")),
            "panel_selection_audit": {"complete": audit_complete, "rule_selected_ids_this_run": historical_rule,
                "matches_approved_panel": historical_rule == list(PANEL10) if audit_complete else None,
                "policy": "Approved panel stays fixed. Changed outcomes on a fresh run are disclosed, never used to silently replace cases."},
            "ordering_classification_counts": dict(Counter(classifications.values())),
            "all_saved_failures": failures, "sources": provenance}


def values_for_table(record, *, markdown=False):
    values, winners = [], []
    for metric_name in ("ari", "ami"):
        pair = [score(record["methods"][method]["selected"], metric_name) for method in METHODS]
        for index, value in enumerate(pair):
            winning = finite(value) and finite(pair[1-index]) and value > pair[1-index]
            text = number(value)
            values.append(f"**{text}**" if winning and markdown else text)
            winners.append(winning)
    return values, winners


def quality_table(records):
    return table(QUALITY_HEADERS,
        [[display_name(row["dataset"]), shape_text(row["shape"]),
          *values_for_table(row, markdown=True)[0]] for row in records])


def make_tables(report):
    by_id = {row["dataset"]: row for row in report["records"]}
    parameters = []
    for row in report["records"]:
        for method in METHODS:
            entry = row["methods"][method]
            chosen = entry["selected"]
            times = [value for value in chosen.get("timing_seconds", []) if finite(value)]
            parameters.append([display_name(row["dataset"]), METHOD_NAMES[method], chosen.get("status", "ok" if chosen else "missing"),
                chosen.get("origin", "not recorded"), f"{chosen['radius']:.10g}" if finite(chosen.get("radius")) else "—",
                chosen.get("radius_index", "—"), chosen.get("min_pts", "—"), chosen.get("merge_scale", "—"),
                entry["successful_grid_settings"], chosen.get("n_clusters", "—"), chosen.get("noise_points", chosen.get("n_noise", "—")),
                f"{chosen['fit_seconds_median']:.6g}" if finite(chosen.get("fit_seconds_median")) else "—", len(times),
                ", ".join(entry["boundary_flags"]) or "none recorded"])
    return {"main_10": quality_table([by_id[name] for name in PANEL10]), "all_20": quality_table(report["records"]),
            "parameters": table(["Dataset", "Method", "Status", "Origin", "Radius", "Radius index", "minPts", "Merge scale", "Grid successes", "Clusters", "Noise", "Median fit s", "Measured fits", "Search boundary"], parameters),
            "datasets": table(["Dataset", "Points × features", "Public source"],
                [[display_name(row["dataset"]), shape_text(row["shape"]),
                  f"[Source]({row['source_url']})" if row["source_url"] else "See manifest"] for row in report["records"]])}


def markdown(report, tables, artifacts):
    audit = report["panel_selection_audit"]
    paragraphs = ["# Experiment 2: original-feature real-data comparison", "",
        "This focused package preserves the approved ten-case panel from the earlier Experiment 2B twenty-dataset screen. It uses standardized original numeric features and includes no UMAP experiment. The complete twenty-dataset appendix remains visible.", "",
        "## Approved ten-dataset table", "",
        "Values are the original-order ARI-selected scores, rounded to four decimals. **Bold** marks the higher full-precision score within each dataset/metric pair; an apparent rounded tie can therefore still have a bold entry. Equal full-precision scores are not bolded. These are not means across permutations or datasets. AMI is evaluated on the ARI-selected assignments and is not independently optimized.", "",
        tables["main_10"], "",
        "The panel is outcome-stratified and illustrative: eight earlier datasets were mandatory anchors; the original frozen rule added the alphabetically first new case favoring Manhattan on both metrics in the original order and all five tested permutations, plus the corresponding Euclidean case if available. No new Euclidean-consistent case qualified in that recorded screen, so the declared alphabetical mixed-case fallback selected Haberman. Breast Tissue was the Manhattan addition. This panel is not an unbiased sample for estimating method win rates.", "",
        "The approved IDs stay fixed on a fresh full run. A diagnostic reapplication of the old rule is recorded in SUMMARY.json; changed outcomes do not silently replace the approved cases.", "",
        f"Selection-rule audit complete for this run: **{audit['complete']}**. Matches the approved IDs: **{audit['matches_approved_panel']}**.", "",
        "## Complete twenty-dataset appendix", "", tables["all_20"], "",
        "## Selected parameters and measured fit timing", "",
        "Each method uses 32 geometry-calibrated radius positions and minPts {0,5,10,20,40,80}, with the recorded one-time matched boundary extension. Selection maximizes full-precision ARI, breaking ties by smaller minPts then smaller radius index. The actual selected radius, minPts and merge scale below are the values used for the displayed scores; the complete saved rows and grid specifications retain full precision.", "",
        tables["parameters"], "",
        "Runtime summaries are medians of the saved measured fits, not arithmetic means. Five measurements follow a warmup where the completed protocol specifies them. Internal method preprocessing is timed; data loading, shared StandardScaler transformation, scoring and tuning are excluded. Imported timings, if present, remain earlier measurements rather than new repetitions.", "",
        "An outer-radius or maximum-minPts winner is a best tested result, not an established global optimum. minPts=0 is a valid setting. The implementations differ in normalization, sorting and compiled/NumPy execution; quality and runtime differences are comparisons of complete methods, not isolated causal effects of the metric.", "",
        "## Data sources and interpretation limits", "", tables["datasets"], "",
        "All scores are label-informed and in-sample; there is no unseen real-data test set. Source labels need not uniquely define natural clusters. The manifests disclose identifier exclusions, repeated observations, preprocessing and original feature types. No statement of statistical significance or universal superiority follows from this exploratory screen. Five row permutations, when available, assess those tested orderings only; they are not independent datasets.", "",
        "## Completion and preservation", "",
        table(["Stage", "Saved", "Current", "Successful", "Origins"],
              [[name, value["saved_records"], value["current_records"], value["successful_records"], json.dumps(value["origins"], sort_keys=True)] for name, value in report["counts"].items()]), "",
        f"All selected scores available: **{report['all_selected_scores_available']}**. Saved failures: **{len(report['all_saved_failures'])}**. Missing/non-successful measurements are never displayed as successful scores.", ""]
    if report["all_saved_failures"]:
        paragraphs.append(table(["Stage", "Dataset", "Method", "Status", "Reason"],
            [[item["stage"], item["record"].get("dataset", item["record"].get("last_progress", {}).get("dataset", "—")),
              item["record"].get("method", item["record"].get("last_progress", {}).get("method", "—")), item["record"].get("status", "failed"),
              item["record"].get("error", item["record"].get("message", item["record"].get("reason", "See saved record")))] for item in report["all_saved_failures"]]))
    if artifacts:
        paragraphs.extend(["", "## Standalone figures", ""])
        for path in artifacts:
            if path.endswith(".png"):
                stem = path[:-4]
                paragraphs.append(f"- [{Path(stem).name}]({path}) · [PDF]({stem}.pdf) · [SVG]({stem}.svg)")
    paragraphs.extend(["", "SUMMARY.json records the selected rows, the fixed-panel audit, all saved failures and hashes of every input file read. This reporter has no dependency on a prior workspace or experiment and never changes input records.", ""])
    return "\n".join(paragraphs)


def render(report, root):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180, "svg.fonttype": "none"})
    destination = root / "plots"
    destination.mkdir(parents=True, exist_ok=True)
    artifacts = []
    def save(fig, name, footer):
        fig.text(.5, .015, footer, ha="center", va="bottom", fontsize=8, color="#444444")
        fig.tight_layout(rect=(0, .09, 1, .95))
        for extension in ("png", "pdf", "svg"):
            path = destination / f"{name}.{extension}"
            fig.savefig(path, bbox_inches="tight")
            artifacts.append(path.relative_to(root).as_posix())
        plt.close(fig)
    by_id = {row["dataset"]: row for row in report["records"]}
    panel = [by_id[name] for name in PANEL10]
    for name, records in (("approved10", panel), ("all20", report["records"])):
        for metric_name in ("ari", "ami"):
            groups = [records[index:index+10] for index in range(0, len(records), 10)]
            fig, axes = plt.subplots(len(groups), 1, figsize=(12, 4.7 if len(groups) == 1 else 7.8), squeeze=False, sharey=True)
            for ax, group in zip(axes[:, 0], groups):
                for offset, method in zip((-.19, .19), METHODS):
                    values = [score(row["methods"][method]["selected"], metric_name) for row in group]
                    positions = np.arange(len(group)) + offset
                    ax.bar(positions, [value if finite(value) else np.nan for value in values], width=.35, color=COLORS[method], label=METHOD_NAMES[method])
                    for position, value in zip(positions, values):
                        if finite(value):
                            ax.annotate(f"{value:.3f}", (position, value), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7)
                        else:
                            ax.text(position, .02, "unavailable", transform=ax.get_xaxis_transform(), ha="center", rotation=90, fontsize=7)
                ax.set_xticks(range(len(group)), [display_name(row["dataset"]) for row in group], rotation=25, ha="right")
                ax.set_ylabel(f"{metric_name.upper()} (higher is better)")
                ax.grid(axis="y", alpha=.2)
                ax.margins(y=.22)
            axes[0, 0].legend(loc="best")
            heading = "Approved ten-dataset panel" if name == "approved10" else "Complete twenty-dataset appendix"
            fig.suptitle(f"Experiment 2: {heading} — {metric_name.upper()}", y=.995)
            save(fig, f"{name}_{metric_name}", "Standardized original features; original-order ARI-selected settings. AMI is evaluated at those same settings.\nThe panel is illustrative; the full screen is retained. No unseen real-data test or metric-only causal claim.")
    fig, ax = plt.subplots(figsize=(11.7, 4.8))
    ax.set_axis_off()
    cell_values = [[display_name(row["dataset"]), shape_text(row["shape"]),
                    *values_for_table(row)[0]] for row in panel]
    drawn = ax.table(cellText=cell_values,
                     colLabels=QUALITY_HEADERS,
                     loc="center", cellLoc="center", colWidths=[.18, .14, .17, .17, .17, .17])
    drawn.auto_set_font_size(False)
    drawn.set_fontsize(10)
    drawn.scale(1, 1.6)
    for (r, c), box in drawn.get_celld().items():
        box.set_linewidth(.4)
        box.set_edgecolor("#cccccc")
        if r == 0:
            box.set_facecolor("#eeeeee")
            box.set_text_props(weight="bold")
        elif c >= 2 and values_for_table(panel[r-1])[1][c-2]:
            box.set_text_props(weight="bold")
    fig.suptitle("Experiment 2: approved original-feature comparison", y=.97)
    save(fig, "approved10_table", "Four decimals; bold marks the higher full-precision score in each dataset/metric pair.\nScores use ARI-selected parameters and are not means across permutations. Full twenty-dataset appendix retained.")
    return artifacts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=RUN)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    report = collect(root)
    tables = make_tables(report)
    artifacts = [] if args.no_plots else render(report, root)
    report["artifacts"] = artifacts
    report["report_source_sha256"] = digest(Path(__file__))
    for item in report["sources"]:
        if digest(root / item["file"]) != item["sha256"]:
            raise RuntimeError("An input changed while reporting; run after fits finish")
    (root / "tables").mkdir(parents=True, exist_ok=True)
    for name, content in tables.items():
        (root / "tables" / f"{name}.md").write_text(content + "\n", encoding="utf-8")
    (root / "SUMMARY.json").write_text(json.dumps(clean_json(report), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (root / "RESULTS.md").write_text(markdown(report, tables, artifacts), encoding="utf-8")
    print(json.dumps({"root": str(root), "approved_panel": list(PANEL10),
                      "all_selected_scores_available": report["all_selected_scores_available"],
                      "plot_files": len(artifacts), "saved_failures": len(report["all_saved_failures"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
