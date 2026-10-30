"""Render saved AMI/silhouette scores only; never fit, rescore or load datasets.

python code/render_scores.py --check
python code/render_scores.py

Reads results/summary.json. Writes plots/{synthetic,landrum},
results/SCORES.md and results/rendering_provenance.json under --root.
Scientific imports are deferred until rendering is explicitly requested.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import textwrap

ROOT = work_root(SOURCE_ROOT)
METHODS = {
    "classix_packed": ("CLASSIX_T (packed)", "red", "x", "--"),
    "dbscan": ("DBSCAN", "blue", "o", "-."),
    "bitbirch": ("BitBIRCH", "green", "v", ":"),
    "bblean": ("BitBIRCH-Lean (reference)", "#9467bd", "s", "--"),
    "butina": ("Taylor-Butina", "#e68613", "D", "-"),
}
SHORT_NAMES = {"classix_packed": "CLASSIX_T", "dbscan": "DBSCAN", "bitbirch": "BitBIRCH",
               "bblean": "Lean", "butina": "Butina"}
SYNTHETIC_STUDIES = (
    ("sample_count", "Scaling with number of data points", "Number of points"),
    ("popcount_separation", "Scaling with separability of data", "Max pop of seed vector"),
    ("cluster_count", "Scaling with number of clusters", "Number of clusters"),
    ("dimension", "Scaling with dimensionality", "Number of dimensions"),
    ("radius", "Sensitivity to Tanimoto radius", "Tanimoto radius r"),
)
LANDRUM_METHODS = ("classix_packed", "bitbirch", "bblean")
METRICS = (
    ("ami", "ami_all", "Adjusted mutual information (all rows)"),
    ("silhouette", "silhouette_mean", "Mean Tanimoto silhouette"),
    ("coverage", "coverage", "Retained points / all input points"),
)
STATUS_NAMES = {"M": "memory-limited fit", "U": "undefined score", "N": "score unavailable",
                "P": "pending", "E": "failed fit"}
FULL_SYNTHETIC_CASES = frozenset(
    [f"sample_count_{n}" for n in range(10000, 90001, 10000)]
    + [f"popcount_separation_{p}" for p in range(50, 501, 50)]
    + [f"cluster_count_{k}" for k in (10, 20, 25, 40, 50, 80, 100)]
    + [f"dimension_{d}" for d in range(500, 5001, 500)]
    + [f"radius_0p{r}" for r in range(2, 8)])
FULL_LANDRUM_CASES = frozenset(f"morgan_min10_landrum_s{i:03d}" for i in range(100))


def summary_contract(summary):
    """Validate the intended experiment, independently of which records exist."""
    scope = summary.get("scope")
    if scope is None:
        # Compatibility for the original combined scratch score report only.
        if summary.get("is_pilot") or "expected_cases" in summary or "expected_methods" in summary:
            raise ValueError("An explicit experiment contract requires scope")
        return {"scope": "both", "is_pilot": False,
                "cases": {"synthetic": FULL_SYNTHETIC_CASES, "landrum": FULL_LANDRUM_CASES},
                "methods": {"synthetic": set(METHODS), "landrum": set(LANDRUM_METHODS)}}
    if scope not in ("synthetic", "landrum") or not isinstance(summary.get("is_pilot"), bool):
        raise ValueError("Scoped summary requires synthetic/landrum scope and boolean is_pilot")
    cases, methods = summary.get("expected_cases"), summary.get("expected_methods")
    for name, values in (("expected_cases", cases), ("expected_methods", methods)):
        if (not isinstance(values, list) or not values or
                any(not isinstance(item, str) or not item for item in values) or len(values) != len(set(values))):
            raise ValueError(name + " must be a nonempty list of unique identifiers")
    full_cases = FULL_SYNTHETIC_CASES if scope == "synthetic" else FULL_LANDRUM_CASES
    full_methods = set(METHODS) if scope == "synthetic" else set(LANDRUM_METHODS)
    if set(methods) != full_methods:
        raise ValueError("Expected methods do not match the complete method roster")
    if summary["is_pilot"]:
        if not set(cases).issubset(full_cases):
            raise ValueError("Pilot cases are not part of the intended full experiment")
        pilot_rows = summary.get("pilot_rows")
        if pilot_rows is not None and (not isinstance(pilot_rows, int) or isinstance(pilot_rows, bool) or pilot_rows < 2):
            raise ValueError("Invalid pilot_rows")
    elif set(cases) != full_cases:
        raise ValueError("A full report must declare every original case; reduced grids require is_pilot=true")
    return {"scope": scope, "is_pilot": summary["is_pilot"],
            "cases": {scope: frozenset(cases)}, "methods": {scope: set(methods)}}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def save(value, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def value_and_status(row, metric):
    """A missing or undefined score never becomes a numeric zero."""
    if row is None:
        return None, "P", "missing setting/method record"
    status = row.get("status")
    if status != "ok":
        mark = "M" if status == "memory_limit" else "P" if status in ("pending", "running") else "E"
        return None, mark, str(status)
    key = dict((name, key) for name, key, _ in METRICS)[metric]
    value = row.get(key)
    if metric == "silhouette":
        quality_status = row.get("silhouette_status")
        if quality_status != "ok":
            mark = "U" if quality_status == "undefined" else "P" if quality_status in ("pending", "running") else "N"
            return None, mark, row.get("silhouette_reason") or str(quality_status)
    if value is None:
        return None, "N", row.get(metric + "_reason") or (key + " unavailable")
    if not finite(value):
        raise ValueError(f"Nonfinite {key}: {row.get('case_id')} / {row.get('method')}")
    return float(value), None, None


def validate_records(summary, allow_partial=False):
    contract = summary_contract(summary)
    records = summary.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("summary.json must contain a nonempty records list")
    seen, case_methods, case_x = set(), {}, {}
    for row in records:
        scope, method = row.get("scope"), row.get("method")
        if scope not in ("synthetic", "landrum") or method not in METHODS:
            raise ValueError("Unknown scope or method")
        if (scope not in contract["cases"] or row.get("case_id") not in contract["cases"][scope]
                or method not in contract["methods"][scope]):
            raise ValueError("Record lies outside the declared experiment contract")
        if scope == "landrum" and method not in LANDRUM_METHODS:
            raise ValueError("Unexpected Landrum method")
        if scope == "synthetic" and row.get("study") not in {item[0] for item in SYNTHETIC_STUDIES}:
            raise ValueError("Unknown synthetic study")
        if not finite(row.get("x")) or not row.get("case_id"):
            raise ValueError("Missing finite x or case ID")
        identity = (scope, row["case_id"], method)
        if identity in seen:
            raise ValueError("Duplicate setting/method record: " + str(identity))
        seen.add(identity)
        case_key = (scope, row["case_id"])
        definition = (row.get("study"), row["x"])
        if case_key in case_x and case_x[case_key] != definition:
            raise ValueError("Inconsistent study/x for one setting")
        case_x[case_key] = definition
        case_methods.setdefault(case_key, set()).add(method)
        for metric, _, _ in METRICS:
            value, _, _ = value_and_status(row, metric)
            if metric == "coverage" and value is not None and not 0 <= value <= 1:
                raise ValueError("Coverage is outside [0,1]")
            if metric == "silhouette" and value is not None and not -1.0000001 <= value <= 1.0000001:
                raise ValueError("Silhouette is outside [-1,1]")
        if row.get("ami_retained") is not None and not finite(row["ami_retained"]):
            raise ValueError("Nonfinite retained-population AMI")
        if row.get("silhouette_status") == "ok" and row.get("status") == "ok":
            if not isinstance(row.get("exact_silhouette"), bool):
                raise ValueError("Successful silhouette must state whether it is exact")
            if not row["exact_silhouette"] and (not finite(row.get("silhouette_se")) or row["silhouette_se"] < 0):
                raise ValueError("Sampled silhouette lacks a valid standard error")
            if not isinstance(row.get("n_queries"), int) or row["n_queries"] <= 0:
                raise ValueError("Successful silhouette lacks query count")
    expected = {(scope, case, method) for scope, cases in contract["cases"].items()
                for case in cases for method in contract["methods"][scope]}
    complete = seen == expected and summary.get("complete_case_results", True) is True
    if not complete and not allow_partial:
        raise ValueError("Missing expected case/method records or unfinished case results; use --allow-partial only for labelled previews")
    return records, complete


def groups(records):
    for study, title, xlabel in SYNTHETIC_STUDIES:
        rows = [row for row in records if row["scope"] == "synthetic" and row["study"] == study]
        if rows:
            yield "synthetic", study, title, xlabel, tuple(METHODS), rows
    rows = [row for row in records if row["scope"] == "landrum"]
    if rows:
        yield "landrum", "similarity", "Landrum / Morgan 1024 / min10 study", "Tanimoto similarity threshold s", LANDRUM_METHODS, rows


def measurement_notes(scope, study, metric, rows, complete):
    notes = [] if complete else ["PARTIAL RESULTS - only saved available records are shown."]
    if metric == "ami":
        notes += ["AMI uses all input rows and saved assignments; retained-population AMI is separate."]
        if scope == "landrum":
            notes += ["Reference: 35 assay-target categories; repeated structures can have different targets."]
    elif scope == "landrum":
        notes += ["CLASSIX: all nonnoise final clusters after minPts=10; BitBIRCH/Lean: size >10."]
    else:
        notes += ["Quality population excludes algorithm noise; coverage reports the retained fraction."]
    if metric == "silhouette":
        valid = [row for row in rows if value_and_status(row, metric)[0] is not None]
        sampled = [row for row in valid if not row["exact_silhouette"]]
        exact = [row for row in valid if row["exact_silhouette"]]
        if sampled:
            q = sorted({row["n_queries"] for row in sampled})
            qtext = str(q[0]) if len(q) == 1 else f"{q[0]}-{q[-1]}"
            notes += [f"Estimated from {qtext} query rows; shading is pointwise +/-1.96 Monte Carlo SE.",
                      "Sampling uncertainty only; no clustering-repeat uncertainty is estimated."]
        if exact:
            notes += ["Some points use all retained rows (exact); no sampling band is drawn for those points."] if sampled else ["All available silhouettes use every retained row (exact; no sampling band)."]
        notes += ["Distances are evaluated against retained reference rows; singleton silhouettes are zero."]
    if scope == "synthetic":
        notes += ["Seed 101; saved partitions checked across repeats; radius study uses 10,000 points."] if study == "radius" else ["Seed 101; one score per setting after checking saved partitions across repeats."]
    else:
        notes += ["CLASSIX radius=1-s, mergeScale=1.35; BitBIRCH/Lean diameter threshold=s."]
    return [line for note in notes for line in textwrap.wrap(note, width=106)]


def finite_marker_indices(values):
    available = [i for i, value in enumerate(values) if math.isfinite(value)]
    if len(values) <= 20:
        return available
    step = max(1, len(values) // 12)
    candidates = set(available[::step])
    for i in available:
        if i == 0 or i == len(values) - 1 or not math.isfinite(values[i-1]) or not math.isfinite(values[i+1]):
            candidates.add(i)
    return sorted(candidates)


def coincident_bb_reference(evaluated):
    """Only annotate equality when both curves have the same defined values."""
    first, second = evaluated.get("bitbirch"), evaluated.get("bblean")
    if not first or not second or len(first) != len(second):
        return False
    defined = 0
    for left, right in zip(first, second):
        if (left[0] is None) != (right[0] is None):
            return False
        if left[0] is not None:
            defined += 1
            if left[0] != right[0]:
                return False
    return defined > 0


def render_figures(summary, root, records, complete):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory
    import numpy as np

    contract = summary_contract(summary)
    is_pilot = contract["is_pilot"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
        "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 10,
        "axes.spines.top": True, "axes.spines.right": True, "axes.grid": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})
    artifacts = []
    for scope, study, title, xlabel, methods, group_rows in groups(records):
        cases = sorted({(row["x"], row["case_id"]) for row in group_rows})
        xs = np.asarray([x for x, _ in cases], dtype=float)
        if len(set(xs)) != len(xs):
            raise ValueError("A plotted study has multiple settings at one x")
        lookup = {(row["case_id"], row["method"]): row for row in group_rows}
        method_rows = {method: [lookup.get((case, method)) for _, case in cases] for method in methods}
        span = float(xs[-1] - xs[0])
        pad = 0.035 * span if span else max(0.02, abs(float(xs[0])) * .04)
        destination = root / "plots" / scope
        destination.mkdir(parents=True, exist_ok=True)
        for metric, _, ylabel in METRICS:
            evaluated = {method: [value_and_status(row, metric) for row in method_rows[method]] for method in methods}
            missing = [(method, values) for method, values in evaluated.items() if any(mark for _, mark, _ in values)]
            used_marks = sorted({mark for _, values in missing for _, mark, _ in values if mark})
            notes = measurement_notes(scope, study, metric, group_rows, complete)
            if is_pilot:
                actual_rows = sorted({row["n_total"] for row in group_rows if isinstance(row.get("n_total"), int)})
                actual = str(actual_rows[0]) if len(actual_rows) == 1 else "/".join(map(str, actual_rows))
                notes = ["PILOT CHECK ONLY - not full experiment results.",
                         f"X positions identify intended full settings; actual pilot input rows: {actual or 'see records'}."] + notes
                notes = [note.replace("radius study uses 10,000 points.", "radius study is a reduced-row pilot.") for note in notes]
                notes = [note.replace("saved partitions checked across repeats", "one pilot fit per setting")
                         .replace("one score per setting after checking saved partitions across repeats", "one pilot fit per setting") for note in notes]
            coincident = scope == "landrum" and coincident_bb_reference(evaluated)
            if coincident:
                notes.append("BitBIRCH and BitBIRCH-Lean values coincide; both curves are drawn without shifting values.")
            note_height = .024 * len(notes)
            key_lines = textwrap.wrap("Status strips: " + "; ".join(f"{mark}={STATUS_NAMES[mark]}" for mark in used_marks), 105) if used_marks else []
            note_top = .020 + .024 * (len(notes) - 1)
            key_top = note_top + .032 + .023 * max(0, len(key_lines) - 1)
            status_bottom = key_top + .032
            bottom = max(.20, .13 + note_height + .027 * len(missing) + .023 * len(key_lines))
            fig, ax = plt.subplots(figsize=(7.4, 5.4))
            fig.subplots_adjust(left=.14, right=.975, top=.90, bottom=bottom)
            if is_pilot:
                ax.text(.5, .55, "PILOT", transform=ax.transAxes, ha="center", va="center",
                        rotation=25, fontsize=54, color="#999999", alpha=.18, zorder=0)
            handles, all_values = {}, []
            bands = []
            # Draw the open reference markers first, retaining the requested legend order.
            draw_methods = sorted(methods, key=lambda method: method != "bblean")
            for method in draw_methods:
                label, color, marker, linestyle = METHODS[method]
                if method == "dbscan" and study == "popcount_separation":
                    linestyle = "--"
                values = evaluated[method]
                ys = np.asarray([value if value is not None else np.nan for value, _, _ in values])
                all_values.extend(ys[np.isfinite(ys)].tolist())
                line_options = {"markerfacecolor": "white", "alpha": .80} if method == "bblean" else {}
                marker_indices = finite_marker_indices(ys)
                if method == "bblean" and len(xs) > 20:
                    # Offset only marker positions, never the plotted coordinates or scores.
                    marker_indices = [i for i in range(4, len(xs), 10) if np.isfinite(ys[i])]
                ax.plot(xs, ys, label=label, color=color, marker=marker, linestyle=linestyle,
                        linewidth=1.3 if method == "bblean" else 1.5, markersize=4.5 if len(xs) > 20 else 5.5,
                        markevery=marker_indices, **line_options)
                handles[method] = Line2D([0], [0], label=label, color=color, marker=marker, linestyle=linestyle,
                                        linewidth=1.5, markersize=5, **line_options)
                if metric == "silhouette":
                    halfwidth = np.asarray([1.96 * row["silhouette_se"] if row is not None and value_and_status(row, metric)[0] is not None
                                            and not row["exact_silhouette"] else np.nan for row in method_rows[method]])
                    sampled = np.isfinite(halfwidth) & np.isfinite(ys)
                    if sampled.any():
                        low, high = ys - halfwidth, ys + halfwidth
                        # Pointwise descriptive Monte Carlo bounds, not a simultaneous interval.
                        ax.fill_between(xs, low, high, where=sampled, color=color, alpha=.12, linewidth=0)
                        # Isolated sampled points still show their SE when there is no neighbouring band segment.
                        if np.count_nonzero(sampled) == 1:
                            ax.errorbar(xs[sampled], ys[sampled], yerr=halfwidth[sampled], fmt="none", color=color, alpha=.5, capsize=2)
                        bands.append({"method": method, "sampled_points": int(sampled.sum())})
                        all_values.extend(low[sampled].tolist() + high[sampled].tolist())
            ax.set_title(("PILOT - " if is_pilot else "") + ("PARTIAL - " if not complete else "") + title)
            ax.set_xlabel(xlabel + (" (intended setting)" if is_pilot else ""))
            ax.set_ylabel(ylabel)
            ax.set_xlim(float(xs[0]) - pad, float(xs[-1]) + pad)
            if len(xs) <= 15:
                ax.set_xticks(xs)
            if study == "sample_count":
                ax.tick_params(axis="x", labelsize=8)
            if metric == "coverage":
                ax.set_ylim(-.025, 1.055)
            elif metric == "ami":
                ax.set_ylim(min(-.05, min(all_values, default=0) - .04), 1.055)
            else:
                ax.set_ylim(min(-.05, min(all_values, default=0) - .05), max(1.05, max(all_values, default=1) + .04))
            ax.legend(handles=[handles[method] for method in methods], loc="best", fontsize=8, frameon=True)
            if not all_values:
                ax.text(.5, .48, "No defined scores", transform=ax.transAxes, ha="center", color="#777777")
            for index, line in enumerate(notes):
                fig.text(.5, note_top - .024 * index, line, ha="center", fontsize=7)
            transform = blended_transform_factory(ax.transData, fig.transFigure)
            status_counts = {}
            for index, (method, values) in enumerate(missing):
                y = status_bottom + .027 * (len(missing) - index - 1)
                color = METHODS[method][1]
                counts = Counter(mark for _, mark, _ in values if mark)
                status_counts[method] = dict(counts)
                label = SHORT_NAMES[method] + " " + "/".join(f"{mark}:{count}" for mark, count in sorted(counts.items()))
                fig.text(.132, y, label, ha="right", va="center", fontsize=6.4, color=color)
                for x, (_, mark, _) in zip(xs, values):
                    if not mark:
                        continue
                    # Dense 100-threshold strips use a narrow tick; the row label identifies U/M/N counts.
                    symbol = mark if len(xs) <= 25 else {"U": "|", "M": "x", "N": ".", "P": ".", "E": "+"}[mark]
                    fig.text(x, y, symbol, transform=transform, ha="center", va="center", fontsize=7, color=color, fontweight="bold")
            for index, line in enumerate(key_lines):
                if len(xs) > 25:
                    line += " (dense strips: U=|, M=x, N/P=.)" if index == len(key_lines) - 1 else ""
                fig.text(.5, key_top - .023 * index, line, ha="center", fontsize=7)
            files = []
            for extension in ("png", "pdf", "svg"):
                path = destination / f"{study}_{metric}.{extension}"
                fig.savefig(path, dpi=300, facecolor="white")
                files.append({"file": path.relative_to(root).as_posix(), "sha256": sha(path)})
            plt.close(fig)
            artifacts.append({"scope": scope, "study": study, "metric": metric, "files": files,
                              "yscale": "linear", "x_limits": [float(xs[0]) - pad, float(xs[-1]) + pad],
                              "status_counts": status_counts, "notes": notes, "sampling_bands": bands,
                              "bitbirch_reference_coincident": coincident, "is_pilot": is_pilot})
    return artifacts


def availability_table(records):
    lines = ["| Scope / study | Method | Records | Fit failures | AMI all available | Silhouette available (exact / sampled) | Undefined / unavailable | Coverage available |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    reasons = []
    for scope, study, _, _, methods, rows in groups(records):
        for method in methods:
            subset = [row for row in rows if row["method"] == method]
            ami = sum(value_and_status(row, "ami")[0] is not None for row in subset)
            coverage = sum(value_and_status(row, "coverage")[0] is not None for row in subset)
            valid = [row for row in subset if value_and_status(row, "silhouette")[0] is not None]
            exact = sum(row["exact_silhouette"] for row in valid)
            missing = Counter(value_and_status(row, "silhouette")[1] for row in subset)
            failed = sum(row["status"] != "ok" for row in subset)
            lines.append(f"| {scope} / {study} | {METHODS[method][0]} | {len(subset)} | {failed} | {ami} | {len(valid)} ({exact} / {len(valid)-exact}) | {missing['U']} / {missing['N']} | {coverage} |")
            counts = Counter(reason for row in subset for _, mark, reason in [value_and_status(row, "silhouette")] if mark)
            if counts:
                reasons.append((scope, study, method, counts))
    lines += ["", "Unavailable silhouette reasons:", ""]
    if not reasons:
        lines.append("None.")
    for scope, study, method, counts in reasons:
        lines.append(f"- {scope}/{study}, {METHODS[method][0]}: " + "; ".join(f"{reason}: {count}" for reason, count in sorted(counts.items())) + ".")
    return lines


def write_report(summary, root, summary_path, records, complete, artifacts):
    counts = Counter((row["scope"], row["status"]) for row in records)
    contract = summary_contract(summary)
    scopes = set(contract["cases"])
    lines = ["# Saved-partition AMI and Tanimoto silhouette", "",
        ("Complete setting coverage." if complete else "**PARTIAL RESULTS: this report does not cover every planned setting.**"), "",
        "This report reads saved score records. It performs no clustering or metric evaluations. Numerical values and scoring provenance remain in [summary.json](summary.json). The main AMI figures use `ami_all`; `ami_retained` remains an explicitly different population in the numerical records.", "",
        "One score is shown per method/setting after the scoring pipeline checks repeated saved partitions. Clustering repeats are not treated as independent silhouette estimates. For sampled silhouettes, shading denotes pointwise mean +/- 1.96 reported Monte Carlo standard errors, not uncertainty across fits or a simultaneous confidence band. Exact full-population silhouettes have no sampling band. The scoring protocol uses a shared seed-101 random priority per dataset and the first min(2000, retained population) eligible query rows. Equal populations share queries; differently filtered populations may use different queries. Each query is compared with every retained reference row. Reported standard errors include the finite-population correction sqrt(1-m/N).", "",
        "Silhouette uses Tanimoto distances in the retained population. Singleton silhouettes use zero; fewer than two retained clusters or a cluster for every retained point is undefined. Undefined scores and failed fits remain disconnected gaps. M denotes a memory-limited fit; U an undefined score; N an unavailable score. Their status strips do not represent numerical observations.", "",
        "AMI compares every saved row with its reference category, including any -1 noise assignments as a predicted category. For Landrum, silhouette and coverage follow the earlier quality populations: all nonnoise final CLASSIX clusters after minPts=10, and only BitBIRCH/Lean clusters of size strictly greater than 10. These populations need not match; compare coverage when interpreting silhouette differences. CLASSIX minPts performs reassignment during fitting and does not guarantee every final cluster reaches 10 points. Synthetic coverage likewise records removal of any algorithm noise; it does not invent assignments for failed fits. Where the saved BitBIRCH and Lean values coincide exactly, the Landrum figure states this explicitly and uses open reference markers at different x positions without shifting any values.", "",
        "The Landrum input has 90,498 pKi>4 records, 35 target categories and 55,513 distinct literal SMILES. It preserves 34,985 repeated-SMILES rows beyond the first occurrence. Of the distinct SMILES, 14,533 appear against multiple targets, accounting for 42,162 rows. Thus these reference categories are not a unique structure-based ground truth. Duplicates are retained in both external-label and internal-distance scores.", "",
        "## Availability", "",
        "; ".join(f"{scope}: {status}={count}" for (scope, status), count in sorted(counts.items())) + ".", ""]
    if "landrum" not in scopes:
        lines = [line for line in lines if not line.startswith("The Landrum input has")]
        lines = [line.replace("For Landrum, silhouette and coverage follow the earlier quality populations: all nonnoise final CLASSIX clusters after minPts=10, and only BitBIRCH/Lean clusters of size strictly greater than 10. These populations need not match; compare coverage when interpreting silhouette differences. CLASSIX minPts performs reassignment during fitting and does not guarantee every final cluster reaches 10 points. ", "").replace(" Where the saved BitBIRCH and Lean values coincide exactly, the Landrum figure states this explicitly and uses open reference markers at different x positions without shifting any values.", "") for line in lines]
    if contract["is_pilot"]:
        lines[0] = "# PILOT ONLY: saved-partition AMI and Tanimoto silhouette"
        lines[2] = "**PILOT CHECK ONLY. These reduced-row results are not the full experiment or publication figures.**"
        lines.insert(4, "Axes identify the selected intended full experiment settings. Actual input rows are recorded per setting; pilot estimates and populations can differ substantially from the full experiment.\n")
        lines = [line.replace("The Landrum input has 90,498", "The full Landrum source (not this reduced-row pilot) has 90,498") for line in lines]
    lines += availability_table(records)
    lines += ["", "## Figures", "", "Every linked PNG also has PDF and SVG versions in the same directory.", ""]
    for scope, study, title, _, _, _ in groups(records):
        links = []
        for metric, _, _ in METRICS:
            links.append(f"[{metric}](../plots/{scope}/{study}_{metric}.png)")
        lines.append(f"- {scope}: {title} - " + ", ".join(links) + ".")
    lines += ["", "## Provenance", "",
        f"- Score summary SHA256: `{sha(summary_path)}`.",
        f"- Renderer SHA256: `{sha(__file__)}`.",
        "- Figure hashes, availability markers and measurement notes: [rendering_provenance.json](rendering_provenance.json).",
        "- Method colours and symbols follow the earlier saved-run renderer; all quality axes are linear. BitBIRCH-Lean remains a separately labelled optimized reference.", ""]
    for key in ("protocol", "scoring_protocol", "sampling", "sources", "provenance"):
        if key in summary:
            lines += [f"Recorded `{key}`:", "", "```json", json.dumps(summary[key], indent=2, allow_nan=False), "```", ""]
    output = root / "results" / "SCORES.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--check", action="store_true", help="Validate record structure without importing plotting libraries or writing output")
    parser.add_argument("--allow-partial", action="store_true", help="Create explicitly labelled previews from incomplete setting coverage")
    args = parser.parse_args()
    root = args.root.resolve()
    summary_path = args.summary.resolve() if args.summary else root / "results" / "summary.json"
    initial_hash = sha(summary_path)
    summary = read(summary_path)
    records, complete = validate_records(summary, args.allow_partial)
    contract = summary_contract(summary)
    pending = sum(row.get("status") in ("pending", "running") or (
        row.get("status") == "ok" and row.get("silhouette_status") in ("pending", "running")) for row in records)
    if args.check:
        print(json.dumps({"records": len(records), "complete_setting_coverage": complete,
                          "scope": contract["scope"], "is_pilot": contract["is_pilot"],
                          "pending_scores": pending, "status": "ok", "writes": 0}))
        return 0
    if pending and not args.allow_partial:
        raise ValueError("Scoring is still pending; wait for completion or explicitly use --allow-partial")
    complete = complete and pending == 0
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    artifacts = render_figures(summary, root, records, complete)
    expected_figures = 3 * sum(1 for _ in groups(records))
    if complete and len(artifacts) != expected_figures:
        raise ValueError("A complete report requires three figures per present study")
    if complete and not contract["is_pilot"] and len(artifacts) != {"synthetic": 15, "landrum": 3, "both": 18}[contract["scope"]]:
        raise ValueError("The full report does not cover every study")
    if sha(summary_path) != initial_hash:
        raise RuntimeError("Score summary changed during rendering")
    report = write_report(summary, root, summary_path, records, complete, artifacts)
    save({"status": "complete" if complete else "partial", "operation": "saved-score rendering only",
          "generated_utc": datetime.now(timezone.utc).isoformat(), "records": len(records),
          "scope": contract["scope"], "is_pilot": contract["is_pilot"],
          "summary_sha256": initial_hash, "renderer_sha256": sha(__file__),
          "fits": 0, "metric_evaluations": 0, "report": report.relative_to(root).as_posix(),
          "figures": artifacts}, root / "results" / "rendering_provenance.json")
    print(json.dumps({"status": "complete" if complete else "partial", "figures": len(artifacts),
                      "files": sum(len(item["files"]) for item in artifacts), "report": str(report)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
