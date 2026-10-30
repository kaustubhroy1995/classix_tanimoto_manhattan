"""Render the paper's four-method view from a verified Experiment 2 run.

Reads generated results only. Does not fit, retune, download, or edit a paper.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics


METHODS = (
    ("classix", "manhattan", "CLASSIX M"),
    ("classix", "euclidean", "CLASSIX E"),
    ("dbscan", "manhattan", "DBSCAN M"),
    ("optics", "manhattan", "OPTICS M"),
)
NAMES = {
    "iris": "Iris", "banknote": "Banknote", "wine": "Wine",
    "breast_cancer": "Breast cancer", "digits": "Digits", "seeds": "Seeds",
    "glass": "Glass", "sonar": "Sonar", "vehicle": "Vehicle",
    "segmentation": "Segmentation", "breast_tissue": "Breast tissue",
    "leaf": "Leaf", "libras": "Libras", "vertebral": "Vertebral",
    "wireless": "Wireless", "parkinsons": "Parkinsons", "coimbra": "Coimbra",
    "raisin": "Raisin", "haberman": "Haberman", "musk": "Musk",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def key(row):
    return row["dataset"], row["algorithm"], row["metric"]


def format_value(value: float, best: float, score: bool) -> str:
    rendered = f"{value:.4f}" if score else f"{value:.3g}"
    if math.isclose(value, best, rel_tol=0, abs_tol=1e-12):
        return r"\textbf{" + rendered + "}"
    return rendered


def render(run_dir: Path, output_dir: Path) -> None:
    policy = read_json(run_dir / "protocol.json")
    datasets = policy["datasets"]
    require(len(datasets) == 20 and len(set(datasets)) == 20,
            "Expected the twenty-dataset Experiment 2 roster")
    require(set(datasets) == set(NAMES), "Dataset roster differs from the paper")
    selected = read_json(run_dir / "results/selected.json")
    lookup = {key(row): row for row in selected}
    require(len(lookup) == len(selected), "Duplicate selected-result keys")
    with (run_dir / "tables/all_20.csv").open(newline="", encoding="utf-8-sig") as stream:
        report_rows = list(csv.DictReader(stream))
    metadata = {key(row): row for row in report_rows}
    require(len(metadata) == len(report_rows), "Duplicate full-report keys")

    quality = [r"\begin{tabular}{llrrrrrrrr}", r"\toprule",
               r"Dataset & $(n,d)$ & \multicolumn{2}{c}{CLASSIX\_M} & \multicolumn{2}{c}{CLASSIX} & \multicolumn{2}{c}{DBSCAN\_M} & \multicolumn{2}{c}{OPTICS\_M} \\",
               r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
               r" & & ARI & AMI & ARI & AMI & ARI & AMI & ARI & AMI \\", r"\midrule"]
    runtime = ["% Median complete fit times in milliseconds; five fits after one warmup.",
               r"\begin{tabular}{lrrrr}", r"\toprule",
               r"Dataset & CLASSIX\_M & CLASSIX & DBSCAN\_M & OPTICS\_M \\", r"\midrule"]
    rows_out = []
    for dataset in datasets:
        rows = []
        dimensions = set()
        for algorithm, metric, _ in METHODS:
            identity = dataset, algorithm, metric
            require(identity in lookup and identity in metadata,
                    f"Missing selected result or report: {identity}")
            row, report = lookup[identity], metadata[identity]
            for field in ("ari", "ami", "fit_seconds_median"):
                value = float(row[field])
                require(math.isfinite(value) and value == float(report[field]),
                        f"Selected/report mismatch: {identity}, {field}")
            times = row["timing_seconds"]
            require(len(times) == 5 and all(math.isfinite(t) and t >= 0 for t in times),
                    f"Expected five finite fit times: {identity}")
            require(statistics.median(times) == row["fit_seconds_median"],
                    f"Timing median differs: {identity}")
            n, d = int(report["n"]), int(report["d"])
            require(n > 0 and d > 0, f"Invalid dimensions: {identity}")
            dimensions.add((n, d))
            rows.append(row)
            rows_out.append(dict(dataset=dataset, algorithm=algorithm, metric=metric,
                                 n=n, d=d, ari=row["ari"], ami=row["ami"],
                                 fit_seconds_median=row["fit_seconds_median"],
                                 fit_ms_median=1000 * row["fit_seconds_median"]))
        require(len(dimensions) == 1, f"Methods use different dimensions: {dataset}")
        n, d = dimensions.pop()
        cells = [NAMES[dataset], f"$({n},{d})$"]
        for row in rows:
            for score in ("ari", "ami"):
                cells.append(format_value(row[score], max(r[score] for r in rows), True))
        quality.append(" & ".join(cells) + r" \\")
        milliseconds = [1000 * row["fit_seconds_median"] for row in rows]
        runtime.append(" & ".join([NAMES[dataset]] + [
            format_value(t, min(milliseconds), False) for t in milliseconds]) + r" \\")
    quality.extend([r"\bottomrule", r"\end{tabular}"])
    runtime.extend([r"\bottomrule", r"\end{tabular}"])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "quality_table.tex").write_text("\n".join(quality) + "\n", encoding="utf-8")
    (output_dir / "runtime_table.tex").write_text("\n".join(runtime) + "\n", encoding="utf-8")
    with (output_dir / "four_method_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_out[0]))
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"Wrote {len(rows_out)} verified method/dataset rows to {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True,
                        help="Completed Experiment 2 density_baselines directory")
    parser.add_argument("--output-dir", type=Path,
                        help="Default: RUN_DIR/tables/manuscript")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    output_dir = args.output_dir.resolve() if args.output_dir else run_dir / "tables/manuscript"
    render(run_dir, output_dir)


if __name__ == "__main__":
    main()
