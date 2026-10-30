from pathlib import Path
import sys
import os
import hashlib
import json
import shutil
import subprocess
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import work_root, read_configuration
HERE = Path(__file__).resolve().parent
RUNTIME_ROOT = work_root(HERE.parent)
CONFIG = HERE.parent / "config.py"
PLAN = read_configuration(CONFIG)
RUN_MODE = os.environ.get("CLASSIX_RUN_MODE", "full")
RUN_DIR = RUNTIME_ROOT / "runs" / ("pilot" if RUN_MODE == "pilot" else "reproduction")
PLOT_DIR = Path(sys.argv[1])
PILOT_CASE = "sample_count_10000"
PILOT_ROWS = int(os.environ.get("CLASSIX_PILOT_ROWS", "256"))
python = sys.executable
assert not list((RUNTIME_ROOT / "runs").rglob("run.lock")), "Wait for all timed fit controllers to exit."
metrics_output = RUNTIME_ROOT / "derived" / ("ami_silhouette_pilot" if RUN_MODE == "pilot" else "ami_silhouette_full")
metrics_command = [str(python), "-B", str(HERE / "saved_scores.py"),
                       "--scope", "synthetic", "--config", str(CONFIG),
                       "--results-dir", str(RUN_DIR / "results"), "--data-dir", str(RUN_DIR / "data"),
                       "--output-dir", str(metrics_output)]
if RUN_MODE == "pilot":
        metrics_command += ["--cases", PILOT_CASE, "--pilot-rows", str(PILOT_ROWS)]
print("Scoring saved assignments:", subprocess.list2cmdline(metrics_command))
subprocess.run(metrics_command, cwd=RUNTIME_ROOT, check=True, creationflags=subprocess.CREATE_NO_WINDOW)
original_report = json.loads((RUN_DIR / "report/summary.json").read_text(encoding="utf-8"))
metric_summary_path = metrics_output / "results/summary.json"
metric_summary = json.loads(metric_summary_path.read_text(encoding="utf-8"))
metric_rendering = json.loads((metrics_output / "results/rendering_provenance.json").read_text(encoding="utf-8"))
metric_verification = json.loads((metrics_output / "verification.json").read_text(encoding="utf-8"))
metric_summary_hash = hashlib.sha256(metric_summary_path.read_bytes()).hexdigest()
assert metric_summary["scope"] == metric_rendering["scope"] == "synthetic"
assert metric_summary["is_pilot"] == metric_rendering["is_pilot"] == (RUN_MODE == "pilot")
assert metric_rendering["status"] == "complete" and metric_verification["status"] == "passed"
assert metric_rendering["summary_sha256"] == metric_verification["summary_sha256"] == metric_summary_hash
assert metric_rendering["renderer_sha256"] == hashlib.sha256((HERE / "render_scores.py").read_bytes()).hexdigest()
assert metric_summary["complete_case_results"] and not any(row["silhouette_status"] == "pending" for row in metric_summary["records"])
expected_metric_cases = [PILOT_CASE] if RUN_MODE == "pilot" else [case["case_id"] for case in PLAN["cases"]]
assert set(metric_summary["expected_cases"]) == set(expected_metric_cases)
assert set(metric_summary["expected_methods"]) == set(PLAN["methods"])
metric_files = [item for figure in metric_rendering["figures"] for item in figure["files"]]
assert len({item["file"] for item in metric_files}) == len(metric_files)
for item in metric_files:
        figure = metrics_output / item["file"]
        assert figure.resolve().is_relative_to((metrics_output / "plots/synthetic").resolve())
        assert hashlib.sha256(figure.read_bytes()).hexdigest() == item["sha256"]
metric_figures = sorted(p for p in (metrics_output / "plots/synthetic").iterdir()
                            if p.is_file() and p.suffix.lower() in {".png", ".pdf", ".svg"})
assert {str(p.relative_to(metrics_output).as_posix()) for p in metric_figures} == {item["file"] for item in metric_files}
assert metric_figures, "No AMI/silhouette/coverage figures were generated."
metric_destination = metrics_output / "plots/synthetic"
if RUN_MODE == "full" and original_report["complete_case_results"] and not original_report["is_pilot"]:
        assert len(metric_figures) == 45, "A complete synthetic sweep must produce 15 figures in three formats."
        metric_destination = PLOT_DIR.resolve()
        metric_destination.mkdir(parents=True, exist_ok=True)
        for source in metric_figures:
            target = metric_destination / source.name
            shutil.copyfile(source, target)
            assert hashlib.sha256(target.read_bytes()).digest() == hashlib.sha256(source.read_bytes()).digest()
print("AMI/silhouette report:", metrics_output)
print("AMI/silhouette/coverage plots:", metric_destination)
