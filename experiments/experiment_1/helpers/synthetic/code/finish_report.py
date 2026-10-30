import json
from pathlib import Path
import sys
from notebook_report import collect_fresh, verify_fresh, figures
root, run_dir, requested = map(Path, sys.argv[1:4])
assert not list((root/"runs").rglob("run.lock")), "Wait for fit controllers to exit."
verification = verify_fresh(root, run_dir)
summary = collect_fresh(root, run_dir)
report = run_dir/"report"
report.mkdir(exist_ok=True)
(report/"verification.json").write_text(json.dumps(verification, indent=2, allow_nan=False)+"\n")
destination = requested if summary["complete_case_results"] and not summary["is_pilot"] else report/"plots"
artifacts = figures(summary, destination)
summary["figures"] = artifacts
(report/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
print(json.dumps({"successful_assignments_checked": verification["successful_assignments_checked"],
                  "terminal_failures": len(verification["terminal_failures"]),
                  "figures": len(artifacts), "plots": str(destination), "report": str(report)}, indent=2))
