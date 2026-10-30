"""Locations of source code and generated numeric experiment files."""
import hashlib
import json
import os
from pathlib import Path

HELPERS = Path(__file__).resolve().parent
REPO = HELPERS.parents[2]
RUN = Path(os.environ.get("CLASSIX_NUMERIC_RUN", HELPERS.parent / "_generated")).resolve()
DENSITY = RUN / "density_baselines"
MANHATTAN = REPO / "src" / "classix_m.py"
PROVENANCE = RUN / "manhattan_source.json"

def ensure_source_provenance():
    record = {
        "compatible_file": "src/classix_m.py",
        "compatible_sha256": hashlib.sha256(MANHATTAN.read_bytes()).hexdigest(),
        "implementation": "Canonical repository CLASSIX_M",
    }
    RUN.mkdir(parents=True, exist_ok=True)
    if PROVENANCE.exists():
        if json.loads(PROVENANCE.read_text(encoding="utf-8")) != record:
            raise RuntimeError("Manhattan source changed; use a fresh generated run directory.")
    else:
        PROVENANCE.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record
