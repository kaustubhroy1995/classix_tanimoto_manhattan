from pathlib import Path
import importlib.metadata as metadata
import json
import platform
import sys
expected = dict(line.split("==") for line in Path(sys.argv[1]).read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))
expected["bblean"] = "0.10.3"
assert sys.version_info[:2] == (3, 12) and sys.maxsize > 2**32 and platform.machine().lower() in {"amd64", "x86_64"}
actual = {name: metadata.version(name) for name in expected}
mismatch = {name: (actual[name], version) for name, version in expected.items() if actual[name] != version}
if mismatch: raise RuntimeError(f"Package versions differ (actual, recorded): {mismatch}")
record = {"python": sys.version, "platform": platform.platform(), "packages": actual}
Path("environment.json").write_text(json.dumps(record, indent=2) + "\n")
print("Pinned Python environment verified.")
