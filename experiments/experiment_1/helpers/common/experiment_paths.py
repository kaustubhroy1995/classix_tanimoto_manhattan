from pathlib import Path
import copy
import os
import runpy

HELPERS = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = HELPERS.parent
REPO = NOTEBOOK_DIR.parents[1]
COMMON = HELPERS / "common"
VENDOR = HELPERS / "vendor/bitbirch"

def work_root(source):
    return Path(os.environ.get("CLASSIX_EXPERIMENT_WORK", NOTEBOOK_DIR / "_generated" / Path(source).name)).resolve()

def read_configuration(path):
    values = runpy.run_path(str(path))
    return copy.deepcopy(values["get_design"]() if "get_design" in values else values["DESIGN"])

def read_repetitions(path):
    return copy.deepcopy(runpy.run_path(str(path)).get("REPETITIONS", {}))
