from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import urllib.request

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
COMMON = HERE / "common"
WHEEL_URL = "https://files.pythonhosted.org/packages/b7/2b/f7f088fe70597d72f3705ee050e71110ca3a56e8bb5c1f954025c6f9fac9/bblean-0.10.3-cp312-cp312-win_amd64.whl"
WHEEL_SHA = "36d10341de071024a0585020d5907d12950c93e251c03f74e92a14eee6a78fd9"
THREADS = {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}

class BinaryExperiment:
    def __init__(self, study, work_dir, mode="full", pilot_rows=256, python=None):
        assert study in {"synthetic", "molecular", "radius"}
        assert mode in {"full", "pilot", "check"}
        self.study, self.mode, self.pilot_rows = study, mode, pilot_rows
        self.source = HERE / study
        self.work = Path(work_dir).resolve()
        self.python = Path(python).resolve() if python else self.work / ".venv/Scripts/python.exe"
        self.external_python = python is not None
        self.work.mkdir(parents=True, exist_ok=True)
        self.env = dict(os.environ, **THREADS, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg",
                        CLASSIX_EXPERIMENT_WORK=str(self.work), CLASSIX_RUN_MODE=mode,
                        CLASSIX_PILOT_ROWS=str(pilot_rows), CLASSIX_SILHOUETTE_DIR=str(self.work / "native"),
                        PYTHONPATH=os.pathsep.join(map(str, [REPO / "src", COMMON, HERE / "vendor/bitbirch"])))
        self.run_dir = self.work / "runs" / ("landrum_fixed" if study == "molecular" else ("pilot" if mode == "pilot" else "reproduction"))

    def command(self, arguments, **kwargs):
        process = subprocess.Popen(list(map(str, arguments)), cwd=self.work, env=self.env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                   encoding="utf-8", errors="replace",
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0, **kwargs)
        try:
            for line in process.stdout:
                print(line, end="", flush=True)
            returncode = process.wait()
        except BaseException:
            process.terminate()
            process.wait()
            raise
        if returncode:
            raise subprocess.CalledProcessError(returncode, process.args)
        return returncode

    def script(self, name, *arguments):
        return self.command([self.python, "-B", self.source / "code" / name, *arguments])

    def setup(self, install=True):
        if os.name != "nt":
            raise RuntimeError("These recorded benchmarks require Windows x64 and Python 3.12.")
        if not self.python.exists() and not self.external_python:
            base = [sys.executable] if sys.version_info[:2] == (3, 12) else ["py", "-3.12"]
            self.command([*base, "-m", "venv", self.work / ".venv"])
        if install:
            self.command([self.python, "-m", "pip", "install", "-r", HERE / "requirements.txt"])
            wheel = self.work / "downloads" / WHEEL_URL.rsplit("/", 1)[-1]
            wheel.parent.mkdir(exist_ok=True)
            if not wheel.exists():
                urllib.request.urlretrieve(WHEEL_URL, wheel)
            assert hashlib.sha256(wheel.read_bytes()).hexdigest() == WHEEL_SHA
            self.command([self.python, "-m", "pip", "install", "--no-deps", wheel])
        self.command([self.python, "-B", COMMON / "check_environment.py", HERE / "requirements.txt"])
        source = REPO / "src/spmv/spsubmatxvec.c"
        library = REPO / "src/spmv/spsubmatxvec.dll"
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        library_hash = hashlib.sha256(library.read_bytes()).hexdigest() if library.exists() else None
        current = False
        for record in (REPO / "_build/native").glob("*/build_metadata.json"):
            info = json.loads(record.read_text(encoding="utf-8-sig"))
            if info.get("source_sha256") == source_hash and info.get("dll_sha256") == library_hash:
                current = True
                break
        if not current:
            self.command(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", REPO / "scripts/build_native.ps1"])
        self.command([self.python, "-c", "from classix_t import CLASSIX_T; from spmv import HARDWARE_POPCOUNT; print('CLASSIX_T ready; POPCNT:', HARDWARE_POPCOUNT)"])
        return self.python

    def prepare_data(self):
        if self.study != "molecular":
            return
        source = REPO / "data/BL_Sets_data_rdkitblog.csv"
        assert hashlib.sha256(source.read_bytes()).hexdigest() == "ddc8dbf5d92b1f8255ef210456c008dd5bf91866bc9f93870592ebcb5868ba1d"
        target = self.work / "data/source" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copyfile(source, target)
        assert target.read_bytes() == source.read_bytes()
        args = ["--dataset", "landrum"]
        if self.mode == "pilot":
            args += ["--pilot-rows", self.pilot_rows]
        self.script("prepare_inputs.py", *args)

    def run(self):
        if self.study == "molecular":
            args = ["--config", self.source / "config.py", "--output-dir", self.run_dir]
            self.script("run_benchmark.py", *args, "--check")
            if self.mode != "check":
                self.script("run_benchmark.py", *args, "--run")
        else:
            args = ["--check" if self.mode == "check" else "--run", "--output-dir", self.run_dir]
            if self.mode == "pilot":
                case = "sample_count_10000" if self.study == "synthetic" else "radius_0.3_n10000_d1024"
                args += ["--cases", case, "--pilot-rows", self.pilot_rows]
            self.script("reproduce.py", *args)
        return self.run_dir

    def plot(self, plot_dir):
        plot_dir = Path(plot_dir).resolve()
        if self.mode == "check":
            return None
        if self.study == "synthetic":
            self.script("finish_report.py", self.work, self.run_dir, plot_dir)
            return self.run_dir / "report"
        if self.study == "radius" and self.mode == "full":
            output = self.work / "derived/radius_sensitivity"
            self.script("score_and_plot.py", "--run-dir", self.run_dir, "--output-dir", output)
            result = json.loads((output / "verification.json").read_text(encoding="utf-8"))
            assert result["status"] == "passed" and result["full_method_settings"] == 30
            assert result["all_methods_complete"] and result["scored_rows_per_fit"] == 10000
            plot_dir.mkdir(exist_ok=True, parents=True)
            for artifact in result["artifacts"]:
                source = Path(artifact["path"])
                assert source.parent.resolve() == (output / "plots").resolve()
                assert hashlib.sha256(source.read_bytes()).hexdigest() == artifact["sha256"]
                shutil.copyfile(source, plot_dir / source.name)
            return output
        return None

    def score(self, plot_dir):
        if self.mode == "check" or (self.study == "molecular" and self.mode == "pilot"):
            return None
        assert not list((self.work / "runs").rglob("run.lock"))
        source = COMMON / "silhouette_kernel.c"
        library = self.work / "native/silhouette_kernel.dll"
        record = self.work / "native/build_metadata.json"
        current = False
        if record.exists() and library.exists():
            info = json.loads(record.read_text(encoding="utf-8-sig"))
            current = (info["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
                       and info["dll_sha256"] == hashlib.sha256(library.read_bytes()).hexdigest())
        if not current:
            self.command(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                          COMMON / "build_silhouette.ps1", "-OutputDirectory", self.work / "native"])
        if self.study == "synthetic":
            self.script("finish_scores.py", Path(plot_dir).resolve())
            return self.work / "derived" / ("ami_silhouette_pilot" if self.mode == "pilot" else "ami_silhouette_full")
        if self.study == "molecular":
            self.script("score_landrum.py", "--root", self.work,
                        "--data-dir", self.work / "data/prepared", "--results-dir", self.run_dir / "results",
                        "--native-code", COMMON / "native_silhouette.py", "--chi-code", self.source / "code/chi_tanimoto.py")
            if self.mode == "full":
                destination = Path(plot_dir).resolve(); destination.mkdir(exist_ok=True, parents=True)
                for source in (self.work / "plots").glob("*"):
                    if source.is_file(): shutil.copyfile(source, destination / source.name)
            return self.work
