"""One fit with explicit pre-fit/post-fit checkpoints and separate validation."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
import argparse
import gc
import json
from pathlib import Path
import sys
import traceback
from time import perf_counter
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = work_root(SOURCE_ROOT)
BASE = ROOT
PACKED = REPO / "src"
from memory_win import snapshot


def send(event, **values):
    print("BENCH_JSON " + json.dumps({"event": event, **values}), flush=True)


def handshake(expected):
    received = sys.stdin.readline().strip()
    if received != expected:
        raise RuntimeError(f"Expected {expected}, received {received}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rows", type=int)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    import numpy as np
    from threadpoolctl import threadpool_limits, threadpool_info
    plan = read_configuration(args.config)
    case = next(c for c in plan["cases"] if c["case_id"] == args.case)
    directory = args.data_dir/case["data_id"]
    manifest = json.loads((directory/"manifest.json").read_text())
    parameters = {}
    library = {}
    if args.variant.startswith("classix"):
        sys.path.insert(0, str(PACKED))
        import classix_t as module
        estimator = module.CLASSIX_T
        kind = "packed_little"
        module.print = lambda *a, **k: None
        module.tqdm = lambda iterable, *a, **k: iterable
        parameters = case["classix"]
    elif args.variant == "bitbirch":
        sys.path.insert(0, str(VENDOR))
        from bitbirch import bitbirch as module
        module.set_merge("diameter")
        estimator = module.BitBirch
        parameters = {k: case["bitbirch"][k] for k in ("threshold", "branching_factor")}
        kind = "int64"
        library = {"module": module.__file__, "merge_criterion": "diameter"}
    elif args.variant == "bblean":
        import bblean
        import bblean.similarity as similarity
        assert similarity._jt_sim_arr_vec_packed.__module__ == "bblean._cpp_similarity"
        estimator = bblean.BitBirch
        parameters = case["bitbirch"]
        kind = "packed_big"
        library = {"version": bblean.__version__, "native_kernel": similarity._jt_sim_arr_vec_packed.__module__,
                   "logical_features": case["features"], "native_features": 8*((case["features"]+7)//8),
                   "zero_padding_bits": (-case["features"]) % 8}
    elif args.variant == "dbscan":
        import sklearn
        from sklearn.cluster import DBSCAN
        estimator = DBSCAN
        parameters = case["dbscan"]
        kind = "int32"
        library = {"version": sklearn.__version__, "algorithm": "auto", "metric": "jaccard",
                   "input_note": "Notebook int32 input; sklearn performs its own boolean conversion inside fit"}
    elif args.variant == "butina":
        import rdkit
        from rdkit import DataStructs
        from rdkit.ML.Cluster import Butina
        parameters = case["butina"]
        kind = "packed_little"
        library = {"version": rdkit.__version__, "implementation": plan["butina_implementation"],
                   "module": Butina.__file__, "logical_features": case["features"],
                   "native_features": 8*((case["features"]+7)//8),
                   "zero_padding_bits": (-case["features"]) % 8}
    else:
        raise ValueError(args.variant)
    with threadpool_limits(limits=1):
        pool_info = threadpool_info()
        gc.collect()
        after_import = snapshot()
        data = np.load(directory/manifest["files"][kind]["file"], allow_pickle=False)
        if args.rows is not None:
            data = data[:args.rows].copy()
        representation_started = perf_counter()
        fingerprints = None
        if args.variant == "butina":
            fingerprints = [DataStructs.CreateFromBinaryText(row.tobytes()) for row in data]
        representation_seconds = perf_counter() - representation_started
        gc.collect()
        before = snapshot()
        send("ready", pid=os.getpid(), variant=args.variant, case=args.case, rows=len(data), features=case["features"],
             input_kind=kind, input_bytes=data.nbytes, parameters=parameters, library=library,
             representation_seconds=representation_seconds,
             threadpools=pool_info, after_import=after_import, before_fit=before)
        handshake("GO")
        started = perf_counter()
        try:
            model = None if args.variant == "butina" else estimator(**parameters)
            if args.variant == "butina":
                distances = []
                for i in range(1, len(fingerprints)):
                    similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[:i])
                    distances.extend([1 - value for value in similarities])
                clusters = Butina.ClusterData(distances, len(fingerprints), parameters["cutoff"],
                                               isDistData=True, reordering=parameters["reordering"])
                fit_end = perf_counter()
                labels = np.full(len(data), -1, dtype=np.int64)
                for cluster_id, members in enumerate(clusters):
                    labels[list(members)] = cluster_id
            elif args.variant == "dbscan":
                model.fit(data)
                fit_end = perf_counter()
                labels = model.labels_
            elif args.variant == "bblean":
                # Native 0.10.3 unpacking requires a byte-aligned dimension.
                # Existing packbits padding is zero and preserves similarity.
                model.fit(data, input_is_packed=True, n_features=8*data.shape[1])
                fit_end = perf_counter()
                labels = model.get_assignments()
            elif args.variant == "bitbirch":
                model.fit(data)
                fit_end = perf_counter()
                labels = model.get_assignments(len(data))
            else:
                model.fit(data)
                fit_end = perf_counter()
                labels = model.labels
            ended = perf_counter()
        except Exception as error:
            send("fit_error", error_type=type(error).__name__, message=str(error),
                 elapsed_seconds=perf_counter()-started, memory=snapshot(), traceback=traceback.format_exc())
            return
        after = snapshot()
        send("fit_done", clustering_seconds=ended-started, construction_and_fit_seconds=fit_end-started,
             assignment_seconds=ended-fit_end, after_fit=after)
        handshake("VALIDATE")
        try:
            from quality import adjusted_rand_index
            from state_checks import array_hash, logical_hash, state_fingerprint
            labels = np.asarray(labels, dtype=np.int64)
            assert labels.shape == (len(data),) and np.all(labels >= (-1 if args.variant == "dbscan" else 0))
            noise_points = int(np.sum(labels == -1))
            ids, counts = np.unique(labels[labels >= 0], return_counts=True)
            truth = np.load(directory/manifest["files"]["truth"]["file"], allow_pickle=False)[:len(data)]
            ari = adjusted_rand_index(truth, labels)
            mapping = {}
            canonical = np.empty_like(labels)
            for i, value in enumerate(labels):
                canonical[i] = mapping.setdefault(int(value), len(mapping))
            args.output.parent.mkdir(parents=True, exist_ok=True)
            np.save(args.output, labels)
            bitorder = "little" if kind == "packed_little" else "big" if kind == "packed_big" else None
            logical = logical_hash(data, case["features"], bitorder)
            if args.rows is None:
                assert logical == manifest["files"]["int32"]["sha256"], "Input changed"
            state_file = None
            if args.variant.startswith("classix"):
                state_file = args.output.with_suffix(".state.json")
                state_file.write_text(json.dumps(state_fingerprint(model, case["features"], args.variant == "classix_packed"), indent=2)+"\n")
            send("validated", ari=ari, clusters=len(ids), noise_points=noise_points,
                 singletons=int(np.sum(counts == 1)), largest_cluster=int(counts.max()) if len(counts) else 0,
                 aggregation_groups=len(model.splist) if args.variant.startswith("classix") else None,
                 adjacency_bytes=model.Adj.nbytes if args.variant.startswith("classix") else None,
                 adjacency_storage="connection bits + sparse minPts markers" if args.variant.startswith("classix") else None,
                 minpts_marker_count=model.Adj.marker_count if args.variant.startswith("classix") else None,
                 minpts_marker_storage_bytes=model.Adj.marker_storage_bytes if args.variant.startswith("classix") else None,
                 state_file=str(state_file) if state_file else None,
                 canonical_partition_sha256=array_hash(canonical), logical_input_sha256=logical,
                 input_unchanged=True if args.rows is None else None, assignments_file=str(args.output))
        except Exception as error:
            send("validation_error", error_type=type(error).__name__, message=str(error), traceback=traceback.format_exc())


if __name__ == "__main__":
    main()
