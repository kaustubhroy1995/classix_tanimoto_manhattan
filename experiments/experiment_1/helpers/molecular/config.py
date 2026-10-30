import os

METHODS = ["classix_packed", "dbscan", "butina", "bitbirch", "bblean"]

def get_design():
    settings = {"mode": os.environ.get("CLASSIX_RUN_MODE", "full"), "pilot_rows": int(os.environ.get("CLASSIX_PILOT_ROWS", "256"))}
    pilot = settings['mode'] == 'pilot'
    case = {
        'case_id': 'landrum_fixed', 'dataset': 'landrum',
        'data_id': 'landrum_pki_gt4_morgan1024',
        'rows': settings['pilot_rows'] if pilot else 90498, 'features': 1024,
        'study': 'fixed_configuration', 'seed': 101, 'methods': list(METHODS),
        'classix': {'radius': 0.3, 'minPts': 100, 'mergeScale': 1.33},
        'dbscan': {'eps': 0.3, 'min_samples': 50, 'metric': 'jaccard', 'algorithm': 'auto', 'n_jobs': 1},
        'butina': {'cutoff': 0.3, 'reordering': False},
        'bitbirch': {'threshold': 9/44, 'branching_factor': 50, 'merge_criterion': 'diameter'},
        'method_repeat_ids': {method: [0] for method in METHODS},
    }
    return {'experiment': '1_landrum', 'stage': 'landrum_fixed', 'seed': 101,
            'methods': list(METHODS), 'cases': [case], 'diagnostic_pilot': pilot,
            'butina_implementation': 'streaming_exact_static_rdkit',
            'dbscan_implementation': 'sklearn_precomputed_exact_packed_radius_graph',
            'fit_timeout_seconds': None, 'job_memory_limit_bytes': 20*1024**3,
            'design': {'scheduled_fits': 5, 'warmups': 0, 'measured_repeats': 1,
                       'parameter_search': False, 'source_order': 'No shuffle or deduplication',
                       'external_scores': 'ARI/AMI on all rows, noise one category',
                       'internal_population': 'All nonnoise clusters, including singletons; no post-size filter',
                       'silhouette': 'Seed 101, at most 2000 retained queries, all retained references, singleton score 0',
                       'runtime': 'Constructor + complete fit + complete labels; preparation and scoring excluded'}}
