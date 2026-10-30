# Experimental settings used for the reported results.

CLASSIX_CONFIG = {'experiment': '2',
 'frozen_before_new_scores_utc': '2026-09-24T17:31:14Z',
 'seed': 101,
 'datasets': ['iris',
              'banknote',
              'wine',
              'breast_cancer',
              'digits',
              'seeds',
              'glass',
              'sonar',
              'vehicle',
              'segmentation',
              'breast_tissue',
              'leaf',
              'libras',
              'vertebral',
              'wireless',
              'parkinsons',
              'coimbra',
              'raisin',
              'haberman',
              'musk'],
 'anchors': ['iris', 'banknote', 'wine', 'breast_cancer', 'digits', 'seeds', 'glass', 'sonar'],
 'preprocessing': 'All original numeric features after shared StandardScaler; identifier columns '
                  'excluded; no UMAP or binarization; dataset-specific exclusions and repeated '
                  'observations disclosed.',
 'radius_count': 32,
 'min_pts': [0, 5, 10, 20, 40, 80],
 'merge_scale': 1.4,
 'radius_low': '0.35 * median fifth-neighbour distance including duplicates, excluding self',
 'radius_high': '0.9 * median strictly-positive pairwise distance',
 'calibration': "Exact distances in each method's metric; geometric physical radii divided by its "
                'internal scalar normalization. Reject zero/nonfinite lower bound and disclose '
                'failure rather than silently changing data or calibration.',
 'extension': 'Once only: if either initial winner is at radius index 0 or 31, add indices -4..-1 '
              'and 32..35 for BOTH; if either initial winner has minPts80, add120,160 for BOTH.',
 'selection': 'Maximize full-precision ARI; ties smaller minPts then radius index. AMI evaluated '
              'at the same setting, never optimized separately.',
 'permutations': 5,
 'permutation_seed': 'SeedSequence([101, dataset index in this fixed roster, replicate 0..4, '
                     '867]); identical row order for both methods; selected parameters fixed, '
                     'labels aligned to original order.',
 'selected_timing_repetitions': 5,
 'panel_selection': {'target_size': 10,
                     'mandatory': 'All eight existing datasets, with no deletions based on '
                                  'outcomes.',
                     'new_manhattan_case': 'First new dataset in alphabetical ID order where the '
                                           'original selected result AND all five paired '
                                           'permutations have strictly positive '
                                           'Manhattan-minus-Euclidean ARI and AMI.',
                     'new_euclidean_case': 'First new dataset in alphabetical ID order where the '
                                           'original selected result AND all five paired '
                                           'permutations have strictly negative '
                                           'Manhattan-minus-Euclidean ARI and AMI.',
                     'fallback': 'Fill missing slots first from alphabetically ordered unselected '
                                 'new mixed/tied/order-sensitive complete cases, then from other '
                                 'alphabetically ordered complete new cases. Never force a missing '
                                 'direction or conceal exclusions.',
                     'interpretation': 'Explicitly outcome-stratified illustrative panel; all20 '
                                       'attempted datasets and descriptive screen counts are shown '
                                       'alongside it. Five permutations assess ordering only, not '
                                       'sampling generalization or statistical significance.'},
 'prior_evidence': 'Standalone reproduction: download and prepare all data and refit all settings; '
                   'no previous run records are imported.',
 'earlier_umap': 'Earlier UMAP and synthetic diagnostics are separate investigations, outside this '
                 'approved original-feature table.',
 'resource_policy': 'Serial fits, BLAS threads1, public small datasets, small unscored readiness '
                    'fits before new full grids. Preserve checkpoints and any failures.',
 'origin': 'Approved real-feature screen previously developed as Experiment 2B.',
 'approved_panel': ['iris',
                    'banknote',
                    'wine',
                    'breast_cancer',
                    'digits',
                    'seeds',
                    'glass',
                    'sonar',
                    'breast_tissue',
                    'haberman'],
 'panel_presentation': 'Keep the author-approved ten cases; omit the permutation column from the '
                       'main table. Retain all twenty datasets and permutation evidence in '
                       'generated supporting outputs.'}

DENSITY_CONFIG = {'experiment': '2C: sklearn density baselines on the original Experiment 2 full-feature screen',
 'frozen_utc': '2026-09-28T16:30:43.939731+00:00',
 'datasets': ['iris',
              'banknote',
              'wine',
              'breast_cancer',
              'digits',
              'seeds',
              'glass',
              'sonar',
              'vehicle',
              'segmentation',
              'breast_tissue',
              'leaf',
              'libras',
              'vertebral',
              'wireless',
              'parkinsons',
              'coimbra',
              'raisin',
              'haberman',
              'musk'],
 'panel': ['iris',
           'banknote',
           'wine',
           'breast_cancer',
           'digits',
           'seeds',
           'glass',
           'sonar',
           'breast_tissue',
           'haberman'],
 'seed': 101,
 'preprocessing': 'All original numeric features after shared StandardScaler; identifier columns '
                  'excluded; no UMAP or binarization; dataset-specific exclusions and repeated '
                  'observations disclosed.',
 'metrics': ['euclidean', 'manhattan'],
 'min_samples': [2, 5, 10, 20, 40, 80],
 'candidate_count': 32,
 'dbscan': {'eps_low': '0.35 * median fifth-other-row neighbour distance',
            'eps_high': '0.9 * median strictly positive pairwise distance',
            'spacing': 'geometric',
            'algorithm': 'auto',
            'n_jobs': 1},
 'optics': {'xi_low': 0.01,
            'xi_high': 0.3,
            'spacing': 'geometric',
            'max_eps': 'infinity',
            'cluster_method': 'xi',
            'min_cluster_size': None,
            'predecessor_correction': True,
            'algorithm': 'auto',
            'n_jobs': 1,
            'memory': None},
 'extension': 'Once per algorithm/dataset, shared by both metrics: if either initial winner has '
              'index0 or31, add indices -4..-1 and32..35 for both; if either winner uses '
              'min_samples80, add120 and160 not exceeding n for both. Reuse original radius/xi '
              'ratio. No further extension.',
 'selection': 'Maximum full-precision all-row ARI; ties smaller min_samples then parameter index. '
              'AMI evaluates the same assignment, never selected separately. Noise -1 is one '
              'predicted category; preserve degenerate outcomes.',
 'optics_search_cache': 'One sklearn OPTICS fit/graph per dataset,metric,min_samples; reuse '
                        'cluster_optics_xi during search. Selected fitting/timing uses fresh full '
                        'OPTICS.fit on raw standardized X with memory=None.',
 'timing': 'One discarded warmup plus5 fresh model.fit repetitions, original standardized float64 '
           'X copied before timing; threadpool_limits1, n_jobs1. Include each method internal '
           'preprocessing/neighbourhood/ordering/merging/extraction in fit; exclude construction, '
           'input copy, external StandardScaler, scoring, calibration, and search. Retiming '
           'existing selected CLASSIX configs without retuning; verify original partitions. '
           'Interleave all6 variants by deterministic shuffled order each repetition.',
 'resource_policy': 'Serial small-data fits on this Windows machine; max n2310,d166; readiness '
                    'checks before full grids; resumable per candidate/graph/dataset, STOP file '
                    'checked between work units; no model parallelism.',
 'versions': {'python': '3.12.10',
              'numpy': '2.2.6',
              'scipy': '1.15.3',
              'scikit-learn': '1.6.1',
              'classixclustering': '1.5.1'}}
