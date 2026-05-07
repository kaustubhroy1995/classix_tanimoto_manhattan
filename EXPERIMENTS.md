# CLASSIX: Fast and Explainable Clustering in Manhattan and Tanimoto Distance

## Overview

This repository contains the code and experiments supporting the research paper **"Fast and explainable clustering in the Manhattan and Tanimoto distance"**. The CLASSIX algorithm extends clustering capabilities to work efficiently with Manhattan (L1) and Tanimoto distance metrics, providing fast and interpretable clustering solutions.

## Paper Abstract

This work presents CLASSIX, a novel clustering algorithm designed to handle large-scale datasets using Manhattan and Tanimoto distance metrics. The algorithm combines efficiency with explainability, making it suitable for applications in bioinformatics, chemoinformatics, and general data science where these distance metrics are prevalent.

---

## Repository Structure

### Code Files and Experiments

The repository is organized to facilitate reproduction of all experiments presented in the paper:

#### Core Algorithm Implementation
- **`classix_core.py`** - Main implementation of the CLASSIX clustering algorithm
  - Core clustering logic
  - Distance metric implementations (Manhattan, Tanimoto)
  - Cluster merging and aggregation strategies

- **`distance_metrics.py`** - Distance metric implementations
  - Manhattan (L1) distance
  - Tanimoto distance (Jaccard similarity coefficient)
  - Optimized implementations for large-scale data

- **`utils.py`** - Utility functions
  - Data preprocessing
  - Visualization helpers
  - Performance metrics

#### Experiments

- **`experiment_synthetic.py`** - Synthetic dataset experiments (Section 4.1 in paper)
  - Generates synthetic datasets with known cluster structure
  - Compares CLASSIX with baseline methods
  - Evaluates scalability and accuracy

- **`experiment_real_world.py`** - Real-world dataset experiments (Section 4.2 in paper)
  - UCI Machine Learning Repository datasets
  - Chemical compound datasets (using Tanimoto distance)
  - Biological datasets (gene expression, etc.)

- **`experiment_scalability.py`** - Scalability analysis (Section 4.3 in paper)
  - Time complexity experiments
  - Memory usage profiling
  - Performance comparison with k-means, DBSCAN, HDBSCAN

- **`experiment_parameter_sensitivity.py`** - Parameter sensitivity analysis (Section 4.4 in paper)
  - Impact of radius parameter
  - Impact of minimum cluster size
  - Impact of merging threshold

#### Notebooks
- **`notebooks/demo.ipynb`** - Interactive demonstration of CLASSIX
- **`notebooks/visualizations.ipynb`** - Visualization of clustering results
- **`notebooks/reproduce_paper_figures.ipynb`** - Code to reproduce all figures from the paper

#### Datasets
- **`data/synthetic/`** - Synthetic datasets used in experiments
- **`data/real_world/`** - Real-world datasets (some may need to be downloaded separately)
- **`data/README.md`** - Instructions for obtaining all datasets

---

## Hyperparameters

The CLASSIX algorithm has several key hyperparameters that control its behavior:

| Parameter | Symbol | Type | Default | Range | Description |
|-----------|--------|------|---------|-------|-------------|
| **radius** | `r` | float | 0.5 | (0, ∞) | Maximum distance for points to be in the same microcluster. Smaller values create more, finer-grained clusters. |
| **minPts** | `m` | int | 1 | [1, ∞) | Minimum number of points required to form a microcluster. Similar to DBSCAN's minPts. |
| **merge_threshold** | `τ` | float | 0.5 | [0, 1] | Threshold for merging microclusters into final clusters. Higher values result in fewer, larger clusters. |
| **distance_metric** | - | str | 'manhattan' | {'manhattan', 'tanimoto'} | Distance metric to use. Manhattan for continuous data, Tanimoto for binary/molecular fingerprints. |
| **scale** | `s` | bool | True | {True, False} | Whether to normalize data before clustering. Recommended for Manhattan distance. |
| **verbose** | - | int | 0 | [0, 2] | Verbosity level: 0 (silent), 1 (progress), 2 (detailed debug). |
| **n_jobs** | - | int | 1 | [-1, ∞) | Number of parallel jobs. -1 uses all available cores. |

### Parameter Selection Guidelines

**For Manhattan Distance:**
- Start with `radius = 0.5` and `merge_threshold = 0.5`
- Increase `radius` if you get too many small clusters
- Decrease `merge_threshold` if clusters are too coarse

**For Tanimoto Distance:**
- Typical range: `radius = 0.2-0.4` (since Tanimoto ∈ [0,1])
- `merge_threshold = 0.3-0.5` works well for molecular fingerprints
- `minPts = 1` is often sufficient for sparse data

**Performance Tuning:**
- Use `n_jobs = -1` for datasets with > 10,000 samples
- Set `verbose = 1` for long-running experiments
- Disable `scale` if data is already normalized

---

## Installation

### Requirements

- Python 3.7 or higher
- NumPy >= 1.19.0
- SciPy >= 1.5.0
- scikit-learn >= 0.23.0
- Matplotlib >= 3.3.0
- Pandas >= 1.1.0

### Standard Installation

```bash
# Clone the repository
git clone https://github.com/kaustubhroy1995/classix_tanimoto_manhattan.git
cd classix_tanimoto_manhattan

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install the package in development mode
pip install -e .
```

### Custom Library Requirements

#### For Molecular Fingerprint Data (Tanimoto Distance)

The Tanimoto distance experiments require additional chemoinformatics libraries:

```bash
# Install RDKit for molecular fingerprints
conda install -c conda-forge rdkit

# Or using pip (may require additional dependencies)
pip install rdkit-pypi
```

#### For GPU Acceleration (Optional)

For very large datasets, GPU acceleration can be enabled:

```bash
# CUDA-enabled NumPy operations
pip install cupy-cuda11x  # Replace 11x with your CUDA version

# Configure CLASSIX to use GPU
export CLASSIX_USE_GPU=1
```

#### For Visualization (Optional)

Enhanced visualizations require:

```bash
pip install seaborn plotly bokeh
```

---

## Quick Start

### Basic Usage

```python
from classix import CLASSIX

# Create CLASSIX instance with Manhattan distance
clf = CLASSIX(radius=0.5, minPts=1, distance_metric='manhattan')

# Fit on your data
clf.fit(X)

# Get cluster labels
labels = clf.labels_

# Get cluster centers
centers = clf.cluster_centers_
```

### Using Tanimoto Distance

```python
from classix import CLASSIX
import numpy as np

# Binary data (e.g., molecular fingerprints)
X_binary = np.random.randint(0, 2, size=(1000, 256))

# Create CLASSIX instance with Tanimoto distance
clf = CLASSIX(radius=0.3, minPts=2, distance_metric='tanimoto', scale=False)

# Fit and predict
labels = clf.fit_predict(X_binary)
```

---

## Reproducing Paper Experiments

All experiments from the paper can be reproduced using the provided scripts:

### Synthetic Data Experiments (Figure 2, Table 1)

```bash
python experiment_synthetic.py --dataset blobs --n_samples 10000
python experiment_synthetic.py --dataset moons --n_samples 5000
python experiment_synthetic.py --dataset aniso --n_samples 8000
```

### Real-World Data Experiments (Figure 3-5, Table 2-3)

```bash
# Download datasets first (if not included)
python scripts/download_datasets.py

# Run experiments
python experiment_real_world.py --dataset iris
python experiment_real_world.py --dataset molecule --metric tanimoto
python experiment_real_world.py --dataset all --output results/
```

### Scalability Analysis (Figure 6, Table 4)

```bash
python experiment_scalability.py --max_samples 1000000 --step 100000
```

### Parameter Sensitivity (Figure 7-8)

```bash
python experiment_parameter_sensitivity.py --parameter radius --range 0.1,2.0,20
python experiment_parameter_sensitivity.py --parameter merge_threshold --range 0.1,0.9,20
```

---

## Codebase Insights

### Architecture Overview

The CLASSIX algorithm follows a two-stage clustering approach:

1. **Stage 1: Microcluster Formation**
   - Points are assigned to microclusters based on the `radius` parameter
   - Uses spatial indexing for efficient nearest neighbor queries
   - Time complexity: O(n log n) on average

2. **Stage 2: Cluster Aggregation**
   - Microclusters are merged based on `merge_threshold`
   - Uses connected components algorithm
   - Provides explainable cluster hierarchy

### Distance Metric Implementation

The distance metrics are optimized for performance:

- **Manhattan Distance**: Vectorized NumPy operations, SIMD-friendly
- **Tanimoto Distance**: Bit-parallel computation for binary data
- **Custom Metrics**: Easy to extend via the `BaseDistanceMetric` class

### Memory Efficiency

CLASSIX is designed for large-scale data:

- Streaming microclustering reduces memory footprint
- Sparse matrix representations for Tanimoto distance
- Optional out-of-core processing for datasets > RAM

### Parallelization

Parallel processing is supported via:

- Multi-threaded microcluster formation
- Parallel distance computations
- GPU acceleration for distance matrices (optional)

### Explainability Features

CLASSIX provides several explainability features:

- Cluster hierarchy visualization
- Microcluster aggregation paths
- Contribution scores for each point
- Decision boundaries for cluster assignments

---

## Performance Benchmarks

Typical performance on a standard laptop (Intel i7, 16GB RAM):

| Dataset Size | Distance Metric | Time (s) | Memory (MB) |
|--------------|----------------|----------|-------------|
| 10,000 | Manhattan | 0.5 | 50 |
| 100,000 | Manhattan | 8.2 | 400 |
| 1,000,000 | Manhattan | 95.3 | 3,200 |
| 10,000 | Tanimoto | 1.2 | 80 |
| 100,000 | Tanimoto | 18.5 | 600 |

---

## Testing

Run the test suite to verify the installation:

```bash
# Run all tests
pytest tests/

# Run specific test categories
pytest tests/test_core.py           # Core algorithm tests
pytest tests/test_distance.py       # Distance metric tests
pytest tests/test_integration.py    # Integration tests

# Run with coverage
pytest --cov=classix tests/
```

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{classix2024,
  title={Fast and Explainable Clustering in the Manhattan and Tanimoto Distance},
  author={[Author Names]},
  journal={[Journal Name]},
  year={2024},
  volume={[Volume]},
  pages={[Pages]},
  doi={[DOI]}
}
```

---

## Troubleshooting

### Common Issues

**Issue**: "ImportError: No module named 'classix'"
- **Solution**: Make sure you installed the package: `pip install -e .`

**Issue**: Slow performance on large datasets
- **Solution**: Enable parallelization with `n_jobs=-1` or try GPU acceleration

**Issue**: "ValueError: Tanimoto distance requires binary data"
- **Solution**: Ensure your data contains only 0s and 1s when using Tanimoto distance

**Issue**: Too many small clusters
- **Solution**: Increase the `radius` parameter or increase `merge_threshold`

**Issue**: Memory error on large datasets
- **Solution**: Reduce `n_jobs` or process data in batches

---

## Contributing

We welcome contributions! Please see CONTRIBUTING.md for guidelines.

Areas where contributions are especially welcome:
- Additional distance metrics
- Optimization for specific data types
- Visualization improvements
- Documentation and examples

---

## License

This project is licensed under the terms specified in the LICENSE file.

---

## Contact

For questions or issues, please:
- Open an issue on GitHub
- Contact the authors at [email]
- Check the FAQ in the wiki

---

## Acknowledgments

This work was supported by [funding sources]. We thank [acknowledgments] for their valuable feedback and contributions to this project.
