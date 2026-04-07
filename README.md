# CLASSIX Tanimoto and Manhattan

Code in support of the article "Fast and explainable clustering in the Manhattan and Tanimoto distance".

## Run In A Fresh Environment

The project is designed to run directly from source in `src/`.

### 1. Clone and enter the repository

```bash
git clone https://github.com/kaustubhroy1995/classix_tanimoto_manhattan.git
cd classix_tanimoto_manhattan
```

### 2. Create and activate a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### 3. Install Python dependencies

```bash
python -m pip install -r requirement.txt
```

### 4. Build the local `spmv` shared library (required for `CLASSIX_T`)

`CLASSIX_T` uses a local C shared object loaded from `src/spmv/lib/spsubmatxvec.so`.

```bash
cd src/spmv
make
mkdir -p lib
cp spsubmatxvec.so lib/spsubmatxvec.so
cd ../..
```

If your system does not have `clang`, edit `src/spmv/Makefile` and replace `clang` with `gcc`.

### 5. Quick smoke test from source

#### CLASSIX_M (Manhattan)

```bash
python - <<'PY'
import sys
sys.path.insert(0, 'src')
import numpy as np
from classix_m import CLASSIX_M

x = np.array([[0.,1.],[1.,0.],[0.8,0.2],[0.2,0.8]])
model = CLASSIX_M(radius=0.5, minPts=1)
model.fit(x)
print('CLASSIX_M OK:', len(model.labels))
PY
```

#### CLASSIX_T (Tanimoto)

```bash
python - <<'PY'
import sys
sys.path.insert(0, 'src')
import numpy as np
from classix_t import CLASSIX_T

x = np.array([[1,0,1],[1,1,0],[0,1,1],[1,0,0]], dtype=np.int32)
model = CLASSIX_T(radius=0.3, minPts=1)
model.fit(x)
print('CLASSIX_T OK:', len(model.labels))
PY
```

## Experiments

The experiment notebooks are in `experiments/`:

1. Synthetic data generation and CLASSIX_T experiments: `classix_t_blobs.ipynb`
2. CLASSIX_T probabilistic analysis simulations: `classix_t_prob.ipynb`
3. Timing with chemdb dataset for Taylor-Butina, DBSCAN and CLASSIX_T: `classix_t_main.ipynb`
4. CLASSIX_M experiments with Iris: `classix_m_expt_iris.ipynb`
5. CLASSIX_M experiments with Banknote: `classix_m_expt_banknote.ipynb`
6. CLASSIX_M experiments with MNIST: `classix_m_expt_mnist.ipynb`

<!# Hyperparameter choices
The hyperparameter ranges for  CLASSIX\_M were $\texttt{radius} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{radius} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{radius} \in (0.01, 0.1)$ in steps of $0.005$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for the MNIST dataset. 

For DBSCAN $\texttt{eps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{eps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{eps} \in (0.1, 1)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for the MNIST dataset. 

For OPTICS $\texttt{maxeps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{maxeps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{maxeps} \in (0.1, 3)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for the MNIST dataset.>
