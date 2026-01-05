# CLASSIX Tanimoto and Manhattan
Code in support of the article "Fast and explainable clustering in the Manhattan and Tanimoto distance"

# Installation
CLASSIX Tanimoto and Manhattan can be installed by cloning this repository and following creating an object of the CLASSIX_T or CLASSIX_M classes.

Optional: Perfomance metrics in CLASSIX Tanimoto require our custom sparse matrix - vector product subroutine, which needs to be installed using the appropriate wheel file for Windows, Mac and Linux. This can be done by using the command - pip install <wheel_file_name>. The extracted package will not contain the source code; only a pre-compiled shared library is provided. This does not require C/C++ compilers to be installed to run.

# Experiments
Our code for experiments can be found in the experiments section.

The following files are used for the following experiments:
1. Synthetic data generation - classix_t_blobs.ipynb
2. Experiments with synthetic data for CLASSIX_T - classix_t_blobs.ipynb
3. CLASSIX_T probabilistic analysis simulations - classix_t_prob.ipynb
4. Timing with chemdb dataset for taylor-butina, DBSCAN and CLASSIX_T - classix_t_main.ipynb
5. CLASSIX_M experiments with Iris dataset - classix_m_iris.ipynb
6. CLASSIX_M experiments with Banknote dataset - classix_m_banknote.ipynb
7. CLASSIX_M experiments with MNIST dataset - classix_m_mnist.ipynb

# Hyperparameter choices
The hyperparameter ranges for  CLASSIX\_M were $\texttt{radius} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{radius} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{radius} \in (0.01, 0.1)$ in steps of $0.005$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for the MNIST dataset. 

For DBSCAN $\texttt{eps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{eps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{eps} \in (0.1, 1)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for the MNIST dataset. 

For OPTICS $\texttt{maxeps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for IRIS; $\texttt{maxeps} \in (0.1, 0.5)$ in steps of $0.025$, $\texttt{minPts} \in (0, 50)$ in steps of $5$ for Banknote; and $\texttt{maxeps} \in (0.1, 3)$ in steps of $0.025$, $\texttt{minsamples} \in (0, 50)$ in steps of $5$ for the MNIST dataset.
