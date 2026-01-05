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
3. Experiments with chemdb dataset with CLASSIX_T - classix_t_main.ipynb
