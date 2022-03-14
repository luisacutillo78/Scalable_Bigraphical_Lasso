# Scalable Bigraphical Lasso code
This repository reproduces the research presented in the paper: 
"Two-way Sparse Network Inference for Count Data", Authors: Li, Sijia; López-García, Martín; Lawrence, Neil D; Cutillo, Luisa, accepted at AISTATS2022.

## Matlab libraries needed: 
We are using  **Matlab 2020b**.
You will need to install the following matlab libraries: ndlutil, rca, glmnet_matlab, L1General. Copies of these are included in this repo.
You will also need the library of Teralasso (Greenewald et al. 2019) for some comparison studies.

## Examples of use
Example on synthetic data: run_synthetic_data.m

Example on real data: main_genes.m

## System requirements
The experiment was run on Intel(R) Core(TM) i5-9500, Windows 10. Please let us know if you encounter any problems either on Windows or other operating system.

## For Figure 1
* run **Figure1.m**.
## For Figure 2
* run **run_syntheticGaussian_Three_method_time.m**.
## For Figure 3, Figure 9 and Figure 10
* run **run_synthetic_data.m**.
## For Figure 4
* run **run_synthetic_data_three_blocks.m**.
## For Figure 5
* run **main_genes.m**.
## For Figure 6
* run **main.R**.


