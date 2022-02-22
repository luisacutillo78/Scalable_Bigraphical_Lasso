# Nonparanormal_bigraphical_lasso code.
This repository reproduces the research presented in the paper: 
"Two-way Sparse Network Inference for Count Data", Authors: Li, Sijia; López-García, Martín; Lawrence, Neil D; Cutillo, Luisa, accepted at AISTATS2022.
## Demo on Real Data
In the file Demo1.m we provide an example of use of our code in case of count data. In Demo2.m we provide an example of use of our code in case of continuous data. 

Data used in Demo1 example are...

Data used in Demo2 example are...

To run the code:
* Clone the repo
* Launch Matlab
* Run Demo.m

## Matlab libraries needed: 
We are using  **Matlab 2020b**.
You will need to install the following matlab libraries: ndlutil, rca, glmnet_matlab, L1General. Copies of these are included in this repo.
You will also need the library of Teralasso (Greenewald et al. 2019) for some comparison studies.

We assume that you'll run these scripts from your current working directory.  Change the `mypath` variable if you do otherwise.

The experiment was run on Windows 10. Please let us know if you encounter any problems either on Windows or other operating system.

## For Figure 1
* run **Rebuttal_Figures.m**.
## For Figure 2
* run **run_syntheticGaussian_Three_method_time.m**.
## For Figure 3, Figure 9 and Figure 10
* run **run_synthetic_data.m**.
## For Figure 4
* run **run_synthetic_data_three_blocks.m**.
## For Figure 6
* run **main_duck_fast_new_colour_ST.m**.
## For Figure 7
* run **main_genes.m**.
## For Figure 8
* run **main.R**.
## For experiment on coil duck pictures
* run **main_duck_fast_new_colour_ST.m**.
