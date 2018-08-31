# Functional Quantile Regression

A method to perform Bayesian functional quantile regression. The proposed method is applicable to data sets given the _N_ by _T_ matrix of response functions _Y_ and _N_ by _p_ design matrix _X_ for functional quantile regression. 

This repository provides 
- a set of MATLAB scripts to implement our proposed Bayesian Functional Quantile Regression model and the naïve Bayesian Quantile Regression; 
- R scripts to implement the bootstrap-based approaches which were compared to our proposed model in the paper; 
- the code to reproduce all the figures in the manuscript and the supplementary
- the code to adjust for block effects from the preprocessed mass spectrometry data and to generate the simulation datasets

File structure: 
- The subfolder “BayesianFQR/” contains a set of MATLAB scripts to implement the Bayesian FQR model with discrete wavelet transform (DWT) on regression coefficient functions and a horseshoe prior on wavelet coefficients. These scripts are called when performing FQR on the real or simulation datasets.
- The subfolder “BayesianQR/” contains a set of MATLAB scripts to implement the naïve Bayesian Quantile Regression. These scripts are called when running this naïve approach on the real or simulation datasets.
- Other scripts or files not described above are either auxiliary functions used for plot or ROC analysis, or various output from running the scripts above.

## Example
Please refer to ```readme_data_and_code_description.docx``` for 
- instructions of how to reproduce results in the paper including both simulations and real data application 
- descriptions of the real data 

## Reproducibility 

All the figures and tables in the paper can be reproduced by running the provided code. Note that to reproduce Table 2 that summarizes the simulation performance for various methods, the posterior (or bootstrap) samples obtained by each approach at each considered quantile level based on each of 100 replicate datasets in each simulation setting are needed. Please refer to ```Reproducibility.md``` for more details about how to fully reproduce results in the paper. 

## Reference

Yusha Liu, Meng Li and Jeffrey S. Morris (2018). Function-on-Scalar Quantile Regression with Application to Mass Spectrometry Proteomics Data. 
