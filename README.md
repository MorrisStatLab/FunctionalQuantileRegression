# Bayesian Function-on-Scalar Quantile Regression

A method to perform Bayesian function-on-scalar quantile regression, i.e., Bayesian FQR. The proposed method is applicable to data sets given the _N_ by _T_ matrix of response functions _Y_ and _N_ by _p_ design matrix _X_ for function-on-scalar quantile regression. 

This repository provides 
- a set of MATLAB scripts to implement our proposed Bayesian FQR model (adjusted or unadjusted), and the naïve pointwise Bayesian quantile regression; 
- R scripts to implement the bootstrap-based approaches which were compared to our proposed model in the paper; 
- the code to reproduce all the figures in the paper;
- the code to adjust for block effects from the preprocessed mass spectrometry data and to generate the simulation datasets.

File structure: 
- The subfolder "BayesianFQR/" contains a set of MATLAB scripts to implement the Bayesian FQR model with discrete wavelet transform (DWT) on regression coefficient functions and a horseshoe prior on wavelet coefficients.
- The subfolder "BayesianFQR_corrected_likelihood/" contains a set of MATLAB scripts to implement the adjusted version of the Bayesian FQR model using sandwich likelihood correction (see Section 2.4 in the paper) with discrete wavelet transform (DWT) on regression coefficient functions and a horseshoe prior on wavelet coefficients.
- The subfolder "BayesianQR/" contains a set of MATLAB scripts to implement the naïve pointwise Bayesian quantile regression using asymmetric Laplace likelihood.
- The subfolder "Code for dwt and idwt/" contains the MATLAB scripts to perform DWT and inverse DWT.
- In the subfolder "realdata/", the file "ProteomicsData.mat" stores the raw and preprocessed pancreatic cancer mass spectrometry data, and the file "model.mat" stores the block effect adjusted mass spectra _Y_, design matrix _X_, spectral locations _x0_ (in Daltons) and DWT specs. This subfolder also includes the MATLAB code to adjust for block effects, and scripts to implement our proposed Bayesian FQR and various alternatives on the mass spectrometry dataset.
- The subfolders "simulations_continuous_x/" contain the scripts to generate simulation datasets and to implement our proposed model and alternative methods.
- Other subfolders not described above contain auxiliary functions that need to be used when running the scripts above.

## Dependency 
- For our proposed method: MATLAB version 2016b or later (wavelet toolbox needed).
- For selected other methods: R version 3.2.2 or later (the package “quantreg” needed).

## Example
Below is an example to run Bayesian FQR on one simulation dataset.

```

% Load in the dataset    
% Note that the input data must be a structure array that at least includes the _N_ by _T_ data matrix _Y_, the _N_ by _p_ design matrix _X_, the sample size _N_, the grid size _T_, the number of regressors _p_, and a structure array _wavespecs_ that includes the DWT and iDWT information.      
% Before running Bayesian FQR, please make sure your input data have the same structure as the following example, and refer to "simulations_continuous_x/simdata2.m" to see how to generate such a structure array in MATLAB.         
load('simulations_continuous_x/data/model1.mat');  

% Add the path of the functions needed to run Bayesian FQR     
addpath('BayesianFQR');

% Specify MCMCspecs, the quantile level desired, and the random seed for reproducibility.        
MCMCspecs.minVC=1e-6;     
MCMCspecs.maxO=1e20;     
MCMCspecs.B=2000;     
MCMCspecs.burnin=2000;     
MCMCspecs.thin=3;     
MCMCspecs.blocksize=2000;     

qt=0.9;       
seed=100;     
  
% Call Bayesian FQR and get posterior samples     
% It takes about 15 minutes to run on a 64-bit operating system with 2 processors and 32GB RAM.       
result=FQR_HS(model, qt, MCMCspecs, seed);     
MCMC_betat=result.MCMC_betat;   

% MCMC_P is a _B_ by _T_ matrix containing _B_ posterior samples for the regression coefficient function coding the group difference.        
MCMC_P=MCMC_betat(:,(model.T+1):end);    

```

## Reproducibility 

All the figures in the paper can be reproduced by running the provided code. Note that to reproduce Table 2 that summarizes the simulation performance for various methods, the posterior (or bootstrap) samples obtained by each approach at each considered quantile level based on the 100 replicate datasets are needed. Please refer to ```Reproducibility.md``` for more details about how to fully reproduce results in the paper. 

## Reference

Yusha Liu, Meng Li and Jeffrey S. Morris (2019). Function-on-Scalar Quantile Regression with Application to Mass Spectrometry Proteomics Data. <https://arxiv.org/abs/1809.00266>
