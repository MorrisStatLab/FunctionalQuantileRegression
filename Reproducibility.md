## Data 

### Description 

Real data: At the University of Texas M.D. Anderson Cancer Center, a study was conducted using matrix assisted laser desorption and ionization, time-of-flight mass spectrometry (MALDI-TOF) to discover potential proteomic markers of pancreatic cancer. In this study, researchers collected the blood serum samples from 139 pancreatic cancer patients and 117 normal controls and ran them on a MALDI-TOF mass spectrometer to produce a mass spectrum for each sample.  These data were preprocessed using the Cromwell pipeline (Coombes, et al. 2005 Proteomics), including baseline correction, normalization, denoising and log_2 transformation, as detailed in Koomen et al. (2005, Clinical Cancer Research).  In our analysis here, we focused on the spectral region [5000D, 8000D] of the original dataset, which includes 1659 observations per spectrum. The 256 samples in this dataset were measured in four different blocks over a span of several months, so we estimated and subtracted the block-specific mean from the preprocessed mass spectra to adjust for the block effects before performing functional quantile regression. The raw dataset, preprocessed dataset, and the block effect adjusted dataset are all provided.

Simulated data: Under the simulation setting described in Section 3 in the paper, we generated 100 replicate datasets. Each replicate dataset includes simulated functions consisting of 301 observations for each of 400 subjects. For space considerations, only 1 replicate simulated dataset is provided. However, all 100 replicates are available upon request, and can also be generated using the MATLAB script we provided.

### How to access data 

Real data: Available in the subfolder ```realdata```. The raw and preprocessed mass spectrometry data are stored as matrices in “.mat” file; the block effect adjusted data matrix _Y_ which is ready for model fitting is stored as a structure array in “.mat” file along with design matrix _X_, spectral locations _x_0_ (in Daltons) and wavelet transform specifications _wavespecs_. See the section below for how to fully reproduce this preprocesing step. 

Simulated data: Available in the subfolder ```simulations_continuous_x```. The simulation datasets are stored both as structure arrays in “.mat” files and matrices in “.txt” files. 

## Reproducibility

### What is to be reproduced    
All the figures in the paper can be reproduced by running the provided code. To reproduce Table 2 that summarizes the simulation performance for various methods, the posterior (or bootstrap) samples obtained by each approach at each considered quantile level based on the 100 replicate datasets are needed. To reproduce Figure S2-5 in the supplement, the posterior (or bootstrap) samples of the regression coefficients obtained by each approach at each quantile level based on the proteomics dataset are needed. For space considerations, we do not provide these posterior (or bootstrap) samples here.  However, they can be generated using the scripts we provided, and are available upon request. We also provide the script to implement the Bayesian FQR model and various alternatives on a given functional dataset and to generate posterior (or bootstrap) samples, and the code to do estimation and inference on functional coefficients using these posterior (or bootstrap) samples.

### How to reproduce analyses     
- Figures: The MATLAB script “figures.m” in the main folder contains the code to reproduce all the figures in the paper. 
- Real data: In the subfolder “realdata/”, the file “ProteomicsData.mat” stores the raw and preprocessed pancreatic cancer mass spectrometry data, the cancer/normal status of each of 256 subjects and the spectral locations of observations (in Daltons); the MATLAB script “realdata_preparation.m” records how to adjust for the block effect of the preprocessed data and generate the file “model.mat” which stores the adjusted mass spectra, design matrix, spectral locations and wavelet transform specifications as a structure array; the MATLAB script “realdata_modelfit.m” and R script “realdata_modelfit_bootstrap.R” record how to implement our proposed Bayesian FQR and various alternative approaches on the adjusted mass spectrometry dataset and to obtain posterior (or bootstrap) samples for the quantile regression coefficient functions. 
- Simulated data: In the subfolder “simulations_continuous_x/”, the MATLAB script “data_generation.m” records how to generate the 100 replicate simulation datasets, and to obtain the ground truth for each regression coefficient function at selected quantile levels; the MATLAB script “data_modelfit.m” and R script “data_modelfit_bootstrap.R” record how to implement our proposed Bayesian FQR and alternative approaches on the simulated datasets and to obtain posterior (or bootstrap) samples; the MATLAB script “beta2_estimation_inference_example.m” illustrates how to summarize the estimation and inference performance for the coefficient function _beta_2_ using posterior (or bootstrap) samples based on a simulated dataset, given the ground truth.


## References
[1] Coombes K R, Tsavachidis S, Morris J S, et al (2005). Improved peak detection and quantification of mass spectrometry data acquired from surface - enhanced laser desorption and ionization by denoising spectra with the undecimated discrete wavelet transform. _Proteomics_, 5(16): 4107-4117.

[2] Koomen J M, Shih L N, Coombes K R, et al (2005). Plasma protein profiling for diagnosis of pancreatic cancer reveals the presence of host response proteins. _Clinical Cancer Research_, 11(3): 1110-1118.

