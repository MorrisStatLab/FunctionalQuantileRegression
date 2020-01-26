clear;

% the quantile levels
qt=[0.1 0.25 0.5 0.75 0.9];

% load in data
load('Y:/submission/realdata/model.mat');

% MCMC specifications 
MCMCspecs.minVC=1e-6;
MCMCspecs.maxO=1e20;

MCMCspecs.B=2000;
MCMCspecs.burnin=5000;
MCMCspecs.thin=5;
MCMCspecs.blocksize=2000;

MCMCspecs.time_update=1000;

%% Bayesian FQR

% add the path of the functions needed to run Bayesian FQR
addpath('Y:/submission/BayesianFQR');

% for space considerations, we only provided the posterior samples for 
% betat at qt=0.9, which is used to plot Figure 4 in the paper. However,
% Users can run the code below to generate the posterior samples for all 
% model parameters at each quantile level. Posterior samples are also 
% available upon request.

% for each quantile level, it takes about 4.5 hours to fit the model on 
% a 64-bit operating system with 2 processors and 32GB RAM

for i=1:5
    result=FQR_HS(model,qt(i),MCMCspecs,100*i);
    MCMC_betat=result.MCMC_betat;
    MCMC_beta=result.MCMC_beta;
    MCMC_lambda2=result.MCMC_lambda2;
    MCMC_tau2=result.MCMC_tau2;
    MCMC_sigma=result.MCMC_sigma;
    MCMC_betat=MCMC_betat(:,(model.T+1):end);
    MCMC_beta=MCMC_beta(:,(model.wavespecs.K+1):end);
    MCMC_lambda2=MCMC_lambda2(:,(model.wavespecs.K+1):end);
    MCMC_tau2=MCMC_tau2(:,(model.wavespecs.J+1):end);
    dlmwrite(sprintf('Y:/submission/realdata/output/HS/MCMC_betat_%d.txt',qt(i)*100),MCMC_betat,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS/MCMC_beta_%d.txt',qt(i)*100),MCMC_beta,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS/MCMC_lambda2_%d.txt',qt(i)*100),MCMC_lambda2,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS/MCMC_tau2_%d.txt',qt(i)*100),MCMC_tau2,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS/MCMC_sigma_%d.txt',qt(i)*100),MCMC_sigma,'delimiter','\t','precision','%12.6e');
end

rmpath('Y:/submission/BayesianFQR');


%% naive Bayesian QR

% add the path of the functions needed to run naive Bayesian FQR
addpath('Y:/submission/BayesianQR');

% for space considerations, we do not provide the posterior samples for 
% model parameters. However, users can run the code below to generate the
% posterior samples for any model parameter at each quantile level. 
% Posterior samples are also available upon request.

% for each quantile level, it takes about 1 hour to fit the model on 
% a 64-bit operating system with 2 processors and 32GB RAM.

for i=1:5
    result=QR(model,qt(i),MCMCspecs,10*i);
    MCMC_betat=result.MCMC_betat;
    MCMC_sigma=result.MCMC_sigma;
    MCMC_betat=MCMC_betat(:,(model.T+1):end);    
    dlmwrite(sprintf('Y:/submission/realdata/output/QR/MCMC_betat_%d.txt',qt(i)*100),MCMC_betat,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/QR/MCMC_sigma_%d.txt',qt(i)*100),MCMC_sigma,'delimiter','\t','precision','%12.6e');
end

rmpath('Y:/submission/BayesianQR');

%% adjusted Bayesian FQR

% add the path of the functions needed to run adjusted Bayesian FQR
addpath('Y:/submission/BayesianFQR_corrected_likelihood');

% for space considerations, we only provided the posterior samples for 
% betat at qt=0.9, which is used to plot Figure 4 in the paper. However,
% Users can run the code below to generate the posterior samples for all 
% model parameters at each quantile level. Posterior samples are also 
% available upon request.

for i=1:5
    result=FQR_corrected_HS(model,ones(model.T,1),qt(i),MCMCspecs,100*i);
    MCMC_betat=result.MCMC_betat;
    MCMC_beta=result.MCMC_beta;
    MCMC_lambda2=result.MCMC_lambda2;
    MCMC_tau2=result.MCMC_tau2;
    MCMC_betat=MCMC_betat(:,(model.T+1):end);
    MCMC_beta=MCMC_beta(:,(model.wavespecs.K+1):end);
    MCMC_lambda2=MCMC_lambda2(:,(model.wavespecs.K+1):end);
    MCMC_tau2=MCMC_tau2(:,(model.wavespecs.J+1):end);
    dlmwrite(sprintf('Y:/submission/realdata/output/HS_corrected/MCMC_betat_%d.txt',qt(i)*100),MCMC_betat,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS_corrected/MCMC_beta_%d.txt',qt(i)*100),MCMC_beta,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS_corrected/MCMC_lambda2_%d.txt',qt(i)*100),MCMC_lambda2,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/realdata/output/HS_corrected/MCMC_tau2_%d.txt',qt(i)*100),MCMC_tau2,'delimiter','\t','precision','%12.6e');
end

rmpath('Y:/submission/BayesianFQR_corrected_likelihood');

%% Bootstrap-based FQR 
% Bootstrap-based FQR approaches are implemented using 'quantreg' package 
% (Koenker, 2017) in R, with post-smoothing on coefficient estimates 
% by splines or wavelet denoising. 

% Please refer to the R code submission/realdata/realdata_modelfit_bootstrap.R.

% for space considerations, we only provided the bootstrap samples for 
% betat obtained by the two-step bootstrap-based FQR with wavelets denoising
% at qt=0.9, which are used to plot Figure 4 in the paper. However, users
% can run the provided scripts to generate the bootstrap samples for
% betat for all bootstrap-based approaches at each quantile level. These
% bootstrap samples are also available upon request.

% Two-step bootstrap-based FQR with wavelets denoising
clear;

qt=[0.1 0.25 0.5 0.75 0.9];

for i1=1:5
    raw=dlmread(sprintf('Y:/submission/realdata/output/freqQR/raw/BS_betat_%d.txt',qt(i1)*100));
    smoothed=NaN(size(raw));
    for i=1:size(raw,1)
        smoothed(i,:)=wden(raw(i,:),'heursure','s','mln',4,'db4');
    end
    dlmwrite(sprintf('Y:/submission/realdata/output/freqQR/wavelets/BS_betat_%d.txt',qt(i1)*100),smoothed,'delimiter','\t','precision','%12.6e');
    clear smoothed;
end

%% Functional mean regression 
% Functional mean regression is performed using the WFMM model by Morris (2006).
% The implementation is done using the MATLAB code by Zhu (2011).

%% References
% Koenker, R. (2017). quantreg: Quantile Regression. R package version 5.33.
% Vienna: R Foundation.

% Morris, J. S., & Carroll, R. J. (2006). Wavelet-based functional mixed models. 
% Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(2), 179-199.

% Zhu, H., Brown, P. J., & Morris, J. S. (2011). 
% Robust, adaptive functional regression in functional mixed model framework. 
% Journal of the American Statistical Association, 106(495), 1167-1179.
