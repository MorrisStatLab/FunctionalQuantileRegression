clear;

% the quantile levels
qt=[0.1 0.2 0.8 0.9];

% MCMC specifications 
MCMCspecs.minVC=1e-6;
MCMCspecs.maxO=1e20;

MCMCspecs.B=2000;
MCMCspecs.burnin=2000;
MCMCspecs.thin=3;
MCMCspecs.blocksize=2000;

MCMCspecs.time_update=1000;

% for space considerations, we only provide 1 replicate simulated dataset,
% and do not provide posterior (or bootstrap) samples for simulated data.  
% However, all simulation datasets and posterior (or bootstrap) samples are
% available upon request, and users can also follow the script below to 
% generate posterior (or bootstrap) samples.

rng(0,'twister');

seed = randi([0 10000],4,100);
seed2 = randi([0 10000],4,100);

%% Bayesian FQR

% add the path of the functions needed to run Bayesian FQR
addpath('Y:/submission/BayesianFQR');

% for each replicate and quantile, it takes about 22 minutes to fit the   
% model on a 64-bit operating system with 2 processors and 256GB RAM 

for i1=1:100
    load(sprintf('Y:/submission/simulations1/data/model%d.mat',i1))
    for i2=1:4
        result=FQR_HS(model,qt(i2),MCMCspecs,seed(i2,i1));
        MCMC_betat=result.MCMC_betat;
        MCMC_P=MCMC_betat(:,(model.T+1):end);
        save(sprintf('Y:/submission/simulations1/output/HS/MCMC_betat_%d_rep_%d.mat',qt(i2)*100,i1),'MCMC_P');
    end
    clear model;
end

rmpath('Y:/submission/BayesianFQR');

%% naive Bayesian QR

% add the path of the functions needed to run naive Bayesian FQR
addpath('Y:/submission/BayesianQR');

% for each replicate and quantile, it takes about 8 minutes to fit the   
% model on a 64-bit operating system with 2 processors and 256GB RAM

for i1=1:100
    load(sprintf('Y:/submission/simulations1/data/model%d.mat',i1))
    for i2=1:4
        result=QR(model,qt(i2),MCMCspecs,seed2(i2,i1));
        MCMC_betat=result.MCMC_betat;
        MCMC_P=MCMC_betat(:,(model.T+1):end);
        save(sprintf('Y:/submission/simulations1/output/QR/MCMC_betat_%d_rep_%d.mat',qt(i2)*100,i1),'MCMC_P');
    end
    clear model;
end

rmpath('Y:/submission/BayesianQR');

%% Bootstrap-based FQR 
% Bootstrap-based FQR approaches are implemented using 'quantreg' package (Koenker, 2017) in R, 
% with post-smoothing on functional coefficient estimates by splines or wavelet denoising. 
% Please refer to the R code submission/simulations1/data_modelfit_bootstrap.R

% Two-step bootstrap-based FQR with wavelets denoising
clear;

qt=[0.1 0.2 0.8 0.9];

for i1=1:100
    for i2=1:4
        raw=dlmread(sprintf('Y:/submission/simulations1/output/freqQR/raw/BS_betat_%d_rep_%d.txt',qt(i2)*100,i1));
        smoothed=NaN(size(raw));
        for i=1:size(raw,1)
            smoothed(i,:)=wden(raw(i,:),'heursure','s','mln',4,'db4');
        end
        dlmwrite(sprintf('Y:/submission/simulations1/output/freqQR/wavelets/BS_betat_%d_rep_%d.txt',qt(i2)*100,i1),smoothed,'delimiter','\t','precision','%12.6e');
        clear smoothed;
    end
end

%% References
% Koenker, R. (2017). quantreg: Quantile Regression. R package version 5.33.
% Vienna: R Foundation.

