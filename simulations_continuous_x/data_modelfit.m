clear;

qt=[0.1 0.2 0.5 0.8 0.9];

% MCMC specifications 
MCMCspecs.minVC=1e-6;
MCMCspecs.maxO=1e20;

MCMCspecs.B=2000;
MCMCspecs.burnin=2000;
MCMCspecs.thin=3;
MCMCspecs.blocksize=2000;

MCMCspecs.time_update=1000;

rng(0,'twister');
seed = randi([0 10000],100,5);
seed2 = randi([10000 20000],100,5);
seed3 = randi([20000 30000],100,5);

%% Bayesian FQR

% add the path of the functions needed to run Bayesian FQR
addpath('Y:/submission/BayesianFQR');

% for each replicate and quantile, it takes about 15 minutes to fit the   
% model on a 64-bit operating system with 2 processors and 32GB RAM 

for i1=1:100
    load(sprintf('Y:/submission/simulations_continuous_x/data/model%d.mat',i1))
    for i2=1:5
        result=FQR_HS(model,qt(i2),MCMCspecs,seed(i1,i2));
        MCMC_betat=result.MCMC_betat;
        save(sprintf('Y:/submission/simulations_continuous_x/output/HS/MCMC_betat_%d_rep_%d.mat',qt(i2)*100,i1),'MCMC_betat');      
    end
    clear model;
end

rmpath('Y:/submission/BayesianFQR');

%% naive Bayesian QR

% add the path of the functions needed to run naive Bayesian FQR
addpath('Y:/submission/BayesianQR');

% for each replicate and quantile, it takes about 8 minutes to fit the   
% model on a 64-bit operating system with 2 processors and 32GB RAM

for i1=1:100
    load(sprintf('Y:/submission/simulations_continuous_x/data/model%d.mat',i1))
    for i2=1:5
        result=QR(model,qt(i2),MCMCspecs,seed2(i1,i2));
        MCMC_betat=result.MCMC_betat;
        save(sprintf('Y:/submission/simulations_continuous_x/output/QR/MCMC_betat_%d_rep_%d.mat',qt(i2)*100,i1),'MCMC_betat');
    end
    clear model;
end

rmpath('Y:/submission/BayesianQR');

%% adjusted Bayesian FQR

% add the path of the functions needed to run adjusted Bayesian FQR
addpath('Y:/submission/BayesianFQR_corrected_likelihood');

% for each replicate and quantile, it takes about 15 minutes to fit the   
% model on a 64-bit operating system with 2 processors and 32GB RAM 

for i1=1:100
    load(sprintf('Y:/submission/simulations_continuous_x/data/model%d.mat',i1))
    for i2=1:5
        result=FQR_corrected_HS(model,ones(model.T,1),qt(i2),MCMCspecs,seed3(i1,i2));
        MCMC_betat=result.MCMC_betat;
        save(sprintf('Y:/submission/simulations_continuous_x/output/HS_corrected/MCMC_betat_%d_rep_%d.mat',qt(i2)*100,i1),'MCMC_betat');      
    end
    clear model;
end

rmpath('Y:/submission/BayesianFQR_corrected_likelihood');