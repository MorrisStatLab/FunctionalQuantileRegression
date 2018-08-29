clear;

% observed grid and groud truth
x0=dlmread('Y:/submission/simulations2/x0.txt');
beta2=dlmread('Y:/submission/simulations2/beta2.txt');

T=size(beta2,2);

% add the path of the functions needed to calculate SimBaS, sensitivity
% and false positive
addpath('Y:/submission/ROCfunctions');

%% For space considerations, we do not provide the posterior (or bootstrap)
% samples based on the simulated datasets, but we illustrate below how to 
% use posterior (and bootstrap) samples to summarize the estimation and
% inferential performance for an approach given the ground truth, 
% including calculation of IMSE and IVar, sensitivity and false positive
% for identifying group difference based on SimBaS and effect size.

% The posterior (and bootstrap) samples obtained by each approach 
% considered in the paper are available upon request, at each quantile 
% level and for all 100 replicate simulation datasets.

% We show how to summarize the estimation and inference performance for  
% qt=0.9, and the procedures are the same for other quantiles.

qt=0.9;
beta=beta2(5,:);

% A magnitude smaller than delta1 in the ground truth is effectively treated
% as 0
delta1=0.011;

% The minimum effect size of practical significance that we want to identify
delta2=0.3;

% Different cut-off values considered for SimBaS
alpha=[0.0005 0.001:0.001:0.005 0.01:0.01:0.99 0.995 0.999 0.9995];

%% naive Bayesian QR
% load in posterior samples of Beta_2(t)
load(sprintf('Y:/submission/simulations2/output/QR/MCMC_betat_%d_rep_1.mat',qt*100));

% sensitivity and false positive for identifying a group difference of at
% least delta2 at various alpha levels for this replicate dataset
[sens1,~,fp1]=simbas_ROC_flag_contiguous(beta,MCMC_P,delta1,delta2,alpha);

% posterior mean estimate of Beta_2(t)
betahat1=mean(MCMC_P);

% IMSE of Beta_2(t) for this replicate
sum((betahat1-beta).^2)
  
% IVar of Beta_2(t) for this replicate
tmp1=bsxfun(@minus,MCMC_P,betahat1);
sum(sum(tmp1.^2))/size(MCMC_P,1)

%% Bayesian FQR with horseshoe prior
% load in posterior samples of Beta_2(t)
load(sprintf('Y:/submission/simulations2/output/HS/MCMC_betat_%d_rep_1.mat',qt*100));

% sensitivity and false positive for identifying a group difference of at
% least delta2 at various alpha levels for this replicate dataset
[sens2,~,fp2]=simbas_ROC_flag_contiguous(beta,MCMC_P,delta1,delta2,alpha);

% posterior mean estimate of Beta_2(t)
betahat2=mean(MCMC_P);

% IMSE of Beta_2(t) for this replicate
sum((betahat2-beta).^2)
    
% IVar of Beta_2(t) for this replicate
tmp2=bsxfun(@minus,MCMC_P,betahat2);
sum(sum(tmp2.^2))/size(MCMC_P,1)

%% Bootstrap-based QR
% load in bootstrap samples of Beta_2(t)
BS_P=dlmread(sprintf('Y:/submission/simulations2/output/freqQR/raw/BS_betat_%d_rep_1.txt',qt*100));

% sensitivity and false positive for identifying a group difference of at
% least delta2 at various alpha levels for this replicate dataset
[sens3,~,fp3]=simbas_ROC_flag_contiguous(beta,BS_P,delta1,delta2,alpha);

% mean estimate of Beta_2(t)
betahat3=mean(BS_P);

% IMSE of Beta_2(t) for this replicate
sum((betahat3-beta).^2)

% IVar of Beta_2(t) for this replicate
tmp3=bsxfun(@minus,BS_P,betahat3);
sum(sum(tmp3.^2))/size(BS_P,1)

%% Bootstrap-based QR with spline smoothing
% load in bootstrap samples of Beta_2(t)
BS_P=dlmread(sprintf('Y:/submission/simulations2/output/freqQR/splines/BS_betat_%d_rep_1.txt',qt*100));

% sensitivity and false positive for identifying a group difference of at
% least delta2 at various alpha levels for this replicate dataset
[sens4,~,fp4]=simbas_ROC_flag_contiguous(beta,BS_P,delta1,delta2,alpha);

% mean estimate of Beta_2(t)
betahat4=mean(BS_P);

% IMSE of Beta_2(t) for this replicate
sum((betahat4-beta).^2)

% IVar of Beta_2(t) for this replicate
tmp4=bsxfun(@minus,BS_P,betahat4);
sum(sum(tmp4.^2))/size(BS_P,1)

%% Bootstrap-based QR with wavelet denoising
% load in bootstrap samples of Beta_2(t)
BS_P=dlmread(sprintf('Y:/submission/simulations2/output/freqQR/wavelets/BS_betat_%d_rep_1.txt',qt*100));

% sensitivity and false positive for identifying a group difference of at
% least delta2 at various alpha levels for this replicate dataset
[sens5,~,fp5]=simbas_ROC_flag_contiguous(beta,BS_P,delta1,delta2,alpha);

% mean estimate of Beta_2(t)
betahat5=mean(BS_P);

% IMSE of Beta_2(t) for this replicate
sum((betahat5-beta).^2)

% IVar of Beta_2(t) for this replicate
tmp5=bsxfun(@minus,BS_P,betahat5);
sum(sum(tmp5.^2))/size(BS_P,1)
