% We summarize the estimation and inference performance for beta2 at qt=0.1
% below, and procedures are the same for beta3 and other quantile levels.
clear;

% load in the ground truth for beta2
beta2=dlmread('Y:/submission/simulations_continuous_x/beta2.txt');
T=size(beta2,2);

% add the path of the functions needed to calculate SimBaS, sensitivity
% and false positive
addpath('Y:/submission/ROCfunctions');

%% For space considerations, we do not provide posterior or bootstrap
% samples based on the simulated datasets, but we illustrate below how to 
% use posterior (or bootstrap) samples to summarize the estimation and
% inferential performance for an approach given the ground truth, 
% including calculation of IMSE, sensitivity and false positive rate for
% identifying functional regions of interest based on SimBaS and effect size.

qt=0.1;
beta=beta2(1,:);

% A magnitude smaller than delta1 is effectively treated as 0
delta1=0.01;

% The minimum effect size of practical significance that we want to detect
delta2=0.3;

% Different cut-off values considered for SimBaS
alpha=[0.0005 0.001:0.001:0.005 0.01:0.01:0.99 0.995 0.999 0.9995];

%% naive Bayesian QR
% load in posterior samples of beta_2(t)
load(sprintf('Y:/submission/simulations_continuous_x/output/QR/MCMC_betat_%d_rep_1.mat',qt*100));
MCMC_P1=MCMC_betat(:,(T+1):2*T);

% sensitivity and false positive rate for detecting regions in beta2 of 
% at least delta2 based on SimBaS at various alpha levels for this replicate
[sens1,~,fp1]=simbas_ROC_flag_contiguous(beta,MCMC_P1,delta1,delta2,alpha);

% posterior mean estimate of beta_2(t)
betahat1=mean(MCMC_P1);

% IMSE of beta_2(t) for this replicate
sum((betahat1-beta).^2)

%% Bayesian FQR with horseshoe prior
% load in posterior samples of beta_2(t)
load(sprintf('Y:/submission/simulations_continuous_x/output/HS/MCMC_betat_%d_rep_1.mat',qt*100));
MCMC_P2=MCMC_betat(:,(T+1):2*T);

% sensitivity and false positive rate for detecting regions in beta2 of 
% at least delta2 based on SimBaS at various alpha levels for this replicate
[sens2,~,fp2]=simbas_ROC_flag_contiguous(beta,MCMC_P2,delta1,delta2,alpha);

% posterior mean estimate of beta_2(t)
betahat2=mean(MCMC_P2);

% IMSE of beta_2(t) for this replicate
sum((betahat2-beta).^2)
    
%% adjusted Bayesian FQR with horseshoe prior
% load in posterior samples of beta_2(t)
load(sprintf('Y:/submission/simulations_continuous_x/output/HS_corrected/MCMC_betat_%d_rep_1.mat',qt*100));
MCMC_P3=MCMC_betat(:,(T+1):2*T);

% sensitivity and false positive rate for detecting regions in beta2 of 
% at least delta2 based on SimBaS at various alpha levels for this replicate
[sens3,~,fp3]=simbas_ROC_flag_contiguous(beta,MCMC_P3,delta1,delta2,alpha);

% posterior mean estimate of beta_2(t)
betahat3=mean(MCMC_P3);

% IMSE of beta_2(t) for this replicate
sum((betahat3-beta).^2)
    