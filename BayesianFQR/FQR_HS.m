function result=FQR_HS(model,qt,MCMCspecs,seed)
%% Set the seed
rng(seed);

%% Parameter initialization
minVC=MCMCspecs.minVC;

Y=model.Y;
T=model.T;
p=model.p;

wavespecs=model.wavespecs;

W0=GetW(model,Y);
[beta_mle,]=sweep_simple_regress(W0,model.n);

% very vague prior for sigma(t), the scale parameter in the data space.
prior_gamma_a=repmat(0.001,T,1);
prior_gamma_b=repmat(0.001,T,1);

% beta refers to the p*K matrix of Beta*(k) in the basis space.
% betat refers to the p*T matrix of Beta(t) in the data space.
beta=beta_mle*wavespecs.W1t;
betat=beta_mle;

% lambda2 is the p*K matrix of local shrinkage parameters in the basis
% space. nu are the auxiliary variables of the same dimension used for half 
% Cauchy.
nu=1./gamrnd(1/2,1,p,wavespecs.K);
lambda2=1./gamrnd(1/2,nu);

% sa2 is the p*1 vector of hyperprior parameters for the global shrinkage 
% parameters tau2 in the basis space.
sa2=ones(p,1);

% tau2 is the p*J matrix of global shrinkage parameters in the basis space.
% gamma are the auxiliary varibles of the same dimension used for half
% Cauchy.
gamma=1./gamrnd(1/2,1,p,wavespecs.J);
tau2=1./gamrnd(1/2,gamma);

%% Other preparations and preallocations
B=MCMCspecs.B;
burnin=MCMCspecs.burnin;
thin=MCMCspecs.thin;
blocksize=MCMCspecs.blocksize;

MCMC_beta=NaN(blocksize,p*wavespecs.K);
MCMC_betat=NaN(blocksize,p*T);
MCMC_lambda2=NaN(blocksize,p*wavespecs.K);
MCMC_tau2=NaN(blocksize,p*wavespecs.J);
MCMC_sa2=NaN(blocksize,p);
MCMC_sigma=NaN(blocksize,T);

%% MCMC
ii=0; % indicator for the samples saved after burnin and thinning.

tic;

for i=1:(B*thin+burnin)
    y=model.Y-model.X*betat;
    
    % update the T*1 scale parameters sigma(t) in the data space.
    sigma=1./gamrnd(prior_gamma_a+model.n,1./(prior_gamma_b+sum(check_function(y,qt))'));
    sigma=max(sigma,minVC);
    
    % update the N*T matrix of latent variables xi in the data space.
    xi=1./invgaussrnd(1./abs(y)/qt/(1-qt),1./repmat(sigma',model.n,1)/qt/(1-qt)/2);
    xi=max(xi,minVC);
    
    % data rescaling that facilitates the update of beta and betat.
    [Vbetans_trans,W_trans,]=transform(qt,xi,sigma,model,wavespecs,1,1);
    
    % update p*K beta (and p*T betat) from the full conditional.
    [beta,betat]=UpdateBetaNoOrthog(beta,betat,Vbetans_trans,W_trans,lambda2,tau2,model,wavespecs,MCMCspecs);
    
    % update the p*K local shrinkage parameters lambda2 and 
    % auxiliary variables nu in the basis space.
    lambda2=UpdateLambda(beta,nu,tau2,wavespecs,MCMCspecs);
    
    nu=1./gamrnd(1,1./(1+1./lambda2));
    
    % update the p*J global shrinkage parameters tau2 and 
    % auxiliary variables gamma in the basis space.
    tau2=UpdateTau(beta,lambda2,gamma,model,wavespecs,MCMCspecs);
    
    gamma=1./gamrnd(1,1./(1./repmat(sa2,1,wavespecs.J)+1./tau2));
    
    % update the p*1 hyperprior parameters for tau2.
    sa2=UpdateSa2(gamma,wavespecs,MCMCspecs);
    
    if (i>burnin)&&(mod(i-burnin,thin)==0)      
        ii=ii+1;
        
        is=mod(ii-1,blocksize)+1; % this is the row number in the block.
        
        MCMC_beta(is,:)=reshape(beta',1,numel(beta)); 
        MCMC_betat(is,:)=reshape(betat',1,numel(betat));
        MCMC_lambda2(is,:)=reshape(lambda2',1,numel(lambda2));
        MCMC_tau2(is,:)=reshape(tau2',1,numel(tau2));
        MCMC_sa2(is,:)=sa2';
        MCMC_sigma(is,:)=sigma';        
    end
        
end
        
time_spend=toc;

%% Posterior mean estimate for beta(t) and sigma(t)
betathat=reshape(mean(MCMC_betat)',T,p)';
sigmahat=mean(MCMC_sigma)';

%% Save results
result.time_spend=time_spend;

result.MCMC_betat=MCMC_betat;
result.MCMC_beta=MCMC_beta;
result.MCMC_lambda2=MCMC_lambda2;
result.MCMC_tau2=MCMC_tau2;
result.MCMC_sa2=MCMC_sa2;
result.MCMC_sigma=MCMC_sigma;

result.betathat=betathat;
result.sigmahat=sigmahat;





