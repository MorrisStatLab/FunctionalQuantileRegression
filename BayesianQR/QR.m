function result=QR(model,qt,MCMCspecs,seed)
%% Set the seed
rng(seed);

%% Parameter initialization
minVC=MCMCspecs.minVC;

Y=model.Y;
T=model.T;
p=model.p;

W0=GetW(model,Y);
[beta_mle,]=sweep_simple_regress(W0,model.n);

% very vague prior for the scale parameter sigma
prior_gamma_a=repmat(0.001,T,1);
prior_gamma_b=repmat(0.001,T,1);

%% Other preparations and preallocations
B=MCMCspecs.B;
burnin=MCMCspecs.burnin;
thin=MCMCspecs.thin;
blocksize=MCMCspecs.blocksize;

MCMC_betat=NaN(blocksize,p*T);
MCMC_sigma=NaN(blocksize,T);

betat=beta_mle;

%% MCMC
ii=0; % indicator for the samples saved after burnin and thinning.

tic;

for i=1:(B*thin+burnin)
    y=Y-model.X*betat;
    
    % update the T*1 scale parameters sigma(t)
    sigma=1./gamrnd(prior_gamma_a+model.n,1./(prior_gamma_b+sum(check_function(y,qt))'));
    sigma=max(sigma,minVC);
    
    % update the N*T matrix of latent variables xi
    xi=1./invgaussrnd(1./abs(y)/qt/(1-qt),1./repmat(sigma',model.n,1)/qt/(1-qt)/2);
    xi=max(xi,minVC);
    
    % data rescaling that facilitates the update of betat
    [Vbetans_trans,W_trans,]=transform(qt,xi,sigma,model,1,1);
    
    % update the p*T quantile regression coefficients
    betat=UpdateBetaNoOrthog(betat,Vbetans_trans,W_trans,model);

    if (i>burnin)&&(mod(i-burnin,thin)==0)      
        ii=ii+1;
        
        is=mod(ii-1,blocksize)+1; % this is the row number in the block.
        
        MCMC_betat(is,:)=reshape(betat',1,numel(betat)); 
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
result.MCMC_sigma=MCMC_sigma;

result.betathat=betathat;
result.sigmahat=sigmahat;
