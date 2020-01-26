function result=FQR_corrected_HS(model,sigma,qt,MCMCspecs,seed)
%% set the seed
rng(seed);

%% general p
T=model.T;
p=model.p;

%% run Bayesian QR to get the posterior mean and covariance matrix for each t
resultQR=QR(model,sigma,qt,MCMCspecs);

betathat=resultQR.betathat;

COV_post=resultQR.postcov;
COV_lik=NaN(p,p,T);
Sn=model.X'*model.X/model.n;

for t=1:T
    COV_lik(:,:,t)=model.n*qt*(1-qt)*COV_post(:,:,t)*Sn*COV_post(:,:,t);
end

Cond1=NaN(p,p-1,T);
Cond2=NaN(p,T);

for i=1:p
    perm=[i, 1:(i-1),(i+1):p];
    for t=1:T
        COV_perm=COV_lik(perm,perm,t);
        Sigma_aa=COV_perm(1,1);
        Sigma_ab=COV_perm(1,2:end);
        Sigma_bb=COV_perm(2:end,2:end);
        Cond1(i,:,t)=Sigma_ab/Sigma_bb;
        Cond2(i,t)=Sigma_aa-(Sigma_ab/Sigma_bb)*(Sigma_ab');
    end
end
                     

%% initialze parameters
wavespecs=model.wavespecs;

beta=betathat*wavespecs.W1t;
betat=betathat;

nu=1./gamrnd(1/2,1,p,wavespecs.K);
lambda2=1./gamrnd(1/2,nu);

sa2=ones(p,1);
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

%% MCMC
ii=0; % indicator for the samples saved after burnin and thinning.

tic;

for i=1:(B*thin+burnin)
    [beta,betat]=UpdateBetaNoOrthog(beta,betat,betathat,Cond1,Cond2,lambda2,tau2,wavespecs,MCMCspecs);
        
    lambda2=UpdateLambda(beta,nu,tau2,wavespecs,MCMCspecs);
    
    nu=1./gamrnd(1,1./(1+1./lambda2));
    
    tau2=UpdateTau(beta,lambda2,gamma,model,wavespecs,MCMCspecs);
    
    gamma=1./gamrnd(1,1./(1./repmat(sa2,1,wavespecs.J)+1./tau2));
    
    sa2=UpdateSa2(gamma,wavespecs,MCMCspecs);
    
    if (i>burnin)&&(mod(i-burnin,thin)==0)      
        ii=ii+1;
        
        is=mod(ii-1,blocksize)+1; % this is the row number in the block.
        
        MCMC_beta(is,:)=reshape(beta',1,numel(beta)); 
        MCMC_betat(is,:)=reshape(betat',1,numel(betat));
        MCMC_lambda2(is,:)=reshape(lambda2',1,numel(lambda2));
        MCMC_tau2(is,:)=reshape(tau2',1,numel(tau2));
        MCMC_sa2(is,:)=sa2';        
    end
        
end
        
time_spend=toc;

%% save results
result.time_spend=time_spend;

result.MCMC_betat=MCMC_betat;
result.MCMC_beta=MCMC_beta;
result.MCMC_lambda2=MCMC_lambda2;
result.MCMC_tau2=MCMC_tau2;
result.MCMC_sa2=MCMC_sa2;






