function result=QR(model,sigma,qt,MCMCspecs)

%% initialze beta
minVC=MCMCspecs.minVC;

Y=model.Y;
T=model.T;
p=model.p;

W0=GetW(model,Y);
[beta_mle,]=sweep_simple_regress(W0,model.n);

%% Other preparations and preallocations
B=MCMCspecs.B;
burnin=MCMCspecs.burnin;
thin=MCMCspecs.thin;
blocksize=MCMCspecs.blocksize;

MCMC_betat=NaN(blocksize,p*T);

betat=beta_mle;

%% MCMC
ii=0; % indicator for the samples saved after burnin and thinning.

tic;

for i=1:(B*thin+burnin)
    y=Y-model.X*betat;
    
    xi=1./invgaussrnd(1./abs(y)/qt/(1-qt),1./repmat(sigma',model.n,1)/qt/(1-qt)/2);
    xi=max(xi,minVC);
    
    [Vbetans_trans,W_trans,]=transform(qt,xi,sigma,model,1,1);
    
    betat=UpdateBetaNoOrthog_QR(betat,Vbetans_trans,W_trans,model);

    if (i>burnin)&&(mod(i-burnin,thin)==0)      
        ii=ii+1;
        
        is=mod(ii-1,blocksize)+1; % this is the row number in the block.
        
        MCMC_betat(is,:)=reshape(betat',1,numel(betat));       
    end
        
end
        
time_spend=toc;


%% post-processing and save results
result.time_spend=time_spend;

result.MCMC_betat=MCMC_betat;

result.betathat=reshape(mean(MCMC_betat)',T,p)';

result.postcov=NaN(p,p,T);

for t=1:T
    tmp=NaN(blocksize,p);
    for l=1:p
        tmp(:,l)=MCMC_betat(:,t+(l-1)*T);
    end
    result.postcov(:,:,t)=cov(tmp);
end


