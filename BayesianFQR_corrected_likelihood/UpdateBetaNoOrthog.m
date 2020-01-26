function [beta,betat]=UpdateBetaNoOrthog(beta,betat,betathat,Cond1,Cond2,lambda2,tau2,wavespecs,MCMCspecs)

%%%%    Updates regression coefficients Beta from its full conditional
%%%%    (based on the corrected Gaussian likelihood) in the wavelet space.
%%%%        * Updates Beta one-at-a-time 
%%%%        * Assumes non-orthogonal design matrix X
%%%%
%%%%
%%%%    Functions needed: UpdateBeta;


% general p
p=size(beta,1);
T=wavespecs.T;

tau2=tau2*uneqkron(wavespecs.Kj)';
    
Betatnsi=NaN(p,T);

for i=1:p
    idx=[1:(i-1),(i+1):p];
    for t=1:T
        Betatnsi(i,t)=betathat(i,t)+Cond1(i,:,t)*(betat(idx,t)-betathat(idx,t));
    end
    [beta(i,:),betat(i,:)]=UpdateBeta(Betatnsi(i,:),Cond2(i,:),lambda2(i,:),tau2(i,:),wavespecs,MCMCspecs);
end


