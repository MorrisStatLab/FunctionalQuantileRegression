function [beta,betat]=UpdateBetaNoOrthog(beta,betat,Vbetat,W,lambda2,tau2,model,wavespecs,MCMCspecs)

%%%%    Update the quantile regression coefficients p*T matrix Beta from 
%%%%    full conditional f(Beta|xi,sigma,lambda2,tau2,Y)  
%%%%    Projection into the basis space and back to the data space are
%%%%    performed in this step.
%%%%
%%%%    Updates Beta one-at-a-time conditional on other covariates
%%%%       
%%%%    Functions needed: UpdateBeta


p=model.p;
T=wavespecs.T;

Btilde=Vbetat.*W.XtD;
Betatnsi=NaN(p,T);

tau2=tau2*uneqkron(wavespecs.Kj)';
    
for i=1:p  % compute beta_i conditional on beta_{-i} here.    
    Bi=NaN(1,T);
    for t=1:T
         Bi(t)=W.XtX(i,:,t)*betat(:,t); % these are the stuff to condition on.
    end 
    Betatnsi(i,:)=Btilde(i,:)+betat(i,:)-Vbetat(i,:).*Bi;
    [beta(i,:),betat(i,:)]=UpdateBeta(Betatnsi(i,:),Vbetat(i,:),lambda2(i,:),tau2(i,:),wavespecs,MCMCspecs);
end 



