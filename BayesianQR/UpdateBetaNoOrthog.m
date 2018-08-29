function betat=UpdateBetaNoOrthog(betat,Vbetat,W,model)

%%%%    Update the quantile regression coefficients p*T matrix Beta from 
%%%%    full conditional f(Beta|xi,sigma,Y)  
%%%%
%%%%    Updates Beta one-at-a-time conditional on other covariates


p=model.p;
T=size(betat,2);

Btilde=Vbetat.*W.XtD;
Betatnsi=NaN(p,T);
    
for i=1:p  % compute beta_i conditional on beta_{-i} here.     
    Bi=NaN(1,T);
    for t=1:T
         Bi(t)=W.XtX(i,:,t)*betat(:,t); % these are the stuff to condition on.
    end 
    Betatnsi(i,:)=Btilde(i,:)+betat(i,:)-Vbetat(i,:).*Bi;
    betat(i,:)=normrnd(Betatnsi(i,:),sqrt(Vbetat(i,:)));
end 



