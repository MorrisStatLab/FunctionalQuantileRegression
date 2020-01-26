function betat=UpdateBetaNoOrthog_QR(betat,Vbetat,W,model)

%%%%    Updates regression coefficients Beta from its full conditional
%%%%    (based on the uncorreted AL likelihood) in the wavelet space.
%%%%        * Updates Beta one-at-a-time 
%%%%        * Assumes non-orthogonal design matrix X
%%%%
%%%%
%%%%    Functions needed: UpdateBeta;


% general p
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



