function [beta,betat]=UpdateBeta(betatns,Vbetatns,lambda2,tau2,wavespecs,MCMCspecs)

W1=wavespecs.W1;
W1t=wavespecs.W1t;

minVC=MCMCspecs.minVC;

betans=betatns*W1t;

Vbetans=W1t'*sparse(diag(Vbetatns))*W1t;

Vbetans=(Vbetans+Vbetans')/2;

scale=max(lambda2.*tau2,minVC^2);

shrink1=Vbetans*sparse(diag(1./scale));

shrink2=shrink1*((eye(wavespecs.K)+shrink1)\eye(wavespecs.K));

Sigma=Vbetans-shrink2*Vbetans;

Sigma=(Sigma+Sigma')/2;

Sigma=Sigma+1e-10*eye(wavespecs.K);

mu=betans'-shrink2*betans';

beta=mvnrnd(mu,Sigma);

betat=beta*W1;
