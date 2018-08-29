function lambda2=UpdateLambda(beta,nu,tau2,wavespecs,MCMCspecs)

%%%% Update the local shrinkage parameters lambda2. Dimension p*K.
tau2=tau2*uneqkron(wavespecs.Kj)';

beta2=beta.^2;

lambda2=1./gamrnd(1,1./(beta2./tau2/2+1./nu));

lambda2=max(lambda2,MCMCspecs.minVC);