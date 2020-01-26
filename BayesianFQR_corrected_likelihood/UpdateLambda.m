function lambda2=UpdateLambda(beta,nu,tau2,wavespecs,MCMCspecs)

tau2=tau2*uneqkron(wavespecs.Kj)';

beta2=beta.^2;

lambda2=1./gamrnd(1,1./(beta2./tau2/2+1./nu));

lambda2=max(lambda2,MCMCspecs.minVC);