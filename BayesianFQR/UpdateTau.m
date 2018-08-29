function tau2=UpdateTau(beta,lambda2,gamma,model,wavespecs,MCMCspecs)

%%%% Update the global shrinkage parameters tau2. Dimension p*J.
sumop=uneqkron(wavespecs.Kj);

alpha=1/2+repmat(wavespecs.Kj,model.p,1)/2;

beta2=beta.^2;

inv_theta=(beta2./lambda2/2)*sumop+(1./gamma);

tau2=1./gamrnd(alpha,1./inv_theta);

tau2=max(tau2,MCMCspecs.minVC);
