function sa2=UpdateSa2(gamma,wavespecs,MCMCspecs)

prior_gamma_a=0.001;
prior_gamma_b=(prior_gamma_a+1)*1;

a=wavespecs.J/2+prior_gamma_a;
inv_theta=sum(1./gamma,2)+prior_gamma_b;

sa2=1./gamrnd(a,1./inv_theta);

sa2=max(sa2,MCMCspecs.minVC);
