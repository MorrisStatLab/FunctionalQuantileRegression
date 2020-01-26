function samp=invgaussrnd(mu,lamb)
% This function generate random numbers from inverse gauss distribution. 
% f(x|mu,lam)=sqrt(lam/(2*pi*x^3))*exp(-lam(x-mu)^2/(2*mu^2*x)),x,mu,lam>0
% Input
%       mu: a vector or matrix or scalar.
%       lam: has the same dimension as mu.
% Output
%       samp: has the same dimension as mu, random samples generated from
%       inverse gaussian using the corresponding mu and lamb.
%
% Reference: Michael, J.R., Schucany, W.R and Haas, R.W. (1976) Generating
% Random Variates Using Transformations With Multiple Roots. Am. Stat.
% 30,88-90

[m,n]=size(mu);
w=mu.*chi2rnd(1,m,n);
c=mu./(2*lamb);
x1=mu+c.*(w-sqrt(w.*(4*lamb+w))); % Note: x1 is always < x2.
iter=0;
while (any(any(x1(:)<=0)))&&(iter<10^6)   % avoid the case that there are x1's that are zero, make up by regenerate x1.
    id0=find(x1<=0);    
    w0=mu(id0).*chi2rnd(1,length(id0),1);
    c0=mu(id0)./(2*lamb(id0));
    x1(id0)=mu(id0)+c0.*(w0-sqrt(w0.*(4*lamb(id0)+w0))); % Note: x1 is always < x2.
    iter=iter+1;
end

p1=mu./(mu+x1);
selx1=unifrnd(0,1,m,n)<p1; % This line evaluate whether select x1 or x2.
samp=selx1.*x1+(1-selx1).*(mu.^2)./x1;

