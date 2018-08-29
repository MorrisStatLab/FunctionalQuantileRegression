function [theta_mle,beta_mle,Var_beta]=regression_mle(W,n,q)
% This function compute the maximum likelihood estimation of the model:
% d_{jk}=XB_{jk}+E_{jk}
% where E(E_{jk})=0, Var(E_{jk})=2*s*_{jk}/q/(1-q).
% Output:
%      theta_mle: the estimated s*_{jk}.
%      beta_mle: the MLE of B.
%      Var_beta: the estimated marginal variance of B_mle.
[beta_mle,theta_mle,XtX_inv]=sweep_simple_regress(W,n);

Var_beta=diag(XtX_inv)*theta_mle'; % The variance of beta_mle
theta_mle=theta_mle*q*(1-q)/2; % The estimate of s*_{jk}