function [B_mle,Sig2_mle,XtX_inv]=sweep_simple_regress(W,n)
% This function estimate the simple regression coefficient using sweep
% operator.
% The simple regression is:
% D=XB+E
% D: n by K, B: p by K, X: n by p, E: n by K. Assume var(E_j)=sig_j^2I, for all
% columns j.
% Input:
%       W: W.XtX, W.XtD, W.DtD
% Output:
%       B_mle=(X'X)^(-1)X'D
%       Sig2_mle_j=1/n(D_j-XB_mlej)^T*(D_j-XB_mlej)
%       XtX_mle=(X'X)^(-1)

p=size(W.XtX,1);
MM=[W.XtX,W.XtD;W.XtD',W.DtD]; 
 for i=1:p
    MM=sweep(MM,i);
 end
 B_mle=MM(1:p,p+1:end);
 Sig2_mle=diag(MM(p+1:end,p+1:end))/n;
 XtX_inv=MM(1:p,1:p);