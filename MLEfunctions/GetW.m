function W=GetW(model,D)

%%%%% W=GetW(model,D)
%%%%%   Compute cross-products matrices
%%%%%           X'X       X'D
%%%%%           D'X       D'D
% Input:    
%
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               D = (n x K) matrix of wavelet coefficients
%
% Output:   W = structure with elements:
%           XtX,XtD,DtD

W.XtX=model.X'*model.X;
W.XtD=model.X'*D;
W.DtD=D'*D;



