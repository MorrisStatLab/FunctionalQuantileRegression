function B=rowprod(v,A)
% This function compute the diag(v)*A, i.e., using each component of v
% multiply each row of A. 
% v: a column vector of length n
% A: a matrix of length n by p.
B=repmat(v,1,size(A,2)).*A;
