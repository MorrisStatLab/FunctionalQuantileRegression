function A=uneqkron(n)

% UNEQKRON Build a vector of ones of dimension n[i], i=1,...,k within identity matrix of size k.
%          Builds a mean or sum operator for unbalanced data.
%          That is, it generalizes kron(eye(k),ones(n,1)) to case where n is unbalanced.
%
%          This is useful for calculating means and SS in ANOVA using
%          matrix form, and also for "expanding" matrices
%
%
% Input n=vector of sample sizes for k groups
%
%

k=length(n);
nsum=sum(n);
ncumsum=cumsum(n);
A=zeros(nsum,k);

A(1:ncumsum(1),1)=ones(n(1),1);

for i=2:k
    A((ncumsum(i-1)+1):ncumsum(i),i)=ones(n(i),1);    
end

%% To get each diagonal block to consist of a n[i] x n[i] matrix of ones,
%% just take A*A'
%%

