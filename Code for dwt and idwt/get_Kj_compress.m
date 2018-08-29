function [Kj,Kj_all]=get_Kj_compress(wavespecs)

%%%% Given wavespecs object with Kj=# coeffs per level and
%%%% keep = vector of 0's and 1's with 1 in coefficients to "keep"
%%%% figure out J and Kj = # levels and # coefficients per level for reduce
%%%% set of coefficients that are kept

J=length(wavespecs.Kj);
Kj_all=wavespecs.Kj;
Kj=zeros(1,J);
temp_first=[1,1+cumsum(Kj_all(1:(end-1)))];
temp_last=cumsum(Kj_all);
for i=1:J;
    Kj(i)=sum(wavespecs.keep(temp_first(i):temp_last(i)));
end;
%%% Also need to add what to do when there are "0 levels" or levels with
%%% very small number of wavelet coefficients which we want to potentially
%%% combine together
%%%
%%% Also, need to think about implications on spike-slab prior since we are
%%% eliminating a lot of the sparsity.