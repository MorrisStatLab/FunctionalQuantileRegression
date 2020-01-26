function pho=check_function(x,q)
% This function calculates the check function pho(x)at quantile q.

[m,n]=size(x);
pho=zeros(m,n);
pho(x>0)=q*x(x>0);
pho(x<0)=(q-1)*x(x<0);


