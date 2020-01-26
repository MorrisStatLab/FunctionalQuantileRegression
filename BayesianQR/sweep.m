function A=sweep(A,k)

%%%%% Sweep over the k^th diagonal element of A.
%%%%%
%%%%% See Seber (1977) or Goodnight (1979 TAS) for description of Sweep
%%%%% operator.
%%%%%
%%%%%

D=A(k,k);
F=A(k,:)/D;
B=A(:,k);
A=A-B*F;
A(:,k)=-B/D;
A(k,:)=F;
A(k,k)=1/D;