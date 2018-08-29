function [Kj]=get_Kj_2d(S);

%%%% Take Matlab wavelet toolbox 3d wavedec object C and output single
%%%% vector of wavelet coefficients;

J=(size(S,1)-2)*3+1;
Kj=zeros(J,1);
temp=S(1:(end-1),1).*S(1:(end-1),2);
K=[1,repmat(3,1,length(temp)-1)]*temp;
Kj=[temp(1),reshape(repmat(temp(2:end),1,3)',1,3*(length(temp)-1))];

