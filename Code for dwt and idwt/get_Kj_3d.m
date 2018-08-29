function [Kj,Kj_all]=get_Kj_3d(C,coeffs_keep);

%%%% Take Matlab wavelet toolbox 3d wavedec object C and output single
%%%% vector of wavelet coefficients;

J=length(C.dec);
Kj_all=zeros(J,1);
Kj=zeros(J,1);
temp=(C.sizes(1:C.level,1).*C.sizes(1:C.level,2).*C.sizes(1:C.level,3));
K=[8,repmat(7,1,C.level-1)]*temp;
Kj_all=[temp(1),reshape(repmat(temp,1,7)',1,7*length(temp))];
temp_first=[1,1+cumsum(Kj_all(1:(end-1)))];
temp_last=cumsum(Kj_all);
Kj=Kj_all;
for (i=1:J);
    Kj(i)=sum(coeffs_keep(temp_first(i):temp_last(i)));
end;
