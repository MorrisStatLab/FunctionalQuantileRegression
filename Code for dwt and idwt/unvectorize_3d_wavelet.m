function Cnew=unvectorize_3d_wavelet(d,C)

%%%% Take Matlab wavelet toolbox 3d wavedec object C and output single
%%%% vector of wavelet coefficients;

J=length(C.dec);
temp=(C.sizes(1:C.level,1).*C.sizes(1:C.level,2).*C.sizes(1:C.level,3));
K=[8,repmat(7,1,C.level-1)]*temp;

Cnew.sizeINI=C.sizeINI;
Cnew.level=C.level;
Cnew.filters=C.filters;
Cnew.mode=C.mode;
Cnew.sizes=C.sizes;

k=C.sizes(1,1)*C.sizes(1,2)*C.sizes(1,3);
Cnew.dec{1}=reshape(d(1:k),C.sizes(1,1),C.sizes(1,2),C.sizes(1,3));
ctr=k+1;
j=1;
for i=1:C.level
    k=C.sizes(i,1)*C.sizes(i,2)*C.sizes(i,3);
    for l=1:7
        j=j+1;
        %d(ctr:(ctr+k-1))=reshape(C.dec{j},k,1);
        Cnew.dec{j}=reshape(d(ctr:(ctr+k-1)),C.sizes(i,1),C.sizes(i,2),C.sizes(i,3));
        ctr=ctr+k;
    end;
end;

    