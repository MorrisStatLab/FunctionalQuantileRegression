function Y=IDWT_rows(D,wavespecs)

% help file needs to be updated
% dwt_rows(infile,nlevels,wavelet): Compute IDWT on rows of infile, using wavelet and 
%                                   number of levels specified.
%
% Input: infile : matrix, with each row containing wavelet coefficients for function on which to perform IDWT (n x k)
%        nlevels: Number of levels of wavelet transform
%       wavespecs: structure with information on wavelet settings,
%                   including:
%           wavelet: Name of wavelet basis to use
%
% Output: D : matrix with each row containing a function (n x T)
%         K : vector of number of coefficients per level (1 x nlevels)
%         T : number of sampled points per function 
%
% Once can change the extension mode  to periodic by using dwtmode('per').
%tic

if wavespecs.compress<1
   temp=zeros(size(D,1),length(wavespecs.keep));
   temp(:,wavespecs.keep==1)=D;
   D=temp;
end;
if wavespecs.ndim==1
    dwtmode(wavespecs.boundary,'nodisp');
    if wavespecs.compress<1
        K=[wavespecs.Kj_all,wavespecs.T];
    else
        K=[wavespecs.Kj,wavespecs.T];
    end;
    Y=NaN(size(D,1),wavespecs.T);
    for i=1:size(D,1)
        y=waverec(D(i,:),K,wavespecs.wavelet);
        Y(i,:)=y;
    end;
elseif wavespecs.ndim==2
    if wavespecs.rectangular==1
        Y=NaN(size(D,1),wavespecs.t(1)*wavespecs.t(2));
        K1=[wavespecs.Kj1,wavespecs.t(1)];
        K2=[wavespecs.Kj2;wavespecs.t(2)];
        D(:,wavespecs.reorder)=D;
        for i=1:size(D,1)
            d=reshape(D(i,:),sum(wavespecs.Kj2),sum(wavespecs.Kj1));
            y2=NaN(wavespecs.t(2),size(d,2));
            dwtmode(wavespecs.boundary(2,:),'nodisp');
            for j=1:size(d,2);
                [y2(:,j)]=waverec(d(:,j),K2,wavespecs.wavelet(2,:));
            end;
            y1=NaN(wavespecs.t(2),wavespecs.t(1));
            dwtmode(wavespecs.boundary(1,:),'nodisp');
            for j=1:size(y2,1);
                y1(j,:)=waverec(y2(j,:),K1,wavespecs.wavelet(1,:));
            end;
            Y(i,:)=reshape(y1,1,wavespecs.t(1)*wavespecs.t(2));
            if (i/100==floor(i/100))
                i,toc
            end;
        end
    else
        Y=NaN(size(D,1),wavespecs.t(1)*wavespecs.t(2));
        [C,S]=wavedec2(reshape(Y(1,:),wavespecs.t(1),wavespecs.t(2)),wavespecs.nlevels,wavespecs.wavelet);
        for i=1:size(D,1)
            y=waverec2(D(i,:),S,wavespecs.wavelet);
            Y(i,:)=reshape(y,1,wavespecs.t(1)*wavespecs.t(2));
        end
    end
else
    C=wavespecs.C;
    T=C.sizeINI(1)*C.sizeINI(2)*C.sizeINI(3);
    Y=NaN(size(D,1),T);
    for i=1:size(D,1)
        y=waverec3(unvectorize_3d_wavelet(D(i,:),C));
        Y(i,:)=reshape(y,1,T);
    end
end




