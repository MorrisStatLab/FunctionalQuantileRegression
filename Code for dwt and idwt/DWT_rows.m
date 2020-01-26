function [D,wavespecs,D_compressed]=DWT_rows(X,wavespecs)

% Compute DWT 
% If ndim=1, then compute DWT on rows of X, using wavelet and number of levels specified.
% If ndim=2, then compute rectangular/square DWT on infile, using wavelet and number of levels specified.
%       If ndim=2, then user must specify wavespecs.t(1) (# rows) and
%       wavespecs.t(2) (# columns) such that wavespecs.t1*wavespecs.t2=T
%       If ndim=3, then user must specify wavepecs.t(1), .t(2), .t(3) such
%                   that t1*t2*t3=T.
% Input: X: matrix, with each row containing a function on which to perform DWT (n x T)
%           if ndim=2, then column stacked vector containing 2d image
%           with wavespecs.t1 rows and wavespecs.t2 columns
%           if ndim=3, then stacked vector stacking depth then column;
%        wavespecs: structure containing wavelet settings, including
%           nlevels: Number of levels of wavelet transform
%                    If ndim=2, then nlevels is vector of length 2
%                    If not specified, per Percival and Walden (page 145),
%                    J = log_2{(T-1)/(L-1)} + 1 where L=length of chosen
%                    wavelet filter.
%           wavelet: Name of wavelet basis to use
%                    Default "db4" Daubechies wavelet with 4 vanishing
%                    moments
%           boundary: a string, indicates which type of Discrete wavelet transform
%                   extension mode to use. Can be: 'sym' (default),'per','zpd'
%                   See help dwtmode.
%                   Note: if ndim=2 and rectangular=1 then wtmode is a list
%                   with 2 elements containing wtmode for rows and columns
%           compress: Level of compression (0<compress<=1) -- keep only
%                   those coefficients in smallest set necessary to maintain 
%                   this proportion of total variability for each of the n
%                   functions (=1 by default meaning no compression)
%           rectangular: If =1 and ndim=2, then rectangular DWT used;
%                       otherwise square wavelet transform used. 
%

% Output: D : matrix with each row containing the wavelet coefficients for 
%               the function (n x K)
%         D_compressed: if compression done, then only the kept set of
%               wavelet coefficients (n x K*)
%         wavespecs: updated wavespecs, including
%               J: number of wavelet levels (nlevels+1)
%               K: total number of wavelet coefficients
%               Kj: vector of length J containing number of compressed
%                   wavelet coefficients at each level
%               Kj_all: if compress=1, then number of total wavelet
%                   coefficients at each level, before compression
%               keep = vector of length K, with indicator "1" if basis kept
%                       and "0" if not kept after compression.

if nargin<2
    wavespecs.wavelet='db4';
    wavespecs.nlevels=floor(log(size(X,2))/log(2));
    wavespecs.ndim=1;
    wavespecs.compress=1;
end;

if isfield(wavespecs,'compress')==0 
    wavespecs.compress=1;
end;
if isfield(wavespecs,'ndim')==0
    wavespecs.ndim=1;
end;
if isfield(wavespecs,'wavelet')==0
    wavespecs.wavelet='db4';
end;
if wavespecs.ndim==2&isfield(wavespecs,'rectangular')==0
    wavespecs.rectangular=1;
end;
if wavespecs.ndim==2&wavespecs.rectangular==1
    wavespecs.wavelet=[wavespecs.wavelet;wavespecs.wavelet];
end;
if wavespecs.ndim==2&isfield(wavespecs,'t')==0
    'For 2D functions must let wavespecs.t=[t1, t2] be a vector of length 2 giving # rows/columns such that ncol(X)=t1*t2'
    return
end;
if wavespecs.ndim==2&size(X,2)~=wavespecs.t(1)*wavespecs.t(2)
    'Number of rows wavespecs.t(1) * Number of columns wavespecs.t(2) is not equal to number of columns of X'
    return
end;
if isfield(wavespecs,'nlevels')==0

    if wavespecs.ndim==1
        L=length(wfilters(wavespecs.wavelet));
        wavespecs.nlevels=floor(log((size(X,2)-1)/(L-1))/log(2));
    elseif wavespecs.ndim==2
        if wavespecs.rectangular==0
            L=length(wfilters(wavespecs.wavelet));
            wavespecs.nlevels=floor(log((min(wavespecs.t(1),wavespecs.t(2))-1)/(L-1))/log(2));
        else
            wavespecs.nlevels=[0,0];
            L1=length(wfilters(wavespecs.wavelet(1,:)));
            L2=length(wfilters(wavespecs.wavelet(2,:)));
            wavespecs.nlevels(1)=floor(log((wavespecs.t(1)-1)/(L1-1))/log(2));
            wavespecs.nlevels(2)=floor(log((wavespecs.t(2)-1)/(L2-1))/log(2));
        end;
    else
        L=length(wfilters(wavespecs.wavelet));
        wavespecs.nlevels=floor(log((min(wavespecs.t)-1)/(L-1))/log(2));
    end;
end;
if isfield(wavespecs,'boundary')==0
    wavespecs.boundary='per';
end;
if wavespecs.ndim==2&wavespecs.rectangular==1&size(wavespecs.boundary,1)==1
    wavespecs.boundary=[wavespecs.boundary;wavespecs.boundary];
end;

tic
if wavespecs.ndim==1
        n=size(X,1);
        x=X(1,:);
        dwtmode(wavespecs.boundary,'nodisp');
        [C,S0]=wavedec(x,wavespecs.nlevels,wavespecs.wavelet);
        D=zeros(n,length(C)); 
        wavespecs.Kj=S0(1:(end-1));
        wavespecs.K=sum(wavespecs.Kj);
        wavespecs.J=length(wavespecs.Kj);
        wavespecs.T=size(X,2);
        D(1,:)=C;
    for i=2:n
        D(i,:)=wavedec(X(i,:),wavespecs.nlevels,wavespecs.wavelet);
    end;
elseif wavespecs.ndim==2
    x=reshape(X(1,:),wavespecs.t(1),wavespecs.t(2));
    n=size(X,1);
    if wavespecs.rectangular==1  
        dwtmode(wavespecs.boundary(1,:),'nodisp');
        [C11,S1]=wavedec(x(1,:),wavespecs.nlevels(1),wavespecs.wavelet(1,:));
        C1=zeros(size(x,1),size(C11,2));
        C1(1,:)=C11;
        for j=2:size(x,1);
            [C1(j,:)]=wavedec(x(j,:),wavespecs.nlevels(1),wavespecs.wavelet(1,:));
        end;
        dwtmode(wavespecs.boundary(2,:),'nodisp');
        [C21,S2]=wavedec(C1(:,1),wavespecs.nlevels(2),wavespecs.wavelet(2,:));
        C2=zeros(size(C21,1),size(C1,2));
        C2(:,1)=C21;
        for j=2:size(C1,2);
            [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels(2),wavespecs.wavelet(2,:));
        end;
        %%%% Should be reordered in blocks of K1 x K2 
        [reorder,wavespecs.Kj]=reorder_rectangular(S1(1:(end-1)),S2(1:(end-1)));
        wavespecs.reorder=reorder;
        wavespecs.Kj1=S1(1:(end-1));
        wavespecs.Kj2=S2(1:(end-1));
        wavespecs.j(1)=length(wavespecs.Kj1);
        wavespecs.j(2)=length(wavespecs.Kj2);
        wavespecs.t(1)=S1(end);
        wavespecs.t(2)=S2(end);
        wavespecs.T=S1(end)*S2(end);
        C=reshape(C2,1,size(C2,1)*size(C2,2));
        D=zeros(n,length(C));  
        D(1,:)=C;
    else
        [C,S]=wavedec2(x,wavespecs.nlevels,wavespecs.wavelet);
        D=zeros(n,length(C)); 
        D(1,:)=C;
        wavespecs.S=S;
        wavespecs.Kj=get_Kj_2d(S);
    end
    for i=2:n
        x=reshape(X(i,:),wavespecs.t(1),wavespecs.t(2));
        if wavespecs.rectangular==1  
            dwtmode(wavespecs.boundary(1,:),'nodisp');
            C1=zeros(size(x,1),size(C2,2));
            for j=1:size(x,1);
                [C1(j,:)]=wavedec(x(j,:),wavespecs.nlevels(1),wavespecs.wavelet(1,:));
            end;
            C2=zeros(size(C2,1),size(C2,2));
            dwtmode(wavespecs.boundary(2,:),'nodisp');
            for j=1:size(C1,2);
                [C2(:,j)]=wavedec(C1(:,j),wavespecs.nlevels(2),wavespecs.wavelet(2,:));
            end;
            %%%% Should be reordered in blocks of K1 x K2 but leave alone for
            %%%% now
            C=reshape(C2,1,size(C2,1)*size(C2,2));
            D(i,:)=C(reorder);
            if (i/100==floor(i/100))
                i,toc
            end;
        else
            D(i,:)=wavedec2(x,wavespecs.nlevels,wavespecs.wavelet);
        end;
    end;
    wavespecs.K=sum(wavespecs.Kj);
else   
    if wavespecs.boundary=='per' 
        wavespecs.boundary='ppd';
    end;
    [C]=wavedec3(reshape(X(1,:),wavespecs.t(1),wavespecs.t(2),wavespecs.t(3)),wavespecs.nlevels,wavespecs.wavelet,wavespecs.boundary);
    temp=(C.sizes(1:C.level,1).*C.sizes(1:C.level,2).*C.sizes(1:C.level,3));
    wavespecs.K=[8,repmat(7,1,C.level-1)]*temp;
    wavespecs.Kj=get_Kj_3d(C,ones(1,wavespecs.K));
    wavespecs.C=C;
    D=NaN(size(X,1),wavespecs.K);
    wavespecs.J=length(wavespecs.Kj);
    D(1,:)=vectorize_3d_wavelet(C);
    for i=1:size(X,1)
        [D(i,:)]=vectorize_3d_wavelet(wavedec3(reshape(X(i,:),wavespecs.t(1),wavespecs.t(2),wavespecs.t(3)),wavespecs.nlevels,wavespecs.wavelet,wavespecs.boundary));
    end
end; %if wavespec.dim

%'Done with Wavelet Transforms ',toc

%% write standalone function for compress so user can try different levels
%% of compression very quickly

if wavespecs.compress<1
    [D,wavespecs]=wavelet_compress(D,wavespecs);
end;
%'Done with wavelet compression',toc




