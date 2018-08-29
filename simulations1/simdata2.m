function model=simdata2(x,loc,scale,sigma,pho,shift,nsize,wavespecs)
%% generate a replicate dataset for the symmetric heavy tailed setting i.
T=size(x,1);

% Two groups
Y1=NaN(nsize,T);
Y2=NaN(nsize,T);

% Gaussian basis
basis=NaN(size(loc,1),T);
for i=1:size(basis,1)
    basis(i,:)=normpdf(x',loc(i),scale(i));
end

% AR(1) Gaussian noise added to the functional observations
corr=NaN(T,T);
for t1=1:T
    for t2=1:T
        corr(t1,t2)=pho^(abs(t2-t1));
    end
end

% Generate the magnitudes (weights) for each Gaussian basis,
% then the observed functions for group 1
for i=1:nsize
    wgt1=normrnd(30,1.5,7,1);
    wgt1(3)=1.75*trnd(2)+30;
    wgt1(6)=normrnd(30,1);
    mean1=sum(bsxfun(@times,basis,wgt1));
    Y1(i,:)=mvnrnd(mean1,sigma^2*corr)+shift;
end

% Generate the magnitudes (weights) for each Gaussian basis,
% then the observed functions for group 2
for i=1:nsize
    wgt2=normrnd(30,1.5,7,1);
    wgt2(3)=normrnd(30,1);
    wgt2(6)=1.75*trnd(2)+30;
    mean2=sum(bsxfun(@times,basis,wgt2));
    Y2(i,:)=mvnrnd(mean2,sigma^2*corr)+shift;
end

% create a struct to store the data
model.Y=vertcat(Y1,Y2);
model.n=size(model.Y,1);
model.T=size(model.Y,2);

model.X=ones(model.n,2);
model.X((nsize+1):end,2)=-1;
model.p=2;

model.x0=x;

% get the DWT and inverse DWT matrix which will be used for model fitting
[D,wavespecs1]=DWT_rows(model.Y,wavespecs);
[~,wavespecs1]=wavelet_compress(D,wavespecs1,0.99975);

[Wt,~]=DWT_rows(eye(model.T),wavespecs);
W1t=Wt(:,wavespecs1.keep==1);

W1=IDWT_rows(eye(wavespecs1.K),wavespecs1);

% save the wavelet transform specifictions, DWT and iDWT matrices
wavespecs1.W1=sparse(W1);
wavespecs1.W1t=sparse(W1t);

wavespecs1.J=sum(wavespecs1.Kj~=0);
wavespecs1.Kj=wavespecs1.Kj(wavespecs1.Kj~=0);

model.wavespecs=wavespecs1;