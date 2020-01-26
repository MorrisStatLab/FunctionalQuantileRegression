function model=simdata2(x,loc,scale,sigma,pho,shift,nsize,wavespecs)
%%
model.n=2*nsize;
model.X=ones(model.n,3);
model.X(1:nsize,2)=-1;
model.X(:,3)=normrnd(0,1,[model.n 1]);
model.p=3;

T=size(x,1);

Y1=NaN(nsize,T);
Y2=NaN(nsize,T);

basis=NaN(size(loc,1),T);
for i=1:size(basis,1)
    basis(i,:)=normpdf(x',loc(i),scale(i));
end

corr=NaN(T,T);
for t1=1:T
    for t2=1:T
        corr(t1,t2)=pho^(abs(t2-t1));
    end
end

for i=1:nsize
    wgt1=normrnd(30,1,4,1);
    wgt1(1)=1.75*trnd(2)+30;
    wgt1(3)=normrnd(30,0.4)+0.5;
    wgt1(4)=wgt1(4)+model.X(i,3);
    mean1=sum(bsxfun(@times,basis,wgt1));
    Y1(i,:)=mvnrnd(mean1,sigma^2*corr)+shift;
end

for i=1:nsize    
    wgt2=normrnd(30,1,4,1);
    wgt2(3)=30+1/gamrnd(1,1/0.35);
    wgt2(4)=wgt2(4)+model.X(i+nsize,3);
    mean2=sum(bsxfun(@times,basis,wgt2));
    Y2(i,:)=mvnrnd(mean2,sigma^2*corr)+shift;
end

model.Y=vertcat(Y1,Y2);
model.T=size(model.Y,2);

model.x0=x;

%% wavelet compression rate=0.9999
[D,wavespecs1]=DWT_rows(model.Y,wavespecs);
[~,wavespecs1]=wavelet_compress(D,wavespecs1,0.9999);

% get the DWT matrix
[Wt,~]=DWT_rows(eye(model.T),wavespecs);
W1t=Wt(:,wavespecs1.keep==1);

W1=IDWT_rows(eye(wavespecs1.K),wavespecs1);

% save the wavespecs to model
wavespecs1.W1=sparse(W1);
wavespecs1.W1t=sparse(W1t);

wavespecs1.J=sum(wavespecs1.Kj~=0);
wavespecs1.Kj=wavespecs1.Kj(wavespecs1.Kj~=0);

model.wavespecs=wavespecs1;