function model=simdata(x,loc,scale,sigma,shift,nsize)
%% generate functional data of enormous sample size for the right skewed
% setting ii to obtain the ground truth for different quantiles
T=size(x,1);

% Two groups
Y1=NaN(nsize,T);
Y2=NaN(nsize,T);

% Gaussian basis
basis=NaN(size(loc,1),T);
for i=1:size(basis,1)
    basis(i,:)=normpdf(x',loc(i),scale(i));
end

% Generate the magnitudes (weights) for each Gaussian basis,
% then the observed functions of enormous sample size for group 1
for i=1:nsize
    wgt1=normrnd(30,1.5,7,1);
    wgt1(2)=30+1/gamrnd(1,1/0.35);
    wgt1(6)=normrnd(30,0.4)+0.5;
    mean1=sum(bsxfun(@times,basis,wgt1));
    Y1(i,:)=normrnd(mean1,sigma)+shift;
end

% Generate the magnitudes (weights) for each Gaussian basis,
% then the observed functions of enormous sample size for group 2
for i=1:nsize
    wgt2=normrnd(30,1.5,7,1);
    wgt2(2)=normrnd(30,0.4)+0.6;
    wgt2(6)=30+1/gamrnd(1,1/0.35);
    mean2=sum(bsxfun(@times,basis,wgt2));
    Y2(i,:)=normrnd(mean2,sigma)+shift;
end

% create a struct to store the data
model.Y=vertcat(Y1,Y2);
model.n=size(model.Y,1);
model.T=size(model.Y,2);

model.X=ones(model.n,2);
model.X((nsize+1):end,2)=-1;
model.p=2;

model.x0=x;