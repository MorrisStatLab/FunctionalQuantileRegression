clear; 

% Assume you are in the parent directory right now. May need to change path if necessary.

% load in the mass spectrometry data that have been preprocessed, including
% baseline correction, normalization, denoising and log2 transformation.
% The preprocessed dataset that we start from covers the spectral region
% [4000D, 19154D] and includes 6654 observations per spectrum.
load('./realdata/ProteomicsData.mat');  

%% block effect adjustment

% add the path of the code needed for block effect adjustment
addpath('./MLEfunctions');

status=zeros(size(Y_6654,1),1);

status(cancer_status==1,1)=1;
status(cancer_status~=1,1)=-1;

model1.Y=Y_6654(block==1,:);
model2.Y=Y_6654(block==2,:);
model3.Y=Y_6654(block==3,:);
model4.Y=Y_6654(block==4,:);

model1.X=horzcat(ones(sum(block==1),1),status(block==1));
model2.X=horzcat(ones(sum(block==2),1),status(block==2));
model3.X=horzcat(ones(sum(block==3),1),status(block==3));
model4.X=horzcat(ones(sum(block==4),1),status(block==4));

[model1.n,model1.p]=size(model1.X);
[model2.n,model2.p]=size(model2.X);
[model3.n,model3.p]=size(model3.X);
[model4.n,model4.p]=size(model4.X);

% Get the mean MLE estimate for each block
W1=GetW(model1,model1.Y);
[beta1_mle,]=sweep_simple_regress(W1,model1.n);

W2=GetW(model2,model2.Y);
[beta2_mle,]=sweep_simple_regress(W2,model2.n);

W3=GetW(model3,model3.Y);
[beta3_mle,]=sweep_simple_regress(W3,model3.n);

W4=GetW(model4,model4.Y);
[beta4_mle,]=sweep_simple_regress(W4,model4.n);


% plots to compare the mean MLE estimate and sample mean for each block
plot(x_6654_new,beta1_mle(1,:),'r-')
hold on
plot(x_6654_new,mean(model1.Y),'k-')
hold off

plot(x_6654_new,beta2_mle(1,:),'r-')
hold on
plot(x_6654_new,mean(model2.Y),'k-')
hold off

plot(x_6654_new,beta3_mle(1,:),'r-')
hold on
plot(x_6654_new,mean(model3.Y),'k-')
hold off

plot(x_6654_new,beta4_mle(1,:),'r-')
hold on
plot(x_6654_new,mean(model4.Y),'k-')
hold off


% calculate a grand mean that is aggregated over 4 blocks
beta_mean=(beta1_mle(1,:)+beta2_mle(1,:)+beta3_mle(1,:)+beta4_mle(1,:))/4;

% subtract the estimated block effect from the mass spectra for each block
diff1=beta1_mle(1,:)-beta_mean;
diff2=beta2_mle(1,:)-beta_mean;
diff3=beta3_mle(1,:)-beta_mean;
diff4=beta4_mle(1,:)-beta_mean;

model1.Ynew=bsxfun(@minus,model1.Y,diff1);
model2.Ynew=bsxfun(@minus,model2.Y,diff2);
model3.Ynew=bsxfun(@minus,model3.Y,diff3);
model4.Ynew=bsxfun(@minus,model4.Y,diff4);

% compare the sample mean of adjusted mass spectra across groups, which
% should be essentially the same
plot(x_6654_new,mean(model1.Ynew),'k-')
hold on
plot(x_6654_new,mean(model2.Ynew),'r-')
plot(x_6654_new,mean(model3.Ynew),'g-')
plot(x_6654_new,mean(model4.Ynew),'b-')
hold off


% aggregate the adjusted mass spectra from each block into one dataset
Y_6654_adjusted=vertcat(model1.Ynew,model2.Ynew,model3.Ynew,model4.Ynew);
X_adjusted=NaN(size(status,1),2);
X_adjusted(:,1)=1;
X_adjusted(:,2)=vertcat(model1.X(:,2),model2.X(:,2),model3.X(:,2),model4.X(:,2));

%% take the subset corresponding to [5000D,8000D] from the block effect 
% adjusted dataset and do basis transform

% add the path of the code for wavelet transform
addpath('./Code for dwt and idwt');

% the subset corresponding to [5000D,8000D]
[n,p]=size(X_adjusted);
T=size(Y_6654_adjusted,2);

idx=(x_6654_new>=5e3)&(x_6654_new<=8e3);

% create a struct to store the data
model.x0=x_6654_new(idx);
model.Y=Y_6654_adjusted(:,idx);
model.X=X_adjusted;

[model.n,model.p]=size(model.X);
model.T=size(model.Y,2);

% perform discrete wavelet transform (DWT) on the functional response Y
wavespecs.wtmode='per';
wavespecs.ndim=1;
wavespecs.boundary='per'; 
wavespecs.wavelet='db4';
wavespecs.compress=1;

[D,wavespecs1]=DWT_rows(model.Y,wavespecs);
[~,wavespecs1]=wavelet_compress(D,wavespecs1,0.99);

% get the DWT and inverse DWT matrix which will be used for model fitting
[Wt,~]=DWT_rows(eye(model.T),wavespecs);
W1t=Wt(:,wavespecs1.keep==1);

W1=IDWT_rows(eye(wavespecs1.K),wavespecs1);

% save the wavelet transform specifictions, DWT and iDWT matrices
wavespecs=wavespecs1;

wavespecs.W1=sparse(W1);
wavespecs.W1t=sparse(W1t);

model.wavespecs=wavespecs;

%% save the struct to mat file and txt files
save('./realdata/model.mat','model');

dlmwrite('./realdata/Y.txt',model.Y,'delimiter','\t','precision','%12.6e');
dlmwrite('./realdata/X.txt',model.X,'delimiter','\t','precision','%12.6e');
dlmwrite('./realdata/x0.txt',model.x0,'delimiter','\t','precision','%12.6e');

