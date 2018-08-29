clear;

% add the path of the code.
addpath('Y:/submission/simulations2'); 
addpath('Y:/submission/Code for dwt and idwt');

%% ground truth
% the discrete grid on which functional data are observed
x=linspace(0,15,501)';
dlmwrite(strcat('Y:/submission/simulations2','/x0.txt'),x,'delimiter','\t','precision','%12.6e');

% Gaussian kernels as described in Table1 in the manuscript
kernels=dlmread(strcat('Y:/submission/simulations2','/kernels.txt'));

sigma=3;

% generate functional data of enormous sample size for the right skewed 
% setting ii to obtain the ground truth for different quantiles
model=simdata(x,kernels(:,1),kernels(:,2),sigma,15,1000000);

quantiles1=quantile(model.Y(1:1000000,:),[0.1 0.2 0.5 0.8 0.9]);
quantiles2=quantile(model.Y(1000001:2000000,:),[0.1 0.2 0.5 0.8 0.9]);

beta1=(quantiles1+quantiles2)/2;
beta2=(quantiles1-quantiles2)/2;

% save groud truth of Beta for different quantiles
dlmwrite(strcat('Y:/submission/simulations2','/beta1.txt'),beta1,'delimiter','\t','precision','%12.6e');
dlmwrite(strcat('Y:/submission/simulations2','/beta2.txt'),beta2,'delimiter','\t','precision','%12.6e');

% plot the ground truth for different quantiles
figure;
plot(x,beta1(1,:),'b-')
hold on
plot(x,beta1(2,:),'y-')
plot(x,beta1(3,:),'k-')
plot(x,beta1(4,:),'g-')
plot(x,beta1(5,:),'m-')
xlim([0 15])
title('true intercept functions for different quantiles')
legend('0.1','0.2','0.5','0.8','0.9')
hold off

figure;
plot(x,beta2(1,:),'b-')
hold on
plot(x,beta2(2,:),'y-')
plot(x,beta2(3,:),'k-')
plot(x,beta2(4,:),'g-')
plot(x,beta2(5,:),'m-')
xlim([0 15])
title('true difference functions for different quantiles')
legend('0.1','0.2','0.5','0.8','0.9')
hold off

%% specify design matrix X
X=ones(400,2);
X(201:end,2)=-1;
dlmwrite('Y:/submission/simulations2/X.txt',X,'delimiter','\t','precision','%12.6e');

%% discrete wavelet transform (DWT) specifications
wavespecs.wtmode='per';
wavespecs.ndim=1;
wavespecs.boundary='per'; 
wavespecs.wavelet='db4';
wavespecs.compress=1;

%% generate 100 replicates and save the data

% for space considerations, we only provide 1 replicate simulation dataset.
% Users can run the code below to generate all replicates.
% All replicates are also available upon request.
for i=1:100
    model=simdata2(x,kernels(:,1),kernels(:,2),sigma,0.8,15,200,wavespecs);
    save(sprintf('Y:/submission/simulations2/data/model%d.mat',i),'model');
    dlmwrite(sprintf('Y:/submission/simulations2/data/model%d.txt',i),model.Y,'delimiter','\t','precision','%12.6e');
end


