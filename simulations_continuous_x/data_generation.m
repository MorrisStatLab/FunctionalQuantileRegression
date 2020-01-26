clear;

% add the path of the code.
addpath('Y:/submission/simulations_continuous_x'); 
addpath('Y:/submission/Code for dwt and idwt');

%% obtain the ground truth
% create the discrete grid on which functional data are observed
x=linspace(0,9,301)';
dlmwrite(strcat('Y:/submission/simulations_continuous_x','/x0.txt'),x,'delimiter','\t','precision','%12.6e');

% Gaussian kernels as described in Table1 in the manuscript
kernels=NaN(4,2);
kernels(:,1)=([1 3.25 5.5 8])';
kernels(:,2)=0.18;

dlmwrite(strcat('Y:/submission/simulations_continuous_x','/kernels.txt'),kernels,'delimiter','\t','precision','%12.6e');

sigma=3;

% generate functional data of enormous sample size 
% to obtain the ground truth of beta1 and beta2 for different quantiles
model=simdata(x,kernels(:,1),kernels(:,2),sigma,15,2000000);

quantiles1=quantile(model.Y(1:2000000,:),[0.1 0.2 0.5 0.8 0.9]);
quantiles2=quantile(model.Y(2000001:end,:),[0.1 0.2 0.5 0.8 0.9]);

beta1=(quantiles1+quantiles2)/2;
beta2=(quantiles2-quantiles1)/2;

% beta3 is available in closed form
beta3=normpdf(x',kernels(4,1),kernels(4,2));

% save groud truth of beta for different quantiles
dlmwrite(strcat('Y:/submission/simulations_continuous_x','/beta1.txt'),beta1,'delimiter','\t','precision','%12.6e');
dlmwrite(strcat('Y:/submission/simulations_continuous_x','/beta2.txt'),beta2,'delimiter','\t','precision','%12.6e');
dlmwrite(strcat('Y:/submission/simulations_continuous_x','/beta3.txt'),beta3,'delimiter','\t','precision','%12.6e');

%% discrete wavelet transform (DWT) specifications
wavespecs.wtmode='per';
wavespecs.ndim=1;
wavespecs.boundary='per'; 
wavespecs.wavelet='db4';
wavespecs.compress=1;

%% generate 100 replicates 

% for space considerations, we only provide 1 replicate simulation dataset.
% Users can run the code below to generate all replicates.
% All replicates are also available upon request.

for i=1:100
    model=simdata2(x,kernels(:,1),kernels(:,2),sigma,0.5,15,200,wavespecs);
    save(sprintf('Y:/submission/simulations_continuous_x/data/model%d.mat',i),'model');
    dlmwrite(sprintf('Y:/submission/simulations_continuous_x/data/Y_rep_%d.txt',i),model.Y,'delimiter','\t','precision','%12.6e');
    dlmwrite(sprintf('Y:/submission/simulations_continuous_x/data/X_rep_%d.txt',i),model.X,'delimiter','\t','precision','%12.6e');
end
