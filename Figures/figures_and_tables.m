clear;

%% Assume you are in the subfolder "Figures/" now. May need to change path if necessary.

%% Figure 1
load('../realdata/ProteomicsData.mat');

% x_12096 and Yraw refer to the spectral locations (in Daltons) and  
% functional responses of the raw (i.e. unpreprocessed) spectrometry   
% dataset that covers the spectral region [4000D, 40000D] and 
% includes 12096 observations per spectrum.

% Y_6654, the starting point of the dataset for which we performed FQR,
% was preprocessed from Yraw including baseline correction, 
% normalization, denoising and log2 transformation.
% x_6654_new and Y_6654 include only part of the spectral region surveyed
% in x_12096 and Yraw, i.e. [4000D, 19154D] with 6654 observations 
% per spectrum.

idx1=(x_12096>=5e3 & x_12096<=8e3);
idx2=(x_6654_new>=5e3 & x_6654_new<=8e3);

g=figure;

subplot(2,2,1)
plot(x_12096(idx1),Yraw(1,idx1),'k-')
xlim([5e3 8e3])
title('(a)')
xlabel('m/z(Daltons)')
ylabel('Intensity')

subplot(2,2,2)
plot(x_6654_new(idx2),Y_6654(1,idx2),'k-')
xlim([5e3 8e3])
title('(b)')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')

subplot(2,2,3)
plot(x_12096(idx1),Yraw(2,idx1),'k-')
xlim([5e3 8e3])
title('(c)')
xlabel('m/z(Daltons)')
ylabel('Intensity')

subplot(2,2,4)
plot(x_6654_new(idx2),Y_6654(2,idx2),'k-')
xlim([5e3 8e3])
title('(d)')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')

savefig(g,'Y:/submission/Figure1.fig')
close(g)

clear;

%% Figure 2
x0=dlmread('../realdata/x0.txt');
X=dlmread('../realdata/X.txt');
Y=dlmread('../realdata/Y.txt');

[n,p]=size(X);
T=size(Y,2);

idx=find(x0>5764 & x0<5765);
y1=Y(X(:,2)==1,idx);
y2=Y(X(:,2)==-1,idx);

g=figure;

subplot(2,2,[1,2])
plot(x0,mean(Y(X(:,2)==1,:))-mean(Y(X(:,2)==-1,:)),'k-')
hold on
plot(x0,quantile(Y(X(:,2)==1,:),0.1)-quantile(Y(X(:,2)==-1,:),0.1),'b-')
plot(x0,quantile(Y(X(:,2)==1,:),0.25)-quantile(Y(X(:,2)==-1,:),0.25),'g-')
plot(x0,quantile(Y(X(:,2)==1,:),0.5)-quantile(Y(X(:,2)==-1,:),0.5),'y-')
plot(x0,quantile(Y(X(:,2)==1,:),0.75)-quantile(Y(X(:,2)==-1,:),0.75),'r-')
plot(x0,quantile(Y(X(:,2)==1,:),0.9)-quantile(Y(X(:,2)==-1,:),0.9),'m-')
xlim([5e3 8e3])
ylim([-1.5 2.5])
plot([5764.1 5764.1],[-2 3],'color',[0.8 0.8 0.8])
title('(a)')
legend('mean','0.1','0.25','0.5','0.75','0.9')
legend('boxoff')
hold off

subplot(2,2,3)
h1=histfit(y1,15,'kernel');
set(h1(1),'EdgeColor','none');
h1(1).FaceColor = [0.8 0.8 0.8];
h1(2).Color = [0.5 0.5 0.5];
hold on
h1a=plot([mean(y1) mean(y1)],[0 22],'k-');
h1b=plot([quantile(y1,0.1) quantile(y1,0.1)],[0 14],'b-');
h1c=plot([quantile(y1,0.25) quantile(y1,0.25)],[0 25],'g-');
h1d=plot([quantile(y1,0.5) quantile(y1,0.5)],[0 27],'y-');
h1e=plot([quantile(y1,0.75) quantile(y1,0.75)],[0 11],'r-');
h1f=plot([quantile(y1,0.9) quantile(y1,0.9)],[0 7],'m-');
xlim([min(Y(:,idx))-0.5 max(Y(:,idx))+0.5])
title('(b)')
legend([h1a h1b h1c h1d h1e h1f],{'mean','0.1','0.25','0.5','0.75','0.9'})
legend('boxoff')
hold off

subplot(2,2,4)
h2=histfit(y2,12,'kernel');
set(h2(1),'EdgeColor','none');
h2(1).FaceColor = [0.8 0.8 0.8];
h2(2).Color = [0.5 0.5 0.5];
hold on
h2a=plot([mean(y2) mean(y2)],[0 33],'k-');
h2b=plot([quantile(y2,0.1) quantile(y2,0.1)],[0 9],'b-');
h2c=plot([quantile(y2,0.25) quantile(y2,0.25)],[0 33],'g-');
h2d=plot([quantile(y2,0.5) quantile(y2,0.5)],[0 33],'y-');
h2e=plot([quantile(y2,0.75) quantile(y2,0.75)],[0 42],'r-');
h2f=plot([quantile(y2,0.9) quantile(y2,0.9)],[0 13],'m-');
xlim([min(Y(:,idx))-0.5 max(Y(:,idx))+0.5])
ylim([0 45])
title('(c)' )
legend([h2a h2b h2c h2d h2e h2f],{'mean','0.1','0.25','0.5','0.75','0.9'})
legend('boxoff')
hold off

savefig(g,'./Figure2.fig')
close(g)

clear;

%% Figure 3
x0=dlmread('../simulations1/x0.txt');

beta2=dlmread('../simulations1/beta2.txt');

g=figure;

subplot(2,2,1)
plot(x0,beta2(1,:),'b-')
hold on
plot(x0,beta2(2,:),'g-')
plot(x0,beta2(3,:),'y-')
plot(x0,beta2(4,:),'r-')
plot(x0,beta2(5,:),'m-')
xlim([0 15])
ylim([-2.7 2.7])
title('(a)')
legend('0.1','0.2','0.5','0.8','0.9')
legend('boxoff')
hold off

x0=dlmread('../simulations2/x0.txt');

beta2=dlmread('../simulations2/beta2.txt');

subplot(2,2,2)
plot(x0,beta2(1,:),'b-')
hold on
plot(x0,beta2(2,:),'g-')
plot(x0,beta2(3,:),'y-')
plot(x0,beta2(4,:),'r-')
plot(x0,beta2(5,:),'m-')
xlim([0 15])
ylim([-2.7 2.7])
title('(b)')
legend('0.1','0.2','0.5','0.8','0.9')
legend('boxoff')
hold off


subplot(2,2,3)

load('../simulations1/data/model1.mat')

plot(model.x0,model.Y(1,:),'r-.')
hold on
plot(model.x0,model.Y(324,:),'k-')
title('(c)')
xlabel('\it t','Interpreter','tex')
ylabel('\it y','Interpreter','tex')
legend('Group 1','Group 2','location','Northeast')
legend('boxoff')
hold off


subplot(2,2,4)

load('../simulations2/data/model1.mat')

plot(model.x0,model.Y(146,:),'r-.')
hold on
plot(model.x0,model.Y(201,:),'k-')
title('(d)')
xlabel('\it t','Interpreter','tex')
ylabel('\it y','Interpreter','tex')
legend('Group 1','Group 2','location','Northeast')
legend('boxoff')
hold off

savefig(g,'./Figure3.fig')
close(g)

clear;

%% Figure 4
addpath('../Plotfunctions');

alpha=0.05;
delta=log2(1.5)/2;

x0=dlmread('../realdata/x0.txt');
Y=dlmread('../realdata/Y.txt');

qt=0.9;

% Bayesian FQR
MCMC_P1=dlmread(sprintf('../realdata/output/HS/MCMC_betat_%d.txt',qt*100));
[SimBaS1,upr1,lwr1]=jointband_simbas(MCMC_P1,alpha);
upr1_pt=mean(MCMC_P1)+1.96*std(MCMC_P1);
lwr1_pt=mean(MCMC_P1)-1.96*std(MCMC_P1);

flag1=(SimBaS1<0.05) & (abs(mean(MCMC_P1))>delta);
flag1=flag_contiguous_sites(flag1);
idx1=find(flag1~=0);

% QR with wavelet denoising
MCMC_P2=dlmread(sprintf('../realdata/output/freqQR/wavelets/BS_betat_%d.txt',qt*100));
[SimBaS2,upr2,lwr2]=jointband_simbas(MCMC_P2,alpha);
upr2_pt=mean(MCMC_P2)+1.96*std(MCMC_P2);
lwr2_pt=mean(MCMC_P2)-1.96*std(MCMC_P2);

flag2=(SimBaS2<0.05) & (abs(mean(MCMC_P2))>delta);
flag2=flag_contiguous_sites(flag2);
idx2=find(flag2~=0);

% GFMM mean regression 
MCMC_P3=dlmread('../realdata/output/GFMM/MCMC_betat.txt');
[SimBaS3,upr3,lwr3]=jointband_simbas(MCMC_P3,alpha);
upr3_pt=mean(MCMC_P3)+1.96*std(MCMC_P3);
lwr3_pt=mean(MCMC_P3)-1.96*std(MCMC_P3);

flag3=(SimBaS3<0.05) & (abs(mean(MCMC_P3))>delta);
flag3=flag_contiguous_sites(flag3);
idx3=find(flag3~=0);

% plot
g=figure;

subplot(3,1,1)
x=[x0',fliplr(x0')];               
y=[lwr1,fliplr(upr1)]; 
ypt=[lwr1_pt,fliplr(upr1_pt)];
h=fill(x,y,[0.8 0.8 0.8],x,ypt,[0.5 0.5 0.5]); 
set(h,'EdgeColor','none');
hold on
plot(x0,mean(MCMC_P1),'k-')
hline = refline([0 delta]);
hline.Color = 'g';
hline = refline([0 -delta]);
hline.Color = 'g';
hline = refline([0 0]);
hline.Color = 'b';
xlim([5e3 8e3])
ylim([-1 1.5])
for i=1:size(idx1,2)
    plot([x0(idx1(i)) x0(idx1(i))],[-1 -0.92],'r-');
end
title('(a)')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,1,2)
x=[x0',fliplr(x0')];               
y=[lwr2,fliplr(upr2)]; 
ypt=[lwr2_pt,fliplr(upr2_pt)];
h=fill(x,y,[0.8 0.8 0.8],x,ypt,[0.5 0.5 0.5]); 
set(h,'EdgeColor','none');
hold on
plot(x0,mean(MCMC_P2),'k-')
hline = refline([0 delta]);
hline.Color = 'g';
hline = refline([0 -delta]);
hline.Color = 'g';
hline = refline([0 0]);
hline.Color = 'b';
xlim([5e3 8e3])
ylim([-1 1.5])
for i=1:size(idx2,2)
    plot([x0(idx2(i)) x0(idx2(i))],[-1 -0.92],'r-');
end
title('(b)')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,1,3)
x=[x0',fliplr(x0')];               
y=[lwr3,fliplr(upr3)]; 
ypt=[lwr3_pt,fliplr(upr3_pt)];
h=fill(x,y,[0.8 0.8 0.8],x,ypt,[0.5 0.5 0.5]); 
set(h,'EdgeColor','none');
hold on
plot(x0,mean(MCMC_P3),'k-')
hline = refline([0 delta]);
hline.Color = 'g';
hline = refline([0 -delta]);
hline.Color = 'g';
hline = refline([0 0]);
hline.Color = 'b';
xlim([5e3 8e3])
ylim([-1 1.5])
for i=1:size(idx3,2)
    plot([x0(idx3(i)) x0(idx3(i))],[-1 -0.92],'r-');
end
title('(c)')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

orient(g,'tall')

savefig(g,'./Figure4.fig')
close(g)

clear;
rmpath('../Plotfunctions');
