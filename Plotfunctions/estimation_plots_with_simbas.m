function g=estimation_plots_with_simbas(x0,qt,alpha,delta)

%% Bayesian QR
MCMC_P1=dlmread(sprintf('Y:/submission/realdata/output/QR/MCMC_betat_%d.txt',qt*100));
[SimBaS1,upr1,lwr1]=jointband_simbas(MCMC_P1,alpha);
upr1_pt=mean(MCMC_P1)+1.96*std(MCMC_P1);
lwr1_pt=mean(MCMC_P1)-1.96*std(MCMC_P1);

flag1=(SimBaS1<0.05) & (abs(mean(MCMC_P1))>delta);
flag1=flag_contiguous_sites(flag1);
idx1=find(flag1~=0);

%% FQR + Horseshoe Prior
MCMC_P2=dlmread(sprintf('Y:/submission/realdata/output/HS/MCMC_betat_%d.txt',qt*100));
[SimBaS2,upr2,lwr2]=jointband_simbas(MCMC_P2,alpha);
upr2_pt=mean(MCMC_P2)+1.96*std(MCMC_P2);
lwr2_pt=mean(MCMC_P2)-1.96*std(MCMC_P2);

flag2=(SimBaS2<0.05) & (abs(mean(MCMC_P2))>delta);
flag2=flag_contiguous_sites(flag2);
idx2=find(flag2~=0);

%% Bootstrap-based QR 
MCMC_P3=dlmread(sprintf('Y:/submission/realdata/output/freqQR/raw/BS_betat_%d.txt',qt*100));
[SimBaS3,upr3,lwr3]=jointband_simbas(MCMC_P3,alpha);
upr3_pt=mean(MCMC_P3)+1.96*std(MCMC_P3);
lwr3_pt=mean(MCMC_P3)-1.96*std(MCMC_P3);

flag3=(SimBaS3<0.05) & (abs(mean(MCMC_P3))>delta);
flag3=flag_contiguous_sites(flag3);
idx3=find(flag3~=0);

%% Bootstrap-based QR with spline smoothing
MCMC_P4=dlmread(sprintf('Y:/submission/realdata/output/freqQR/splines/BS_betat_%d.txt',qt*100));
[SimBaS4,upr4,lwr4]=jointband_simbas(MCMC_P4,alpha);
upr4_pt=mean(MCMC_P4)+1.96*std(MCMC_P4);
lwr4_pt=mean(MCMC_P4)-1.96*std(MCMC_P4);

flag4=(SimBaS4<0.05) & (abs(mean(MCMC_P4))>delta);
flag4=flag_contiguous_sites(flag4);
idx4=find(flag4~=0);

%% Bootstrap-based QR with wavelet denoising
MCMC_P5=dlmread(sprintf('Y:/submission/realdata/output/freqQR/wavelets/BS_betat_%d.txt',qt*100));
[SimBaS5,upr5,lwr5]=jointband_simbas(MCMC_P5,alpha);
upr5_pt=mean(MCMC_P5)+1.96*std(MCMC_P5);
lwr5_pt=mean(MCMC_P5)-1.96*std(MCMC_P5);

flag5=(SimBaS5<0.05) & (abs(mean(MCMC_P5))>delta);
flag5=flag_contiguous_sites(flag5);
idx5=find(flag5~=0);

%% plot
g=figure;

subplot(3,2,1)
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
title('Bayesian QR')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,2,3)
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
title('Bayesian FQR with wavelet transform and horseshoe prior')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,2,2)
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
title('Bootstrap-based QR')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,2,6)
x=[x0',fliplr(x0')];               
y=[lwr4,fliplr(upr4)]; 
ypt=[lwr4_pt,fliplr(upr4_pt)];
h=fill(x,y,[0.8 0.8 0.8],x,ypt,[0.5 0.5 0.5]); 
set(h,'EdgeColor','none');
hold on
plot(x0,mean(MCMC_P4),'k-')
hline = refline([0 delta]);
hline.Color = 'g';
hline = refline([0 -delta]);
hline.Color = 'g';
hline = refline([0 0]);
hline.Color = 'b';
xlim([5e3 8e3])
ylim([-1 1.5])
for i=1:size(idx4,2)
    plot([x0(idx4(i)) x0(idx4(i))],[-1 -0.92],'r-');
end
title('Bootstrap-based QR with spline smoothing')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

subplot(3,2,4)
x=[x0',fliplr(x0')];               
y=[lwr5,fliplr(upr5)]; 
ypt=[lwr5_pt,fliplr(upr5_pt)];
h=fill(x,y,[0.8 0.8 0.8],x,ypt,[0.5 0.5 0.5]); 
set(h,'EdgeColor','none');
hold on
plot(x0,mean(MCMC_P5),'k-')
hline = refline([0 delta]);
hline.Color = 'g';
hline = refline([0 -delta]);
hline.Color = 'g';
hline = refline([0 0]);
hline.Color = 'b';
xlim([5e3 8e3])
ylim([-1 1.5])
for i=1:size(idx5,2)
    plot([x0(idx5(i)) x0(idx5(i))],[-1 -0.92],'r-');
end
title('Bootstrap-based QR with wavelet denoising')
xlabel('m/z(Daltons)')
ylabel('Intensity (log_2 scale)')
hold off

orient(g,'tall')


