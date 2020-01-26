clear;

%% Figure S1
load('Y:/submission/realdata/ProteomicsData.mat');

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

Yraw_cancer=Yraw(cancer_status==1,:);
Yraw_normal=Yraw(cancer_status~=1,:);

Ylog_cancer=Y_6654(cancer_status==1,:);
Ylog_normal=Y_6654(cancer_status~=1,:);

idx1=(x_12096>=5e3 & x_12096<=8e3);
idx2=(x_6654_new>=5e3 & x_6654_new<=8e3);

g=figure;

for i=1:20
    subplot(5,4,i)
    plot(x_12096(idx1),Yraw_normal(i,idx1),'k-')
    xlim([5e3 8e3])
    xlabel('m/z(Daltons)')
    ylabel('Intensity')
end

savefig(g,'Y:/submission/FigureS1_raw_normal.fig')
close(g)


g=figure;

for i=1:20
    subplot(5,4,i)
    plot(x_6654_new(idx2),Ylog_normal(i,idx2),'k-')
    xlim([5e3 8e3])
    ylim([-2.05 7.5])
    xlabel('m/z(Daltons)')
    ylabel('Intensity')
end

savefig(g,'Y:/submission/FigureS1_log_normal.fig')
close(g)


g=figure;

for i=1:20
    subplot(5,4,i)
    plot(x_12096(idx1),Yraw_cancer(i,idx1),'k-')
    xlim([5e3 8e3])
    xlabel('m/z(Daltons)')
    ylabel('Intensity')
end

savefig(g,'Y:/submission/FigureS1_raw_cancer.fig')
close(g)


g=figure;

for i=1:20
    subplot(5,4,i)
    plot(x_6654_new(idx2),Ylog_cancer(i,idx2),'k-')
    xlim([5e3 8e3])
    ylim([-2.05 7.5])
    xlabel('m/z(Daltons)')
    ylabel('Intensity')
end

savefig(g,'Y:/submission/FigureS1_log_cancer.fig')
close(g)

clear;
