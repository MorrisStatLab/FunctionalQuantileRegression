
%% Figure S1, S2 and S3 are plotted in R. Please refer to the R code file
% supplementary_figures_1_2_3.R.
% To reproduce Figure S1-S3, posterior samples for betat and sigma 
% obtained by Bayesian FQR based on the proteomics dataset at each quantile 
% level are needed. We do not provide them due to space considerations, 
% but users can refer to submission/realdata/realdata_modelfit.m to 
% generate them. They are also available upon request.

%% Supplementary Figure S4
% To reproduce Figure S4, posterior (or bootstrap) samples for betat 
% obtained by each approach based on the proteomics dataset at each
% quantile are needed. We do not provide them due to space considerations,
% but users can refer to submission/realdata/realdata_modelfit.m and 
% realdata_modelfit_bootstrap.R to generate the posterior (or bootstrap) samples. 
% They are also available upon request.

addpath('Y:/submission/Plotfunctions');

alpha=0.05;
delta=log2(1.5)/2;

x0=dlmread('Y:/submission/realdata/x0.txt');

qt=[0.1 0.25 0.5 0.75 0.9];

for i=1:5
    g=estimation_plots_with_simbas(x0,qt(i),alpha,delta);
    savefig(g,sprintf('Y:/submission/FigureS4_qt%d.fig',qt(i)*100))
    close(g)
end

rmpath('Y:/submission/Plotfunctions');