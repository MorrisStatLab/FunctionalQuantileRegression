function [sens,spec,fp]=simbas_ROC_flag_contiguous(beta,MCMC_P,delta1,delta2,alpha)

% calculate SimBaS at each t based on posterior (or bootstrap) samples
SimBaS=jointband_simbas2(MCMC_P);

T=size(beta,2);

% assign positive/negative based on ground truth
% the locations with a magnitude between delta1 and delta2 are removed from
% ROC analysis
truth=NaN(1,T);

truth(abs(beta)<=delta1)=0;
truth((abs(beta)>delta1)&(abs(beta)<delta2))=-1;
truth(abs(beta)>=delta2)=1;

% calculate sensitivity and specificity
sens=NaN(size(alpha));
spec=NaN(size(alpha));

for i=1:size(alpha,2)
    flag=(SimBaS<=alpha(i));
    flag=flag_contiguous_sites(flag);
    sens(i)=sum(flag==1 & truth==1)/sum(truth==1);
    spec(i)=sum(flag==0 & truth==0)/sum(truth==0);
end

fp=1-spec;
