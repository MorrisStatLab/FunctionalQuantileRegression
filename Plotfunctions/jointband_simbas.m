function [SimBaS,upper_CI, lower_CI] = jointband_simbas(MCMC_P,alpha)
  % function to calculate (1-alpha) simultaneous credible band and 
  % simultaneous band scores (SimBaS) based on MCMC samples contained in 
  % a B*T matrix 'MCMC_P'. For details on the simultaneous credible
  % band, refer to Crainiceanu et al 2007. At each location, 
  % the SimBaS is defined as the minimum significance level at which the 
  % simultaneous credible band excludes zero. 
   
   % Inputs
   %    'MCMC_P' - a B by T matrix containing MCMC samples. 
   %        B - number of MCMC iterations (MCMCspecs.B)
   %        T - number of function samples (number of columns in Y)
   %    'alpha' - a row vector containing significance levels at which 
   %              simultaneous credible bands are to be returned. 
   
   % Outputs
   %   'SimBaS' - simultaneous band scores. SimBaS are truncated at 
   %              0.001 for reducing computational burden.
   %   'upper_CI' - a length(alpha) by T matrix containing the upper bounds 
   %        of the simultaneous credible bands. The first row corresponds to 
   %        the first level in alpha.       
   %   'lower_CI' - a length(alpha) by T matrix containing the lower bounds 
   %        of the simultaneous credible bands. 
    
   [B,T] = size(MCMC_P);  
   sd_P = NaN(1,T);
   for i=1:T
       sd_P(i) = std(MCMC_P(:,i));
   end
   mean_P = mean(MCMC_P);
   
   z_P = NaN(1,B);
   for j=1:B
       z_P(j) = max(abs((MCMC_P(j,:)-mean_P)./sd_P));
   end
   
   levels = [0.99:-0.01:0.01,0.009:-0.001:0.001]; 
   cb1 = quantile(z_P,1-levels);
   SimBaS = 1*ones(1,T);
   for k=1:length(levels)
       temp_ind = (mean_P - cb1(k)*sd_P).*(mean_P + cb1(k)*sd_P);
       SimBaS(temp_ind>0) = levels(k);
   end
   
   m = length(alpha);
   upper_CI = zeros(m,T);
   lower_CI = zeros(m,T);
   cb2 = quantile(z_P,1-alpha);
   for j=1:m
       upper_CI(j,:) = mean_P + cb2(j)*sd_P;
       lower_CI(j,:) = mean_P - cb2(j)*sd_P;
   end    
end

