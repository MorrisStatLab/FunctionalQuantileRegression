function SimBaS = jointband_simbas2(MCMC_P)
  % function to calculate 
  % simultaneous band scores (SimBaS) based on MCMC samples contained in 
  % a B*T matrix 'MCMC_P'. For details on the simultaneous credible
  % band, refer to Crainiceanu et al 2007. At each location, 
  % the SimBaS is defined as the minimum significance level at which the 
  % simultaneous credible band excludes zero. 
   
   % Inputs
   %    'MCMC_P' - a B by T matrix containing MCMC samples. 
   %        B - number of MCMC iterations (MCMCspecs.B)
   %        T - number of function samples (number of columns in Y)
   
   % Outputs
   %   'SimBaS' - simultaneous band scores.  
    
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
   
   levels = 0.9995:-0.0005:0.0005; 
   cb1 = quantile(z_P,1-levels);
   SimBaS = 1*ones(1,T);
   for k=1:length(levels)
       temp_ind = (mean_P - cb1(k)*sd_P).*(mean_P + cb1(k)*sd_P);
       SimBaS(temp_ind>0) = levels(k);
   end
   

