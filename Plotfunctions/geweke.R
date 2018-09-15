geweke <- function(qt)
{
   library(coda)

   pvalue <- list(NULL)


   # pvalues from Geweke test on posterior samples of Beta_2(t) at all t
   MCMC_betat <- read.table(paste0("realdata/output/HS/MCMC_betat_",qt*100,".txt"),header=FALSE,sep='\t')

   betat <- rep(NA,ncol(MCMC_betat))

   for(t in 1:ncol(MCMC_betat))
   {
      betat[t] <- geweke.diag(MCMC_betat[,t],frac1=0.45,frac2=0.45)[[1]]
   }

   pvalue[[1]] <- unlist(lapply(betat,function(x){2*pnorm(-abs(x))}))


   # pvalues from Geweke test on posterior samples of sigma(t) at all t
   MCMC_sigma <- read.table(paste0("realdata/output/HS/MCMC_sigma_",qt*100,".txt"),header=FALSE,sep='\t')

   sigma <- rep(NA,ncol(MCMC_sigma))

   for(t in 1:ncol(MCMC_sigma))
   {
      sigma[t] <- geweke.diag(MCMC_sigma[,t],frac1=0.45,frac2=0.45)[[1]]
   }

   pvalue[[2]] <- unlist(lapply(sigma,function(x){2*pnorm(-abs(x))}))


   names(pvalue) <- c("betat","sigma")

   save(pvalue,file=paste0("realdata/output/geweke_pvalues_Bayesian_FQR_qt",qt*100,".Rdata"))


   # calculate the proportion of pvalues below 0.05 and 0.1
   pvalues <- c(pvalue[[1]],pvalue[[2]])

   prop <- c(mean(abs(pvalues)<0.05),mean(abs(pvalues)<0.1))

   return(prop)
}
