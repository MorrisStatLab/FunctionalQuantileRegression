traceplot <- function(qt,idx)
{

   # trace plots of posterior samples of Beta_2(t) at selected t
   MCMC_betat <- read.table(paste0("realdata/output/HS/MCMC_betat_",qt*100,".txt"),header=FALSE,sep="\t")

   B <- nrow(MCMC_betat)

   png(paste0("traceplot_betat_",qt*100,".png"),width = 960, height = 960)

   par(mfrow=c(4,4))

   par(mar=c(2,2.5,1.5,1))

   for(i in 1:16)
   {
      plot(1:B,MCMC_betat[,idx[i]],type="l",xaxt="n",xlab="iteration",ylab="",main=paste0("l=",idx[i]))
      axis(1, at=c(0,2000),labels=c(0,2000))
   }

   dev.off()


   # trace plots of posterior samples of sigma(t) at selected t
   MCMC_sigma <- read.table(paste0("realdata/output/HS/MCMC_sigma_",qt*100,".txt"),header=FALSE,sep="\t")

   png(paste0("traceplot_sigma_",qt*100,".png"),width = 960, height = 960)

   par(mfrow=c(4,4))

   par(mar=c(2,2.5,1.5,1))

   for(i in 1:16)
   {
      plot(1:B,MCMC_sigma[,idx[i]],type="l",xaxt="n",xlab="iteration",ylab="",main=paste0("l=",idx[i]))
      axis(1, at=c(0,2000),labels=c(0,2000))
   }

   dev.off()
}
