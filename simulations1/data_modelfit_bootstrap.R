setwd("Y:/submission/simulations1")

X <- read.table("X.txt",header=FALSE,sep="\t")

time <- read.table("x0.txt",header=FALSE,sep="\t")

time <- as.numeric(time[,1])

qt <- c(0.1,0.2,0.5,0.8,0.9)

# Bootstrap-based FQR relies on the "quantreg" package, so please install it before running the code

source("Y:/submission/freqQR.R")


setwd("data")

# for space considerations, we do not provide the bootstrap samples. Users can run the script below to generate them, and they are also available upon request.
# for each replicate, it takes about 70 minutes to fit the model (for all quantiles) on a 64-bit operating system with 2 processors and 256GB RAM

for(i in 1:100)
{
   Y <- read.table(paste0("model",i,".txt"),header=FALSE,sep="\t")   
   result <- freqQR(Y=Y,X=X,time=time,qt=qt,B=2000)

   write.table(result$coef1,paste0("Y:/submission/simulations1/output/freqQR/raw/BS_betat_10_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef2,paste0("Y:/submission/simulations1/output/freqQR/raw/BS_betat_20_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef4,paste0("Y:/submission/simulations1/output/freqQR/raw/BS_betat_80_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef5,paste0("Y:/submission/simulations1/output/freqQR/raw/BS_betat_90_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

   write.table(result$coef1.ss,paste0("Y:/submission/simulations1/output/freqQR/splines/BS_betat_10_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef2.ss,paste0("Y:/submission/simulations1/output/freqQR/splines/BS_betat_20_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef4.ss,paste0("Y:/submission/simulations1/output/freqQR/splines/BS_betat_80_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$coef5.ss,paste0("Y:/submission/simulations1/output/freqQR/splines/BS_betat_90_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}
