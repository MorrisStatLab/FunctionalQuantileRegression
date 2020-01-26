# Bootstrap-based FQR relies on the "quantreg" package, so please install it before running the code
source("Y:/submission/simulations_continuous_x/freqQR.R")

setwd("Y:/submission/simulations_continuous_x")

qt <- c(0.1,0.2,0.5,0.8,0.9)

 
# Users can run the script below to generate the bootstrap samples of beta2 and beta3 based on naive pointwise QR, 
# which can be post-smoothed to obtain the bootstrap samples of beta2 and beta3 by spline or wavelet smoothing, etc.
# For each replicate, it takes about 60 minutes to fit the model (for all 5 quantiles) on a 64-bit operating system with 2 processors and 32GB RAM.

for(i in 1:100)
{ 
   X <- read.table(paste0("data/X_rep_",i,".txt"),header=FALSE,sep="\t")
   Y <- read.table(paste0("data/Y_rep_",i,".txt"),header=FALSE,sep="\t")   
   result <- freqQR(Y=Y,X=X,qt=qt,B=2000)

   write.table(result$beta2.1,paste0("output/freqQR/raw/beta2_10_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta2.2,paste0("output/freqQR/raw/beta2_20_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta2.3,paste0("output/freqQR/raw/beta2_50_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta2.4,paste0("output/freqQR/raw/beta2_80_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta2.5,paste0("output/freqQR/raw/beta2_90_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

   write.table(result$beta3.1,paste0("output/freqQR/raw/beta3_10_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta3.2,paste0("output/freqQR/raw/beta3_20_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta3.3,paste0("output/freqQR/raw/beta3_50_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta3.4,paste0("output/freqQR/raw/beta3_80_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
   write.table(result$beta3.5,paste0("output/freqQR/raw/beta3_90_rep_",i,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

   

