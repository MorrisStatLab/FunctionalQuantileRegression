# call the function "freqQR" to perform FQR on the real data
# Bootstrap-based FQR relies on the "quantreg" package, so please install it before running the code
source("Y:/submission/realdata/freqQR.R")

setwd("Y:/submission/realdata")

X <- read.table("X.txt",header=FALSE,sep="\t")

Y <- read.table("Y.txt",header=FALSE,sep="\t")   

time <- read.table("x0.txt",header=FALSE,sep="\t")

time <- as.numeric(time[,1])

qt <- c(0.1,0.25,0.5,0.75,0.9)


# for space considerations, we do not provide the bootstrap samples. Users can run the script below to generate them, and they are also available upon request.
# it takes about 12.5 hours to fit the model (for all 5 quantile levels) on a 64-bit operating system with 2 processors and 32GB RAM

result <- freqQR(Y=Y,X=X,time=time,qt=qt,B=2000)


# save the bootstrap samples obtained using naive Frequentist FQR 
write.table(result$coef1,"Y:/submission/realdata/output/freqQR/raw/BS_betat_10.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef2,"Y:/submission/realdata/output/freqQR/raw/BS_betat_25.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef3,"Y:/submission/realdata/output/freqQR/raw/BS_betat_50.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef4,"Y:/submission/realdata/output/freqQR/raw/BS_betat_75.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef5,"Y:/submission/realdata/output/freqQR/raw/BS_betat_90.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

# save the bootstrap samples obtained using two-step Frequentist FQR with spline smoothing
write.table(result$coef1.ss,"Y:/submission/realdata/output/freqQR/splines/BS_betat_10.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef2.ss,"Y:/submission/realdata/output/freqQR/splines/BS_betat_25.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef3.ss,"Y:/submission/realdata/output/freqQR/splines/BS_betat_50.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef4.ss,"Y:/submission/realdata/output/freqQR/splines/BS_betat_75.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(result$coef5.ss,"Y:/submission/realdata/output/freqQR/splines/BS_betat_90.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

