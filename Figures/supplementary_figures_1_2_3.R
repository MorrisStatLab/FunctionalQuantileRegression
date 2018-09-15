# To reproduce Figure S1-S3, posterior samples for betat and sigma obtained by
# Bayesian FQR based on the proteomics dataset at each quantile are needed. We do not
# provide them due to space considerations, but they are available upon request,
# and users can also refer to submission/realdata/realdata_modelfit.m to generate them.


# Figure S1 and S2

setwd("./")

source("../Plotfunctions/traceplot.R")

idx <- seq(100,1600,100)

qt <- c(0.1,0.25,0.5,0.75,0.9)

for(i in 1:5)
{
   traceplot(qt[i],idx)
}


# Figure S3

source("../Plotfunctions/geweke.R")

# relies on the "coda" package to calculate Geweke test statistics
geweke(0.1)
[1] 0.0575648 0.1051838

geweke(0.25)
[1] 0.05846896 0.11332128

geweke(0.5)
[1] 0.06118143 0.11512960

geweke(0.75)
[1] 0.05304400 0.09975889

geweke(0.9)
[1] 0.05937312 0.11573237


load("../realdata/output/geweke_pvalues_Bayesian_FQR_qt10.RData")
pvalues1 <- c(pvalue[[1]],pvalue[[2]])

load("../realdata/output/geweke_pvalues_Bayesian_FQR_qt25.RData")
pvalues2 <- c(pvalue[[1]],pvalue[[2]])

load("../realdata/output/geweke_pvalues_Bayesian_FQR_qt50.RData")
pvalues3 <- c(pvalue[[1]],pvalue[[2]])

load("../realdata/output/geweke_pvalues_Bayesian_FQR_qt75.RData")
pvalues4 <- c(pvalue[[1]],pvalue[[2]])

load("../realdata/output/geweke_pvalues_Bayesian_FQR_qt90.RData")
pvalues5 <- c(pvalue[[1]],pvalue[[2]])


png("./FigureS3.png",width = 960, height = 960)
par(mfrow=c(3,2))
hist(pvalues1,ylim=c(0,1.25),freq=FALSE,cex.lab=2,xlab="Geweke p-value",cex.lab = 1.8,main="(a)", cex.main=2.5)
hist(pvalues2,ylim=c(0,1.25),freq=FALSE,cex.lab=2,xlab="Geweke p-value",cex.lab = 1.8,main="(b)", cex.main=2.5)
hist(pvalues3,ylim=c(0,1.25),freq=FALSE,cex.lab=2,xlab="Geweke p-value",cex.lab = 1.8,main="(c)", cex.main=2.5)
hist(pvalues4,ylim=c(0,1.25),freq=FALSE,cex.lab=2,xlab="Geweke p-value",cex.lab = 1.8,main="(d)", cex.main=2.5)
hist(pvalues5,ylim=c(0,1.25),freq=FALSE,cex.lab=2,xlab="Geweke p-value",cex.lab = 1.8,main="(e)", cex.main=2.5)
dev.off()
