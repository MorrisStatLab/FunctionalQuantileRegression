freqQR <- function(Y,X,time,qt,B)
{
   # the function to run naive bootstrap-based FQR and generate B bootstrap samples
   # the function to run two-step bootstrap-based FQR with spline smoothing and generate B bootstrap samples
    
   # Bootstrap-based FQR relies on the "quantreg" package, so please install it before running the code
   # qt is a vector of length 5 specifying the quantile levels used in our work

   library(quantreg)

   set.seed(123)

   coef1 <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef2 <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef3 <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef4 <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef5 <- matrix(NA,nrow=B,ncol=ncol(Y))

   coef1.ss <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef2.ss <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef3.ss <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef4.ss <- matrix(NA,nrow=B,ncol=ncol(Y))
   coef5.ss <- matrix(NA,nrow=B,ncol=ncol(Y))   

   for(i in 1:B)
   {
      ind <- sample(1:nrow(Y),replace=TRUE)
      Y0 <- Y[ind,]
      X0 <- X[ind,]
      
      for(t in 1:ncol(Y0))
      {
         tmp.data <- data.frame(y=Y0[,t],x2=X0[,2])
         suppressWarnings(tmp.coef <- coef(rq(y~x2,tau=qt,data=tmp.data,ci=FALSE)))
         coef1[i,t] <- tmp.coef[2,1]
         coef2[i,t] <- tmp.coef[2,2]
         coef3[i,t] <- tmp.coef[2,3]
         coef4[i,t] <- tmp.coef[2,4]
         coef5[i,t] <- tmp.coef[2,5]
      }

      fit1 <- smooth.spline(x=time,y=coef1[i,])
      fit2 <- smooth.spline(x=time,y=coef2[i,])
      fit3 <- smooth.spline(x=time,y=coef3[i,])
      fit4 <- smooth.spline(x=time,y=coef4[i,])
      fit5 <- smooth.spline(x=time,y=coef5[i,])

      coef1.ss[i,] <- as.numeric(fit1$y)
      coef2.ss[i,] <- as.numeric(fit2$y)
      coef3.ss[i,] <- as.numeric(fit3$y)
      coef4.ss[i,] <- as.numeric(fit4$y)
      coef5.ss[i,] <- as.numeric(fit5$y)     
   }

   return(list(coef1=coef1,coef2=coef2,coef3=coef3,coef4=coef4,coef5=coef5,
               coef1.ss=coef1.ss,coef2.ss=coef2.ss,coef3.ss=coef3.ss,coef4.ss=coef4.ss,coef5.ss=coef5.ss))
}



      

      