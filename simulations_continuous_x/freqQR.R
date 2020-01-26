freqQR <- function(Y,X,qt,B)
{
   library(quantreg)

   set.seed(12345)

   beta2.1 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta2.2 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta2.3 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta2.4 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta2.5 <- matrix(NA,nrow=B,ncol=ncol(Y))

   beta3.1 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta3.2 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta3.3 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta3.4 <- matrix(NA,nrow=B,ncol=ncol(Y))
   beta3.5 <- matrix(NA,nrow=B,ncol=ncol(Y))
  

   for(i in 1:B)
   {
      ind <- sample(1:nrow(Y),replace=TRUE)
      Y0 <- Y[ind,]
      X0 <- X[ind,]
      
      for(t in 1:ncol(Y0))
      {
         tmp.data <- data.frame(y=Y0[,t],x2=X0[,2],x3=X0[,3])

         suppressWarnings(tmp.coef <- coef(rq(y~x2+x3,tau=qt,data=tmp.data,ci=FALSE)))

         beta2.1[i,t] <- tmp.coef[2,1]
         beta2.2[i,t] <- tmp.coef[2,2]
         beta2.3[i,t] <- tmp.coef[2,3]
         beta2.4[i,t] <- tmp.coef[2,4]
         beta2.5[i,t] <- tmp.coef[2,5]

         beta3.1[i,t] <- tmp.coef[3,1]
         beta3.2[i,t] <- tmp.coef[3,2]
         beta3.3[i,t] <- tmp.coef[3,3]
         beta3.4[i,t] <- tmp.coef[3,4]
         beta3.5[i,t] <- tmp.coef[3,5]
      }
    
   }

   return(list(beta2.1=beta2.1,beta2.2=beta2.2,beta2.3=beta2.3,beta2.4=beta2.4,beta2.5=beta2.5,
               beta3.1=beta3.1,beta3.2=beta3.2,beta3.3=beta3.3,beta3.4=beta3.4,beta3.5=beta3.5))
}



      

      