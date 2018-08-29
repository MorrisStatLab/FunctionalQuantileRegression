setwd("Y:/submission/simulations2")

library(FDboost)

X <- read.table("X.txt",header=FALSE,sep="\t")

x0 <- read.table("x0.txt",header=FALSE,sep="\t")

x0 <- as.numeric(x0[,1])


# replicate 1
Y <- read.table("data/model1.txt",header=FALSE,sep="\t")
data <- list(Y = as.matrix(Y), group = as.factor(X[,2]), time = x0, group1 = as.factor(1*(X[,2]==1)) )


# Fit FDboost to simulated functional data using default specs (except for larger df for the base learner modeling "time" variable)
# Get an identically 0 estimate of the group effect
start_time <- Sys.time()
model <- FDboost(Y ~ 1 + bolsc(group,df=1) %A0% bbs(time,df=20),
                 timeformula=~bbs(time,df=20),
                 numInt="equal",family=QuantReg(tau=0.9),
                 offset=NULL,data=data,control=boost_control(mstop=100,nu=0.2))
end_time <- Sys.time()

end_time-start_time
Time difference of 6.8212 secs


start_time <- Sys.time()
app <- applyFolds(model,folds=cv(rep(1,length(unique(model$id))),B = 5),grid=seq(1,2001,20))
end_time <- Sys.time()

end_time-start_time
Time difference of 8.560407 mins

mstop(app)
[1] 2001

plot(app,main="5-fold cross validation")


start_time <- Sys.time()
model <- model[mstop(app)]
end_time <- Sys.time()

end_time-start_time
Time difference of 1.524513 mins

# This fitting process including parameter tuning takes about 10 minutes,
# on a 64-bit operating system with 2 processors and 256GB RAM.

summary(coef(model))
        Length Class  Mode
offset  6      -none- list
smterms 2      -none- list

coef(model)$smterms
$`bols(ONEx) %A0% bbs(time)`
$`bols(ONEx) %A0% bbs(time)`$x
 [1]  0.0000000  0.3846154  0.7692308  1.1538462  1.5384615  1.9230769
 [7]  2.3076923  2.6923077  3.0769231  3.4615385  3.8461538  4.2307692
[13]  4.6153846  5.0000000  5.3846154  5.7692308  6.1538462  6.5384615
[19]  6.9230769  7.3076923  7.6923077  8.0769231  8.4615385  8.8461538
[25]  9.2307692  9.6153846 10.0000000 10.3846154 10.7692308 11.1538462
[31] 11.5384615 11.9230769 12.3076923 12.6923077 13.0769231 13.4615385
[37] 13.8461538 14.2307692 14.6153846 15.0000000

$`bols(ONEx) %A0% bbs(time)`$xlab
[1] "time"

$`bols(ONEx) %A0% bbs(time)`$xlim
[1]  0 15

$`bols(ONEx) %A0% bbs(time)`$value
 [1] 19.1699605 17.2901320 21.4683637 39.1921451 34.7447650  3.5155290
 [7]  6.1445536 43.6376023 41.8373477  5.0812016 -0.3714026 23.4088268
[13] 39.2458751 39.5293571 29.2240967 17.8464774 14.4352398 23.6826259
[19] 39.1401468 36.6958010 15.4051081  3.3141528 10.2810434 28.2834160
[25] 36.0476861 13.9258175 -3.1238392 21.0289418 43.9271945 29.7768462
[31] 11.6371233 12.4475323 17.9128675 20.5291348 25.0586285 28.2007711
[37] 20.8511414  8.6001821  3.4048207  4.3596563

$`bols(ONEx) %A0% bbs(time)`$dim
[1] 1

$`bols(ONEx) %A0% bbs(time)`$main
[1] "bols(ONEx) %A0% bbs(time)"


$`bolsc(group) %A0% bbs(time)`
$`bolsc(group) %A0% bbs(time)`$x
[1] -1 1 
Levels: -1 1

$`bolsc(group) %A0% bbs(time)`$y
 [1]  0.0000000  0.3846154  0.7692308  1.1538462  1.5384615  1.9230769
 [7]  2.3076923  2.6923077  3.0769231  3.4615385  3.8461538  4.2307692
[13]  4.6153846  5.0000000  5.3846154  5.7692308  6.1538462  6.5384615
[19]  6.9230769  7.3076923  7.6923077  8.0769231  8.4615385  8.8461538
[25]  9.2307692  9.6153846 10.0000000 10.3846154 10.7692308 11.1538462
[31] 11.5384615 11.9230769 12.3076923 12.6923077 13.0769231 13.4615385
[37] 13.8461538 14.2307692 14.6153846 15.0000000

$`bolsc(group) %A0% bbs(time)`$xlab
[1] "group"

$`bolsc(group) %A0% bbs(time)`$ylab
[1] "time"

$`bolsc(group) %A0% bbs(time)`$ylim
[1]  0 15

$`bolsc(group) %A0% bbs(time)`$xlim
[1] NA NA

$`bolsc(group) %A0% bbs(time)`$z
NULL

$`bolsc(group) %A0% bbs(time)`$zlab
[1] NA

$`bolsc(group) %A0% bbs(time)`$vecStand
NULL

$`bolsc(group) %A0% bbs(time)`$value
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
[1,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
[2,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
     [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
[1,]     0     0     0     0     0     0     0     0     0     0     0     0
[2,]     0     0     0     0     0     0     0     0     0     0     0     0
     [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
[1,]     0     0     0     0     0     0     0     0     0     0     0     0
[2,]     0     0     0     0     0     0     0     0     0     0     0     0
     [,39] [,40]
[1,]     0     0
[2,]     0     0

$`bolsc(group) %A0% bbs(time)`$dim
[1] 2

$`bolsc(group) %A0% bbs(time)`$main
[1] "bolsc(group) %A0% bbs(time)"

plot(model,main="Intercept function for tau=0.9\nnumber of boosting iterations=2001")

model.coef <- coef(model)$smterms

save(model.coef,file="output/FDboost/fitmodel1.RData")




# Fit FDboost to simulated functional data using more suitable, manually chosen specs
# This time we are able to get reasonable but quite noisy point estimate of the group effect
data <- list(Y = as.matrix(Y), group = as.factor(X[,2]), time = x0, group1 = as.factor(1*(X[,2]==1)) )

# just use a dummy for the effect of group1
# use a large and spiky basis for these spiky data without the inappropriate smoothness penalty
# use a reasonable offset instead of the (oversmoothed) default

start_time <- Sys.time()
model <- FDboost(Y ~ 1 + bols(group1, intercept=FALSE), 
                 timeformula=~bbs(time, degree=1, knots=150, differences=0, lambda=0),
                 numInt = "equal", family=QuantReg(tau=0.9), offset = apply(data$Y, 2, median),
                 data = data, control = boost_control(mstop = 100, nu = 0.2))
end_time <- Sys.time()

end_time-start_time
Time difference of 7.114039 secs


start_time <- Sys.time()
app <- applyFolds(model,folds=cv(rep(1,length(unique(model$id))),B = 5),grid=seq(1,1001,20))
model <- model[mstop(app)]
end_time <- Sys.time()

end_time-start_time
Time difference of 5.296196 mins

mstop(app)
[1] 1001

plot(app,main="5-fold cross validation")

# This fitting process including parameter tuning takes about 5 minutes,
# on a 64-bit operating system with 2 processors and 256GB RAM.


coef(model)$smterms
$`bols(ONEx) %A0% bbs(time)`
$`bols(ONEx) %A0% bbs(time)`$x
 [1]  0.0000000  0.3846154  0.7692308  1.1538462  1.5384615  1.9230769
 [7]  2.3076923  2.6923077  3.0769231  3.4615385  3.8461538  4.2307692
[13]  4.6153846  5.0000000  5.3846154  5.7692308  6.1538462  6.5384615
[19]  6.9230769  7.3076923  7.6923077  8.0769231  8.4615385  8.8461538
[25]  9.2307692  9.6153846 10.0000000 10.3846154 10.7692308 11.1538462
[31] 11.5384615 11.9230769 12.3076923 12.6923077 13.0769231 13.4615385
[37] 13.8461538 14.2307692 14.6153846 15.0000000

$`bols(ONEx) %A0% bbs(time)`$xlab
[1] "time"

$`bols(ONEx) %A0% bbs(time)`$xlim
[1]  0 15

$`bols(ONEx) %A0% bbs(time)`$value
 [1] 3.296209 3.193010 3.932183 5.733566 4.119804 3.950184 3.665792 3.752361
 [9] 3.863966 3.920670 3.596951 3.919132 5.494921 3.908634 3.975414 3.669607
[17] 4.086677 4.044263 4.162523 6.002289 3.618991 3.457433 3.610892 4.458655
[25] 4.611020 3.665233 3.597145 5.227501 4.735743 3.783176 3.925507 3.943706
[33] 3.181930 3.918631 4.822072 4.429694 3.669730 3.821208 3.812960 3.530105

$`bols(ONEx) %A0% bbs(time)`$dim
[1] 1

$`bols(ONEx) %A0% bbs(time)`$main
[1] "bols(ONEx) %A0% bbs(time)"


$`bols(group1) %O% bbs(time)`
$`bols(group1) %O% bbs(time)`$x
[1] 0 1
Levels: 0 1

$`bols(group1) %O% bbs(time)`$y
 [1]  0.0000000  0.3846154  0.7692308  1.1538462  1.5384615  1.9230769
 [7]  2.3076923  2.6923077  3.0769231  3.4615385  3.8461538  4.2307692
[13]  4.6153846  5.0000000  5.3846154  5.7692308  6.1538462  6.5384615
[19]  6.9230769  7.3076923  7.6923077  8.0769231  8.4615385  8.8461538
[25]  9.2307692  9.6153846 10.0000000 10.3846154 10.7692308 11.1538462
[31] 11.5384615 11.9230769 12.3076923 12.6923077 13.0769231 13.4615385
[37] 13.8461538 14.2307692 14.6153846 15.0000000

$`bols(group1) %O% bbs(time)`$xlab
[1] "group1"

$`bols(group1) %O% bbs(time)`$ylab
[1] "time"

$`bols(group1) %O% bbs(time)`$ylim
[1]  0 15

$`bols(group1) %O% bbs(time)`$xlim
[1] NA NA

$`bols(group1) %O% bbs(time)`$z
NULL

$`bols(group1) %O% bbs(time)`$zlab
[1] NA

$`bols(group1) %O% bbs(time)`$vecStand
NULL

$`bols(group1) %O% bbs(time)`$value
         [,1]      [,2]       [,3]       [,4]       [,5]       [,6]
[1,] 0.000000 0.0000000  0.0000000  0.0000000 0.00000000  0.0000000
[2,] 0.531365 0.8965906 -0.3249032 -0.4699844 0.05964217 -0.1990593
            [,7]      [,8]     [,9]       [,10]     [,11]     [,12]
[1,]  0.00000000 0.0000000 0.000000  0.00000000 0.0000000 0.0000000
[2,] -0.08466183 0.9674014 2.459784 -0.06537297 0.4884009 0.1143211
         [,13]     [,14]      [,15]     [,16]      [,17]      [,18]
[1,] 0.0000000 0.0000000  0.0000000 0.0000000  0.0000000  0.0000000
[2,] 0.1558577 0.1692969 -0.2782685 0.4827008 -0.5689437 -0.4171744
            [,19]     [,20]     [,21]     [,22]       [,23]     [,24]
[1,]  0.000000000 0.0000000 0.0000000 0.0000000  0.00000000 0.0000000
[2,] -0.006957003 0.3815852 0.2398819 0.3688478 -0.07461796 0.6719259
          [,25]     [,26]     [,27]     [,28]      [,29]      [,30]
[1,]  0.0000000 0.0000000 0.0000000  0.000000  0.0000000 0.00000000
[2,] -0.3728626 0.5112674 0.5198682 -1.998216 -0.9697207 0.02260837
          [,31]     [,32]      [,33]      [,34]       [,35]      [,36]
[1,]  0.0000000  0.000000 0.00000000  0.0000000  0.00000000  0.0000000
[2,] -0.2026788 -0.061247 0.09443297 -0.0937098 -0.01684742 -0.6737562
          [,37]    [,38]      [,39]     [,40]
[1,]  0.0000000  0.00000  0.0000000  0.000000
[2,] -0.4266212 -0.11017 -0.5267375 -0.157575

$`bols(group1) %O% bbs(time)`$dim
[1] 2

$`bols(group1) %O% bbs(time)`$main
[1] "bols(group1) %O% bbs(time)"

plot(model,main="Intercept function for tau=0.9\nnumber of boosting iterations=1001")





# for space considerations, the replicate datasets 2-10 are not provided but are available upon request.
for (i in 2:10)
{
   Y <- read.table(paste0("data/model",i,".txt"),header=FALSE,sep="\t")

   data <- list(Y=as.matrix(Y),group=as.factor(X[,2]),time=x0)

   model <- FDboost(Y ~ 1 +
                    bolsc(group,df=1) %A0% bbs(time,df=20),
                    timeformula=~bbs(time,df=20),
                    numInt="equal",family=QuantReg(tau=0.9),
                    offset=NULL,data=data,control=boost_control(mstop=100,nu=0.2))

   app <- applyFolds(model,folds=cv(rep(1,length(unique(model$id))),B = 5),grid=seq(1,2001,20))

   model <- model[mstop(app)]

   model.coef <- coef(model)$smterms

   save(model.coef,file=paste0("output/FDboost/fitmodel",i,".RData"))
}


