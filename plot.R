.libPaths("/home/grad/rondai/RLibrary")
setwd("/home/grad/rondai/vsv")

suppressMessages({
  
  if (!require(doParallel)) install.packages('doParallel', repos =
                                               'http://cran.r-project.org')
  if (!require(LaplacesDemon)) install.packages('LaplacesDemon', repos =
                                                  'http://cran.r-project.org')
  if (!require(ncvreg)) install.packages('ncvreg', repos =
                                           'http://cran.r-project.org')
  
  if (!require(rqPen)) install.packages('rqPen', repos =
                                          'http://cran.r-project.org')
  
  
  
  library(doParallel)
  library(LaplacesDemon)
  library(ncvreg)
  library(rqPen)
  library(MASS)
  
})


source("vsv_functions.R")

registerDoParallel(detectCores())
##########################
# Basic setting
s = 200
n = 200
n1 = 10*n
# sigma=sqrt(3)
tau=seq(1/10, 9/10, 1/10)
k = length(tau)

p = 6 # model size
vt = rep(0, p)
vt[c(1, 2, 5)] = 1
ind = which(vt != 0)
# vt[c(1, 2, 5)] = c(3, 1.5, 2)
# vt[c(1,3,5,8,10,13,16)]=c(2,1.5,0.8,1,1.75,0.75,0.3)
cut.off=seq(ceiling(0.5 * k), k-1, 1)
u=length(cut.off)
##############
lambda.qr=exp(seq(-5,5,length=100)) # lambda set for qr
nlambda = 100

resample.method = 2
# 0 is no resampling. 1 is resampling n with replication. 2 is resample n/2 without replication.
##############


##################
# Covariance matrix of the covariates
covm=matrix(0, p, p)
for(i in 1 : p)
{
  for (j in 1 : p) covm[i, j]=0.5 ^ (abs(i - j))
}
###################

# out=numeric(16)

distribution = 1
res = foreach(i = 1 : s, .combine = "rbind") %dopar%
  {
    
    set.seed(i)
    
    # vt[c(1, 2, 5, 8)] = runif(length(ind), 0.05, 2)
    
    ##############################
    # Training data
    x = mvrnorm(n, rep(0, p), covm)
    ##############################
    # Various error distributions
    if(distribution == 1) e = rcauchy(n)
    if(distribution == 2) e = rnormm(n, c(0.5, 0.5), c(-5, 5), c(1, 1))
    ##############################
    y = (x %*% vt + e)[,1]
    ##############################
    
    ##############################
    # Validation data
    x.v = mvrnorm(n1, rep(0, p), covm)
    ##############################
    # Various error distributions
    if(distribution == 1) e.v = rcauchy(n1)
    if(distribution == 2) e.v = rnormm(n1, c(0.5, 0.5), c(-5, 5), c(1, 1))
    ##############################
    y.v = (x.v %*% vt + e.v)[, 1]
    ##############################
    
    ##############################
    # estimator with resampling and multipe quantiles, least squares and least absolute deviation
    resample.mq = matrix(0, nrow=k, ncol=p)
    resample.ls = matrix(0, nrow=k, ncol=p)
    resample.lad = matrix(0, nrow=k, ncol=p)
    for(j in 1 : k)
    {
      
      set.seed(j)
      
      if(resample.method == 0) resample = 1 : n
      if(resample.method == 1) resample = sample(1 : n, replace = T)
      if(resample.method == 2) resample = sample(x = 1 : n, size = ceiling(n / 2), replace = F)
      x1 = x[resample, ]
      y1 = y[resample]
      
      fit.mq = pqr(y, x, tau = tau[j], lambda.qr, y.v, x.v)
      resample.mq[j, ] = fit.mq[-1]
      
      fit.ls = lsr(y1, x1, nlambda, y.v, x.v)
      resample.ls[j, ] = fit.ls[-1]
      
      
      if(j != ceiling(k / 2)) 
      {
        
        fit.lad = pqr(y1, x1, tau = 0.5, lambda.qr, y.v, x.v)
        resample.lad[j, ] = fit.lad[-1]
        
      } else {
        
        resample.lad[j, ] = fit.mq[-1]
        
      }
      
    }
    ##############################
    
    
    # ##############################
    # # estimator with least squares
    # est.ls = lsr(y, x, nlambda = nlambda, bict = bict)
    # est.ls = est.ls[-1]
    # ##############################
    # 
    # ##############################
    # # estimator with least absolute deviation
    # est.lad = pqr(list(y), list(x), tau = 0.5, lambda.qr, bict = bict) $ coefficients
    # est.lad = est.lad[[1]]
    # est.lad = est.ls[-1]
    # ##############################
    
    count.mq = count(resample.mq)
    count.ls = count(resample.ls)
    count.lad = count(resample.lad)
    
    rbind(count.mq, count.ls, count.lad)
    
  }


res1 = matrix(0, s, p)
res2 = matrix(0, s, p)
res3 = matrix(0, s, p)

for(i in 1 : s)
{
  res1[i, ] = res[3 * (i - 1) + 1, ]
  res2[i, ] = res[3 * (i - 1) + 2, ]
  res3[i, ] = res[3 * (i - 1) + 3, ]
}

distribution = 2
res = foreach(i = 1 : s, .combine = "rbind") %dopar%
  {
    
    set.seed(i)
    
    # vt[c(1, 2, 5, 8)] = runif(length(ind), 0.05, 2)
    
    ##############################
    # Training data
    x = mvrnorm(n, rep(0, p), covm)
    ##############################
    # Various error distributions
    if(distribution == 1) e = rcauchy(n)
    if(distribution == 2) e = rnormm(n, c(0.5, 0.5), c(-5, 5), c(1, 1))
    ##############################
    y = (x %*% vt + e)[,1]
    ##############################
    
    ##############################
    # Validation data
    x.v = mvrnorm(n1, rep(0, p), covm)
    ##############################
    # Various error distributions
    if(distribution == 1) e.v = rcauchy(n1)
    if(distribution == 2) e.v = rnormm(n1, c(0.5, 0.5), c(-5, 5), c(1, 1))
    ##############################
    y.v = (x.v %*% vt + e.v)[, 1]
    ##############################
    
    ##############################
    # estimator with resampling and multipe quantiles, least squares and least absolute deviation
    resample.mq = matrix(0, nrow=k, ncol=p)
    resample.ls = matrix(0, nrow=k, ncol=p)
    resample.lad = matrix(0, nrow=k, ncol=p)
    for(j in 1 : k)
    {
      
      set.seed(j)
      
      if(resample.method == 0) resample = 1 : n
      if(resample.method == 1) resample = sample(1 : n, replace = T)
      if(resample.method == 2) resample = sample(x = 1 : n, size = ceiling(n / 2), replace = F)
      x1 = x[resample, ]
      y1 = y[resample]
      
      fit.mq = pqr(y, x, tau = tau[j], lambda.qr, y.v, x.v)
      resample.mq[j, ] = fit.mq[-1]
      
      fit.ls = lsr(y1, x1, nlambda, y.v, x.v)
      resample.ls[j, ] = fit.ls[-1]
      
      
      if(j != ceiling(k / 2)) 
      {
        
        fit.lad = pqr(y1, x1, tau = 0.5, lambda.qr, y.v, x.v)
        resample.lad[j, ] = fit.lad[-1]
        
      } else {
        
        resample.lad[j, ] = fit.mq[-1]
        
      }
      
    }
    ##############################
    
    
    # ##############################
    # # estimator with least squares
    # est.ls = lsr(y, x, nlambda = nlambda, bict = bict)
    # est.ls = est.ls[-1]
    # ##############################
    # 
    # ##############################
    # # estimator with least absolute deviation
    # est.lad = pqr(list(y), list(x), tau = 0.5, lambda.qr, bict = bict) $ coefficients
    # est.lad = est.lad[[1]]
    # est.lad = est.ls[-1]
    # ##############################
    
    count.mq = count(resample.mq)
    count.ls = count(resample.ls)
    count.lad = count(resample.lad)
    
    rbind(count.mq, count.ls, count.lad)
    
  }


res4 = matrix(0, s, p)
res5 = matrix(0, s, p)
res6 = matrix(0, s, p)

for(i in 1 : s)
{
  res4[i, ] = res[3 * (i - 1) + 1, ]
  res5[i, ] = res[3 * (i - 1) + 2, ]
  res6[i, ] = res[3 * (i - 1) + 3, ]
}

lab=expression(X[1],X[2],X[3],X[4],X[5],X[6])
bcol=rep(0,p)
bcol[ind]="red"
bcol[-ind] = "NA"

# if(resample.method == 0) name = "count0.eps"
# if(resample.method == 1) name = "count1.eps"
# if(resample.method == 2) name = "count2.eps"

name = "Cauchy_LMN.eps"

setEPS()
postscript(name)
par(mfrow=c(3,2))
boxplot(res1,main="Multiple Quantile Loss Functions",names=lab,ylab="Votes",cex.lab=1.2,cex.main=1,
        cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")

boxplot(res4,main="Multiple Quantile Loss Functions",names=lab,ylab="Votes",cex.lab=1.2,cex.main=1,
        cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")

boxplot(res2,main="Quadratic Loss Function with Resampling",names=lab,ylab="Votes",cex.lab=1.2,
        cex.main=1,cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")

boxplot(res5,main="Quadratic Loss Function with Resampling",names=lab,ylab="Votes",cex.lab=1.2,
        cex.main=1,cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")

boxplot(res3,main="Absolute Loss Function with Resampling",names=lab,ylab="Votes",cex.lab=1.2,
        cex.main=1,cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")

boxplot(res6,main="Absolute Loss Function with Resampling",names=lab,ylab="Votes",cex.lab=1.2,
        cex.main=1,cex.axis=1.2,xaxt="n",col=bcol)
axis(1,cex.axis=1.2,at=1:p,labels=lab)
axis(1,cex.axis=1.2,at=ind,labels=lab[ind],col.axis="red")
dev.off()

par(mfrow=c(1,1))





