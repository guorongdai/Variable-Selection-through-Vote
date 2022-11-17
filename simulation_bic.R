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
  
  if (!require(writexl)) install.packages('writexl', repos =
                                            'http://cran.r-project.org')
  
  
  
  library(doParallel)
  library(LaplacesDemon)
  library(ncvreg)
  library(rqPen)
  library(writexl)
  library(MASS)
  
})


source("vsv_functions.R")

registerDoParallel(detectCores())


simulation.no = 1
if(simulation.no == 1) distribution.set = c(1, 2)
if(simulation.no == 2) distribution.set = c(3, 4)
if(simulation.no == 3) distribution.set = c(5, 6)
zz = length(distribution.set)

##########################
# Basic setting
s = 200
n = 200
t = 6
tau=seq(1/10, 9/10, 1/10)
k = length(tau)

p = 400 # model size

#####################################
# setting parameters
vt = rep(0, p)
vt[c(1, 3, 5, 8, 10, 13, 16)]=c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.5)
############################ 

ind = which(vt != 0)
#####################################

cut.off=seq(ceiling(0.5 * k), k, 1)
u=length(cut.off)
##############
lambda.qr=exp(seq(-5,5,length=100)) # lambda set for qr
nlambda = 100
##############


##################
# Covariance matrix of the covariates
covm=matrix(0, p, p)
for(i in 1 : p)
{
  for (j in 1 : p) covm[i, j]=0.5 ^ (abs(i - j))
}
###################

table.ave = matrix(0, 3 * u + 2, 12)
table.std = matrix(0, 3 * u + 2, 12)


res = foreach(i = 1 : (s * zz)) %dopar%
  {
    
    distribution = distribution.set[ceiling(i / s)]
    
    set.seed(i)
    
    ##############################
    # Data
    x = mvrnorm(n, rep(0, p), covm)
    ##############################
    # Various error distributions
    if(distribution == 1) e = rnorm(n, 0, sqrt(3))
    if(distribution == 2) e = rt(n, 2)
    if(distribution == 3) e = rlaplace(n)
    if(distribution == 4) e = rnormm(n, c(0.5, 0.5), c(-2, 2), c(1, 1))
    if(distribution == 5) e = rnormm(n, c(0.1, 0.9), c(0, 0), c(5, 1))
    if(distribution == 6) 
    {
      e = rcauchy(n)
      t = 20
    }
    ##############################
    y = (x %*% vt + e)[,1]
    ##############################
    
    ##############################
    # estimator with resampling and multipe quantiles, least squares and least absolute deviation
    mq = matrix(0, nrow=k, ncol=p)
    resample.ls = matrix(0, nrow=k, ncol=p)
    resample.lad = matrix(0, nrow=k, ncol=p)
    for(j in 1 : k)
    {
      
      set.seed(j)
      
      fit.mq = pqr.bic(y, x, tau[j], lambda.qr, t)
      mq[j, ] = fit.mq[-1]
      
    }
    ##############################
    
    
    ##############################
    # estimator with least squares
    est.ls = numeric(p)
    ##############################
    
    ##############################
    # estimator with least absolute deviation
    est.lad = numeric(p)
    ##############################
    
    count.mq = count(mq)
    count.ls = count(resample.ls)
    count.lad = count(resample.lad)
    
    count.matrix = rbind(count.mq, count.ls, count.lad)
    n.method = nrow(count.matrix)
    
    
    outcome = matrix(0, n.method * u + 2, 2)
    
    for(j in 1 : u)
    {
      
      threshold = cut.off[j]
      
      selection = (count.matrix >= threshold)
      
      fn = apply(selection, 1, function(a) sum(a[ind] == F)) # false negative
      fp = apply(selection, 1, function(a) sum(a[-ind] == T)) # false positive
      
      outcome[(j - 1) * n.method + (1 : n.method), 1] = fn
      outcome[(j - 1) * n.method + (1 : n.method), 2] = fp
      
    }
    
    ls.lad = rbind(est.ls, est.lad)
    fn.extra = apply(ls.lad, 1, function(a) sum(a[ind] == 0))
    fp.extra = apply(ls.lad, 1, function(a) sum(a[-ind] != 0))
    
    outcome[n.method * u + (1 : 2), 1] = fn.extra
    outcome[n.method * u + (1 : 2), 2] = fp.extra
    
    outcome
    
  }

for(jj in 1 : zz)
{
  
  res.tt = res[(jj - 1) * s + (1 : s)]
  
  tt = distribution.set[jj]
  
  table.ave[, (tt - 1) * 2 + (1 : 2)] = apply(simplify2array(res.tt), 1 : 2, mean)
  table.std[, (tt - 1) * 2 + (1 : 2)] = apply(simplify2array(res.tt), 1 : 2, sd)
  
}

table.ave
table.std

final.table = list("mean" = data.frame(table.ave), "sd" = data.frame(table.std))

location = paste("/home/grad/rondai/vsv/simulation_", simulation.no, ".xlsx", sep = "")
write_xlsx(final.table, path = location)

registerDoSEQ()