.libPaths("/home/grad/rondai/RLibrary")
suppressMessages({
  
  if (!require(ncvreg)) install.packages('glmnet', repos =
                                           'https://cran.revolutionanalytics.com/')
  
  if (!require(doParallel)) install.packages('doParallel', repos =
                                               'http://cran.r-project.org')
  
  
  library(ncvreg)
  
  library(doParallel)
  
})

###########################
# bic-related methods
lsr.bic = function(y, x, nlambda = 100, t = 1, penalty = "SCAD")
{
  
  n = nrow(x)
  p = ncol(x)
  
  mod = ncvreg(x, y, family = "gaussian", penalty = penalty, nlambda = nlambda) 
  
  fit = mod $ beta
  
  loss = mod $ loss
  card = count(fit[-1, ])
  lsbic = bic(loss, card, t, n, p)
  
  which.lambda = which.min(lsbic)
  
  return(fit[, which.lambda])
  
  
}

pqr.bic = function(y, x, tau, lambdaset, t = 1, penalty="SCAD", 
                   maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                   a=3.7)
{
  
  n = nrow(x)
  p = ncol(x)
  
  fit = di.qr(y, x, tau = tau, lambdaset, penalty=penalty, maxin=maxin, 
              maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
  
  m = ncol(fit)
  residual = matrix(rep(y, m), n, m) - cbind(rep(1, n), x) %*% fit
  loss = residual * (tau - (residual < 0))
  loss = colSums(loss)
  card = count(fit[-1, ])
  qbic = bic(loss, card, t, n, p)
  
  which.lambda = which.min(qbic)
  
  return(fit[, which.lambda])
  
}

bic = function(loss, card, t, n, p) log(loss) + card * log(log(n)) * log(p) / (n * t)
###########################

###########################
# cv-realted methods
lsr.cv = function(y, x, nlambda = 100, seed, nfold = 5, penalty="SCAD")
{
  
  mod = cv.ncvreg(X = x, y = y, family = "gaussian", penalty = penalty, nlambda = nlambda,
                  seed = seed, nfolds = nfold) 
  res = mod $ fit $ beta[, mod $ min]
  
  return(res)
  
}


pqr.cv = function(y, x, tau, lambdaset, seed, nfold = 5, penalty="SCAD", 
                  maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08, a=3.7)
{
  
  dyn.load("cqpen.so")
  
  
  n = length(y)
  
  set.seed(seed)
  perm = sample(1 : n, n)
  x = x[perm, ]
  y = y[perm]
  
  n1 = floor(n / nfold)
  n2 = n - (nfold - 1) * n1
  
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype = as.integer(0)  
  } else if (penalty == "MCP"){
    pentype = as.integer(1)
  } else{
    pentype = as.integer(2)
  }
  
  p = as.integer( ncol(x) )
  n = as.integer(n)
  n1 = as.integer(n1)
  n2 = as.integer(n2)
  nfold = as.integer(nfold)
  ll = as.integer(length(lambdaset))
  
  residuals=matrix(0, n * ll, 1)
  
  
  betaset = matrix(as.double(0), ll, p)
  residuals = matrix(y, n, 1)
  
  y         = as.double(y)
  xdouble     = matrix(as.double(unlist(x)),n,p)
  tau       = as.double(tau)
  a         = as.double(a)
  eps       = as.double(eps)
  maxin     = as.integer(maxin)
  maxout=as.integer(maxout)
  lambdaset    = as.double(lambdaset)
  loss = as.double(numeric(ll))
  
  intvalset = matrix(as.double(0), ll, 1) 
  
  out=.Fortran("cvpqr",xdouble,y,loss,betaset,intvalset,residuals,
               n,p,tau,a,eps,maxin,maxout,lambdaset,pentype,ll,nfold,n1,n2)
  
  
  loss = out[[3]]
  
  
  lambda = lambdaset[which.min(loss)]
  
  fit = di.qr(y, x, tau = tau, lambda, penalty=penalty, maxin=maxin, 
              maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
  
  return(fit[, 1])
  
}

# pqr.cv = function(y, x, tau, lambdaset, seed, nfold = 10, penalty="SCAD", 
#                   maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08, a=3.7)
# {
#   
#   dyn.load("cqpen.so")
#   
#   
#   n = length(y)
#   
#   set.seed(seed)
#   perm = sample(1 : n, n)
#   x = x[perm, ]
#   y = y[perm]
#   
#   n1 = floor(n / nfold)
#   n2 = n - (nfold - 1) * n1
#   
#   y = list(y)
#   x = list(x)
#   
#   # Set penalty function
#   if(penalty == "SCAD"){
#     pentype = as.integer(0)  
#   } else if (penalty == "MCP"){
#     pentype = as.integer(1)
#   } else{
#     pentype = as.integer(2)
#   }
#   
#   K = as.integer(length(y))
#   p = as.integer( ncol(x[[1]]) )
#   m = as.integer(length(tau))
#   n = as.integer(n)
#   n1 = as.integer(n1)
#   n2 = as.integer(n2)
#   nfold = as.integer(nfold)
#   ll = as.integer(length(lambdaset))
#   
#   residuals=matrix(0, n * ll, m * K)
#   
#   
#   betaset = matrix(as.double(0), ll * nfold, p * K * m)
#   
#   for(i in 1 : K)
#   {
#     
#     residuals[, (m * (i - 1) + 1) : (m * i)] = matrix(as.double(rep(y[[i]], m)), n, m)
#     
#   }
#   
#   y         = lapply(y,as.double)
#   xdouble     = matrix(as.double(unlist(x)),n,p*K)
#   tau       = as.double(tau)
#   a         = as.double(a)
#   eps       = as.double(eps)
#   maxin     = as.integer(maxin)
#   maxout=as.integer(maxout)
#   lambdaset    = as.double(lambdaset)
#   
#   intvalset = matrix(as.double(0), ll * nfold, m * K) 
#   
#   out=.Fortran("cvpqr",xdouble,betaset,intvalset,residuals,
#                n,p,K,tau,m,a,eps,maxin,maxout,lambdaset,pentype,ll,nfold,n1,n2)
#   
#   
#   betaset = out[[2]]
#   intvalset = out[[3]]
#   
#   betaset[abs(betaset)<coef.cutoff]=0
#   
#   coefficients = rbind(intvalset[, 1], t(betaset))
#   
#   
#   x = x[[1]]
#   y = y[[1]]
#   
#   loss = numeric(ll)
#   
#   for(i in 1 : nfold)
#   {
#     
#     if(i < nfold)
#     {
#       
#       ind = (i - 1) * n1 + (1 : n1)
#       
#       x.v = x[ind, ]
#       y.v = y[ind]
#       x.t = x[-ind, ]
#       y.t = y[-ind]
#       
#       coef.v = coefficients[, (i - 1) * ll + (1 : ll)]
#       
#       error = matrix(rep(y.v, ll), n1, ll) - cbind(rep(1, n1), x.v) %*% coef.v
#       
#       loss = loss + colSums(error * (tau - (error < 0)))
#       
#       
#     } else {
#       
#       ind = ((nfold - 1) * n1 + 1) : n
#       
#       x.v = x[ind, ]
#       y.v = y[ind]
#       x.t = x[-ind, ]
#       y.t = y[-ind]
#       
#       coef.v = coefficients[, (i - 1) * ll + (1 : ll)]
#       
#       error = matrix(rep(y.v, ll), n2, ll) - cbind(rep(1, n2), x.v) %*% coef.v
#       
#       loss = loss + colSums(error * (tau - (error < 0)))
#       
#     }
#     
#   }
#   
#   lambda = lambdaset[which.min(loss)]
#   
#   fit = di.qr(y, x, tau = tau, lambda, penalty=penalty, maxin=maxin, 
#               maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
#   
#   return(fit[, 1])
#   
# }

# pqr.cv = function(y, x, tau, lambdaset, seed, nfold = 5, penalty="SCAD", parallel = F,
#                   maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08, a=3.7)
# {
#   
#   n = length(y)
#   
#   set.seed(seed)
#   perm = sample(1 : n, n)
#   x = x[perm, ]
#   y = y[perm]
#   
#   n.cv = floor(n / nfold)
#   
#   if(parallel == T)
#   {
#     error = foreach(i = 1 : nfold, .combine = "rbind") %dopar%
#       {
#         
#         if(i < nfold)
#         {
#           
#           ind = (i - 1) * n.cv + (1 : n.cv)
#           
#           x.v = x[ind, ]
#           y.v = y[ind]
#           x.t = x[-ind, ]
#           y.t = y[-ind]
#           
#           loss = pqr.error(y = y.t, x = x.t, tau = tau, lambdaset = lambdaset, y.v = y.v, 
#                            x.v = x.v, penalty = penalty, maxin=maxin, maxout=maxout, eps = eps, 
#                            coef.cutoff = coef.cutoff, a = a)
#           
#         } else {
#           
#           ind = ((nfold - 1) * n.cv + 1) : n
#           
#           x.v = x[ind, ]
#           y.v = y[ind]
#           x.t = x[-ind, ]
#           y.t = y[-ind]
#           
#           loss = pqr.error(y = y.t, x = x.t, tau = tau, lambdaset = lambdaset, y.v = y.v, 
#                            x.v = x.v, penalty = penalty, maxin=maxin, maxout=maxout, eps = eps, 
#                            coef.cutoff = coef.cutoff, a = a)
#           
#         }
#         
#         
#         loss
#         
#       }
#     
#   } else {
#     
#     error = foreach(i = 1 : nfold, .combine = "rbind") %do%
#       {
#         
#         if(i < nfold)
#         {
#           
#           ind = (i - 1) * n.cv + (1 : n.cv)
#           
#           x.v = x[ind, ]
#           y.v = y[ind]
#           x.t = x[-ind, ]
#           y.t = y[-ind]
#           
#           loss = pqr.error(y = y.t, x = x.t, tau = tau, lambdaset = lambdaset, y.v = y.v, 
#                            x.v = x.v, penalty = penalty, maxin=maxin, maxout=maxout, eps = eps, 
#                            coef.cutoff = coef.cutoff, a = a)
#           
#         } else {
#           
#           ind = ((nfold - 1) * n.cv + 1) : n
#           
#           x.v = x[ind, ]
#           y.v = y[ind]
#           x.t = x[-ind, ]
#           y.t = y[-ind]
#           
#           loss = pqr.error(y = y.t, x = x.t, tau = tau, lambdaset = lambdaset, y.v = y.v, 
#                            x.v = x.v, penalty = penalty, maxin=maxin, maxout=maxout, eps = eps, 
#                            coef.cutoff = coef.cutoff, a = a)
#           
#         }
#         
#         
#         loss
#         
#       }
#     
#   }
#   
#   error.sum = colSums(error)
#   lambda = lambdaset[which.min(error.sum)]
#   
#   fit = di.qr(y, x, tau = tau, lambda, penalty=penalty, maxin=maxin, 
#               maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
#   
#   return(fit[, 1])
#   
# }
# 
# 
# pqr.error = function(y, x, tau, lambdaset, y.v, x.v, penalty="SCAD", 
#                      maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
#                      a=3.7)
# {
#   
#   fit = di.qr(y, x, tau = tau, lambdaset, penalty=penalty, maxin=maxin, 
#               maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
#   
#   m = ncol(fit)
#   n1 = length(y.v)
#   residual = matrix(rep(y.v, m), n1, m) - cbind(rep(1, n1), x.v) %*% fit
#   loss = residual * (tau - (residual < 0))
#   bic = colSums(loss)
#   
#   return(bic)
#   
# }
###########################

########################### 
lsr = function(y, x, nlambda = 100, y.v, x.v, penalty = "SCAD")
{
  
  n = nrow(x)
  p = ncol(x)
  
  mod = ncvreg(x, y, family = "gaussian", penalty = penalty, nlambda = nlambda) 
  
  fit = mod $ beta
  
  # residual = mod $ loss
  # card = apply(fit[-1, ], 2, function(a) sum(a != 0))
  # bic = bic(residual, n, p, card, bict)
  
  # bic = log(sum(residual)) + apply(fit[-1, ], 2, function(a) sum(a != 0)) * log(log(n)) * log(p) / n
  
  m = ncol(fit)
  n1 = length(y.v)
  residual = matrix(rep(y.v, m), n1, m) - cbind(rep(1, n1), x.v) %*% fit
  loss = residual ^ 2
  bic = colSums(loss)
  
  which.lambda = which.min(bic)
  
  return(fit[, which.lambda])
  
  
}
########################### 

########################### 
pqr = function(y, x, tau, lambdaset, y.v, x.v, penalty="SCAD", 
               maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
               a=3.7)
{
  
  fit = di.qr(y, x, tau = tau, lambdaset, penalty=penalty, maxin=maxin, 
              maxout=maxout, eps = eps, coef.cutoff=coef.cutoff, a=a) 
  
  m = ncol(fit)
  n1 = length(y.v)
  residual = matrix(rep(y.v, m), n1, m) - cbind(rep(1, n1), x.v) %*% fit
  loss = residual * (tau - (residual < 0))
  bic = colSums(loss)
  
  which.lambda = which.min(bic)
  
  return(fit[, which.lambda])
  
}


di.qr=function(y, x, tau, lambdaset, penalty="SCAD", maxin=100, maxout=20, 
               eps = 1e-05, coef.cutoff=1e-08, a=3.7)
{
  
  dyn.load("cqpen.so")
  
  y = list(y)
  x = list(x)
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype = as.integer(0)  
  } else if (penalty == "MCP"){
    pentype = as.integer(1)
  } else{
    pentype = as.integer(2)
  }
  
  K = as.integer(length(y))
  p = as.integer( ncol(x[[1]]) )
  m = as.integer(length(tau))
  n = as.integer(length(y[[1]]))
  ll = as.integer(length(lambdaset))
  
  residualsx=matrix(0, n * ll, m * K)
  
  
  betax = matrix(as.double(0), ll, p * K * m)
  
  for(i in 1 : K)
  {
    
    residualsx[, (m * (i - 1) + 1) : (m * i)] = matrix(as.double(rep(y[[i]], ll * m)), n * ll, m)
    
  }
  
  y         = lapply(y,as.double)
  xdouble     = matrix(as.double(unlist(x)),n,p*K)
  tau       = as.double(tau)
  a         = as.double(a)
  eps       = as.double(eps)
  maxin     = as.integer(maxin)
  maxout=as.integer(maxout)
  lambdaset    = as.double(lambdaset)
  
  intvalx = matrix(as.double(0), ll, m * K) 
  
  
  out=.Fortran("pmqmlamint",xdouble,betax,intvalx,residualsx,
               n,p,K,tau,m,a,eps,maxin,maxout,lambdaset,pentype,ll)
  
  
  
  betaset = out[[2]]
  intset = out[[3]]
  
  betaset[abs(betaset)<coef.cutoff]=0
  
  coefficients = rbind(intset[, 1], t(betaset))
  
  return(coefficients)
  
}


count = function(x)
{
  res = apply(x, 2, function(a) sum(a != 0))
  return(res)
}

