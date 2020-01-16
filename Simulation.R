#install.packages("glmnet",repos="https://cran.revolutionanalytics.com/")
#install.packages("rqPen",repos="https://cran.revolutionanalytics.com/")
#install.packages("doParallel",repos='https://cran.revolutionanalytics.com/')
#install.packages("LaplacesDemon",repos='https://cran.revolutionanalytics.com/')
#install.packages("cqrReg",repos='https://cran.revolutionanalytics.com/')
#install.packages("ncvreg",repos="https://cran.revolutionanalytics.com/")
#install.packages("quantreg",repos="https://cran.revolutionanalytics.com/")

suppressMessages({
  library(MASS)
  library(glmnet)
  library(rqPen)
  library(doParallel)
  library(LaplacesDemon)
  library(cqrReg)
  library(ncvreg)
  library(quantreg)
})


##############################################################
# Functions

##################################
# composite quantile regression with scad
cqr = function(y, x, tau=.5, lambda, penalty="SCAD", 
                         initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                         a=3.7, ...)
{
  #cleanInputs(y, x, lambda, initial_beta, penalty, a)
  dyn.load("cqpen.so")
  require(rqPen)
  
  # Set penalty function
  if(penalty == "SCAD"){
    pentype = as.integer(0)  
  } else if (penalty == "MCP"){
    pentype = as.integer(1)
  } else{
    pentype = as.integer(2)
  }
  
  if( is.null(initial_beta) ){
    if(length(tau)>1)
    {
      initial_beta = mapply(LASSO.fit,y=y,x=x,tau=0.5,lambda=lambda,
                            coef.cutoff=coef.cutoff,intercept=1,SIMPLIFY=F)
    } else
    {
      initial_beta = mapply(LASSO.fit,y=y,x=x,tau=tau,lambda=lambda,
                            coef.cutoff=coef.cutoff,intercept=1,SIMPLIFY=F)
    }
  }
  
  beta=lapply(initial_beta,function(x)x[-1])
  intval = lapply(y,quantile,probs=tau)
  residuals = mapply(myf11,y=y,x=x,beta=beta,intval=intval,SIMPLIFY=F)
  # intval = initial_beta[1]
  
  K         = as.integer(length(y))
  n         = as.integer(length(y[[1]]))
  y         = lapply(y,as.double)
  p         = as.integer( ncol(x[[1]]) )
  xdouble     = matrix(as.double(unlist(x)),n,p*K)
  tau       = as.double(tau)
  m         = as.integer(length(tau))
  a         = as.double(a)
  eps       = as.double(eps)
  maxin     = as.integer(maxin)
  maxout=as.integer(maxout)
  betadouble=as.double(unlist(beta))
  intvaldouble=as.double(unlist(intval))
  residualsdouble=matrix(as.double(unlist(residuals)),n,m*K)
  lambda    = as.double(lambda)
  
  # groupl1 = rep(0, p)
  
  out=.Fortran("multicqpen",xdouble,betadouble,intvaldouble,residualsdouble,
               n,p,K,tau,m,a,eps,maxin,maxout,lambda,pentype)
  
  out[[2]][abs(out[[2]])<coef.cutoff]=0
  
  coefficients=list()
  
  for (i in 1:K)
  {
    beta[[i]]=out[[2]][((i-1)*p+1):(i*p)]
    intval[[i]]=out[[3]][((i-1)*m+1):(i*m)]
    coefficients[[i]]=c(intval[[i]],beta[[i]])
  }
  
  return(coefficients)
}

myf11=function(y,x,beta,intval) sapply(intval,function(y,x,beta,intval) y-x%*%beta-intval,
                                       y=y,x=x,beta=beta)

##################################

# kernel density
kern=function(x.new,x,b)
{
  n=length(x.new)
  y=numeric(n)
  for(i in 1:n) y[i]=mean(dnorm((x.new[i]-x)/b))/b
  
  return(y)
}

# weighted multiple quantile estimator
wmq=function(theta,n,p,x,y,tau,rinv,inc)
{
  theta0=apply(as.matrix(theta),2,mean)
  #theta0=lm(y~x)$coefficients[-1]
  he=(y-x%*%as.matrix(theta0))[,1]
  hbeta=quantile(he,tau,type=6)
  band=0.9*(n^(-1/5))*min(sd(he),IQR(he,type=6)/1.34)
  hden=kern(x.new=hbeta,x=he,b=band)
  
  f=diag(hden)
  h=f%*%rinv%*%f
  w.theta=(h%*%rep(1,k))[,1]/sum(h)
  owtheta=rep(0,p)
  owtheta[inc]=w.theta%*%as.matrix(theta) 
  
  return(list("owtheta"=owtheta,"hden"=hden))
}

wmq1=function(theta,n,p,x,y,tau,inc)
{
  k=length(tau)
  
  theta0=apply(as.matrix(theta),2,mean)
  #theta0=lm(y~x)$coefficients[-1]
  he=(y-x%*%as.matrix(theta0))[,1]
  hbeta=quantile(he,tau,type=6)
  band=0.9*(n^(-1/5))*min(sd(he),IQR(he,type=6)/1.34)
  hden=kern(x.new=hbeta,x=he,b=band)
  
  
  hh=matrix(0,k,k)
  for(j in 1:k)
  {
    for(m in 1:k) hh[(j),(m)]=(min(tau[j],tau[m])-tau[j]*tau[m])/(hden[j]*hden[m])
  }
  h=solve(hh)
  
  w.theta=as.vector(h%*%rep(1,k))/sum(h)
  owtheta=rep(0,p)
  owtheta[inc]=w.theta%*%as.matrix(theta) 
  
  return(list("owtheta"=owtheta,"hden"=hden))
}
##############################################################

registerDoParallel(56)
##########################
# Basic setting
s=200
n=200
n1=10*n
sigma=sqrt(3)
tau=seq(1/10,9/10,1/10)
beta=qnorm(tau,0,sigma)
k=length(tau)
mid=ceiling(k/2)

r=matrix(0,k,k)
for(j in 1:k)
{
  for(m in 1:k) r[j,m]=(min(tau[j],tau[m])-tau[j]*tau[m])
}
rinv=solve(r)

p=12
#p=400 # model size
vt=rep(0,p)
vt[c(1,2,5)]=c(3,1.5,2)
q=sum(vt!=0)
ind=which(vt!=0)
vt.ind=vt[ind]
cut.off=seq(ceiling(0.5*k),k-1,1)
u=length(cut.off)
##############
# lambda set for p=12
lambda.cqr=exp(seq(-5,5,length=100)) # lambda set for cqr
lambda.qr=exp(seq(-5,5,length=100)) # lambda set for qr 
lambda.ls=exp(seq(-5,5,length=100)) # lambda set for ls
lcqr=length(lambda.cqr)
lqr=length(lambda.qr)
lls=length(lambda.ls)
##############

##############
# lambda set for p=400
#lambda.cqr=exp(seq(-2,0.5,length=20)) # lambda set for cqr
#lambda.qr=exp(seq(-4,-1,length=20)) # lambda set for qr 
#lambda.ls=exp(seq(-4,0,length=20)) # lambda set for ls
#lcqr=length(lambda.cqr)
#lqr=length(lambda.qr)
#lls=length(lambda.ls)
##############
##########################

##################
# Covariance matrix of the covariates
covm=matrix(0,p,p)
for(i in 1:p)
{
  for (j in 1:p) covm[i,j]=0.5^(abs(i-j))
}
###################

out=numeric(16)

error=foreach(i=1:s,.errorhandling="remove") %dopar%
{
  set.seed(i)
  
  ##############################
  # Training data
  x=mvrnorm(n,rep(0,p),covm)
  x.oracle=x[,ind]
  ##############################
  # Various error distributions
  e=rnorm(n,0,sqrt(3))
  #e=rt(n,2)
  #e=rbeta(n,1,3)
  #e=rlap(n)
  #e=rgamma(n,1,1)
  #e=rnormm(n,c(0.5,0.5),c(-2,2),c(1,1))
  #e=rnormm(n,c(0.5,0.5),c(0,0),c(1,0.5^3))
  #e=runif(n,-1,1)
  ##############################
  y=(x%*%vt+e)[,1]
  ##############################
  
  ##############################
  # Validation data
  x.v=mvrnorm(n1,rep(0,p),covm)
  ##############################
  # Various error distributions
  e.v=rnorm(n1,0,sqrt(3))
  #e.v=rt(n1,2)
  #e.v=rbeta(n1,1,3)
  #e.v=rlap(n1)
  #e.v=rgamma(n1,1,1)
  #e.v=rnormm(n1,c(0.5,0.5),c(-2,2),c(1,1))
  #e.v=rnormm(n1,c(0.5,0.5),c(0,0),c(1,0.5^3))
  #e.v=runif(n1,-1,1)
  ##############################
  y.v=(x.v%*%vt+e.v)[,1]
  ##############################
  
  ##############################################################
  # least square regression
  fitls.m=matrix(0,lls,p)
  pels=numeric(lls)
  for(m in 1:lls)
  {
    #fitls=pls(x,y,lambda=lambda.ls[m],intercept=T,penalty="scad",initial="lasso")
    fitls=ncvreg(x,y,family="gaussian",lambda=lambda.ls[m],penalty="SCAD")$beta[,1]
    
    fitls.m[m,]=fitls[-1]
    resid.ls=y.v-fitls[1]-x.v%*%fitls[-1]
    
    #if(length(fitls)==p)
    #{
      #fitls.m[m,]=fitls
      #resid.ls=y.v-x.v%*%fitls
    #} else
    #{
      #fitls.m[m,]=fitls[-1]
      #resid.ls=y.v-fitls[1]-x.v%*%fitls[-1]
    #}
    
    pels[m]=sum(resid.ls^2)
  }
  hvtls=fitls.m[which.min(pels),] # ls estimator
  hvtls.oracle=rep(0,p)
  #hvtls.oracle[ind]=lm(y~x.oracle)$coefficients[-1]
  hvtls.oracle[ind]=ncvreg(x.oracle,y,family="gaussian",lambda=0,penalty="SCAD")$beta[-1,1]
  ##############################################################
  
  ##############################################################
  # composite quantile regression
  fitcqr.m=matrix(0,lcqr,p)
  pecqr=numeric(lcqr)
  resid.cqr=numeric(k)
  for(m in 1:lcqr)
  {
    fitcqr=cqr(list(y),list(x),tau=tau,lambda=lambda.cqr[m],penalty="SCAD")[[1]]
    
    fitcqr.m[m,]=fitcqr[-(1:k)]
    resid.cqr=sapply(1:k,function(z) y.v-fitcqr[z]-x.v%*%fitcqr[-(1:k)])
    
    for(j in 1:k) pecqr[m]=pecqr[m]+sum(resid.cqr[,j]*(tau[j]-(resid.cqr[,j]<0)))
    
  }
  hvtcqr=fitcqr.m[which.min(pecqr),] # cqr estimator
  hvtcqr.oracle=rep(0,p)
  hvtcqr.oracle[ind]=cqr(list(y),list(x.oracle),tau=tau,lambda=0,penalty="SCAD")[[1]][-(1:k)]
  ##############################################################
  
  ##############################################################
  # multiple quantile regressions
  hvtq=matrix(0,nrow=k,ncol=p)
  for(j in 1:k)
  {
    fit.m=matrix(0,lqr,p+1)
    pe=numeric(lqr)
    for(m in 1:lqr)
    {
      #fit=pqr(x,y,tau=tau[j],lambda=lambda.qr[m],intercept=T,penalty="scad",initial="lasso")
      #fit=QICD(y,x,tau=tau[j],lambda=lambda.qr[m],intercept=T,penalty="SCAD")
      fit=cqr(list(y),list(x),tau=tau[j],lambda=lambda.qr[m])[[1]]
      #fit=suppressWarnings(cv.rq.pen(x,y,tau=tau[j],lambda=lambda.qr[m],
      #                               intercept=T,penalty="SCAD"))$models[[1]]$coefficients
      fit.m[m,]=fit
      
      
      resid=y.v-fit[1]-x.v%*%fit[-1]
      pe[m]=sum(resid*(tau[j]-(resid<0)))
    }
    hvtq[j,]=fit.m[which.min(pe),-1]
    
  }
  # hvtq is the matrix of the single quantile estimators.
  
  hvtq.oracle=matrix(0,k,p)
  for(w in 1:k) hvtq.oracle[w,ind]=cqr(list(y),list(x.oracle),tau=tau[w],lambda=0)[[1]][-1]
  #hvtq.oracle[,ind]=t(rq(y~x.oracle,tau=tau)$coefficients)[,-1]
  ##############################################################
  
  ##############################################################
  # vote
  count=apply(hvtq,2,function(x) sum(x!=0))
  ##############################################################
  
  ##############################################################
  # choose the threshold kappa by validation data
  owqhvt.m=matrix(0,nrow=u,ncol=p)
  pe.cv.w=numeric(u)
  
  for(j in 1:u)
  {
    kappa=cut.off[j]
    if (all(count<kappa))
    {
      hvtq.cv=t(suppressWarnings(rq(y~1,tau=tau))$coefficients) # naive estimator
      tvt=apply(hvtq.cv,2,mean)
      he=y
      hbeta=quantile(he,tau,type=6)
      band=0.9*(n^(-1/5))*min(sd(he),IQR(he)/1.34)
      hden=kern(x.new=hbeta,x=he,b=band)
      
      f=diag(hden)
      h=f%*%rinv%*%f
      w.hvtq=(h%*%rep(1,k))[,1]/sum(h)
      v=rep(0,p)
      owqhvt.m[j,]=v
      
      w.ob=(rinv%*%hden)[,1]
      
      resid.cv=apply(hvtq.cv,1,function(z) y.v-z[1])
      pe.cv=sapply(1:k,function(z) sum(resid.cv[,z]*(tau[z]-(resid.cv[,z]<0))))
      pe.cv.w[j]=sum(pe.cv*w.ob) 
    } else
    {
      inc=which(!(count<kappa))
      
      x1=as.matrix(x[,inc])
      x.v1=as.matrix(x.v[,inc])
      hvtq.cv=matrix(0,nrow=k,ncol=length(inc)+1)
      
      for(w in 1:k) hvtq.cv[w,]=cqr(list(y),list(x1),tau=tau[w],lambda=0)[[1]]
      #hvtq.cv=t(rq(y~x1,tau=tau)$coefficients) # estimators with threshold kappa
      v=wmq(hvtq.cv[,-1],n,p,x1,y,tau,rinv,inc)
      #v=wmq1(hvtq.cv[,-1],n,p,x1,y,tau,inc)
      owqhvt.m[j,]=v$owtheta
      
      hden=v$hden
      w.ob=(rinv%*%hden)[,1]
      
      resid.cv=apply(hvtq.cv,1,function(z) y.v-x.v1%*%as.matrix(z[-1])-z[1])
      pe.cv=sapply(1:k,function(z) sum(resid.cv[,z]*(tau[z]-(resid.cv[,z]<0))))
      pe.cv.w[j]=sum(pe.cv*w.ob) 
    }
  }
  
  owqhvt=owqhvt.m[which.min(pe.cv.w),] # optimal multiple quantiles estimator
  
  ##############################################################
  ##############################################################
  # oracle optimal multiple quantiles estimator
  owqhvt.oracle=wmq(hvtq.oracle[,ind],n,p,x.oracle,y,tau,rinv,ind)$owtheta
  #owqhvt.oracle=wmq1(hvtq.oracle[,ind],n,p,x.oracle,y,tau,ind)$owtheta
  
  ##############################################################
  
  out[1]=sum(owqhvt[ind]!=0)
  out[2]=sum(owqhvt[-ind]!=0)
  out[3]=sum((owqhvt-vt)^2)
  out[4]=sum((owqhvt.oracle-vt)^2)
  
  out[5]=sum(hvtq[mid,ind]!=0)
  out[6]=sum(hvtq[mid,-ind]!=0)
  out[7]=sum((hvtq[mid,]-vt)^2)
  out[8]=sum((hvtq.oracle[mid,]-vt)^2)
  
  out[9]=sum(hvtls[ind]!=0)
  out[10]=sum(hvtls[-ind]!=0)
  out[11]=sum((hvtls-vt)^2)
  out[12]=sum((hvtls.oracle-vt)^2)
  
  out[13]=sum(hvtcqr[ind]!=0)
  out[14]=sum(hvtcqr[-ind]!=0)
  out[15]=sum((hvtcqr-vt)^2)
  out[16]=sum((hvtcqr.oracle-vt)^2)
  
  out
}


res=t(matrix(unlist(error),nrow=16,ncol=length(error)))
colnames(res)=c("NC of OWQ","NIC of OWQ","RE of OWQ","RE of OWQ-Oracle",
                "NC of LAD","NIC of LAD","RE of LAD","RE of LAD-Oracle",
                "NC of LS","NIC of LS","RE of LS","RE of LS-Oracle",
                "NC of CQR","NIC of CQR","RE of CQR","RE of CQR-Oracle")


print("#########################Mean#########################")
ave=apply(res,2,mean)
baseline.ave=ave[3]
for(i in 0:3) 
{
  ave[3+i*4]=ave[3+i*4]/baseline.ave
  ave[4+i*4]=ave[4+i*4]/baseline.ave
}
ave

print("#########################SD#########################")
apply(res,2,sd)
