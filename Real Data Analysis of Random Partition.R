#install.packages("glmnet",repos="https://cran.revolutionanalytics.com/")
#install.packages("rqPen",repos="https://cran.revolutionanalytics.com/")
#install.packages("doParallel",repos='https://cran.revolutionanalytics.com/')
#install.packages("cqrReg",repos='https://cran.revolutionanalytics.com/')
#install.packages("ncvreg",repos="https://cran.revolutionanalytics.com/")
#install.packages("quantreg",repos="https://cran.revolutionanalytics.com/")


suppressMessages({
  library(MASS)
  library(rqPen)
  library(doParallel)
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

registerDoParallel(detectCores())

##################
# data
load("eye.RData") # data
##################

##########################
# Basic setting
s=50
n=80
tau=seq(1/10,9/10,1/10)
k=length(tau)
cut.off=seq(ceiling(0.5*k)+2,k,1) # threshold
u=length(cut.off) 
p=ncol(x.raw)

r=matrix(0,k,k)
for(j in 1:k)
{
  for(m in 1:k) r[j,m]=(min(tau[j],tau[m])-tau[j]*tau[m])
}
rinv=solve(r)
##############


##############
# lambda set 
lambda.cqr=exp(seq(-0.8,0,length=10)) # lambda set for cqr
lambda.qr=exp(seq(-5.4,-2.7,length=50)) # lambda set for owqr 
lambda.ls=exp(seq(-5.3,-4.5,length=10)) # lambda set for ls
lambda.lad=exp(seq(-3.3,-2.5,length=10)) # lambda set for lad
lcqr=length(lambda.cqr)
lqr=length(lambda.qr)
lls=length(lambda.ls)
llad=length(lambda.lad)
##############
##########################

ve=numeric(12)

res=foreach(uu=1:s,.errorhandling="remove") %dopar%
  {
    set.seed(77*uu)
    train=sample(1:length(y.raw),n)
    
    ############
    # training data
    xx.raw=x.raw[train,]
    yy.raw=y.raw[train]
    ############
    
    ############
    # validation data
    xx.v=x.raw[-train,]
    yy.v=y.raw[-train]
    ############
    
    
    
    ##################
    # cross validation
    cv=5 # number of folds
    n1=length(yy.raw)/cv # number observations in each fold
    set.seed(77*uu)
    ind.cv=sample(1:n) # rearrange the observations for cross-validation
    ##################
    
    
    ##############
    # least square
    error1=foreach(i=1:lls,.errorhandling="remove") %do%
      {
        pels=0
        for(j in 1:cv)
        {
          ind.v=ind.cv[((j-1)*n1+1):(j*n1)]
          x=xx.raw[-ind.v,]
          y=yy.raw[-ind.v]
          x.v=xx.raw[ind.v,]
          y.v=yy.raw[ind.v]
          
          fitls=ncvreg(x,y,family="gaussian",lambda=lambda.ls[i],penalty="SCAD")$beta[,1]
          resid.ls=y.v-fitls[1]-x.v%*%fitls[-1]
          
          pels=pels+sum(resid.ls^2)
        }
        pels
      }
    hvtls=ncvreg(xx.raw,yy.raw,family="gaussian",
                 lambda=lambda.ls[which.min(error1)],penalty="SCAD")$beta[,1]
    perror.ls=yy.v-hvtls[1]-xx.v%*%hvtls[-1]
    ##############
    
    ##############
    # composite quantile
    error2=foreach(i=1:lcqr,.errorhandling="remove") %do%
      {
        pecqr=0
        for(j in 1:cv)
        {
          ind.v=ind.cv[((j-1)*n1+1):(j*n1)]
          x=xx.raw[-ind.v,]
          y=yy.raw[-ind.v]
          x.v=xx.raw[ind.v,]
          y.v=yy.raw[ind.v]
          
          fitcqr=cqr(list(y),list(x),tau=tau,lambda=lambda.cqr[i],penalty="SCAD")[[1]]
          resid.cqr=sapply(1:k,function(z) y.v-fitcqr[z]-x.v%*%fitcqr[-(1:k)])
          
          for(m in 1:k) pecqr=pecqr+sum(resid.cqr[,m]*(tau[m]-(resid.cqr[,m]<0)))
        }
        pecqr
      }
    hvtcqr=cqr(list(yy.raw),list(xx.raw),tau=tau,
               lambda=lambda.cqr[which.min(error2)],penalty="SCAD")[[1]]
    perror.cqr=sapply(1:k,function(z) yy.v-hvtcqr[z]-xx.v%*%hvtcqr[-(1:k)])
    ##############
    
    ##############
    # least absolute deviation
    error4=foreach(m=1:llad,.errorhandling="pass")%do%
      {
        pe=0
        
        for(l in 1:cv)
        {
          ind.v=ind.cv[((l-1)*n1+1):(l*n1)]
          x=xx.raw[-ind.v,]
          y=yy.raw[-ind.v]
          x.v=xx.raw[ind.v,]
          y.v=yy.raw[ind.v]
          
          fit=cqr(list(y),list(x),tau=0.5,lambda=lambda.lad[m])[[1]]
          resid=y.v-fit[1]-x.v%*%fit[-1]
          pe=pe+sum(resid*(0.5-(resid<0)))
        }
        pe
      }
    
    for (tt in 1:llad)
    {
      if(!is.numeric(error4[[tt]])) error4[[tt]]=100
    }
    
    hvtlad=cqr(list(yy.raw),list(xx.raw),tau=0.5,
               lambda=lambda.lad[which.min(error4)],penalty="SCAD")[[1]]
    perror.lad=(yy.v-hvtlad[1]-xx.v%*%hvtlad[-1])
    ##############
    
    
    ##############
    # multiple quantiles variable selection
    hvtq=matrix(0,nrow=k,ncol=p+1)
    for(j in 1:k)
    {
      error3=foreach(m=1:lqr,.errorhandling="pass")%do%
        {
          pe=0
          
          for(l in 1:cv)
          {
            ind.v=ind.cv[((l-1)*n1+1):(l*n1)]
            x=xx.raw[-ind.v,]
            y=yy.raw[-ind.v]
            x.v=xx.raw[ind.v,]
            y.v=yy.raw[ind.v]
            
            fit=cqr(list(y),list(x),tau=tau[j],lambda=lambda.qr[m])[[1]]
            resid=y.v-fit[1]-x.v%*%fit[-1]
            pe=pe+sum(resid*(tau[j]-(resid<0)))
          }
          pe
        }
      
      for (tt in 1:lqr)
      {
        if(!is.numeric(error3[[tt]])) error3[[tt]]=100
      }
      
      hvtq[j,]=cqr(list(yy.raw),list(xx.raw),tau=tau[j],
                   lambda=lambda.qr[which.min(error3)],penalty="SCAD")[[1]]
      
    }
    hvtq=hvtq[,-1]
    ##############
    
    
    
    ##############################################################
    # vote
    count=apply(hvtq,2,function(x) sum(x!=0))
    for(i in 1:u) 
    {
      if(!(sum(!(count<cut.off[i]))<n)) cut.off[i]=0
    }
    if(any(cut.off==0)) cut.off=cut.off[-which(cut.off==0)]  
    u=length(cut.off)
    ##############################################################
    
    ##############################################################
    # choose the threshold kappa by validation data
    owqhvt.m=matrix(0,nrow=u,ncol=p)
    
    
    pe.vote=foreach(j=1:u)%do%
      {
        pe.cv.w=0
        pe.cv.lad=0
        pe.cv.cqr=0
        
        kappa=cut.off[j]
        for (l in 1:cv)
        {
          ind.v=ind.cv[((l-1)*n1+1):(l*n1)]
          x=xx.raw[-ind.v,]
          y=yy.raw[-ind.v]
          x.v=xx.raw[ind.v,]
          y.v=yy.raw[ind.v]
          
          if (all(count<kappa))
          {
            pe.cv.w=pe.cv.w+sum((y.v-mean(y))^2) 
            pe.cv.lad=pe.cv.lad+sum(abs(y.v-median(y)))
            for (ii in 1:k) 
            {
              qq=quantile(y,tau[ii],type=6)
              pe.cv.cqr=pe.cv.cqr+sum((y.v-qq)*(tau[ii]-(y.v<qq)))
            }
              
          } else
          {
            inc=which(!(count<kappa))
            
            x1=as.matrix(x[,inc])
            x.v1=as.matrix(x.v[,inc])
            
            lsv.cv=lm(y~x1)$coefficients
            ladv.cv=rq(y~x1)$coefficients
            cqrv.cv=cqr(list(y),list(x1),tau,lambda=0)[[1]]
            
            pe.cv.w=pe.cv.w+sum((y.v-lsv.cv[1]-x.v1%*%lsv.cv[-1])^2) 
            pe.cv.lad=pe.cv.lad+sum(abs(y.v-ladv.cv[1]-x.v1%*%ladv.cv[-1])) 
            for(ii in 1:k) 
            {
              res.cqrv=y.v-cqrv.cv[ii]-x.v1%*%cqrv.cv[-(1:k)]
              pe.cv.cqr=pe.cv.cqr+sum(res.cqrv*(tau[ii]-(res.cqrv<0)))
            }
              
          }
        }
        c(pe.cv.w,pe.cv.lad,pe.cv.cqr)
      }
    
    pe.vote.m=matrix(unlist(pe.vote),nrow=3,ncol=u)
    kappa=cut.off[which.min(pe.vote.m[1,])]
    kappa.lad=cut.off[which.min(pe.vote.m[2,])]
    kappa.cqr=cut.off[which.min(pe.vote.m[3,])]
    
    if (all(count<kappa))
    {
      perror.lsv=yy.v-mean(yy.raw)
    } else
    {
      inc=which(!(count<kappa))
      
      x1=as.matrix(xx.raw[,inc])
      ls.v=lm(yy.raw~x1)$coefficients
      perror.lsv=yy.v-ls.v[1]-as.matrix(xx.v[,inc])%*%ls.v[-1]
    }
    
    if (all(count<kappa.lad))
    {
      perror.ladv=yy.v-median(yy.raw)
    } else
    {
      inc=which(!(count<kappa.lad))
      
      x1=as.matrix(xx.raw[,inc])
      lad.v=rq(yy.raw~x1)$coefficients
      perror.ladv=(yy.v-lad.v[1]-as.matrix(xx.v[,inc])%*%lad.v[-1])
    }
    
    if (all(count<kappa.cqr))
    {
      perror.cqrv=sapply(1:k,function(z) yy.v-quantile(yy.raw,tau[z],type=6))
    } else
    {
      inc=which(!(count<kappa.cqr))
      
      x1=as.matrix(xx.raw[,inc])
      cqr.v=cqr(list(yy.raw),list(x1),tau,lambda=0)[[1]]
      perror.cqrv=sapply(1:k,function(z) yy.v-cqr.v[z]-as.matrix(xx.v[,inc])%*%cqr.v[-(1:k)])
    }
    
    
    ##############################################################
    
    ve[1]=sum(!(count<kappa))
    ve[2]=sum(perror.lsv^2)
    
    ve[3]=(sum(hvtls!=0))-1
    ve[4]=sum(perror.ls^2)
    
    ve[5]=sum(!(count<kappa))
    ve[6]=sum(abs(perror.ladv))
    
    ve[7]=(sum(hvtlad!=0))-1
    ve[8]=sum(abs(perror.lad))
    
    ve[9]=sum(!(count<kappa))
    ve[10]=0
    for(ii in 1:k) ve[10]=ve[10]+sum(perror.cqrv[,ii]*(tau[ii]-(perror.cqrv[,ii]<0)))
    
    ve[11]=(sum(hvtcqr!=0))-k
    ve[12]=0
    for(ii in 1:k) ve[12]=ve[12]+sum(perror.cqr[,ii]*(tau[ii]-(perror.cqr[,ii]<0)))
    
    ve
  }

res.m=t(matrix(unlist(res),nrow=length(ve),ncol=length(res)))
colnames(res.m)=c("selection of lsv","pe of lsv",
                  "selection of ls","pe of ls",
                  "selection of ladv","pe of ladv",
                  "selection of lad","pe of lad",
                  "selection of cqrv","pe of cqrv",
                  "selection of cqr","pe of cqr")

print("#########################Mean#########################")
apply(res.m,2,mean)

print("#########################SD#########################")
apply(res.m,2,sd)


