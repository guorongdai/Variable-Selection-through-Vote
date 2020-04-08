# install.packages("Biobase")
# install.packages("GEOquery")

library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE5680", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

dat=t(exprs(gset))
x0=dat[,c(1:21230,21232:31042)]
y.raw=dat[,21231]
thr=quantile(dat,0.25)

ind=numeric(0)
for(i in 1:ncol(x0)) 
{
  xx=x0[,i]
  if((max(xx)>thr) & ((max(xx)-min(xx))>1)) ind=c(ind,i)
}
x1=x0[,ind]

corr=numeric(0)
for(i in 1:ncol(x1)) corr=c(corr,abs(cor(x1[,i],y)))
x.raw=x1[,order(corr,decreasing=T)[1:300]]

save(x.raw,y.raw,file="eye.RData")

