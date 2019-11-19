rm(list = ls())
library('PIE')
library('glmnet')
set.seed(123)
p=10000
n=400;
Omega=matrix(0,nrow=p,ncol=p);
Omega[6,6]=1
Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
beta<-rep(0,p);
X=matrix(rnorm(n*p),nrow =n);
Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
## PIE basd on the response
system.time(yOmega<-PIE(X,Y))
##Show the estimation
1+tril(yOmega)@i 
rep(1:p,diff(tril(yOmega)@p))
tril(yOmega)@x