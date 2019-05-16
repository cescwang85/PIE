#' Main function for the Penalized Interaction Estimation (PIE) where the tuning parameter is chosen by BIC
#' @param X covariates data matrix of dimension n*p after centering.
#' @param Y response variable with length n.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min.ratio smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. Default is 0.1.
#' @param err the precision used to stop the convergence. Default is 1e-4.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A sparse asymmetric p*p matrix for the interactions.
#' @examples
#' rm(list = ls())
#'library('PIE')
#'library('glmnet')
#'set.seed(99)
#'p=100
#'n=200;
#'Omega=matrix(0,nrow=p,ncol=p);
#'Omega[6,6]=1
#'Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
#'beta<-rep(0,p);beta[c(1,6,10)]=1;
#'X=matrix(rnorm(n*p),nrow =n);
#'Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
##Centering
## PIE basd on the response
#'yOme<-PIE(X,Y)
##PIE based on the residual
#'hbeta<-as.vector(coef(cv.glmnet(X,Y,nfolds =5),s="lambda.min"))[-1];  
#'rOme<-PIE(X,Y-X%*%hbeta)
#'yOme[1:10,1:10]
#'rOme[1:10,1:10]
PIE<-function(X,Y,lambda=NULL,nlambda=50,lambda.min.ratio=0.1,err=1e-4,maxIter=1000,rho=1)
{
  X<-as.matrix(X);
  X<-scale(X,center =TRUE,scale =FALSE);
  Y<-as.vector(Y);
  n=nrow(X);
  p=ncol(X);
  if(is.null(lambda)){
    An<-t(X)%*%diag(Y-mean(Y))%*%X/n;
    lmax<-max(max(abs(An)));
    lambda<-exp(log(lmax)+seq(log(lambda.min.ratio),0,length.out =nlambda))}
  nlambda=length(lambda);
  obj<-hOmega(X,Y,lambda =lambda,err=err,maxIter=maxIter,rho=rho);
  bic<-rep(0,nlambda);
  for (k in 1:nlambda)
  {Ome<-obj$Omega[[k]];
    if ((nnzero(tril(Ome))>0)&&(nnzero(tril(Ome))<(n/log(n)))){
      Ome1<-tril(Ome);
      idx<-1+Ome1@i
      idy<-rep(1:p,diff(Ome1@p))
      Z=matrix(0,nrow=n,ncol=length(idx));
      for (i in 1:n){
        Xi=X[i,];
        Z[i,]=Xi[idx]*Xi[idy];}
      RSS=deviance(lm(Y~Z))/n;}
    else {RSS=var(Y)}
  bic[k]=n*log(RSS)+nnzero(tril(Ome))*log(n);
  }
id<-max(which(bic==min(bic)));
Omega<-obj$Omega[[id]];
tOmega<-matrix(0,p,p);
if (nnzero(Omega)>0){
idx<-1+tril(Omega)@i
idy<-rep(1:p,diff(tril(Omega)@p))
Z=matrix(0,nrow=n,ncol=length(idx));
for (i in 1:n){
  Xi=X[i,];
  Z[i,]=Xi[idx]*Xi[idy];}  
tOmega[idx+p*(idy-1)]=lm(Y~Z)$coefficients[-1];}
return(Matrix((tOmega+t(tOmega))/2, sparse = TRUE))
}
