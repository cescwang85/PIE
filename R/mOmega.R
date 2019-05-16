#' Matrix version solution path function for the Penalized Interaction Estimation (PIE)
#' @param Sn the sample covariance matrix.
#' @param An the sample kernel matrix.
#' @param lambda user supplied tuning parameter; 
#' @param err the precision used to control the convergence. Default is 1e-4.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list of sparse asymmetric p*p matrices, lambda and also the number of iterations.
mOmega<-function(Sn,An,lambda=NULL,nlambda=50,lambda.min.ratio=0.01,err=10^(-4),maxIter=10^3,rho=1)
{p=ncol(An);
if(is.null(lambda)){
  lmax<-max(max(abs(An)));
  lambda<-exp(log(lmax)+seq(log(lambda.min.ratio),0,length.out =nlambda))}
lambda<-sort(lambda);
nlambda=length(lambda);
obj<-eigen(Sn);
U<-obj$vectors;
la<-obj$values;
la=la[la>10^(-7)];
U<-U[,1:length(la)];
NN<-length(lambda);
Rn=rep(0,NN);
Omega_all=NULL;
Z0=matrix(0,p,p);
V0=matrix(0,p,p);
for (k in 1:NN){
  un<-unit(la,U,An,lam=lambda[k],Z0=Z0,V0=V0,err=err,maxIter=maxIter,rho=rho)
  Omega_all<-c(Omega_all,list(as(un[[1]], "sparseMatrix")))
  Rn[k]=un[[3]];
  Z0=un[[1]];
  V0=un[[2]];
}
return(list(Omega=Omega_all,lambda=lambda,niter=Rn))
}