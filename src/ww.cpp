#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List unit(arma::vec eigd,arma::mat U,arma::mat An,double lambda,Rcpp::Nullable<arma::mat> Z0=R_NilValue, 
                Rcpp::Nullable<arma::mat> V0 = R_NilValue,double err=1e-4,int maxIter=1e3,double rho=1){
  int ee;
  int p=An.n_rows;
  arma::mat B=(2*eigd*eigd.t())/(2*eigd*eigd.t()+rho);
  arma::mat Z=arma::zeros(p,p);
  arma::mat V=arma::zeros(p,p);
  arma::mat Z1;
  arma::mat X;  
  arma::mat L;
  if (Z0.isNotNull()) {Z=Rcpp::as<arma::mat>(Z0);}
  if (V0.isNotNull()) {V=Rcpp::as<arma::mat>(V0);}
  int i=0;
  while (((i<maxIter)&&(ee<1))|(i==0))
  {Z1=Z;
    L=An/rho+Z-V;
    X=L-U*(B%(U.t()*L*U))*U.t();
    Z=(X+V>lambda/rho)%(X+V-lambda/rho)+(X+V<-(lambda/rho))%(X+V+lambda/rho);
    V=V+X-Z;
    ee=(rho*norm(Z-Z1,"fro")/(p+rho*10*norm(V,"fro"))<err)*(norm(X-Z,"fro")/(p+10*(norm(Z,"fro")+norm(X,"fro"))/2)<err);
    i=i+1;
  }
  return Rcpp::List::create(Rcpp::Named("Z") =Z,
                            Rcpp::Named("V") =V,
                            Rcpp::Named("i") =i);
}
