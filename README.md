# PIE
Penalized Interaction Estimation(PIE) for high dimensional quadratic regression

# Reference 
Wang, Cheng, Binyan Jiang, and Liping Zhu. "[Penalized Interaction Estimation for Ultrahigh Dimensional Quadratic Regression](https://arxiv.org/abs/1901.07147)." arXiv preprint arXiv:1901.07147 (2019).

## Getting Started
These instructions will give you a toy example for implementing the package.

### Prerequisites
What things you need to install the software and how to install them

```
install.packages("MASS")
install.packages("Rcpp")
install.packages("Matrix")
install.packages("devtools")
```
### Install PIE

```
library("devtools")
devtools::install_github("cescwang85/PIE")
```

### Toy example 

```
#!/usr/local/bin/Rscript
rm(list = ls())
library('PIE')
library('glmnet')
set.seed(99)
p=100
n=200;
Omega=matrix(0,nrow=p,ncol=p);
Omega[6,6]=1
Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
beta<-rep(0,p);beta[c(1,6,10)]=1;
X=matrix(rnorm(n*p),nrow =n);
Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
##Centering
## PIE basd on the response
yOme<-PIE(X,Y)
##PIE based on the residual
hbeta<-as.vector(coef(cv.glmnet(X,Y,nfolds =5),s="lambda.min"))[-1];  
rOme<-PIE(X,Y-X%*%hbeta)

##Show the estimation
yOme[1:10,1:10]
rOme[1:10,1:10]

```
Overall, the computation complexity for each iteration of the algorithm is linear in both the sample size(n) and the number of parameters (p^2).  Welcome any comments or suggesiongs.


