#Author: Sam Lee
#Date: 05.10.2024

#This script is written following "Cluster Robust Inference in Linear Models with Many Covariates" (Cattaneo, Gong, Jansson, & Newey, 2022)

#I define a function that runs a least squares regression and returns estimates for the covariates of interest, X

#Following the paper, the linear regression of interest can be defined as follows:
# Y = X%*%beta + Gamma%*%W + U, where W is a set of non-vanishing covariates

#Used for the DGP
source("clustermanycov.largedata.main.fun.R")

crve.fit <- function(Y, X, W, cluster, beta=0){
  N = length(Y); G = length(unique(cluster)); K = dim(W)[2]; d = dim(as.matrix(X))[2]; n_j = table(cluster)
  n_cumu = vector(); n_cumu[1] = 0; n_cumu_sqr = vector(); n_cumu_sqr[1] = 0
  
  for (i in 1:G){
    n_cumu[i+1] = sum(n_j[1:i])
    n_cumu_sqr[i+1] = sum(n_j[1:i]^{2})
  }
  
  qr.W = qr(W); M = -tcrossprod(qr.Q(qr.W)); diag(M) = 1 + diag(M)
  
  MX = M%*%X;
  XMX = crossprod(MX);
  XMX_inv = solve(XMX);
  XMY = crossprod(MX,Y);
  
  beta.hat = XMX_inv%*%XMY
  
  u.hat = (M%*%Y - MX%*%beta.hat)
  
  #Variance object of covariates of interest
  Sigma_kappa = matrix(0,d^2,G)
  
  #Define Kappa Matrix
  test = try(chol2inv(chol(star(cluster,M))), TRUE)
  if (is.matrix(test)) {
    kappa = test
  } else {
    kappa = pinv(star(dgp$cluster,M))
  }
  
  for (g in 1:G){
    V.hat_g = MX[cluster==g,]
    gamma_g = kronecker(V.hat_g,V.hat_g)
    u.hat_g = u.hat[cluster==g]
    
    if (is.matrix(kappa)){
      var_cr = matrix(0,n_j[g]^{2},G)
      for (h in 1:G){
        kappa_gh = kappa[(n_cumu_sqr[g]+1):n_cumu_sqr[g+1],(n_cumu_sqr[h]+1):n_cumu_sqr[h+1]]
        var_cr[,h] = kappa_gh%*%kronecker(u.hat[cluster==h],u.hat[cluster==h])
      }
      var_g_cr = rowSums(var_cr)
      Sigma_kappa[,g] = t(gamma_g)%*%var_g_cr
    } else {Variance_kappa = NA}
  }
  variance_kappa = rowSums(Sigma_kappa)
  if (is.matrix(kappa)){
    Meat_CR = resize (variance_kappa, d, d)
    V_CR  = (XMX_inv)%*%Meat_CR%*%(XMX_inv)
    T_CR  = (beta.hat - beta) / sqrt(diag(V_CR))
  } else {V_CR = NA; T_CR=NA;}
  
  return(list(
      beta.hat = beta.hat,
      V.CR = V_CR,
      T.CR = T_CR
    )
  )
}

##Example usages
N = 600
q = 150
n = 3
n_prime = 5
nx = 15
qx = floor(N/nx)
Jx = 8
Ju = 2
rho = .5
l = .7
beta = rep(1, 10)
## Nuisance covariates
K.grid = 1:3;
K.grid.LM = c(150,200,250);
dgp = gen.data(q=q,qx=qx,n=n,n_prime=n_prime,nx=nx,Jx=Jx,Ju=Ju,rho=rho,l=l,4,beta=beta, K.i=1, K=K.grid.LM[1]); W = dgp$W
crve.fit(dgp$Y,dgp$X,dgp$W,dgp$cluster, beta=beta)
