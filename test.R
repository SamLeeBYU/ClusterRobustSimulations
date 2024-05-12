library(devtools)
devtools::install_github("SamLeeBYU/ClusterRobustSimulations")

library(CRVE)

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
