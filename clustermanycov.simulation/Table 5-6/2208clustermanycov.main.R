#########################################################################################
## R CODE FOR CLUSTER MANY COVARIATES
## Cluster Robust Inference in Linear Regression Models with Many Covariates
## DATE: 14-August-2022
## Main Simulation Code
#########################################################################################

rm(list=ls(all=TRUE))

#load packages
require("ramify")
require("expm")
require("pracma")

#Set path and load function
setwd("clustermanycov.simulation/Table 5-6")
source("2208clustermanycov.main.fun.R")

# Set seed
set.seed(666)

#########################################################################################
## MONTECARLO SETUP
#########################################################################################
## Number of simulations (S)
S = 1000
## Sample size (N)
N = 600
## Cluster size (n) for homogeneous Cluster size
q = 150
n = 3
n_prime = 5
nx = 15
qx = floor(N/nx)
Jx = 8
Ju = 2
rho = .5
l = .7
## d = dim(beta)
d = 1; beta = 1
## Nuisance covariates
K.grid = 1:5;
K.grid.LM = c(0,25,50,75,100);
models = 5:6
## Load population scaling factors
popval = read.csv("main.popval.csv", row.names=1)
## Models to run:we have four designs, independent and dependent regressors
## homo and hete cluster sizes so in total we have models=1:4.


#########################################################################################
# Run Monte Carlo Experiment
#########################################################################################
# m = 1
for (m in models){
  ## Output Tables
  col.names1 = c("N","d","G","beta","m","s","K","rank","ratio")
  col.names2 = c("beta.hat","V_LZ","V_HC1","V_HC2","V_HC3","V_LOO","V_CR","Neg",
                 "T_LZ","T_HC1","T_HC2","T_HC3","T_LOO","T_CR")
  col.names=c(col.names1,col.names2)
  out = matrix(NA, nrow=S*length(K.grid), ncol=length(col.names), dimnames=list(NULL,col.names))
  
  
  message("Simulations began for Model ", m,". Time: ", Sys.time()); showevery=.01
  row=1;
  
  # s=1; K.i=3;
  for (s in 1:S) {
    if (max(s==seq(0,S,S*showevery))==1) {message("Simulations Completed: ",s," of ",S," (",round(s/S*100),"%) - ", Sys.time())}
    for (K.i in K.grid) {
      if ( m == 1 | m == 3) {dgp = gen.data(q=q,qx=qx,n=(n+n_prime)/2,n_prime=(n+n_prime)/2,nx=nx,Jx=Jx,Ju=Ju,rho=rho,l=l,m,beta=beta ,K.i=K.i, K=K.grid.LM[K.i]); W = dgp$W}
      if ( m == 2 | m == 4 | m == 5 | m == 6) {dgp = gen.data(q=q,qx=qx,n=n,n_prime=n_prime,nx=nx,Jx=Jx,Ju=Ju,rho=rho,l=l,m,beta=beta ,K.i=K.i, K=K.grid.LM[K.i]); W = dgp$W}
      
      K = ncol(W); qr.W = qr(W); if (qr.W$rank==K){
        M = -tcrossprod(qr.Q(qr.W)); diag(M) = 1 + diag(M);
        Mn = 1-min(diag(M))
        if (m == 1| m == 3) {ratio = Mn*((n+n_prime)/2)^{2}}
        if (m == 2| m == 4 | m == 5 | m == 6) {ratio = Mn*(max(n,n_prime))^{2}}
        test = try(chol2inv(chol(star(dgp$cluster,M))), TRUE)
        if (is.matrix(test)) {
          kappa = test
        } else {
          kappa = pinv(star(dgp$cluster,M))
        }
        out[row,1:length(col.names)] = c(N,d,q,beta,m,s,K,qr.W$rank,ratio,lsfit(dgp$Y,dgp$X,dgp$cluster,M,kappa,beta,K));
      }
      
      row=row+1;
    }
  }
  ## Save final table
  filename = paste0("output/output_m",m,".csv");
  if (file.exists(filename)) write.table(out, file=filename, append=TRUE, col.names=FALSE, row.names=FALSE, sep=",") else write.table(out, file=filename, row.names=FALSE, sep=",")
}