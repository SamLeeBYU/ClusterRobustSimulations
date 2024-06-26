#########################################################################################
## R CODE FOR CLUSTER MANY COVARIATES
## Cluster Robust Inference in Linear Regression Models with Many Covariates
## DATE: 14-August-2022
## FUNCTIONS
#########################################################################################

rm(list=ls(all=TRUE))

require("ramify")
require("expm")
require("pracma")
path="clustermanycov.simulation/Table 5-6"

rmixnorm = function(n,m1,m2,s1,s2,alpha) {I = runif(n)<alpha; rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2));}

#--------------------------------------------
#Setup sample size for cluster design
# n=n_prime represents homogeneous cluster size while different values represent 
# heterogeneous cases. N represents total sample size and q represents the number
# of clusters.
cluster.size<-function(q,n,n_prime){
  half.1= floor(q/2);
  half.2= q-floor(q/2);
  if (n==n_prime){
    N=n*q;
    n_j = c(rep(n,q));
  } else {
    n_j = c(rep(n,each=half.1),rep(n_prime,each=half.2));    
    N = sum(n_j);
  }
  return(N)
}

#Setup function to build dependent vectors and errors based on random factors#
# Ng represents the cluster size for each structure, J is the number of 
# unobservable random factors rho is the correlation coefficients and l is the 
# selected parameter. The difference is that covariates mainly follow uniform 
# distribution while the error term evolves from a standard normal. Further, 
# although we use J to represent both random factors, the dimension of J is 
# allowed to be different between errors and regressors.
reg.dep<-function(Ng,J,rho,l){
  lambda = rep(0,J)
  lambda[1] = runif(1,min=-1,max=1)
  for ( i in 2:J ) {
    lambda[i] = rho*lambda[i-1]+sqrt(1-rho^2)*runif(1,min=-1,max=1)
  }
  L = matrix(0,Ng,J)
  for ( i in 1: Ng){
    for (j in 1: J){
      L[i,j] = l*ifelse(j==floor((i-1)*J/Ng)+1,1,0)
    }
  }
  w = L%*%lambda + sqrt(1-l^2)*runif(Ng,min=-1,max=1)
  return(w)
}

err.dep<-function(Ng,J,rho,l){
  lambda = rep(0,J)
  lambda[1] = rnorm(1)
  for ( i in 2:J ) {
    lambda[i] = rho*lambda[i-1]+sqrt(1-rho^2)*rnorm(1)
  }
  L = matrix(0,Ng,J)
  for ( i in 1: Ng){
    for (j in 1: J){
      L[i,j] = l*ifelse(j==floor((i-1)*J/Ng)+1,1,0)
    }
  }
  w = L%*%lambda
  return(w)
}

err.dep.asy<-function(Ng,J,rho,l){
  lambda = rep(0,J)
  lambda[1] = rmixnorm(1,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2)
  for ( i in 2:J ) {
    lambda[i] = rho*lambda[i-1]+sqrt(1-rho^2)*rmixnorm(1,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2)
  }
  L = matrix(0,Ng,J)
  for ( i in 1: Ng){
    for (j in 1: J){
      L[i,j] = l*ifelse(j==floor((i-1)*J/Ng)+1,1,0)
    }
  }
  w = L%*%lambda
  return(w)
}

err.dep.bimo<-function(Ng,J,rho,l){
  lambda = rep(0,J)
  lambda[1] = rmixnorm(1,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2)
  for ( i in 2:J ) {
    lambda[i] = rho*lambda[i-1]+sqrt(1-rho^2)*rmixnorm(1,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2)
  }
  L = matrix(0,Ng,J)
  for ( i in 1: Ng){
    for (j in 1: J){
      L[i,j] = l*ifelse(j==floor((i-1)*J/Ng)+1,1,0)
    }
  }
  w = L%*%lambda
  return(w)
}

fixed.effects = function(N,T,G) {
  ## Mapping: n = N*T | K = (N+G)/n | G <= N-1 | T >= 3
  GS = floor(N/G)
  out = matrix(0,N*T,N+G)
  ## Unit fixed effects
  out[,1] = rep.int(1,N*T)
  if (N>1) for (j in 2:N) out[((j-1)*T+1):((j-1)*T+T),j] = 1
  
  ## Group fixed effects
  if (G>0) for (j in 1:G) out[seq.int(1+(j-1)*GS*T,(j<G)*j*GS*T+(j==G)*N*T,by=T),N+j] = 1
  return(out)
}

#-------------------------------------------------------------------
gen.data <-function(q,qx,n,n_prime,nx,Jx,Ju,rho,l,m,beta,K.i,K,varkappa.v=NULL,varkappa.u=NULL){
  #-------------------------------------------------------------------
  # Explain: Generates data for the simulations
  #-------------------------------------------------------------------
  # INPUTS: - q: the number of clusters
  #         - qx: the number of cluster for regressors
  #         - n and n_prime: the reference sample size within cluster
  #         - nx: the reference sample size for regressors
  #         - Jx and Ju: Random Vector dimensions for regressors and error terms
  #         - rho: correlation coefficient
  #         - l: weighted parameters in D.G.P
  #         - K: the number of covariates
  #         - m: sets the dgp(1,2 for independent covariates with homo and hete
  #                           cluster size;3,4 for dependent covariates with 
  #                           homo and hete cluster size)
  #------------------------------------------------------------------
  # RETURNS: df = (cluster, unit, Yij, Zij)
  #   Clust Unit Yij       Zij      
  #1    1    1  
  #2    1    2 
  #3    1    3  
  #4    1    4  
  #5    2    1 
  #6    2    2  
  #------------------------------------------------------------------
  #-------------------------------------------------------------------
  #Useful Functions
  if(is.null(varkappa.v)) varkappa.v = popval[is.na(popval[,"m"])==FALSE & popval[,"m"]==m & popval[,"K.i"]==K.i , "varkappa.v"];
  if(is.null(varkappa.u)) varkappa.u = popval[is.na(popval[,"m"])==FALSE & popval[,"m"]==m & popval[,"K.i"]==K.i , "varkappa.u"];
  trunc = function(x){c=2; return(-c*(x < -c) + x*(abs(x) <= c) + c*(x > c))};
  # Setup sample sizes per cluster
  half.1= floor(q/2);
  half.2= q-floor(q/2);
  
  if (n==n_prime){
    N=n*q;
    n_j = c(rep(n,q));
  } else {
    n_j = c(rep(n,each=half.1),rep(n_prime,each=half.2));    
    N = sum(n_j);
  }
  nx_j = c(rep(nx,qx))
  
  # Cluster indicator (n>0 for n_j=n_j' - n=0 for n_j != n_j')
  if (n==n_prime){
    cluster = rep(c(1:q),each=n);
    unit    = rep(c(1:n),q);
  } else {
    cluster  = rep(c(1:q),n_j)
    unit.1    = rep(c(1:n),half.1);
    unit.2    = rep(c(1:n_prime),half.2);
    unit    = c(unit.1,unit.2);
  }
  unitx  = rep(c(1:nx),qx)
  
  # Generate regressors
  if ( m == 1 | m == 2 ){
    W = matrix(1,N,1+K); if (K>0) {
      W[,2:(K+1)] = matrix(runif(N*K,min=-1,max=1),N);
    }
  }
  
  if ( m == 3 | m == 4){
    W = matrix(1,N,1+K); if (K>0) {
      for (t in 2:(K+1)){
        W[1:nx_j[1],t] = reg.dep(nx_j[1],Jx,rho,l)
        for (i in 2:qx) {
          W[(sum(nx_j[1:(i-1)])+1):sum(nx_j[1:i]),t] = reg.dep(nx_j[i],Jx,rho,l)
        }
      }
    }
  }
  
  if (m == 5 | m == 6 ){
    W = matrix(1,N,(1+K*1.6)); if (K>0) {
      W = W[,1:(1.6*K)]
      W[,(K+2):(1.6*K)]=fixed.effects(N,1,0.6*K)[,(N+2):(N+0.6*K)]
      for (t in 2:(K+1)){
        W[1:nx_j[1],t] = reg.dep(nx_j[1],Jx,rho,l)
        for (i in 2:qx) {
          W[(sum(nx_j[1:(i-1)])+1):sum(nx_j[1:i]),t] = reg.dep(nx_j[i],Jx,rho,l)
        }
      }
    }
  }
  
  if ( m == 1 | m == 2 | m == 3 | m == 4){
  Wg = W%*%rep.int(1,1+K)
  };
  
  if ( m == 5 | m == 6){
    if (K>0){
      Wg = W%*%rep.int(1,1.6*K)
    }
    else {
      Wg = W%*%rep.int(1,1+1.6*K)
    }
  };
  
  sigma.x.fn = function(Wg) return((1+Wg^2));
  sigma.u.fn = function(Xb,Wg) return((trunc(Xb)+Wg)^2);
  
  v = matrix(rnorm(N),N,1)
  V.HE = sqrt(varkappa.v * sigma.x.fn(Wg)) * v
  X = V.HE
  
  #Generate Error terms
  
  Ud = matrix(0,N,1)
  
  if ( m == 1 | m == 2 | m == 3 | m == 4){
  
    Ud[1:n_j[1]] = err.dep(n_j[1],Ju,rho,l)
    for (i in 2:q) {
      Ud[(sum(n_j[1:(i-1)])+1):sum(n_j[1:i])] = err.dep(n_j[i],Ju,rho,l)
    }
    U = Ud + sqrt(1-l^2)*sqrt(varkappa.u*sigma.u.fn(X,Wg))*matrix(rnorm(N),N,1)
  }
  
  if ( m == 5 ){
    Ud[1:n_j[1]] = err.dep.asy(n_j[1],Ju,rho,l)
    for (i in 2:q) {
      Ud[(sum(n_j[1:(i-1)])+1):sum(n_j[1:i])] = err.dep.asy(n_j[i],Ju,rho,l)
    }
    U = Ud + sqrt(1-l^2)*sqrt(varkappa.u*sigma.u.fn(X,Wg))*matrix(rnorm(N),N,1)
  }
  
  if ( m == 6 ){
    Ud[1:n_j[1]] = err.dep.bimo(n_j[1],Ju,rho,l)
    for (i in 2:q) {
      Ud[(sum(n_j[1:(i-1)])+1):sum(n_j[1:i])] = err.dep.bimo(n_j[i],Ju,rho,l)
    }
    U = Ud + sqrt(1-l^2)*sqrt(varkappa.u*sigma.u.fn(X,Wg))*matrix(rnorm(N),N,1)
  }
  
  ;
  
  # Generate Outcome
  
  Y = beta * X + U
  # Return the data frame
  list  = list(cluster = cluster,unit = unit, unitx = unitx, Y = Y, X = X, V.HE = V.HE, W = W, U = U);
  return(list);
}

## DGPs: Population Scaling
#########################################################################################
## setwd("D:/Research Project/ClusterManyCov/Simulation and Application 2208/Simulation/"); 
## source("2208clustermanycov.main.fun.R"); gen.table(N=,models=1:4,K.grid=1:5)

gen.table = function(q,n,n_prime,models,K.grid){
  ## NOTE: this function uses environment variables n (sample size) and K.grid (grid for K)
  dimnames = list(NULL,c("m","N","n","n_prime","K.i","K","varkappa.v","varkappa.u"))
  popval = matrix(NA,nrow=(length(models)*length(K.grid)),ncol=length(dimnames[[2]]),dimnames=dimnames)
  
  I = 6000 ; row = 1 ; nx.i = 15 ; Jx = 8 ; Ju = 2 ; rho = 0.5 ; l = 0.7
  
  for (m in models){
    message("\nComputing Scaling Constants for Model ",m,".")
    for (K.i in K.grid){
      
      if (m == 1 | m == 3 ) {
        n.i = (n+n_prime)/2
        N = cluster.size(I,n.i,n.i)
        qx.i = N/nx.i
        N_prime = cluster.size(q,n.i,n.i)
      }
      
      if (m == 2 | m == 4 | m == 5 | m == 6) {
        N = cluster.size(I,n,n_prime)
        qx.i = N/nx.i
        N_prime = cluster.size(q,n,n_prime)
      } 

      K.grid.LM = c(0,25,50,75,100);
      K = K.grid.LM[K.i]
      
      if (m == 1 | m == 3 ) {
        dgp = gen.data(q=I,qx=qx.i,n=n.i,n_prime=n.i,nx=nx.i,Jx=Jx,Ju=Ju,rho=rho,l=l,m=m,K.i=K.i,K=K,varkappa.v=1,varkappa.u=1 , beta=1)
        sigma2.v = as.numeric(crossprod(dgp$V.HE[dgp$unitx==1])/qx.i);
      
        dgp = gen.data(q=I,qx=qx.i,n=n.i,n_prime=n.i,nx=nx.i,Jx=Jx,Ju=Ju,rho=rho,l=l,m=m,K.i=K.i,K=K,varkappa.v=1,varkappa.u=1 , beta=1)
        sigma2.u = as.numeric(crossprod(dgp$U[dgp$unit==1])/I);
      
      popval[row,] = c(m,N_prime,n.i,n.i,K.i,K+1,1/sigma2.v,1/sigma2.u)
      
      row=row+1
      } 
      if (m == 2 | m == 4 | m == 5 | m == 6) {
        dgp = gen.data(q=I,qx=qx.i,n=n,n_prime=n_prime,nx=nx.i,Jx=Jx,Ju=Ju,rho=rho,l=l,m=m,K.i=K.i,K=K,varkappa.v=1,varkappa.u=1 , beta=1)
        sigma2.v = as.numeric(crossprod(dgp$V.HE[dgp$unitx==1])/qx.i);
        
        dgp = gen.data(q=I,qx=qx.i,n=n,n_prime=n_prime,nx=nx.i,Jx=Jx,Ju=Ju,rho=rho,l=l,m=m,K.i=K.i,K=K,varkappa.v=1,varkappa.u=1 , beta=1)
        sigma2.u = as.numeric(crossprod(dgp$U[dgp$unit==1])/I);
        
        popval[row,] = c(m,N_prime,n,n_prime,K.i,K+1,1/sigma2.v,1/sigma2.u)
        
        row=row+1
      } 
    }
    write.csv(popval, file=paste0(path,"main.popval.csv"))
  }
}

#########################################################################################
## LSFIT
#########################################################################################
lsfit = function(Y,X,cluster,M,kappa,beta,K) {
  N = length(Y)  ; G = length(unique(cluster)) ; d = dim(as.matrix(X))[2] ; n_j = table(cluster)
  n_cumu = vector() ; n_cumu[1] = 0 ; n_cumu_sqr = vector() ; n_cumu_sqr[1]=0
  
  for (i in 1:G){
    n_cumu[i+1] = sum(n_j[1:i])
    n_cumu_sqr[i+1] = sum(n_j[1:i]^{2})
  }
  
  MX = M%*%X;
  XMX = crossprod(MX);
  XMX_inv = solve(XMX);
  XMY = crossprod(MX,Y);
  
  beta.hat = XMX_inv%*%XMY
  
  u.hat = (M%*%Y - MX%*%beta.hat)
  u2.hat = u.hat^2
  
  Sigma = matrix(0,d^2,G)
  Sigma_HC2 = matrix(0,d^2,G)
  Sigma_HC3 = matrix(0,d^2,G)
  Sigma_LOO = matrix(0,d^2,G)
  Sigma_kappa = matrix(0,d^2,G)
  
  for (g in 1:G){
    V.hat_g = MX[cluster==g]
    gamma_g = kronecker(V.hat_g,V.hat_g)
    u.hat_g = u.hat[cluster==g]
    Y_g = Y[cluster==g]
    M_g = M[(n_cumu[g]+1):n_cumu[g+1],(n_cumu[g]+1):n_cumu[g+1]]
    
    if (det(M_g)>10^{-7}){
      M_g_inv = solve(M_g)
    } else {
      M_g_inv = pinv(M_g)
    }
    
    #Start building different version of clustered error
    #LZ-standard error
    var_g = kronecker(u.hat_g,u.hat_g)
    Sigma[,g] = t(gamma_g)%*%var_g
    
    #HC-2
    var_g_HC2= (kronecker(u.hat_g,M_g_inv%*%u.hat_g)+kronecker(M_g_inv%*%u.hat_g,u.hat_g))/2
    Sigma_HC2[,g] = t(gamma_g)%*%var_g_HC2
    
    #HC-3
    var_g_HC3 = kronecker(M_g_inv%*%u.hat_g,M_g_inv%*%u.hat_g)
    Sigma_HC3[,g] = t(gamma_g)%*%var_g_HC3
    
    
    #LOO
    var_g_LOO = (kronecker(Y_g,M_g_inv%*%u.hat_g)+kronecker(M_g_inv%*%u.hat_g,Y_g))/2
    Sigma_LOO[,g] = t(gamma_g)%*%var_g_LOO
    
    
    #Given kappa
    if (is.matrix(kappa)){
      var_cr = matrix(0,n_j[g]^{2},G)
      for (h in 1: G){
        kappa_gh = kappa[(n_cumu_sqr[g]+1):n_cumu_sqr[g+1],(n_cumu_sqr[h]+1):n_cumu_sqr[h+1]]
        var_cr[,h]   = kappa_gh%*%kronecker(u.hat[cluster==h],u.hat[cluster==h])
      }
      var_g_cr = rowSums(var_cr)
      Sigma_kappa[,g] = t(gamma_g)%*%var_g_cr
    } else {Variance_kappa = NA}
    
  }
  
  variance = rowSums(Sigma)
  variance_HC2 = rowSums(Sigma_HC2)
  variance_HC3 = rowSums(Sigma_HC3)
  variance_LOO = rowSums(Sigma_LOO)
  variance_kappa = rowSums(Sigma_kappa)
  
  neg_LOO = ifelse(variance_LOO<0,1,0)
  if (variance_LOO<0){
    variance_LOO = rowSums(ifelse(Sigma_LOO>0,Sigma_LOO,0))
  } 
  
  #Liang-Zeger Estimator
  Meat_LZ = resize (variance, d, d) 
  V_LZ  = (XMX_inv)%*%Meat_LZ%*%(XMX_inv)
  T_LZ  = (beta.hat - beta) / sqrt(V_LZ)
  
  #HC1 
  V_HC1 = V_LZ * G/(G-1) * (N-1)/(N-K)
  T_HC1 = (beta.hat - beta) / sqrt(V_HC1)
  
  #HC2
  Meat_HC2 = resize (variance_HC2, d, d)
  V_HC2  = (XMX_inv)%*%Meat_HC2%*%(XMX_inv)
  T_HC2  = (beta.hat - beta) / sqrt(V_HC2)
  
  #HC3
  Meat_HC3 = resize (variance_HC3, d, d)
  V_HC3  = (XMX_inv)%*%Meat_HC3%*%(XMX_inv)
  T_HC3  = (beta.hat - beta) / sqrt(V_HC3)
  
  #LOO
  Meat_LOO = resize (variance_LOO, d, d)
  V_LOO  = (XMX_inv)%*%Meat_LOO%*%(XMX_inv)
  T_LOO  = (beta.hat - beta) / sqrt(V_LOO)
  
  #CR
  if (is.matrix(kappa)){
    Meat_CR = resize (variance_kappa, d, d)
    V_CR  = (XMX_inv)%*%Meat_CR%*%(XMX_inv)
    T_CR  = (beta.hat - beta) / sqrt(V_CR)
  } else {V_CR = NA; T_CR=NA;}
  
  return(c(beta.hat,V_LZ,V_HC1,V_HC2,V_HC3,V_LOO,V_CR,neg_LOO,
           T_LZ,T_HC1,T_HC2,T_HC3,T_LOO,T_CR))
}

#Calculation of the kappa matrix using M
star<-function(cluster, M){
  G = length(unique(cluster))
  n_cumu = vector()
  n_cumu[1]=0
  n_j = table(cluster)
  n_cumu_sqr = vector() ; n_cumu_sqr[1]=0
  
  for (i in 1:G){
    n_cumu[i+1] = sum(n_j[1:i])
    n_cumu_sqr[i+1] = sum(n_j[1:i]^{2})
  }
  n_tau = n_cumu_sqr[G+1]
  M_star = matrix(0,n_tau,n_tau)
  
  for (i in 1:G){
    for (j in 1:G){
      M_star[(n_cumu_sqr[i]+1):n_cumu_sqr[i+1],(n_cumu_sqr[j]+1):n_cumu_sqr[j+1]] = kronecker(M[(n_cumu[i]+1):n_cumu[i+1],(n_cumu[j]+1):n_cumu[j+1]],M[(n_cumu[i]+1):n_cumu[i+1],(n_cumu[j]+1):n_cumu[j+1]])
    }
  }
  return(M_star)
}
