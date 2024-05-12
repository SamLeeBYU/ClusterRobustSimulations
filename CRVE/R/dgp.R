#Used for resizing matrices and arrays
require("ramify")

require("expm")

#Used for pinv Pseudo-inverse (Moore-Penrose generalized inverse) function
require("pracma")

popval = read.csv("data/popval.csv", row.names=1)

#path="D:/Research Project/ClusterManyCov/Simulation and Application 2208/Simulation/"

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

#reg.dep(15, 8, 0.5, 0.7)

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

#err.dep(15, 8, 0.5, 0.7)

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

  #I modified the DGP to include a generating process for a multi-dimensioned beta vector (with dimension d) - Sam Lee
  d = length(beta)
  Wg = W%*%matrix(1, ncol=d, nrow=1+K)


  sigma.x.fn = function(Wg) return((1+Wg^2));
  sigma.u.fn = function(Xb,Wg) return((trunc(Xb)+Wg)^2);

  v = matrix(rnorm(N*d),N,d)
  V.HE = sqrt(varkappa.v * sigma.x.fn(Wg)) * v
  X = V.HE

  #Generate Error terms
  Ud = matrix(0,N,d)

  Ud[1:n_j[1],] = err.dep(n_j[1],Ju,rho,l)
  for (i in 2:q) {
    Ud[(sum(n_j[1:(i-1)])+1):sum(n_j[1:i]),] = err.dep(n_j[i],Ju,rho,l)
  }
  U = Ud + sqrt(1-l^2)*sqrt(varkappa.u*sigma.u.fn(X,Wg))*matrix(rnorm(N*d),N,d);

  # Generate Outcome
  Y = X%*%beta + (1/d)*U%*%rep.int(1,d)

  # Return the data frame
  list  = list(cluster = cluster,unit = unit, unitx = unitx, Y = Y, X = X, V.HE = V.HE, W = W, U = U);
  return(list);
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
