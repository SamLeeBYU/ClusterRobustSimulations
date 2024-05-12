---
Title: OLS Regression with Cluster-Robust Variance Estimates
Author: Sam Lee
---

# Introduction

In this project I create an R function that runs a least squares regression after the following regression model:

$y_i = x_i'\beta + w_i'\gamma + u_i$, where $X$ is a set of covariates to obtain unbiased and consitent estimators for, and $W$ is a set of non-vanishing covariates.

This project follows directly after the work of Matias D. Cattaneo, Aibo Gong, Michael Jansson, and Whitney K. Newey. The function I create is simply a compilation and modificaiton of the R code they already written following their paper, "Cluster Robust Inference in Linear Models with Many Covariates" (2022).

In this particular regression model, we are interested in obtaining consistent estimates for $\beta$. We consider a situation were observations are correlated between clusters. Utilizing methods discussed by the authors mentioned previously, we account for this and create a new CRVE (cluster-robust variance estimator) to appropriately adjust for the standard errors of the relevant covariates.

# CRVE Package

The function that I wrote to estimate $\beta$ and the respective cluster-robust variance matrix has been compiled into an R package on this repository. The source code can be found [here](https://github.com/SamLeeBYU/ClusterRobustSimulations/tree/main/CRVE).

## Installation

```
library(devtools)
devtools::install_github("SamLeeBYU/ClusterRobustSimulations/CRVE")

library(CRVE)
```

Utilizing the data generating process designed by Cattaneo, Gong, Jansson, and Whitney, we can perform a direct application of the `crve.fit` function:

```
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
```

In the source code, I modified the DGP to allow for a beta-vector with dimension $d$, where $d \geq 0$. The original authors model the simulation assuming $d=1$. The hyperparameters specified above are parameters to generate correlation within each cluster (heterogeneous cluster sizes in this case). The `crve.fit` function will work on any specificied $Y$, $X$, and $W$ matrices regardless if there is heterogeneous correlation between the clusters.

## Further explanation of the DGP Hyperparameters

The following notes are directly from the authors' original code documentation:

- q: the number of clusters
- qx: the number of cluster for regressors
- n and n_prime: the reference sample size within cluster
- nx: the reference sample size for regressors
- Jx and Ju: Random Vector dimensions for regressors and error terms
- rho: correlation coefficien
- l: weighted parameters in D.G.P
- K: the number of covariates
- m: sets the dgp(1,2 for independent covariates with homogeneous and heterogeneous cluster size; 3,4 for dependent covariates with homogeneous and heterogeneous cluster size)

## More on the `crve.fit` Function

This function estimates the covariates of interest (for the covariate matrix X) along with estimates for the (nuisance) covariate matrix W. This function also returns the cluster-robust variance matrix (V.CR) and the appropriate t-statistic (T.CR) for a given coefficient vector (assumed 0).

### Arguments

- $Y$: An $N\times 1$ vector/matrix; The response variable.
- $X$: An $N\times d$ matrix representing the covariates of interest. The dimension $d$ is the number of covariates (dimension of beta).
- $W$: An $N\times K$ matrix representing the nuisance covariates (e.g. fixed effects). The regression function estimates these but doesn't return them explicitly.
- `cluster`: A one-dimension vector of length $N$ representing the indices of the clusters for each observation. Each unique number in this vector represents a unique cluster.
- `beta`: A $d \times 1$ vector to evaluate the test-statistic of the relevant coefficients. Useful in Monte Carlo simulations.

### Return Values

`beta.hat`: A $d\times 1$ vector of the estimates of the relevant covariates pertaining to the $X$ covariate matrix.
`V.CR`: The dxd cluster-robust variance matrix.
`T.CR`: A $d\times 1$ vector of the t-statistic evaluated relative to the 'true' beta vector passed in (assumed 0).

