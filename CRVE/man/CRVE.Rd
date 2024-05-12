\name{crve.fit}
\alias{CRVE}
\title{
Ordinary Least-Squares Regression with Cluster-Robust Variance Estimation
}
\description{
This function estimates the covariates of interest (for the covariate matrix X) along with estimates for the (nuisance) covariate matrix W. This function also returns the cluster-robust variance matrix (V.CR) and the appropriate t-statistic (T.CR) for a given coefficient vector (assumed 0).
}
\usage{
crve.fit(Y, X, W, cluster, beta=0)
}
\arguments{
  \item{Y}{An Nx1 vector/matrix; The response variable.}
  \item{X}{An Nxd matrix representing the covariates of interest. The dimension d is the number of covariates (dimension of beta).}
  \item{W}{An NxK matrix representing the nuisance covariates (e.g. fixed effects). The regression function estimates these but doesn't return them explicitly.}
  \item{cluster}{A one-dimension vector of length N representing the indices of the clusters for each observation. Each unique number in this vector represents a unique cluster.}
  \item{beta}{A dx1 vector to evaluate the test-statistic of the relevant coefficients. Useful in Monte Carlo simulations.}
}
\value{
Returns a list with three elements:
\item{beta.hat}{A dx1 vector of the estimates of the relevant covariates pertaining to the X covariate matrix.}
\item{V.CR}{The dxd cluster-robust variance matrix.}
\item{T.CR}{A dx1 vector of the t-statistic evaluated relative to the 'true' beta vector passed in (assumed 0).}
}
\references{
Original R code from simulation and cluster-robust variance theory is accredited to the work of Matias D. Cattaneo, Aibo Gong, Michael Jansson, and Whitney K. Newey following their paper, "Cluster Robust Inference in Linear Models with Many Covariates" (2022).
}
\author{
Sam Lee
}

\seealso{
\link{https://github.com/SamLeeBYU/ClusterRobustSimulations/test.R}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
