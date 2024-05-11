---
Title: OLS Regression with Cluster-Robust Variance Estimates
Author: Sam Lee
---

# Introduction

In this project I create an R function that runs a least squares regression after the following regression model:

$y_i = x_i'\beta + w_i'\gamma + u_i$, where $X$ is a set of covariates to obtain unbiased and consitent estimators for, and $W$ is a set of non-vanishing covariates.

This project follows directly after the work of Matias D. Cattaneo, Aibo Gong, Michael Jansson, and Whitney K. Newey. The function I create is simply a compilation and modificaiton of the R code they have already written following their paper, "Cluster Robust Inference in Linear Models with Many Covariates" (2022).

In this particular regression model, we are interested in obtaining consistent estimates for $\beta$. We consider a situation were observations are correlate between clusters. Utilizing methods discussed by the authors mentioned previously, we account for this and create a new CRVE (cluster-robust variance estimator) to appropriately adjust for the standard errors of the relevant covariates.