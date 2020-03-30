# infinitefactor
Bayesian infinite factor modelling

This package was developed for modular construction of Bayesian factor samplers using a variety of data models and priors. Bayesian factor models are key tools for performing reproducible and robust dimension reduction and modeling complex relationships among sets of regressors. Accounting for structure with factor models allows both efficiencies in estimation, and inference on that structure. Popular shrinkage priors (MGSP, Dirichlet-Laplace, others) have been formulated for the factor loadings matrix that do not require prespecified covariate grouping or order-dependence. This package aims to implement those priors for use with multiple data models in an easy, straightforward, and fast environment for application, as well as encourage further model and prior development through use of modular and exchangeable parameter updates.

## Augmentation

The `augment` sampler input argument can be used to modify sampler behavior for non-standard applications. This argument takes an R expression to be evaluated every sampler iteration after the standard parameter updates and before samples are stored. A straightforward application of this argument is to change a sampling parameter from a single value to a random variable. This can be done to place further hyperpriors on shrinkage parameters, or to perform data augmentation. Refer to the sample source code for parameter names.

### Missingness

A common form of data augmentation is imputation due to missingness. In the included factor models we have simple representations of the posterior predictive distributions of the factorized data matrix X and the interaction response y. Here we provide example code to sample missing entries of X under a Missing at Random (MAR) assumption. 

```R
# for data matrix "data" with NA missing values

completeX = function(X, Xmiss, lambda, eta, ps){
  noise = t(replicate(nrow(X), rnorm(ncol(X), 0, sqrt(1/ps))))
  X[is.na(Xmiss)] = (tcrossprod(eta, lambda) +noise)[is.na(Xmiss)]
  return(X)}

missing = expression({X = completeX(X, data, lambda, eta, ps)})

X = data
X[is.na(X)] = rnorm(sum(is.na(X)))
                    
sample = linearMGSP(X, 10000, 5000, augment = missing)
```

Similarly we can define an augmentation for missingness due to a (lower) limit of detection. Here we allow a different LoD for every element of X, corresponding to test and batch LoD variablility, and those limits are encoded in the matrix `LOD` where dim(`LOD`)=dim(`X`). In this case the appropriate posterior predictive is the same as above, but truncated at the limit of detection. We use the `truncnorm` package for truncated normal sampling.

```R
# for data matrix "data" with NA missing values

lodX = function(X, Xmiss, lambda, eta, ps){
  vars = matrix(1/ps, ncol(X), nrow(X), byrow=T)[is.na(Xmiss)]
  nmiss = sum(is.na(Xmiss))
  means = tcrossprod(eta, lambda)[is.na(Xmiss)]
  X[is.na(Xmiss)] = rtruncnorm(nmiss, b = LOD[is.na(X)], means,sqrt(vars))
  return(X)}

missing = expression({X = lodX(X, data, lambda, eta, ps)})

X = data

X[is.na(X)] = LOD[is.na(X)]

sample = linearMGSP(X, 10000, 5000, augment = lodX)
```

This code can be easily altered for upper limits of detection (see `rtruncnorm` documentation).



