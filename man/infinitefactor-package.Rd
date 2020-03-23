\name{infinitefactor-package}
\alias{infinitefactor-package}
\alias{infinitefactor}
\docType{package}
\title{
\packageTitle{infinitefactor}
}
\description{
\packageDescription{infinitefactor}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{infinitefactor}
\packageIndices{infinitefactor}
Perform sampling with the linearMGSP() and linearDL() functions for linear factor models, or interactionMGSP() and interactionDL() functions for factor regression models including 2-way interactions. See jointRot() or msf() for postprocessing.
}
\author{
\packageAuthor{infinitefactor}

Maintainer: \packageMaintainer{infinitefactor}
}
\references{
Bhattacharya, Anirban, and David B. Dunson. "Sparse Bayesian infinite factor models." Biometrika (2011): 291-306.

Bhattacharya, Anirban, et al. "Dirichlet-Laplace priors for optimal shrinkage." Journal of the American Statistical Association 110.512 (2015): 1479-1490.

Ferrari, Federico, and David B. Dunson. "Bayesian Factor Analysis for Inference on Interactions." arXiv preprint arXiv:1904.11603 (2019).
}
\keyword{ package }
\examples{
k0 = 5
p = 20
n = 100

lambda = matrix(rnorm(p*k0, 0, 0.01), ncol = k0)
lambda[sample.int(p, 40, replace = TRUE) +
         p*(sample.int(k0, 40, replace = TRUE)-1)] = rnorm(40, 0, 1)
lambda[1:7, 1] = rnorm(7, 2, 0.5)
lambda[8:14, 2] = rnorm(7, -2, 0.5)
lambda[15:20, 3] = rnorm(6, 2, 0.5)
lambda[,4] = rnorm(p, 0, 0.5)
lambda[,5] = rnorm(p, 0, 0.5)
plotmat(varimax(lambda)[[1]])

X = matrix(rnorm(n*k0),n,k0)\%*\%t(lambda) + matrix(rnorm(n*p), n, p)

out = linearMGSP(X = X, nrun = 1000, burn = 500, adapt = FALSE)

aligned = jointRot(out$lambdaSamps, out$etaSamps)

plotmat(lmean(aligned$lambda))
}
