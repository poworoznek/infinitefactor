\name{linearMGSP}
\alias{linearMGSP}

\title{
Sample Bayesian linear infinite factor models with the Multiplicative Gamma Shrinkage Prior
}
\description{
Perform Bayesian factor analysis by sampling the posterior distribution of parameters in a factor model with the Multiplicative Gamma Shrinkage Prior of Bhattacharya and Dunson
}
\usage{
linearMGSP(X, nrun, burn, thin = 1, prop = 1, 
epsilon = 1e-3, kinit = NULL, adapt = TRUE, 
output = c("covMean", "covSamples", "factSamples",
"sigSamples", "numFactors"), verbose = TRUE, 
dump = FALSE, filename = "samps.Rds", buffer = 10000,
augment = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{Data matrix (n x p)}
\item{nrun}{number of iterations}
\item{burn}{burn-in period}
\item{thin}{thinning interval}
\item{prop}{proportion of elements in each column less than epsilon in magnitude cutoff}
\item{epsilon}{tolerance}
\item{kinit}{initial value for the number of factors}
\item{adapt}{logical. Whether or not to adapt number of factors across sampling}
\item{output}{output type, a vector including some of: c("covMean", "covSamples", "factSamples", "sigSamples", "numFactors")}
\item{verbose}{logical. Show progress bar?}
\item{dump}{logical. Save output object during sampling?}
\item{filename}{if dump, filename for output}
\item{buffer}{if dump, frequency of saving}
\item{augment}{additional sampling steps as an expression}
}

\value{
some of:
\item{covMean}{X covariance posterior mean}
\item{omegaSamps}{X covariance posterior samples}
\item{lambdaSamps}{Posterior factor loadings samples (rotationally ambiguous)}
\item{etaSamps}{Posterior factor samples (rotationally ambiguous)}
\item{sigmaSamps}{Posterior marginal variance samples (see notation in Bhattacharya and Dunson (2011))}
\item{numFacts}{Number of factors for each iteration}
}
\references{
Bhattacharya, Anirban, and David B. Dunson. "Sparse Bayesian infinite factor models." Biometrika (2011): 291-306.
}
\author{
Evan Poworoznek
}

\seealso{
\code{\link{linearDL}}
}
\examples{
k0 = 5
p = 20
n = 50

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

out = linearMGSP(X = X, nrun = 1000, burn = 500)
}
