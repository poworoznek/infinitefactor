\name{interactionDL}
\alias{interactionDL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor regression model with interactions using the Dirichlet-Laplace shrinkage prior
}
\description{
Perform a regression of y onto X and all 2 way interactions in X using the latent factor model introduced in Ferrari and Dunson (2020). This version uses the Dirichlet-Laplace shrinkage prior as in the original paper.
}
\usage{
interactionDL(y, X, nrun, burn = 0, thin = 1, 
              delta_rw = 0.0526749, a = 1/2, k = NULL, 
              output = c("covMean", "covSamples", "factSamples", 
              "sigSamples", "coefSamples","errSamples"), 
              verbose = TRUE, dump = FALSE, filename = "samps.Rds", 
              buffer = 10000, adapt = "burn", augment = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{response vector.}
\item{X}{predictor matrix (n x p).}
\item{nrun}{number of iterations.}
\item{burn}{burn-in period.}
\item{thin}{thinning interval.}
\item{delta_rw}{metropolis-hastings proposal variance.}
\item{a}{shrinkage hyperparameter.}
\item{k}{number of factors.}
\item{output}{output type, a vector including some of: c("covMean", "covSamples", "factSamples", "sigSamples", "coefSamples", "numFactors", "errSamples").}
\item{verbose}{logical. Show progress bar?}
\item{dump}{logical. Save samples to a file during sampling?}
\item{filename}{if dump: filename to address list of posterior samples}
\item{buffer}{if dump: how often to save samples}
\item{adapt}{logical or "burn". Adapt proposal variance in metropolis hastings step? if "burn", will adapt during burn in and not after.}
\item{augment}{additional sampling steps as an expression}
}
\value{
some of:
\item{covMean}{X covariance posterior mean}
\item{omegaSamps}{X covariance posterior samples}
\item{lambdaSamps}{Posterior factor loadings samples (rotationally ambiguous)}
\item{etaSamps}{Posterior factor samples (rotationally ambiguous)}
\item{sigmaSamps}{Posterior marginal variance samples (see notation in Bhattacharya and Dunson (2011))}
\item{phiSamps}{Posterior main effect coefficient samples in factor form (rotationally ambiguous)}
\item{PsiSamps}{Posterior interaction effect coefficient samples in factor form (rotationally ambiguous)}
\item{interceptSamps}{Posterior induced intercept samples}
\item{mainEffectSamps}{Posterior induced main effect coefficient samples}
\item{interactionSamps}{Posterior induced interaction coefficient samples}
\item{ssySamps}{Posterior irreducible error samples}
}
\references{
Ferrari, Federico, and David B. Dunson. "Bayesian Factor Analysis for Inference on Interactions." arXiv preprint arXiv:1904.11603 (2019).
}
\author{
Evan Poworoznek

Federico Ferrari}

\seealso{
\code{\link{interactionMGSP}}
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

beta_true = numeric(p); beta_true[c(1,3,6,8,10,11)] =c(1,1,0.5,-1,-2,-0.5)
Omega_true = matrix(0,p,p)
Omega_true[1,2] = 1; Omega_true[5,2] = -1; Omega_true[10,8] = 1; 
Omega_true[11,5] = -2; Omega_true[1,1] = 0.5; 
Omega_true[2,3] = 0.5; 
Omega_true = Omega_true + t(Omega_true)
y = X\%*\%beta_true + diag(X\%*\%Omega_true\%*\%t(X)) +  rnorm(n,0.5)

intdl = interactionDL(y, X, 1000, 500, k = 5)
}

