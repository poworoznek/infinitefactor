\name{jointRot}
\alias{jointRot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resolve rotational ambiguity in samples of factor loadings and factors jointly
}
\description{
Performs the varimax rotation on the factor loadings samples and column-based matching to resolve resultant sign and label switching. Rotates the factors along with the loadings to induce identifiability jointly. Note this method will only work on lists of factors and factor loadings that share the same constant number of factors (k) across all samples, and will likely crash the session if this is not the case.
}
\usage{
jointRot(lambda, eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{list of factor loadings samples}
  \item{eta}{list of factor samples}
}

\value{
\item{lambda}{rotationally aligned factor loadings samples}
\item{eta}{rotationally aligned factor samples}
%% ...
}
\references{
coming soon...
}
\author{
Evan Poworoznek
}

\seealso{
\code{\link{msf}}
}
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

