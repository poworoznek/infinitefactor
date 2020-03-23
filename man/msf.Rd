\name{msf}
\alias{msf}
\alias{msfOUT}
\alias{aplr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resolve label and sign switching in random matrix samples
}
\description{
The msf() function performs column-based matching of a matrix to a pivot to resolve rotational ambiguity remaining after the application of an orthogonalisation procedure on a list of Bayesian matrix samples. The msfOUT() and aplr() functions perform this same matching but instead of returning aligned samples as does msf(), msfOUT outputs the list of permutations and sign switches needed for alignment and aplr outputs a list of matrices permuted and re-signed by msfOUT() output. msfOUT() and aplr() are used in jointRot(). These functions are written in C++ and may crash the R session if passed inappropriate input.
}
\usage{
msf(lambda, pivot)

msfOUT(lambda, pivot)

aplr(matr, perm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{matrix to be aligned, named for a factor loadings matrix as in the Bhattacharya and Dunson 2011 notation}
  \item{pivot}{matrix to align with which to align lambda}
  \item{matr}{a matrix to apply permutations to}
  \item{perm}{a (possibly signed) permutation order for the matr matrix}
}
\details{
see the examples for suggested usage of msf and jointRot() for suggested usage of msfOUT() and aplr().
}

\author{
Evan Poworoznek
}


\seealso{
\code{\link{jointRot}}
}
\examples{
lambda = diag(10)[,sample(10)] + 0.001
pivot = diag(10)
msf(lambda, pivot)

# fast implementation for a list of samples
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

vari = lapply(out$lambdaSamps, varimax)
loads = lapply(vari, `[[`, 1)

norms = sapply(loads, norm, "2")
pivot = loads[order(norms)][[250]]

aligned = lapply(loads, msf, pivot)
plotmat(summat(aligned))
}
