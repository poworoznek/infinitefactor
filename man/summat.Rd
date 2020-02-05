\name{summat}
\alias{summat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarise a matrix from posterior samples
}
\description{
Provide a summary matrix from a list of matrix-valued parameter samples, returning the mean value for each element with 0 not included in its quantile-based posterior credible interval, and 0 for each element for which 0 is included in its posterior CI.
}
\usage{
summat(list, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list}{list of matrix valued parameter samples of the same dimensions}
  \item{alpha}{type I error probability}
}

\value{
a matrix
}

\author{
Evan Poworoznek
}

\seealso{
\code{\link{lmean}}
}
\examples{
list = replicate(1000, matrix(rnorm(100), ncol = 10) + 
                 10*diag(10), simplify = FALSE)
lmean(list)
summat(list)
plotmat(summat(list))
}
