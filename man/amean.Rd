\name{amean}
\alias{amean}
\title{
Average over the third index of an array
}
\description{
Convenience function to compute matrix sample means when samples are stored as a 3rd order array. Sampling index should be the third mode.
}
\usage{
amean(ar)
}
\arguments{
  \item{ar}{a 3rd order array
}
}
\value{
matrix of dimension dim(ar)[-3]
}

\author{
Evan Poworoznek
}

\seealso{
\code{\link{lmean}}
}
\examples{
ar = array(rnorm(10000), dim = c(10, 10, 100))
amean(ar)
}
