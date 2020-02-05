\name{lmean}
\alias{lmean}
\title{
Average elements of a list
}
\description{
Convenience function to compute sample means when samples are stored as a list. List elements should be compatible with addition and scalar division (e.g. must share the same dimensions).
}
\usage{
lmean(list)
}
\arguments{
  \item{list}{a list of parameter samples}
}
\value{
same type as a single element of the input list
}

\author{
Evan Poworoznek
}

\seealso{
\code{\link{amean}}
}
\examples{
l = replicate(100, rnorm(10), simplify = FALSE)
lmean(l)
}
