\name{plotmat}
\alias{plotmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a matrix
}
\description{
Plot an image of a matrix using ggplot2
}
\usage{
plotmat(mat, color = "green", title = NULL, args = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{Matrix to plot}
  \item{color}{Color scheme: "green", "red", or "wes"}
  \item{title}{optional plot title}
  \item{args}{optional additional ggplot arguments}
}

\value{
sends image to active graphics device or outputs a ggplot object
}

\author{
Evan Poworoznek
}
\note{
Uses reshape2::melt which may be aliased with reshape::melt
}

\examples{
mat = diag(1:9 - 5)
plotmat(mat)
}

