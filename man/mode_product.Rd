\name{mode_product}
\alias{mode_product}
\title{Multiply a tensor along a selected mode}
\usage{
mode_product(x, M, mode)
}
\arguments{
\item{x}{Finite numeric array of order at least two.}
\item{M}{Finite numeric matrix; its column count matches the selected mode.}
\item{mode}{Integer axis index from 1 to \code{length(dim(x))}.}
}
\value{Numeric array after mode-n contraction.}
\description{Contracts a tensor axis with a matrix along the designated mode.}
\examples{
t <- array(1:24, dim = c(2, 3, 4))
m <- matrix(1:6, nrow = 3, ncol = 2)
mode_product(t, m, mode = 1)
}
