\name{embodied_transfer}
\alias{embodied_transfer}
\title{Compute Leontief-based final-demand attribution}
\usage{
embodied_transfer(A, y, intensity, regions = NULL, sectors = NULL,
  max_nodes = 1000L)
}
\arguments{
\item{A}{Nonnegative square technical-coefficient matrix.}
\item{y}{Nonnegative final demand vector.}
\item{intensity}{Nonnegative direct environmental intensity vector.}
\item{regions}{Optional region labels, paired with \code{sectors}.}
\item{sectors}{Optional sector labels, paired with \code{regions}.}
\item{max_nodes}{Safety threshold for dense inverse dimension.}
}
\value{List with pairwise attribution \code{matrix}, optional \code{tensor}, \code{total}, and \code{output}.}
\description{Evaluates the embodied impact attribution \code{F[i, j] = intensity[i] * L[i, j] * y[j]}
with \code{L = (I - A)^(-1)}, allocating upstream production impacts to final consumption demand.}
\examples{
e <- example_mrio()
embodied_transfer(e$A, e$y, e$intensity)$total
}
