\name{path_order_contributions}
\alias{path_order_contributions}
\title{Contributions by Leontief production tier order}
\usage{
path_order_contributions(A, y, intensity, max_order = 3L)
}
\arguments{
\item{A}{Nonnegative technical-coefficient matrix.}
\item{y}{Nonnegative final demand vector.}
\item{intensity}{Nonnegative direct environmental intensity vector.}
\item{max_order}{Maximum Leontief power order to evaluate.}
}
\value{Named numeric vector of contributions by power order.}
\description{Calculates tier contributions \code{intensity' A^k y} via the Neumann series expansion.}
\examples{
e <- example_mrio()
path_order_contributions(e$A, e$y, e$intensity, max_order = 3)
}
