\name{leontief_output}
\alias{leontief_output}
\title{Solve a Leontief input-output balance system}
\usage{leontief_output(A, y)}
\arguments{
\item{A}{Nonnegative numeric square technical coefficients.}
\item{y}{Nonnegative numeric final demand vector.}
}
\value{Numeric gross-output vector.}
\description{Solves \code{(I - A) x = y} without explicitly computing the inverse.}
\examples{
e <- example_mrio()
leontief_output(e$A, e$y)
}
