\name{as_mrio_matrix}
\alias{as_mrio_matrix}
\title{Recover a matrix from a four-mode MRIO tensor}
\usage{as_mrio_matrix(x)}
\arguments{\item{x}{Four-mode object from \code{as_mrio_tensor()}.}}
\value{Numeric matrix in region-major, sector-minor order.}
\description{Inverse of the lossless four-mode index transformation.}
\examples{
e <- example_mrio()
t <- as_mrio_tensor(e$A, e$regions, e$sectors)
as_mrio_matrix(t)
}
