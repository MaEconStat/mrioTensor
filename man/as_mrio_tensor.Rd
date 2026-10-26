\name{as_mrio_tensor}
\alias{as_mrio_tensor}
\title{Convert a region-major MRIO matrix to an explicit four-mode tensor}
\usage{as_mrio_tensor(x, regions, sectors)}
\arguments{
\item{x}{Finite square numeric matrix in region-major, sector-minor order.}
\item{regions}{Unique region labels in matrix order.}
\item{sectors}{Unique sector labels within each region.}
}
\value{An \code{mrio_tensor} array with region, sector, region, sector axes.}
\description{Re-index a composite-node matrix without changing values or accounting identities.}
\examples{
e <- example_mrio()
t <- as_mrio_tensor(e$A, e$regions, e$sectors)
dim(t)
}
