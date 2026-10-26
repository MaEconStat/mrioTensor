\name{mrio_hosvd}
\alias{mrio_hosvd}
\title{Higher-order singular value decomposition (HOSVD)}
\usage{
mrio_hosvd(x, ranks, max_entries = 1e+06)
}
\arguments{
\item{x}{Finite numeric array of order at least two.}
\item{ranks}{Integer target rank per mode.}
\item{max_entries}{Safety ceiling for number of elements in \code{x}.}
}
\value{List containing \code{factors}, \code{core}, \code{reconstruction}, and \code{relative_error}.}
\description{Computes truncated mode-wise SVDs and projects onto the Tucker core array.}
\examples{
t <- array(rnorm(16), dim = c(2, 2, 2, 2))
fit <- mrio_hosvd(t, ranks = c(2, 2, 2, 2))
fit$relative_error
}
