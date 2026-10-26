\name{aggregate_mrio}
\alias{aggregate_mrio}
\title{Aggregate an MRIO tensor across spatial or sectoral dimensions}
\usage{
aggregate_mrio(x, by = c("region", "sector"))
}
\arguments{
\item{x}{Four-mode \code{mrio_tensor}.}
\item{by}{Either \code{"region"} or \code{"sector"}.}
}
\value{A square numeric matrix of aggregated marginal flows.}
\description{Evaluates marginal sums of an MRIO tensor over spatial or industrial modes.}
\examples{
e <- example_mrio()
t <- as_mrio_tensor(e$A, e$regions, e$sectors)
aggregate_mrio(t, by = "region")
aggregate_mrio(t, by = "sector")
}
