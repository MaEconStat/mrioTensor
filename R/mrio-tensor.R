#' Convert a region-major MRIO matrix to an explicit four-mode tensor
#'
#' Reshapes a 2D composite region-sector matrix into a 4-mode tensor with axes
#' \code{(origin_region, origin_sector, destination_region, destination_sector)}.
#' Indexing follows the standard region-major, sector-minor layout:
#' \code{i = (region - 1) * S + sector}, where \code{S} is the number of sectors
#' per region.
#'
#' @param x Finite square numeric matrix, ordered region first, sector second.
#' @param regions Unique region names in matrix order.
#' @param sectors Unique sector names in matrix order, shared across regions.
#' @return An \code{mrio_tensor} array with dimensions (origin_region, origin_sector, destination_region, destination_sector).
#' @export
as_mrio_tensor <- function(x, regions, sectors) {
  .assert_matrix(x, "x", square = TRUE)
  .assert_axes(regions, sectors, nrow(x))
  R <- length(regions)
  S <- length(sectors)
  # R stores arrays column-major, so the sector indices run fastest.
  ans <- aperm(array(as.vector(x), dim = c(S, R, S, R)), c(2L, 1L, 4L, 3L))
  dimnames(ans) <- list(origin_region = regions, origin_sector = sectors,
                        destination_region = regions,
                        destination_sector = sectors)
  class(ans) <- c("mrio_tensor", "array")
  ans
}

#' Recover the original two-dimensional MRIO matrix
#'
#' Re-flattens a 4-mode \code{mrio_tensor} back into its canonical 2D matrix
#' representation in region-major, sector-minor order.
#'
#' @param x Four-mode object created by \code{as_mrio_tensor()}.
#' @return Numeric matrix in region-major, sector-minor order.
#' @export
as_mrio_matrix <- function(x) {
  .assert_mrio_tensor(x)
  R <- dim(x)[1L]
  S <- dim(x)[2L]
  labels <- dimnames(x)
  # Build node labels in the documented region-major order.
  src <- paste(rep(labels[[1L]], each = S), rep(labels[[2L]], times = R),
               sep = "/")
  dst <- paste(rep(labels[[3L]], each = S), rep(labels[[4L]], times = R),
               sep = "/")
  matrix(as.vector(aperm(unclass(x), c(2L, 1L, 4L, 3L))),
         nrow = R * S, ncol = R * S,
         dimnames = list(src, dst))
}

#' Aggregate an MRIO tensor across spatial or sectoral dimensions
#'
#' Evaluates mode-wise marginal sums of an MRIO tensor over regional or
#' industrial dimensions, yielding either an aggregated inter-regional flow
#' matrix or an inter-sectoral technical flow matrix.
#'
#' @param x Four-mode MRIO tensor.
#' @param by Either \code{"region"} (marginalize out sectors) or \code{"sector"} (marginalize out regions).
#' @return A square numeric matrix of aggregated marginal flows.
#' @export
aggregate_mrio <- function(x, by = c("region", "sector")) {
  .assert_mrio_tensor(x)
  by <- match.arg(by)
  keep <- if (by == "region") c(1L, 3L) else c(2L, 4L)
  ans <- apply(x, keep, sum)
  d <- dim(x)[keep]
  matrix(as.vector(ans), nrow = d[1L], ncol = d[2L],
         dimnames = dimnames(x)[keep])
}

#' @export
print.mrio_tensor <- function(x, ...) {
  d <- dim(x)
  cat("<mrio_tensor> ", d[1L], " regions x ", d[2L],
      " sectors x ", d[3L], " regions x ", d[4L],
      " sectors\n", sep = "")
  invisible(x)
}
