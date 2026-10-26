#' Generate a small fully synthetic MRIO demonstration
#'
#' Returns a four-node, two-region, two-sector example. All values are made
#' up for teaching and testing and are not derived from any real MRIO database.
#'
#' @return List containing `A`, `y`, `intensity`, `regions`, `sectors`.
#' @export
example_mrio <- function() {
  A <- matrix(c(
    0.16, 0.06, 0.02, 0.01,
    0.04, 0.18, 0.01, 0.03,
    0.03, 0.02, 0.12, 0.05,
    0.01, 0.04, 0.06, 0.15
  ), nrow = 4L, byrow = TRUE)
  regions <- c("North", "South")
  sectors <- c("Agriculture", "Manufacturing")
  nodes <- paste(rep(regions, each = 2L), rep(sectors, times = 2L),
                 sep = "/")
  dimnames(A) <- list(nodes, nodes)
  list(A = A, y = c(40, 70, 30, 55),
       intensity = c(0.30, 0.65, 0.25, 0.50),
       regions = regions, sectors = sectors)
}
