#' Solve the Leontief production balance system
#'
#' For a technical coefficient matrix \code{A} and final demand vector \code{y},
#' calculates gross output \code{x} satisfying \code{(I - A)x = y}. Solves the
#' linear system directly via LU decomposition without forming the explicit inverse.
#'
#' @param A Nonnegative numeric \code{N x N} technical coefficient matrix.
#' @param y Nonnegative final demand vector of length \code{N}.
#' @return Numeric gross output vector \code{x}.
#' @export
leontief_output <- function(A, y) {
  .assert_technology(A)
  n <- nrow(A)
  .assert_nonnegative(y, n, "y")
  x <- .safe_solve(diag(n) - A, y)
  as.vector(x)
}

#' Compute pairwise Leontief embodied environmental transfer matrix
#'
#' Evaluates the embodied impact attribution \code{F[i, j] = intensity[i] * L[i, j] * y[j]},
#' where \code{L = (I - A)^(-1)}. Allocates upstream direct environmental footprints
#' from producing sectors to final consumption demand categories.
#'
#' @param A Nonnegative square technical-coefficient matrix.
#' @param y Nonnegative final demand vector (\code{N} sectors).
#' @param intensity Nonnegative direct environmental intensity per unit output.
#' @param regions Optional region labels; must be supplied together with \code{sectors}.
#' @param sectors Optional sector labels; must be supplied together with \code{regions}.
#' @param max_nodes Safety threshold for dense inverse calculation (default 1000).
#' @return List with \code{matrix} (pairwise attribution), optional \code{tensor} (if labels supplied), \code{total}, and \code{output}.
#' @export
embodied_transfer <- function(A, y, intensity, regions = NULL, sectors = NULL,
                              max_nodes = 1000L) {
  .assert_technology(A)
  n <- nrow(A)
  .assert_nonnegative(y, n, "y")
  .assert_nonnegative(intensity, n, "intensity")
  if (!is.numeric(max_nodes) || length(max_nodes) != 1L ||
      !is.finite(max_nodes) || max_nodes < 1 || n > max_nodes) {
    stop("Dense inverse safety limit reached or max_nodes invalid. ",
         "Use leontief_output() for totals; increase max_nodes explicitly ",
         "only with sufficient memory.", call. = FALSE)
  }
  if (xor(is.null(regions), is.null(sectors))) {
    stop("Supply regions and sectors together, or omit both.", call. = FALSE)
  }
  if (!is.null(regions)) .assert_axes(regions, sectors, n)
  L <- .safe_solve(diag(n) - A)
  flows <- sweep(sweep(L, 1L, intensity, `*`), 2L, y, `*`)
  output <- as.vector(L %*% y)
  tensor <- if (!is.null(regions)) {
    as_mrio_tensor(flows, regions, sectors)
  } else NULL
  list(matrix = flows, tensor = tensor, total = sum(flows), output = output)
}

#' Contributions by Leontief production tier order
#'
#' Evaluates the Neumann series expansion \code{c[k] = intensity' A^k y} for
#' production tier \code{k = 0, ..., max_order}, capturing direct final demand
#' requirements (\code{k=0}), first-tier intermediate supply (\code{k=1}), and
#' higher-order upstream production cascades.
#'
#' @param A Nonnegative technical-coefficient matrix.
#' @param y Nonnegative final demand vector.
#' @param intensity Nonnegative direct environmental intensity vector.
#' @param max_order Nonnegative integer indicating the maximum power of \code{A} to evaluate.
#' @return Numeric vector of length \code{max_order + 1} with names indicating tier orders.
#' @export
path_order_contributions <- function(A, y, intensity, max_order = 3L) {
  .assert_technology(A)
  n <- nrow(A)
  .assert_nonnegative(y, n, "y")
  .assert_nonnegative(intensity, n, "intensity")
  if (!is.numeric(max_order) || length(max_order) != 1L ||
      is.na(max_order) || !is.finite(max_order) ||
      max_order != as.integer(max_order) || max_order < 0L) {
    stop("max_order must be a non-negative integer.", call. = FALSE)
  }
  max_order <- as.integer(max_order)
  contributions <- numeric(max_order + 1L)
  v <- y
  contributions[1L] <- sum(intensity * v)
  if (max_order > 0L) {
    for (k in seq_len(max_order)) {
      v <- as.vector(A %*% v)
      contributions[k + 1L] <- sum(intensity * v)
    }
  }
  names(contributions) <- paste0("order_", 0L:max_order)
  contributions
}
