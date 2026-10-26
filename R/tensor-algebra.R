#' Mode-n matrix product of a tensor
#'
#' Evaluates the tensor-matrix contraction along axis \code{mode}. For tensor \code{x}
#' of dimensions \code{(d_1, ..., d_K)} and matrix \code{M} of dimensions \code{(J, d_mode)},
#' contracts the \code{mode}-th axis such that the resulting array has dimensions
#' \code{(d_1, ..., J, ..., d_K)}.
#'
#' @param x Finite numeric array of order at least two.
#' @param M Finite numeric matrix whose column count matches \code{dim(x)[mode]}.
#' @param mode Integer axis index, from 1 to \code{length(dim(x))}.
#' @return Numeric array after mode-n contraction.
#' @export
mode_product <- function(x, M, mode) {
  if (!is.array(x) || !is.numeric(x) || length(dim(x)) < 2L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric array of order >= 2.", call. = FALSE)
  }
  .assert_matrix(M, "M")
  d <- dim(x)
  if (!is.numeric(mode) || length(mode) != 1L || is.na(mode) ||
      !is.finite(mode) || mode != as.integer(mode) ||
      mode < 1L || mode > length(d)) {
    stop("mode must be an integer index of a tensor axis.", call. = FALSE)
  }
  mode <- as.integer(mode)
  if (ncol(M) != d[mode]) {
    stop("ncol(M) must match the length of the selected tensor mode.",
         call. = FALSE)
  }
  perm <- c(mode, setdiff(seq_along(d), mode))
  unfolded <- matrix(aperm(x, perm), nrow = d[mode])
  result_permuted <- array(M %*% unfolded,
                           dim = c(nrow(M), d[-mode]))
  result <- aperm(result_permuted, order(perm))
  labels <- dimnames(x)
  if (!is.null(labels)) {
    labels[mode] <- list(rownames(M))
    dimnames(result) <- labels
  }
  result
}

#' Truncated Higher-Order Singular Value Decomposition (HOSVD)
#'
#' Computes mode-wise principal left singular vectors of tensor unfoldings
#' and projects the tensor onto the corresponding Tucker core array.
#'
#' @param x Finite numeric array of order at least two.
#' @param ranks Integer vector specifying the target truncation rank per mode.
#' @param max_entries Safety ceiling for the number of elements in \code{x} (default 1e6).
#' @return List containing \code{factors} (list of mode singular matrices), \code{core} (Tucker core tensor), \code{reconstruction}, and \code{relative_error}.
#' @export
mrio_hosvd <- function(x, ranks, max_entries = 1e6) {
  if (!is.array(x) || !is.numeric(x) || length(dim(x)) < 2L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric array of order >= 2.", call. = FALSE)
  }
  d <- dim(x)
  if (!is.numeric(max_entries) || length(max_entries) != 1L ||
      is.na(max_entries) || !is.finite(max_entries) || max_entries < 1) {
    stop("max_entries must be a finite positive number.", call. = FALSE)
  }
  if (length(x) > max_entries) {
    stop("HOSVD preview limit exceeded; use smaller data or a scalable solver.",
         call. = FALSE)
  }
  max_rank <- vapply(seq_along(d), function(m) {
    min(d[m], prod(d[-m]))
  }, numeric(1L))
  if (!is.numeric(ranks) || length(ranks) != length(d) ||
      any(is.na(ranks)) || any(!is.finite(ranks)) ||
      any(ranks != as.integer(ranks)) ||
      any(ranks < 1L) || any(ranks > max_rank)) {
    stop("ranks must be positive integers within valid mode dimensions.",
         call. = FALSE)
  }
  ranks <- as.integer(ranks)
  factors <- vector("list", length(d))
  for (m in seq_along(d)) {
    perm <- c(m, setdiff(seq_along(d), m))
    unfolded <- matrix(aperm(x, perm), nrow = d[m])
    s <- svd(unfolded, nu = ranks[m], nv = 0L)
    factors[[m]] <- s$u[, seq_len(ranks[m]), drop = FALSE]
  }
  core <- x
  for (m in seq_along(d)) {
    core <- mode_product(core, t(factors[[m]]), mode = m)
  }
  recon <- core
  for (m in seq_along(d)) {
    recon <- mode_product(recon, factors[[m]], mode = m)
  }
  norm_diff <- sqrt(sum((x - recon)^2))
  norm_orig <- sqrt(sum(x^2))
  rel_err <- if (norm_orig == 0) 0 else norm_diff / norm_orig
  names(factors) <- paste0("mode_", seq_along(d))
  list(factors = factors, core = core, reconstruction = recon,
       relative_error = rel_err)
}
