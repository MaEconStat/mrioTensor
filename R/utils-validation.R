.assert_matrix <- function(x, name, square = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) ||
      anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.", call. = FALSE)
  }
  if (square && nrow(x) != ncol(x)) {
    stop(name, " must be square.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_axes <- function(regions, sectors, n) {
  check <- function(x, name) {
    if (!is.character(x) || length(x) == 0L || anyNA(x) ||
        any(!nzchar(x)) || anyDuplicated(x)) {
      stop(name, " must contain unique non-empty strings.", call. = FALSE)
    }
  }
  check(regions, "regions")
  check(sectors, "sectors")
  if (length(regions) * length(sectors) != n) {
    stop("length(regions) * length(sectors) must equal nrow(x).",
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_mrio_tensor <- function(x) {
  if (!inherits(x, "mrio_tensor") || !is.array(x) ||
      length(dim(x)) != 4L || !is.numeric(x) ||
      anyNA(x) || any(!is.finite(x)) ||
      dim(x)[1L] != dim(x)[3L] || dim(x)[2L] != dim(x)[4L]) {
    stop("x must be a finite four-mode mrio_tensor with matching axes.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_technology <- function(A) {
  .assert_matrix(A, "A", square = TRUE)
  if (any(A < 0)) {
    stop("A must be a nonnegative technical-coefficient matrix.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_nonnegative <- function(x, n, name) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != n ||
      anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    stop(name, " must be a finite, nonnegative numeric vector of length ",
         n, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.safe_solve <- function(B, rhs = NULL) {
  tryCatch({
    if (is.null(rhs)) solve(B) else solve(B, rhs)
  }, error = function(e) {
    stop("The Leontief system could not be solved; check I - A for ",
         "singularity or poor conditioning.", call. = FALSE)
  })
}
