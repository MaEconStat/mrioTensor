test_that("mode multiplication by identity preserves arrays", {
  x <- array(seq_len(24), dim = c(2, 3, 4))
  for (i in seq_along(dim(x))) {
    expect_equal(mode_product(x, diag(dim(x)[i]), i), x)
  }
})

test_that("mode multiplication sums the selected mode correctly", {
  x <- array(seq_len(24), dim = c(2, 3, 4))
  out <- mode_product(x, matrix(1, nrow = 1, ncol = 3), 2)
  expect_equal(dim(out), c(2L, 1L, 4L))
  expect_equal(drop(out), apply(x, c(1, 3), sum))
})

test_that("full-rank HOSVD reconstructs the small input", {
  x <- array(seq_len(16) / 17, dim = c(2, 2, 2, 2))
  fit <- mrio_hosvd(x, ranks = c(2, 2, 2, 2))
  expect_equal(fit$reconstruction, x, tolerance = 1e-10)
  expect_lt(fit$relative_error, 1e-10)
  expect_equal(dim(fit$core), c(2L, 2L, 2L, 2L))
})

test_that("HOSVD refuses large dense inputs by default", {
  x <- array(1, dim = c(5, 5, 5))
  expect_error(mrio_hosvd(x, c(2, 2, 2), max_entries = 100), "limit")
})
