test_that("Leontief output solves the balance equations", {
  e <- example_mrio()
  x <- leontief_output(e$A, e$y)
  expect_equal(as.vector((diag(4) - e$A) %*% x), e$y,
               tolerance = 1e-10)
})

test_that("pairwise emissions aggregate to total footprint", {
  e <- example_mrio()
  result <- embodied_transfer(e$A, e$y, e$intensity,
                              e$regions, e$sectors)
  expect_equal(result$total,
               sum(e$intensity * leontief_output(e$A, e$y)),
               tolerance = 1e-10)
  expect_equal(unname(as_mrio_matrix(result$tensor)),
               unname(result$matrix))
})

test_that("zero technical coefficients yield only same-sector effects", {
  A <- matrix(0, 2, 2)
  y <- c(4, 5)
  f <- c(0.5, 2)
  result <- embodied_transfer(A, y, f)
  expect_equal(unname(result$matrix), diag(f * y))
  expect_equal(leontief_output(A, y), y)
})

test_that("successive input-output path orders follow A^k", {
  A <- matrix(c(0, 0.2, 0, 0), 2, 2)
  o <- path_order_contributions(A, c(10, 20), c(2, 3), max_order = 3)
  expect_equal(o$contribution, c(80, 6, 0, 0))
  expect_equal(tail(o$cumulative, 1), 86)
})

test_that("nonnegative and safety preconditions are enforced", {
  e <- example_mrio()
  bad <- e$A
  bad[1, 1] <- -0.1
  expect_error(leontief_output(bad, e$y), "nonnegative")
  expect_error(embodied_transfer(e$A, e$y, e$intensity, max_nodes = 3),
               "limit")
  expect_error(path_order_contributions(e$A, e$y, e$intensity,
                                        max_order = -1), "nonnegative integer")
})
