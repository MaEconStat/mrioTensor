test_that("matrix to tensor indexing is exact and reversible", {
  m <- matrix(seq_len(36), 6, 6)
  t <- as_mrio_tensor(m, c("A", "B"), c("s1", "s2", "s3"))
  expect_s3_class(t, "mrio_tensor")
  expect_equal(dim(t), c(2L, 3L, 2L, 3L))
  expect_equal(t[2, 1, 1, 3], m[4, 3])
  expect_equal(unname(as_mrio_matrix(t)), unname(m))
})

test_that("region and sector marginals match explicit sums", {
  e <- example_mrio()
  t <- as_mrio_tensor(e$A, e$regions, e$sectors)
  reg <- aggregate_mrio(t, by = "region")
  sec <- aggregate_mrio(t, by = "sector")
  expect_equal(reg[1, 2], sum(e$A[1:2, 3:4]))
  expect_equal(sec[1, 2], sum(e$A[c(1, 3), c(2, 4)]))
  expect_equal(sum(reg), sum(e$A))
  expect_equal(sum(sec), sum(e$A))
})

test_that("invalid axes and non-square inputs are rejected", {
  expect_error(as_mrio_tensor(matrix(1:6, 2, 3), "A", "a"), "square")
  expect_error(as_mrio_tensor(diag(2), "A", c("x", "y")), "must equal")
  expect_error(as_mrio_tensor(diag(2), "A", c("x", "x")), "unique")
})
