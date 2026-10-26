# Manual smoke checks after R CMD INSTALL .
# All values are synthetic.
library(mrioTensor)
x <- example_mrio()
T <- as_mrio_tensor(x$A, x$regions, x$sectors)
stopifnot(isTRUE(all.equal(unname(as_mrio_matrix(T)), unname(x$A))))
cat("Round trip and dimensions OK: ", paste(dim(T), collapse = " x "), "\n")
cat("Gross output: ", paste(round(leontief_output(x$A, x$y), 3),
                             collapse = ", "), "\n")
res <- embodied_transfer(x$A, x$y, x$intensity,
                         x$regions, x$sectors)
stopifnot(isTRUE(all.equal(res$total,
      sum(x$intensity * leontief_output(x$A, x$y)))))
cat("Attributed footprint: ", round(res$total, 4), "\n")
fit <- mrio_hosvd(T, ranks = c(2, 2, 2, 2))
cat("Full-rank HOSVD error: ", signif(fit$relative_error, 3), "\n")
