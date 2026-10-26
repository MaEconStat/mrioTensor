# Contributing

This is a small research-software preview. Corrections to documentation,
index consistency, input validation and synthetic test cases are welcome.
Please open an issue describing a minimal reproducible example first.
Do not upload licensed MRIO databases, unpublished manuscripts or identifiable
research datasets to issues or pull requests.

## Local quality checks

```bash
R CMD INSTALL .
Rscript scripts/quickcheck.R
R CMD build .
R CMD check --no-manual mrioTensor_0.1.0.tar.gz
```

Run unit tests after installing `testthat`:

```r
install.packages("testthat")
testthat::test_local(".")
```
