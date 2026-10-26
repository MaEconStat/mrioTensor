# mrioTensor 0.1.0

- Initial release of core MRIO matrix to 4-mode tensor transformation interfaces (`as_mrio_tensor()`, `as_mrio_matrix()`).
- Added mode-wise marginal aggregations (`aggregate_mrio()`) and general mode-n tensor contractions (`mode_product()`).
- Implemented standard Leontief gross output solve (`leontief_output()`) and embodied environmental intensity attribution (`embodied_transfer()`).
- Added production tier Neumann-series decomposition (`path_order_contributions()`) and truncated Higher-Order SVD (`mrio_hosvd()`).
- Added synthetic benchmark data fixture (`example_mrio()`), testthat test suites, and GitHub Actions CI workflow.
