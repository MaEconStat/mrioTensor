# mrioTensor: Tensor-Aware Multi-Regional Input-Output Analysis

[![R-CMD-check](https://github.com/MaEconStat/mrioTensor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MaEconStat/mrioTensor/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

`mrioTensor` provides base-R tools for reshaping Multi-Regional Input-Output (MRIO) tables into semantically structured fourth-order tensors, evaluating multilinear contractions, and performing Leontief-based output and embodied environmental flow accounting.

---

## Motivation

In global economic and environmental accounting (EEIOA), multi-regional input-output data are conventionally stored and manipulated as flattened 2D block matrices $Z \in \mathbb{R}^{N \times N}$, where each index $i = (r - 1)S + s$ conflates geographic location (region $r$) with production technology (sector $s$).

While standard 2D operations suffice for basic matrix inversion, flattening obscures the multi-way structural properties of global value chains. In particular:
1. **Structural semantics**: Geographic trade linkages (region $\to$ region) and industrial input-output linkages (sector $\to$ sector) are entangled in a single composite index.
2. **Multilinear operations**: Marginalizing, contracting, or low-rank approximating intermediate flows across specific dimensions (e.g., aggregating out sectors to obtain pure inter-regional trade, or aggregating out regions to evaluate global technology matrices) requires clumsy matrix slicing.
3. **Higher-order representations**: Generalizing to multiple environmental extensions or multi-year panels naturally yields higher-order arrays (5-th or 6-th order tensors) that benefit from multilinear algebra (such as Higher-Order SVD and Tucker decomposition).

`mrioTensor` bridges standard Leontief input-output analysis and multilinear array operations in base R.

---

## Core Capabilities

| Function | Description | Mathematical Operation |
| :--- | :--- | :--- |
| `as_mrio_tensor()` | Reshape 2D composite matrix to 4D tensor | $\mathcal{T}_{r, s, r', s'} \leftarrow Z_{(r-1)S+s, (r'-1)S+s'}$ |
| `as_mrio_matrix()` | Flatten 4D tensor back to 2D composite matrix | $Z_{(r-1)S+s, (r'-1)S+s'} \leftarrow \mathcal{T}_{r, s, r', s'}$ |
| `aggregate_mrio()` | Marginalize over regional or sectoral dimensions | Inter-regional $\sum_{s, s'} \mathcal{T}$ or Inter-sectoral $\sum_{r, r'} \mathcal{T}$ |
| `mode_product()` | Mode-$n$ matrix-tensor contraction | $\mathcal{Y} = \mathcal{X} \times_n M$ |
| `leontief_output()` | Solve gross output balance without explicit matrix inversion | $(I - A) x = y$ |
| `embodied_transfer()` | Pairwise embodied environmental intensity attribution | $F_{ij} = f_i L_{ij} y_j$, where $L = (I - A)^{-1}$ |
| `path_order_contributions()` | Evaluate Neumann series tier contributions | $c_k = f^\top A^k y$ for tier $k = 0, 1, \dots$ |
| `mrio_hosvd()` | Truncated Higher-Order SVD (Tucker representation) | $\mathcal{T} \approx \mathcal{G} \times_1 U_1 \times_2 U_2 \times_3 U_3 \times_4 U_4$ |
| `example_mrio()` | Synthetic 2-region, 2-sector benchmark data | Verification and testing fixture |

---

## Installation

`mrioTensor` requires **R >= 4.1.0** and relies exclusively on base R for core calculations (zero external runtime dependencies).

```r
# Install development version from GitHub:
# install.packages("remotes")
remotes::install_github("MaEconStat/mrioTensor")
```

Or install from a local source clone:

```bash
git clone https://github.com/MaEconStat/mrioTensor.git
R CMD INSTALL mrioTensor
```

---

## Quickstart

```r
library(mrioTensor)

# Load synthetic 2-region (North, South) x 2-sector (Agriculture, Manufacturing) system
e <- example_mrio()

# 1. Exact 2D -> 4D tensor transformation
T_arr <- as_mrio_tensor(e$A, regions = e$regions, sectors = e$sectors)
dim(T_arr)
# [1] 2 2 2 2
# Dimensions: (origin_region, origin_sector, destination_region, destination_sector)

# Round-trip verification:
stopifnot(isTRUE(all.equal(unname(as_mrio_matrix(T_arr)), unname(e$A))))

# 2. Mode-wise marginal aggregations
reg_trade  <- aggregate_mrio(T_arr, by = "region")  # Inter-regional flow matrix (2 x 2)
tech_links <- aggregate_mrio(T_arr, by = "sector")  # Inter-sectoral technical matrix (2 x 2)

# 3. Leontief gross output solve: (I - A) x = y
x <- leontief_output(e$A, e$y)

# 4. Embodied environmental intensity accounting
# Allocates upstream emissions from supply sectors to final consumption categories
account <- embodied_transfer(
  A = e$A, y = e$y, intensity = e$intensity,
  regions = e$regions, sectors = e$sectors
)

# Total footprint identity: sum(F) == f' x
stopifnot(isTRUE(all.equal(account$total, sum(e$intensity * x))))

# 5. Neumann-series production tier decomposition: c_k = f' A^k y
path_order_contributions(e$A, e$y, e$intensity, max_order = 4)

# 6. Truncated Higher-Order SVD (Tucker decomposition)
hosvd_fit <- mrio_hosvd(T_arr, ranks = c(2, 2, 2, 2))
cat("Relative approximation error:", round(hosvd_fit$relative_error, 4), "\n")
```

---

## Mathematical Formulation

### 1. Tensor Indexing Layout
Let $R$ denote the number of regions and $S$ the number of sectors per region, yielding $N = RS$ composite nodes. Under the canonical **region-major, sector-minor** ordering:
$$i = (r - 1)S + s, \quad j = (r' - 1)S + s'$$
The 4-mode tensor $\mathcal{T} \in \mathbb{R}^{R \times S \times R \times S}$ maps bijectively to the flattened matrix $Z \in \mathbb{R}^{N \times N}$ via:
$$\mathcal{T}_{r, s, r', s'} = Z_{(r-1)S+s, \, (r'-1)S+s'}$$

### 2. Leontief Balance and Output Solves
Given technical coefficient matrix $A = Z \hat{x}^{-1}$ and final demand vector $y$, the standard Leontief balance:
$$x = A x + y \implies (I - A)x = y$$
is solved directly via LU factorization (`solve(I - A, y)`), avoiding dense matrix inversion.

When the spectral radius $\rho(A) < 1$, $(I - A)^{-1}$ admits the convergent Neumann series:
$$L = (I - A)^{-1} = \sum_{k=0}^{\infty} A^k$$
Tier contributions $c_k = f^\top A^k y$ are evaluated via sequential matrix-vector contractions $v_{k} = A v_{k-1}$ ($v_0 = y$), maintaining $O(k N^2)$ complexity rather than $O(k N^3)$.

### 3. Embodied Intensity Allocation
For direct environmental intensity vector $f \in \mathbb{R}^N$ (e.g., carbon emissions per unit of gross output), the pairwise attribution matrix $F \in \mathbb{R}^{N \times N}$ is:
$$F = \hat{f} (I - A)^{-1} \hat{y}$$
where $\hat{\cdot}$ denotes diagonalized vectors. Row sums give sector-level production emissions, column sums give footprint attributed to final consumption, and $\mathbf{1}^\top F \mathbf{1} = f^\top x$.

---

## Computational Notes

- **Memory scaling**: The 4-mode array stores $R^2 S^2 = N^2$ double-precision values. For medium-sized tables (e.g., $R=189, S=26 \implies N=4,914$), one dense array occupies $\approx 192\text{ MB}$.
- **Dense inverses**: For large-scale MRIO systems, users should rely on `leontief_output()` rather than computing the full dense inverse $L = (I - A)^{-1}$. `embodied_transfer()` includes a configurable safety threshold (`max_nodes = 1000L`).

---

## Related Projects

- **[LN-MMSBM](https://github.com/MaEconStat/ln-mmsbm-python)**: Python research prototype for continuous-weight mixture multilayer stochastic block models and regularized tensor decomposition (TWIST) on directed economic flow networks.

---

## License

MIT License. Copyright (c) 2026 Zhen Ma (MaEconStat).
