# Mathematical notes and precise interpretation

## 1. Index layout and lossless tensorization

With `R` regions and `S` sectors, each composite node is `i=(r-1)S+s`.
For a square input matrix `Z in R^(RS x RS)`, define

`T[r,s,r2,s2] = Z[(r-1)S+s, (r2-1)S+s2]`.

This is an index re-expression; it preserves every element. R uses column-major
array storage, so the implementation reshapes to `(S,R,S,R)` before applying
permutation `(2,1,4,3)`. The inverse applies the reverse transformation.
Matrices with different node ordering must be reordered by the caller first.

## 2. Input-output accounting identities

Given an intermediate transaction matrix `Z` and a positive output vector `x`,
the familiar direct requirements matrix uses `A_ij = Z_ij / x_j` (purchaser
output in the denominator). This package deliberately starts *from* `A`, `y`
and `f` to avoid embedding data-cleaning and untested accounting decisions.

`x = A x + y` implies `(I - A)x = y`.
The software solves this system without materializing its full inverse.
If `rho(A) < 1`, the inverse admits the convergent Neumann representation:

`L = (I - A)^(-1) = I + A + A^2 + ...`.

Neither an invertible `(I-A)` nor nonnegative `A` alone ensures a productive
input-output system. Feasibility and units remain the analyst's responsibility.

## 3. Embodied-intensity attribution

Let `f_i` be the direct environmental intensity of producer sector `i`, and
`y_j` final demand for the product of sector `j`.

`F_ij = f_i L_ij y_j`

Then `sum_ij F_ij = f' L y`. Row sums attribute emissions to producers; column
sums attribute upstream embodied impacts to final-demand products.

**Semantic boundary:** This decomposition is not the intermediate transaction
matrix `Z`, nor a complete bilateral transfer/trade flow tracing method. Users
must not relabel its cells as physical shipments between sectors.

## 4. Aggregation versus identification

Summing over origin/destination sectors leaves region-by-region marginals;
summing over origin/destination regions leaves sector-by-sector marginals.
These are descriptive projections. They do not identify a technology effect
or geography effect; those require additional modeling assumptions.

## 5. Mode product and HOSVD

For a K-mode tensor `T` with dimensions `d1 x ... x dK`, its mode-m product
with a `J x dm` matrix `M` has dimensions `(d1,...,J,...,dK)` and is defined by
summing over the old mode-m index. A truncated HOSVD computes leading left
singular vectors of each mode unfolding and projects to a smaller Tucker core.

These subspaces are basis-dependent up to rotations/signs, generally describe
approximation geometry rather than causal economic forces, and do not remove
sampling or model uncertainty.

## 6. Order totals are not full structural path analysis

`c_k = f' A^k y` counts aggregate contribution from supply-chain order `k`.
It does *not* list or rank individual paths and does not eliminate combinatorial
complexity when the analyst needs all individual routes. The sequence is
computed via successive matrix-vector multiplication.

## 7. Scalability boundary

Dense tensorization stores `O((RS)^2)` elements. A full Leontief inverse uses
`O((RS)^2)` storage and is generally more costly to compute than a single
right-hand-side solve. The package makes no speed or storage advantage claim
relative to established, optimized matrix implementations.
