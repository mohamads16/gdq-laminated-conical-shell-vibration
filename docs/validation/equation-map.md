# Equation Map to Shu (1996)

This document maps the MATLAB implementation to the paper:

Shu, C. (1996). "Free vibration analysis of composite laminated conical shells by generalized differential quadrature." *Journal of Sound and Vibration*, 194(4), 587-604.

The goal is traceability, not a full independent derivation. Coefficients copied from the paper appendices are marked as transcribed where they have not been independently rederived in this repository.

## Geometry

Paper:

- Figure 1 defines the conical shell geometry.
- Equation (9): `R(x) = R1 + x sin(alpha)`.
- `L = (R2 - R1) / sin(alpha)`.

Implementation:

- [build_conical_shell_gdq_model.m](../../src/build_conical_shell_gdq_model.m)
  - reads `L`, `R1`, `alpha`, and `Nx` from `params`
  - computes `R = R1 + x*sin(alpha)`
  - stores the large-edge radius as `mdl.geom.R2 = R1 + L*sin(alpha)`

## GDQ Nodes and Weighting Matrices

Paper:

- Equations (5)-(8): GDQ derivative approximation and recurrence for higher-order weighting coefficients.
- Equation (37): `x_i = (L/2)[1 - cos(pi(i-1)/(N-1))]`, so `x_1 = 0` and `x_N = L`.

Implementation:

- [gdq_chebyshev_lobatto.m](../../src/gdq_chebyshev_lobatto.m)
  - builds ascending Chebyshev-Lobatto nodes on `[xL, xR]` to match equation (37)
  - constructs the first derivative matrix using barycentric Chebyshev-Lobatto weights
  - constructs higher derivative matrices using Shu-style recurrence

Validation:

- [test_gdq_derivatives.m](../../tests/test_gdq_derivatives.m) verifies low-order polynomial derivatives to roundoff and checks node ordering.

## Governing Equation Coefficients

Paper:

- Equations (23a)-(23c): coupled shell equations after substituting resultants.
- Equations (24a)-(24c): modal displacement form with circumferential wave number `n`.
- Equations (25a)-(25c): separated ODE system in `U(x)`, `V(x)`, and `W(x)`.
- Appendix 1: definitions of coefficients `S110` through `S334`.

Implementation:

- [build_conical_shell_gdq_model.m](../../src/build_conical_shell_gdq_model.m)
  - extracts `A11`, `A12`, `A22`, `A66`, `B11`, `B12`, `B22`, `B66`, `D11`, `D12`, `D22`, and `D66`
  - computes `S110` through `S334`
  - stores them under `mdl.S`
- [solve_conical_shell_modes.m](../../src/solve_conical_shell_modes.m)
  - forms the block operators for the `U`, `V`, and `W` equations using `diag(S)*Dk`

Status:

- The Appendix 1 coefficient expressions are transcribed from Shu (1996).
- They are benchmarked indirectly against selected paper tables.
- They have not been independently rederived term-by-term in this repository.

## Boundary Conditions

Paper:

- Equation (26): simply supported end conditions `V=0`, `W=0`, `Nx=0`, `Mx=0`.
- Equations (27a)-(27c): simply supported conditions rewritten in terms of `U`, `V`, `W`, and derivatives using coefficients `a1...a5`, `b1...b5`.
- Equation (28): clamped end conditions `U=0`, `V=0`, `W=0`, `W' = 0`.
- Appendix 2: definitions of `a1...a5` and `b1...b5`.
- Equations (30)-(31): GDQ-discretized boundary rows at small and large edges.

Implementation:

- [build_conical_shell_gdq_model.m](../../src/build_conical_shell_gdq_model.m)
  - computes Appendix 2 coefficients for `BC.small` at `x=0`/`R1`
  - computes Appendix 2 coefficients for `BC.large` at `x=L`/`R2`
- [solve_conical_shell_modes.m](../../src/solve_conical_shell_modes.m)
  - `build_boundary_rows` assembles the eight boundary equations
  - `add_SS_rows` implements the simply supported rows
  - `add_C_rows` implements the clamped rows

Status:

- The Appendix 2 coefficient expressions are transcribed from Shu (1996).
- Edge ordering is explicitly tied to equation (37), so `bcSmall` applies at `x=0` and `bcLarge` applies at `x=L`.

## Condensed Eigenproblem

Paper:

- Equation (32): interior unknown vector `XI`.
- Equation (33): boundary unknown vector `XB`.
- Equation (34): interior equation partition.
- Equation (35): boundary equation partition.
- Equation (36): condensed eigenproblem.

Implementation:

- [solve_conical_shell_modes.m](../../src/solve_conical_shell_modes.m)
  - uses interior rows `U(2:N-1)`, `V(2:N-1)`, and `W(3:N-2)`
  - uses boundary unknowns `U(1)`, `U(N)`, `V(1)`, `V(N)`, `W(1)`, `W(2)`, `W(N-1)`, `W(N)`
  - computes `T = -(B_B \ B_I)` so `X_B = T X_I`
  - forms `Kred = -(A_II + A_IB*T)` to solve for positive `mu = rho*h*omega^2`

Diagnostics:

- `sol.meta.boundaryBlockRcond` stores `rcond(B_B)`.
- A warning is emitted when the boundary block reciprocal condition estimate is below `1e-12`.
- `sol.meta.eigenFilter` records eigenvalue filtering tolerances and rejected eigenvalues.

## Frequency Normalization

Paper:

- Equation (38): nondimensional frequency parameter `lambda = R2 * sqrt(rho*h/A11) * omega`.

Implementation:

- [solve_conical_shell_modes.m](../../src/solve_conical_shell_modes.m)
  - solves for `mu = rho*h*omega^2`
  - computes `omega = sqrt(mu/(rho*h))`
  - computes `lambda = R2 * sqrt(mu/A11)`

## Cross-Ply Stiffness Formulas

Paper:

- Equation (39): material ratios for the antisymmetric cross-ply examples.
- Equation (40): simplified `A`, `B`, and `D` stiffness coefficients for antisymmetric cross-ply shells.
- Text following the tables states that the maximum-coupling two-layer case has `B11 = -B22 = (Q11-Q22)h^2/8`.

Implementation:

- [laminate_shu_crossply.m](../../src/laminate_shu_crossply.m)
  - implements the equation (39) material ratios
  - implements equation (40)-style `A` and `D` terms
  - uses `B11 = (h^2/(4*Nlayers))*(Q11-Q22)` and `B22 = -B11`
  - sets `B = 0` for `Nlayers = Inf`

Validation:

- [validate_table2_crossply.m](../../examples/validate_table2_crossply.m) benchmarks the Table 2 convergence cases.

## Laminate ABD Helper Scope

Paper:

- Equations (20)-(22) define laminate `A`, `B`, and `D` terms.
- The implemented governing equations use the `11`, `12`, `22`, and `66` ABD terms.

Implementation:

- [laminate_abd_clt.m](../../src/laminate_abd_clt.m)
  - computes CLT `A`, `B`, and `D` matrices
  - rejects non-negligible `16/26` terms because the current solver does not include the fully anisotropic governing-equation terms
- [laminate_from_unidirectional_plies.m](../../src/laminate_from_unidirectional_plies.m)
  - convenience wrapper around `laminate_abd_clt`

Status:

- Isotropic, cross-ply, and specially orthotropic-style stiffnesses are supported.
- Fully anisotropic angle-ply laminates are intentionally rejected rather than silently truncated.
