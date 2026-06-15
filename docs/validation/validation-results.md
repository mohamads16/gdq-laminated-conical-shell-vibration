# Final Validation Report

Validation was run after the project structure cleanup, API rename, boundary-label fix, and Shu cross-ply coupling correction.

## Environment

- MATLAB: `25.2.0.2998904 (R2025b)`
- Working directory: repository root
- The local paper file `shu1996.pdf` is ignored by Git and is not read by any example or test.

## Commands Run

```matlab
run('examples/validate_table1_isotropic.m')
run('examples/validate_table2_crossply.m')
run('examples/convergence_isotropic_ssss.m')
set(groot,'DefaultFigureVisible','off'); run('examples/plot_isotropic_mode_shapes.m'); close all
results = runtests('tests'); disp(table(results)); assertSuccess(results)
```

## Table 1 Isotropic Validation

Reference: Shu (1996), Table 1, isotropic conical shells with `h/R2 = 0.01`, `nu = 0.3`, `alpha = 30 deg`, and `L sin(alpha)/R2 = 0.75`.

The validation script uses `Nx = 99` to report the converged behavior of this implementation. Maximum absolute errors by boundary condition:

| Boundary condition | Max absolute error |
| --- | ---: |
| `SS-SS` | `0.00179234` |
| `SS-C` | `7.20528e-05` |
| `C-SS` | `9.44814e-05` |
| `C-C` | `0.000640566` |

The mixed-boundary labels now match the physical small/large edge labels in the paper: `SS-C` applies simply supported conditions at `x=0`/`R1` and clamped conditions at `x=L`/`R2`.

## Table 2 Cross-Ply Validation

Reference: Shu (1996), Table 2, antisymmetric cross-ply conical shells with `h/R2 = 0.01`, `alpha = 30 deg`, `L sin(alpha)/R2 = 0.25`, `n = 0`, and `SS-SS` boundaries.

Maximum absolute errors:

| Case | Max absolute error |
| --- | ---: |
| Without coupling, `Nl = Inf` | `0.000291399` |
| Maximum coupling, `Nl = 2` | `0.00595844` |

For the commonly cited converged row, `Nl = 2`, `N = 13`, the corrected code gives `lambda = 0.17786` versus the paper value `0.1799`, with absolute error `0.0020395`.

## Test Result Summary

Final MATLAB unit test run:

| Test | Result |
| --- | --- |
| `test_gdq_derivatives/testPolynomialDerivatives` | passed |
| `test_gdq_derivatives/testNodeOrdering` | passed |
| `test_laminate_support/testCrossPlyLaminateIsAccepted` | passed |
| `test_laminate_support/testAnglePlyLaminateIsRejected` | passed |
| `test_paper_validation/testTable1IsotropicAtHighResolution` | passed |
| `test_paper_validation/testTable2CrossplyConvergence` | passed |
| `test_paper_validation/testKnownNx13IsotropicMismatch` | passed |
| `test_solver_metadata/testEigenFilterAndBoundaryDiagnosticsExist` | passed |

Overall: 8 passed, 0 failed, 0 incomplete.

## Unresolved Nx=13 Discrepancy

Shu Table 1 states that 13 grid points were used. For the isotropic `SS-SS`, `n=0` case, this implementation gives:

| `Nx` | `lambda` | Paper Table 1 |
| ---: | ---: | ---: |
| 13 | `0.18016` | `0.1493` |
| 17 | `0.15324` | `0.1493` |
| 21 | `0.14948` | `0.1493` |
| 99 | `0.15059` | `0.1493` |

Diagnostics showed that the full eigenspectrum agrees with the selected eigenvalue, and the GDQ derivative matrices differentiate low-order polynomials to roundoff. No paper-supported formula change was identified, so this discrepancy is documented rather than tuned away.
