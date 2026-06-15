# Known Issues

## Shu Table 1 Nx=13 Isotropic SS-SS n=0 Discrepancy

Shu (1996) Table 1 states that 13 grid points were used for the isotropic conical shell validation cases. For the `SS-SS`, `n=0` case, this implementation gives:

```text
Nx = 13: lambda = 0.18016
```

The paper table reports:

```text
lambda = 0.1493
```

Higher grid resolutions converge near the paper value. The current validation notes record:

```text
Nx = 17: lambda = 0.15324
Nx = 21: lambda = 0.14948
Nx = 99: lambda = 0.15059
```

Diagnostics so far indicate that the selected eigenvalue is not an `eigs` targeting artifact, and GDQ derivative matrix tests pass for low-order polynomial derivatives. No paper-supported formula change has been identified.

This discrepancy is documented and regression-tested in `tests/test_paper_validation.m`; it has not been tuned away.

## Appendix 1 Coefficients

The Appendix 1 `S_ijk` coefficient expressions are transcribed into the implementation and benchmarked against selected Shu (1996) validation tables. They have not been independently rederived coefficient-by-coefficient in this repository.

