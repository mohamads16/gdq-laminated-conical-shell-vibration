# Numerical Fixes and Evidence

## Boundary Label Fix

Old behavior:

- GDQ nodes were ordered from `x=L` to `x=0`.
- `bcSmall` was applied to the large-radius edge and `bcLarge` to the small-radius edge.
- The original high-resolution `SS-C` output matched Shu Table 1 `C-SS`, and original `C-SS` matched Shu Table 1 `SS-C`.

Evidence:

- Shu equation (37) defines `x_i = (L/2)[1 - cos(pi(i-1)/(N-1))]`, so `x_1=0` and `x_N=L`.
- Geometry equation (9) gives `R(x)=R1+x sin(alpha)`, so `x=0` is the small edge.

New behavior:

- `gdq_chebyshev_lobatto` returns ascending nodes from `xL` to `xR`.
- `bcSmall` is applied at `x=0`, and `bcLarge` is applied at `x=L`.

Validation:

```text
After fix, Nx=99:
SS-C 0.97486 0.74509 0.45060 0.29766 0.25218 0.26397 0.29343 0.33489 0.38692 0.44849
C-SS 0.60171 0.71847 0.37629 0.23236 0.19770 0.21654 0.25241 0.29769 0.35157 0.41371
```

These now align with Shu Table 1 labels.

## Cross-Ply Coupling Factor Fix

Old formula:

```matlab
B11 = 2 * (h^2 / (4 * Nlayers)) * (Q11 - Q22);
B22 = -B11;
```

Old behavior:

```text
Nl=2, h/R2=0.01, N=13: lambda = 0.28535
Paper Table 2:             lambda = 0.1799
```

Evidence:

- Shu equation (40) gives the cross-ply stiffness simplification.
- The text immediately after Table 5 states that for maximum coupling, `B11 = -B22 = (Q11-Q22) h^2/8`.
- For `Nlayers=2`, that corresponds to `B11 = (h^2/(4*Nlayers))*(Q11-Q22)`.

New formula:

```matlab
B11 = (h^2 / (4 * Nlayers)) * (Q11 - Q22);
B22 = -B11;
```

Validation:

```text
Nl=2, h/R2=0.01, N=13: lambda = 0.17786
Paper Table 2:             lambda = 0.1799
Remaining discrepancy:     -0.00204
```

## Minimum Grid Count

Old behavior:

- `build_conical_shell_gdq_model` rejected `Nx < 7`.
- Shu Table 2 includes `N=5`.

New behavior:

- The guard is `Nx >= 5`, because the highest derivative order is four.
- `N=5` now runs and reproduces the coarse Table 2 row qualitatively.

Validation:

```text
Table 2 after fix:
Nl=Inf N5=0.31578 N7=0.22398 N9=0.19636 N11=0.19739 N13=0.19778 N21=0.19779
Nl=2   N5=0.32146 N7=0.22097 N9=0.17876 N11=0.17660 N13=0.17786 N21=0.17792
```

## Unresolved Nx=13 Isotropic SS-SS n=0 Discrepancy

Observation:

```text
Shu Table 1, SS-SS, n=0: 0.1493
This code, Nx=13:          0.18016
This code, Nx=21:          0.14948
This code, Nx=99:          0.15059
```

Diagnostics:

- Full eigenvalue solve gives the same first eigenvalue as `eigs`, so this is not an eigensolver target-selection issue.
- GDQ derivative matrices differentiate low-order polynomials to roundoff.
- No paper-supported formula change was identified.

Decision:

- No tuning or unsupported numerical change was made.
- The discrepancy is documented and covered by the convergence example.

