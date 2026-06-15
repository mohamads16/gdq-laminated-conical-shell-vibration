# GDQ Laminated Conical Shell Vibration

MATLAB implementation of a generalized differential quadrature (GDQ) method for free vibration analysis of isotropic and selected laminated conical shells, based on Shu (1996).

The code solves the separated conical-shell equations for circumferential wave number `n`, simply supported (`SS`) and clamped (`C`) edge combinations, and reports the nondimensional frequency parameter

```text
lambda = R2 * sqrt(rho*h/A11) * omega
```

## Quick Start

Requires MATLAB R2025b or a compatible recent MATLAB release.

```matlab
addpath('src')

R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;

lam = laminate_isotropic(1, 0.3, 1, h);
params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',3, ...
                'rho',lam.rho, 'h',lam.h, ...
                'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',33);

mdl = build_conical_shell_gdq_model(params);
sol = solve_conical_shell_modes(mdl, 'SS', 'SS', 6);
disp(sol.lambda)
```

## Examples

Run from the repository root:

```matlab
run('examples/validate_table1_isotropic.m')
run('examples/validate_table2_crossply.m')
run('examples/convergence_isotropic_ssss.m')
run('examples/plot_isotropic_mode_shapes.m')
```

## Tests

```matlab
results = runtests('tests');
table(results)
```

The tests check GDQ derivative weights, edge-label behavior, solver metadata, laminate coupling guardrails, and selected paper-table reproductions.

## Repository Layout

```text
src/                 Core MATLAB functions
examples/            Runnable examples and paper-table validation scripts
tests/               MATLAB unit tests
docs/validation/     Baseline and validation notes
archive/             Broken, duplicate, or legacy original files
```

## Validation Status

Implemented and benchmarked against selected Shu (1996) tables:

- Shu Appendix 1 coefficient assembly for equations (25a)-(25c), transcribed from the paper
- Shu Appendix 2 simply supported boundary coefficients, transcribed from the paper
- Clamped and simply supported boundary combinations
- Boundary condensation eigenproblem from equations (32)-(36)
- Table 1 isotropic conical shell validation at high grid resolution
- Table 2 antisymmetric cross-ply convergence validation

## Laminate Scope

The current solver uses the ABD terms that appear in Shu's governing equations: `11`, `12`, `22`, and `66`. It supports isotropic, cross-ply, and specially orthotropic-style stiffness data where `16/26` coupling terms are negligible.

`laminate_abd_clt` and `laminate_from_unidirectional_plies` reject laminates with non-negligible `A16/A26`, `B16/B26`, or `D16/D26` terms. Fully anisotropic angle-ply laminates need additional governing-equation terms and are intentionally not silently truncated.

Known limitation:

- Shu Table 1 states that 13 grid points were used, but this implementation gives `lambda=0.18016` for the isotropic `SS-SS`, `n=0`, `Nx=13` case and converges to about `0.15059` for larger `Nx`, close to the table value `0.1493`. Derivative matrix and eigenspectrum diagnostics did not identify a confirmed code bug. This is documented in `docs/validation/numerical-fixes.md`.
- Appendix 1 coefficient formulas are transcribed into code but not independently rederived in this repository; see `docs/validation/equation-map.md`.

## Paper PDF

The local file `shu1996.pdf` is ignored by Git. Do not publish it in this repository unless redistribution rights are confirmed. The paper citation is included in `CITATION.cff`.

## Citation

If you use this code, cite both this repository and:

Shu, C. (1996). Free vibration analysis of composite laminated conical shells by generalized differential quadrature. *Journal of Sound and Vibration*, 194(4), 587-604.

## License

MIT License. See `LICENSE`.
