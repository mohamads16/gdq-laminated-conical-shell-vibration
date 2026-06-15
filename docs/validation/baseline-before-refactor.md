# Baseline Before Refactor

These outputs were captured from the original flat MATLAB folder before moving files or changing numerical behavior.

## Isotropic Table 1-Style Outputs

### Nx = 13

```text
SS-SS 0.18016 0.69275 0.36005 0.21236 0.19080 0.21607 0.25238 0.29765 0.35151 0.41359
SS-C  0.59886 0.71810 0.37719 0.23304 0.19788 0.21651 0.25238 0.29765 0.35151 0.41359
C-SS  0.98032 0.74552 0.45113 0.29816 0.25229 0.26388 0.29331 0.33475 0.38677 0.44836
C-C   0.99293 0.75261 0.46606 0.31564 0.26009 0.26458 0.29332 0.33475 0.38677 0.44836
```

### Nx = 99

```text
SS-SS 0.15059 0.69453 0.35930 0.21170 0.19066 0.21608 0.25241 0.29769 0.35157 0.41371
SS-C  0.60171 0.71847 0.37629 0.23236 0.19770 0.21654 0.25241 0.29769 0.35157 0.41371
C-SS  0.97486 0.74509 0.45060 0.29766 0.25218 0.26397 0.29343 0.33489 0.38692 0.44849
C-C   0.98851 0.75254 0.46592 0.31535 0.25996 0.26466 0.29344 0.33489 0.38692 0.44849
```

The original mixed-boundary labels were physically swapped relative to Shu Table 1.

## Cross-Ply Table 2-Style Outputs

```text
layers=Inf N7=0.22398 N9=0.19636 N11=0.19739 N13=0.19778 N21=0.19779
layers=2   N7=0.21468 N9=0.28132 N11=0.28559 N13=0.28535 N21=0.28535
```

The original `Nl=2` maximum-coupling case was far from the paper values.

## Broken/Suspicious Original Files

- `Draft.m`: failed with undefined `geom.L`, after printing a scratch row-index vector.
- `step2.asv`: MATLAB autosave duplicate of `step2.m`.
- `edit1.m`: debug plotting script.
- `untitled3.m` and `untitled10.m`: duplicate validation/plot scripts.
- `laminate_from_ud`: local function hidden inside `laminate_isotropic.m`, not externally callable.
- `bary_weights_safe`: unused local helper.

