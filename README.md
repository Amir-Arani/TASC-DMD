# TASC-DMD

MATLAB implementation of **Time-Augmented, Space-Contracted Dynamic Mode Decomposition (TASC-DMD)**.

> **Paper:** *[TASC DMD Paper — title and citation to be added after publication]*

---

## Repository Structure

```
TASC-DMD/
├── src/
│   ├── DMD_exact.m                  Core DMD solver with time-delay embedding
│   ├── fbDMD.m                      Forward-Backward DMD (noise-robust)
│   └── nls_rhs.m                    ODE RHS for the NLS benchmark
├── examples/
│   ├── generate_NLS_data.m          Generate NLS benchmark data  [paper]
│   ├── Schrodinger_equ_TASCDMD.m    NLS benchmark — reproduces paper results  [paper]
│   ├── generate_AdvDiff_data.m      Generate advection-diffusion data  [NOT in paper]
│   └── AdvDiff_TASCDMD.m            Advection-diffusion benchmark  [NOT in paper]
├── docs/
│   └── index.html                   GitHub Pages project site
└── README.md
```

---

## Quick Start

```matlab
addpath('src');
cd examples

% --- NLS benchmark (reproduces paper results) ---
run generate_NLS_data.m
run Schrodinger_equ_TASCDMD.m

% --- Advection-Diffusion (NOT in paper — additional example) ---
run generate_AdvDiff_data.m
run AdvDiff_TASCDMD.m
```

---

## Benchmarks

### NLS — Nonlinear Schrödinger Equation *(paper)*
Two-soliton solution. TASC-DMD vs Exact DMD vs fbDMD under 20% noise.

### Advection-Diffusion *(NOT in paper — additional example)*
Linear PDE `u_t + U0*u_x = nu*u_xx` with analytically known eigenvalues:

```
omega_n = -nu*k_n^2 - i*U0*k_n
```

These lie on the parabola `Re(ω) = -(ν/U₀²)·Im(ω)²`, which is overlaid on the
eigenvalue plot as ground truth. Tested at **50% noise**. Cv reported for all methods.

---

## Key Parameters

| Parameter | Symbol | Description | NLS | AdvDiff |
|-----------|--------|-------------|-----|---------|
| `d`  | *d*  | Time-delay levels (Hankel depth) | 11 | 11 |
| `r0` | *r₀* | Spatial basis rank per delay     | 11 | 11 |
| `r1` | *r₁* | Final DMD rank                   | 9  | 9  |
| `noise_level` | — | Gaussian noise fraction     | 0.20 | 0.50 |

---

## Citation

> *[To be populated after publication]*

---

## License

MIT License. See `LICENSE` for details.
