# TASC-DMD

MATLAB implementation of **Time-Augmented, Space-Contracted Dynamic Mode Decomposition (TASC-DMD)**.

> **Paper:** *[TASC DMD Paper — title and citation to be added after publication]*

---

## Repository Structure

```
TASC-DMD/
├── src/
│   ├── DMD_exact.m          Core DMD solver with time-delay embedding
│   ├── fbDMD.m              Forward-Backward DMD (noise-robust)
│   └── nls_rhs.m            ODE RHS for the NLS benchmark equation
├── examples/
│   ├── generate_NLS_data.m         Generate NLS benchmark data (paper)
│   ├── Schrodinger_equ_TASCDMD.m   Main NLS benchmark (reproduces paper results)
│   ├── generate_AdvDiff_data.m     Generate advection-diffusion data (paper)
│   ├── AdvDiff_TASCDMD.m           Advection-diffusion benchmark (paper)
│   ├── generate_Burgers_data.m     Generate inviscid Burgers data (NOT in paper)
│   └── Burgers_TASCDMD.m           Burgers benchmark (NOT in paper)
├── docs/
│   └── index.html           GitHub Pages project site
└── README.md
```

---

## Quick Start

```matlab
% Add src/ to path
addpath('src');
cd examples

% --- NLS benchmark (reproduces paper results) ---
run generate_NLS_data.m          % saves NLS_data.mat
run Schrodinger_equ_TASCDMD.m    % comparison figures

% --- Advection-Diffusion benchmark (known eigenvalues, paper) ---
run generate_AdvDiff_data.m      % saves AdvDiff_data.mat
run AdvDiff_TASCDMD.m            % eigenvalue parabola comparison

% --- Inviscid Burgers (NOT in paper — additional example) ---
run generate_Burgers_data.m      % saves Burgers_data.mat
run Burgers_TASCDMD.m            % reconstruction + Cv comparison
```

---

## Benchmarks

### NLS — Nonlinear Schrödinger Equation *(paper)*
Two-soliton solution. Tests eigenvalue recovery under 20% noise.

### Advection-Diffusion *(paper)*
Linear PDE with analytically known eigenvalues lying on a parabola
`Re(ω) = −(ν/U₀²)·Im(ω)²`. Analytical curve overlaid on DMD spectrum.
Tested at 50% noise.

### Inviscid Burgers *(NOT in paper — additional example)*
Small-amplitude traveling waves with near-imaginary eigenvalues.
Demonstrates TASC-DMD generality on a nonlinear hyperbolic PDE.
Tested at 50% noise. Cv reported for all methods.

---

## Key Parameters

| Parameter | Symbol | Description | Default |
|-----------|--------|-------------|---------|
| `d`  | *d*  | Time-delay levels (Hankel depth) | 11 |
| `r0` | *r₀* | Spatial basis rank per delay | 11 |
| `r1` | *r₁* | Final DMD rank | 9 |
| `noise_level` | — | Gaussian noise fraction | 0.20–0.50 |

---

## Citation

> *[To be populated after publication]*

---

## License

MIT License. See `LICENSE` for details.
