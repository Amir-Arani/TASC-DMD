# TASC-DMD: Time-Augmented, Space-Contracted Dynamic Mode Decomposition

**Paper:** *A Novel Spatiotemporal Decomposition and Identification of Sparse Equations for Human Brain Deformation*

A.H.G. Arani¹², A.A. Alshareef³, D.L. Pham⁴, R.J. Okamoto¹, P.V. Bayly¹

¹ Mechanical Engineering and Materials Science, Washington University in St. Louis  
² Department of Neurosurgery, Washington University School of Medicine  
³ Department of Mechanical Engineering, University of South Carolina  
⁴ Department of Radiology and Radiological Sciences, Uniformed Services University  

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022a%2B-blue)](https://www.mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

---

## Overview

TASC-DMD is a novel dynamic mode decomposition framework for discovering low-dimensional models of complex dynamical systems from spatiotemporal data. It combines:

- **Time-Augmented embedding** — a 3D Hankel-based delay embedding that captures slow-time dynamics
- **Space-Contracted projection** — a sequential SVD reduction that compresses high-dimensional spatial fields into a compact latent space
- **DMD in the latent space** — efficient spectral decomposition with improved noise robustness

The algorithm is validated on the nonlinear Schrödinger equation and applied to characterize 4D strain fields from tagged MRI brain data. A companion SINDy step (TASC-SINDy) discovers sparse governing equations directly in TASC coordinates.

---

## Repository Structure

```
TASC-DMD/
├── src/
│   ├── DMD_exact.m       — Exact (or standard) DMD with optional time-delay embedding
│   ├── fbDMD.m           — Forward-Backward DMD with optional time-delay embedding
│   └── nls_rhs.m         — ODE RHS for the nonlinear Schrödinger equation (Fourier space)
│
├── examples/
│   ├── generate_NLS_data.m        — Solve the NLS equation and save data
│   └── Schrodinger_equ_TASCDMD.m  — Benchmark TASC-DMD vs Exact DMD vs fbDMD
│
├── docs/
│   └── index.html         — GitHub Pages project website
│
└── README.md
```

---

## Requirements

| Requirement        | Version        |
|--------------------|----------------|
| MATLAB             | R2022a or later |
| Signal Processing Toolbox | (optional, for `normalize`) |
| Statistics and Machine Learning Toolbox | (optional) |

> **Note:** The `normalize` function used in the Cv calculation requires MATLAB R2020b or later.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/TASC-DMD.git
cd TASC-DMD
```

### 2. Generate the benchmark data

In MATLAB:
```matlab
cd examples
run generate_NLS_data.m
% Saves NLS_data.mat to the examples/ folder
```

### 3. Run the TASC-DMD benchmark

```matlab
run Schrodinger_equ_TASCDMD.m
```

This reproduces the eigenvalue spectrum comparison and field reconstructions from Figure 2 of the paper.

---

## Function Reference

### `DMD_exact(X, dt, p, tde, r, var)`

Exact or standard DMD with optional time-delay embedding.

| Parameter | Description |
|-----------|-------------|
| `X`       | Data matrix `[n × m]` |
| `dt`      | Time step |
| `p`       | Temporal subsampling factor |
| `tde`     | Number of time delays (1, 2, or 3) |
| `r`       | Rank truncation |
| `var`     | `0` = exact DMD, `1` = standard DMD |

Returns: `[nt2, Phi, lambda, omega, time_dynamics, X_dmd, Cv]`

---

### `fbDMD(X, dt, p, tde, r)`

Forward-Backward DMD — reduces bias from sensor noise by combining forward and backward linear operators via their geometric mean.

Same input/output interface as `DMD_exact` (without the `var` argument).

---

### `nls_rhs(t, psit, dummy, k)`

ODE right-hand side for the NLS equation in Fourier space. Pass to `ode45`.

```matlab
[t, psit] = ode45(@(tt, pp) nls_rhs(tt, pp, [], k), t_span, psit0, opts);
```

---

## Key Parameters (TASC-DMD)

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `d`       | *d*    | Number of time-delay levels (Hankel depth) |
| `r0`      | *r₀*   | Rank of spatial basis per delay level (space-contraction) |
| `r1`      | *r₁*   | Rank for the final DMD step (mode truncation) |

For the NLS benchmark: `d = 11`, `r0 = 11`, `r1 = 9`.

---

## Citation

If you use this code, please cite:

```bibtex
@article{arani2025tascdmd,
  title   = {A Novel Spatiotemporal Decomposition and Identification of
             Sparse Equations for Human Brain Deformation},
  author  = {Arani, A.H.G. and Alshareef, A.A. and Pham, D.L. and
             Okamoto, R.J. and Bayly, P.V.},
  journal = {(journal name)},
  year    = {2025}
}
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
