# Unified State Function — MCMC Cosmological Analysis

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)
[![Cobaya 3.6](https://img.shields.io/badge/cobaya-3.6.2-orange.svg)](https://cobaya.readthedocs.io/)
[![GetDist 1.7](https://img.shields.io/badge/getdist-1.7.7-green.svg)](https://getdist.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: Ready for Publication](https://img.shields.io/badge/status-ready%20for%20publication-brightgreen.svg)]()
[![Buchalter Prize 2026](https://img.shields.io/badge/submission-Buchalter%20Cosmology%20Prize%202026-purple.svg)]()

> **Public code, MCMC chains, and Cobaya configuration files for the paper:**  
> *"Unified State Function: A Quantum-Geometric Framework for Gravitation, Cosmology, and Particle Physics"*  
> Submitted to the **Buchalter Cosmology Prize 2026**.

This repository provides full reproducibility of the Bayesian MCMC analysis comparing the Standard Model (ΛCDM) against the Unified State Function (USF) framework, using a joint dataset of DESI DR2 BAO, Pantheon+, and SH0ES measurements.

---

## Table of Contents

1. [Physics Background](#1-physics-background)
2. [Key Results](#2-key-results)
3. [Datasets](#3-datasets)
4. [Repository Structure](#4-repository-structure)
5. [Installation](#5-installation)
6. [Data Setup](#6-data-setup)
7. [How to Run](#7-how-to-run)
8. [Output Files](#8-output-files)
9. [Authors & Contact](#9-authors--contact)
10. [Citation](#10-citation)
11. [License](#11-license)

---

## 1. Physics Background

### The Hubble Tension

The **Hubble Tension** is one of the most significant open problems in modern cosmology: a persistent ~5σ discrepancy between two independent measurements of the present-day expansion rate of the universe:

| Probe | H₀ (km/s/Mpc) | Method |
|---|---|---|
| CMB (Planck 2018, ΛCDM) | 67.4 ± 0.5 | Early-universe inference |
| SH0ES (Riess et al. 2023) | 73.0 ± 1.04 | Local distance ladder |

Standard ΛCDM cannot reconcile these values without invoking modifications that are either ad hoc or poorly motivated theoretically.

### The Unified State Function Framework

The USF framework introduces a **quantum-geometric correction** to the Friedmann equation, motivated by the low-energy phenomenology of Loop Quantum Gravity and String Theory. Rather than introducing a new dark energy component, the USF posits that the large-scale geometry of spacetime carries an intrinsic, redshift-dependent correction Ω_USF(z) that becomes significant near a characteristic transition epoch.

The modified dimensionless Hubble parameter takes the form:

$$E^2(z) = \Omega_m (1+z)^3 + \Omega_r (1+z)^4 + \Omega_\Lambda + \Omega_\text{USF}(z)$$

where the geometric correction is:

$$\Omega_\text{USF}(z) = \frac{\alpha_q \cdot \sigma(z;\, z_\text{trans},\, \delta) \cdot (1+z)^3}{1 + (1+z)^2 / z_\text{trans}^2}$$

and $\sigma$ is a logistic suppression function centred on the transition redshift $z_\text{trans}$ with fixed width $\delta = 0.25$. The flat-universe closure condition $\Omega_m + \Omega_r + \Omega_\Lambda + \Omega_\text{USF}(0) = 1$ is enforced exactly, guaranteeing $E(0) = 1$ by construction.

### Model Parameters

| Parameter | Symbol | Description |
|---|---|---|
| Hubble constant | H₀ | Present expansion rate (km/s/Mpc) |
| Matter density | Ω_m | Present-day matter fraction |
| Geometric amplitude | α_q | Amplitude of the quantum-geometric correction |
| Transition redshift | z_trans | Redshift at the onset of the geometric correction |

The ΛCDM control model is recovered exactly by fixing α_q = 0.

---

## 2. Key Results

| Observable | USF | ΛCDM |
|---|---|---|
| H₀ (km/s/Mpc) | **72.8 ± 0.9** | 68.3 ± 0.6 |
| Ω_m | 0.305 ± 0.012 | 0.312 ± 0.010 |
| α_q | **0.015 ± 0.003** | — (fixed to 0) |
| z_trans | 1.8 ± 0.3 | — (fixed) |

**Statistical comparison:**

| Metric | Value | Interpretation |
|---|---|---|
| ΔAIC (USF − ΛCDM) | **−3.0** | Positive evidence in favour of USF |
| α_q significance | **> 5σ** | Geometric correction detected |
| Hubble Tension | **Resolved** | USF H₀ consistent with SH0ES |

The USF model resolves the Hubble Tension without introducing dark energy as an ad hoc field, instead attributing the effective H₀ shift to the quantum-geometric correction becoming active at z ≈ 1.8.

---

## 3. Datasets

The MCMC chains are conditioned on three complementary observational datasets:

### DESI DR2 — Baryon Acoustic Oscillations

The Dark Energy Spectroscopic Instrument Data Release 2 provides BAO measurements across seven effective redshift bins. The likelihood evaluates three geometric observables:

- **D_M / r_d** — transverse comoving distance ratio
- **D_H / r_d** — Hubble distance ratio  
- **D_V / r_d** — angle-averaged distance ratio

The full covariance matrix is read directly from the official DESI release files (`desi_gaussian_bao_ALL_GCcomb_mean.txt` and `desi_gaussian_bao_ALL_GCcomb_cov.txt`).

### Pantheon+ — Type Ia Supernovae

The Pantheon+ compilation of 1701 Type Ia supernova light curves (Scolnic et al. 2022) is handled via Cobaya's native `sn.pantheonplus` likelihood. Native marginalisation over the supernova absolute magnitude M_B is employed, avoiding the need to sample M_B explicitly and reducing the parameter space to the four USF cosmological parameters.

> **Note:** The full 1701 × 1701 covariance matrix inversion performed by the Pantheon+ likelihood is the dominant computational cost of the pipeline.

### SH0ES — Local H₀ Prior

A custom Gaussian likelihood encodes the SH0ES 2023 local distance-ladder constraint:

$$\ln \mathcal{L}_\text{SH0ES} = -\frac{1}{2} \left(\frac{H_0 - 73.0}{1.04}\right)^2$$

This is implemented as a standalone Cobaya likelihood class (`SH0ESLikelihood`) rather than a hard prior, allowing it to be toggled independently in future analyses.

---

## 4. Repository Structure

```
usf-desi-cosmology/
│
├── cobaya_configs/             # Cobaya YAML run configurations
│   ├── feu_bao_sn_sh0es.yaml   # USF model: DESI DR2 + Pantheon+ + SH0ES
│   └── lcdm_bao_sn_sh0es.yaml  # ΛCDM control: same data, α_q = 0 fixed
│
├── theory_module/              # Custom Cobaya theory and likelihood classes
│   ├── feu_background.py       # USF theory: modified Friedmann E(z), distances, r_s
│   ├── sh0es_likelihood.py     # Gaussian H₀ prior from SH0ES 2023
│   └── bao_desi_likelihood.py  # DESI DR2 BAO geometry observables
│
├── plot_scripts/               # GetDist visualisation scripts
│   ├── plot_triangle.py        # Full triangle plot + H₀ posterior comparison
│   ├── plot_desi_bao.py        # BAO distance ratios vs. DESI data points
│   └── plot_hubble_diag.py     # H(z) Hubble diagram with residuals
│
├── chains/                     # Output directory for Cobaya MCMC files
│   └── (generated at runtime — excluded from version control)
│
├── results/                    # Output directory for generated PDF figures
│   └── (generated at runtime)
│
├── cobaya_packages/            # Cobaya external data packages
│   └── data/bao_data/desi_bao_dr2/   # DESI DR2 mean and covariance files
│
├── usf_constants.py            # Shared physical constants (c, Ω_b h², Ω_γ h², z_drag, δ)
├── run_pipeline.sh             # Convenience script: clean → MCMC → plots
├── requirements.txt            # Pinned Python dependencies
└── LICENSE                     # MIT License
```

### Module Descriptions

**`usf_constants.py`** — Single source of truth for all fixed numerical values used across the theory module, likelihoods, and plot scripts. Importing from this file ensures that the MCMC engine and the visualisation layer use identical constants.

**`theory_module/feu_background.py`** — Implements the `FEU` Cobaya theory class. Provides `get_H(z)`, `get_comoving_distance(z)`, `get_angular_diameter_distance(z)`, and `get_r_s()`. Enforces the flat-universe closure condition at every MCMC step.

**`theory_module/bao_desi_likelihood.py`** — Reads the official DESI DR2 data files at initialisation, precomputes the inverse covariance matrix, and evaluates the multivariate Gaussian log-likelihood at each sampler step.

**`theory_module/sh0es_likelihood.py`** — A minimal Cobaya likelihood that penalises H₀ values away from the SH0ES central value with a Gaussian weight.

**`cobaya_configs/lcdm_bao_sn_sh0es.yaml`** — The ΛCDM control run uses the same `FEU` theory class as USF, but with `alpha_q` fixed to zero. This guarantees identical numerical integration paths and eliminates systematic differences between the two model evaluations.

---

## 5. Installation

### Prerequisites

- Python 3.12
- A Unix-like shell (macOS or Linux)

### Setup

```bash
# Clone the repository
git clone https://github.com/PantaleonSystems/usf-desi-cosmology.git
cd usf-desi-cosmology

# Create and activate a virtual environment (recommended)
python3.12 -m venv .venv
source .venv/bin/activate

# Install all pinned dependencies
pip install -r requirements.txt
```

The `requirements.txt` pins all dependency versions for full reproducibility. The key packages are:

| Package | Version | Role |
|---|---|---|
| `cobaya` | 3.6.2 | MCMC sampler and theory framework |
| `getdist` | 1.7.7 | Posterior analysis and plotting |
| `numpy` | 2.4.4 | Numerical arrays |
| `scipy` | 1.17.1 | Numerical integration (`quad`) |
| `matplotlib` | 3.10.9 | Figure rendering |

---

## 6. Data Setup

> **This step is mandatory.** The `cobaya_packages/` directory is committed to the repository as an empty scaffold — its `data/` and `code/` subdirectories are excluded from version control (see `cobaya_packages/.gitignore`). After cloning, you must populate them before any MCMC run can execute.

The pipeline requires two external datasets: the **Pantheon+** supernova compilation and the **DESI DR2** BAO measurements.

---

### 6.1 — Pantheon+ (automatic via `cobaya-install`)

Cobaya's built-in `sn.pantheonplus` likelihood knows how to fetch its own data. Run the following command once from the repository root:

```bash
cobaya-install cobaya_configs/feu_bao_sn_sh0es.yaml --packages-path ./cobaya_packages
```

Cobaya will resolve all external requirements declared in the YAML file — including the Pantheon+SH0ES catalogue and its full 1701 × 1701 covariance matrix — and download them into `cobaya_packages/data/sn_data/PantheonPlus/`. An internet connection is required.

Expected output inside `cobaya_packages/data/sn_data/PantheonPlus/` after installation:

```
Pantheon+SH0ES.dat
Pantheon+SH0ES_STAT+SYS.cov
config.dataset
```

> **Alternatively**, `cobaya-run` will also trigger the download automatically on its first invocation if the data is absent and Cobaya prompts you to install missing packages. Either workflow is valid.

---

### 6.2 — DESI DR2 BAO (automatic via `cobaya-install`)

The DESI DR2 BAO data files are distributed through the official [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data) repository on GitHub and can be installed automatically by Cobaya. Run the following command once from the repository root:

```bash
cobaya-install bao.desi_2024_bao_all -p ./cobaya_packages
```

Cobaya will fetch the mean vector and covariance matrix files from the `CobayaSampler/bao_data` registry and place them in `cobaya_packages/data/bao_data/`. An internet connection is required.

**Files installed** (mean and covariance for each tracer bin):

| File | Description |
|---|---|
| `desi_gaussian_bao_ALL_GCcomb_mean.txt` | Combined BAO mean vector (all tracers) |
| `desi_gaussian_bao_ALL_GCcomb_cov.txt` | Full covariance matrix |
| `desi_gaussian_bao_BGS_BRIGHT-21.35_GCcomb_*.txt` | BGS tracer (z ≈ 0.30) |
| `desi_gaussian_bao_LRG_GCcomb_z0.4-0.6_*.txt` | LRG bin 1 (z ≈ 0.51) |
| `desi_gaussian_bao_LRG_GCcomb_z0.6-0.8_*.txt` | LRG bin 2 (z ≈ 0.71) |
| `desi_gaussian_bao_LRG+ELG_LOPnotqso_GCcomb_*.txt` | LRG+ELG overlap bin |
| `desi_gaussian_bao_ELG_LOPnotqso_GCcomb_z1.1-1.6_*.txt` | ELG tracer |
| `desi_gaussian_bao_QSO_GCcomb_*.txt` | QSO tracer (z ≈ 1.49) |
| `desi_gaussian_bao_Lya_GCcomb_*.txt` | Lyα forest (z ≈ 2.33) |

The `BAODESILikelihood` class reads only `desi_gaussian_bao_ALL_GCcomb_mean.txt` and `desi_gaussian_bao_ALL_GCcomb_cov.txt` at runtime; the per-tracer files are available for reference and cross-checks.

---

### 6.3 — Verify the layout

After completing both steps, your `cobaya_packages/` directory should look like:

```
cobaya_packages/
├── data/
│   ├── bao_data/
│   │   └── desi_bao_dr2/
│   │       ├── desi_gaussian_bao_ALL_GCcomb_mean.txt   ← required
│   │       ├── desi_gaussian_bao_ALL_GCcomb_cov.txt    ← required
│   │       └── (per-tracer files...)
│   └── sn_data/
│       └── PantheonPlus/
│           ├── Pantheon+SH0ES.dat                      ← required
│           ├── Pantheon+SH0ES_STAT+SYS.cov             ← required
│           └── config.dataset
└── code/
    └── (Cobaya external code packages, if any)
```

You are now ready to run the pipeline.

---

## 7. How to Run

> All commands should be run from the repository root with the virtual environment active.

### Option A — Automated Pipeline (Recommended)

The `run_pipeline.sh` script orchestrates the full workflow: clean old chains, run both MCMC jobs in sequence, and generate all figures.

```bash
# Full pipeline: clean + sample + plot
bash run_pipeline.sh

# Flags for partial execution:
bash run_pipeline.sh --clean        # Only delete old chains and results, then exit
bash run_pipeline.sh --sample-only  # Run MCMC only (skip plots)
bash run_pipeline.sh --plot-only    # Generate plots from existing chains (skip MCMC)
```

### Option B — Manual Step-by-Step

#### Step 1 — Clear Previous Chains

```bash
rm -rf chains/*
```

This step is mandatory before restarting a run. Cobaya's Metropolis–Hastings sampler detects existing chain files and resumes from them; stale or aborted chains from a previous session will corrupt the new run.

#### Step 2 — Run the USF Model

```bash
cobaya-run cobaya_configs/feu_bao_sn_sh0es.yaml
```

The sampler runs until the Gelman–Rubin convergence criterion R−1 < 0.01 is satisfied across all four free parameters. The dominant cost is the per-step inversion of the Pantheon+ 1701 × 1701 covariance matrix. **Expected wall time: 2–6 hours** on a modern laptop CPU.

Chain files are written incrementally to `chains/feu_bao_sn_sh0es.*`.

#### Step 3 — Run the ΛCDM Control Model

```bash
cobaya-run cobaya_configs/lcdm_bao_sn_sh0es.yaml
```

This run samples only H₀ and Ω_m (α_q and z_trans are fixed), so convergence is faster. Chain files are written to `chains/lcdm_bao_sn_sh0es.*`.

#### Step 4 — Generate Plots

```bash
# Full triangle plot with H₀ posterior comparison
python plot_scripts/plot_triangle.py

# BAO distance ratios overlaid on DESI DR2 data points
python plot_scripts/plot_desi_bao.py

# H(z) Hubble diagram with model curves and residuals
python plot_scripts/plot_hubble_diag.py
```

All figures are saved as PDF files in the `results/` directory.

### Convergence Diagnostics

Cobaya writes a `.log` file alongside each chain set. The final line reports the Gelman–Rubin R−1 values for each parameter. A run is considered converged when all R−1 values fall below 0.01. Do not use chains from unconverged runs for scientific conclusions.

---

## 8. Output Files

After a successful run, the following files are produced:

```
chains/
├── feu_bao_sn_sh0es.1.txt      # USF MCMC samples (chain 1)
├── feu_bao_sn_sh0es.2.txt      # USF MCMC samples (chain 2)
├── feu_bao_sn_sh0es.input.yaml # Copy of the input configuration
├── feu_bao_sn_sh0es.log        # Convergence diagnostics log
├── lcdm_bao_sn_sh0es.1.txt     # ΛCDM MCMC samples
└── ...

results/
├── triangle_usf_lcdm.pdf       # Parameter triangle plot
├── bao_observables.pdf         # BAO distance ratios
└── hubble_diagram.pdf          # H(z) diagram
```

---

## 9. Authors & Contact

| Name | Affiliation | Contact |
|---|---|---|
| **Efrain Marcelo Pulgar Pantaleon** | Graduate Program in Electrical Engineering and Computer Science, UFRN | [![GitHub](https://img.shields.io/badge/GitHub-efrainmpp1-181717?logo=github)](https://github.com/efrainmpp1) · [efrain.pulgar.110@ufrn.edu.br](mailto:efrain.pulgar.110@ufrn.edu.br) |
| **Efrain Pantaleon Matamoros** | School of Science and Technology, UFRN | [efrain.pantaleon@ufrn.br](mailto:efrain.pantaleon@ufrn.edu.br) |

For questions, reproducibility issues, or collaboration proposals, please open a [GitHub Issue](../../issues) or contact the authors directly.

---

## 10. Citation

If you use this code or the results of this analysis in your research, please cite:

```bibtex
@article{PantaleonSystems2026_USF,
  author       = {Pantaleon Systems},
  title        = {{Unified State Function: A Quantum-Geometric Framework
                   for Gravitation, Cosmology, and Particle Physics}},
  journal      = {Buchalter Cosmology Prize 2026},
  year         = {2026},
  note         = {Code and chains available at
                  \url{https://github.com/PantaleonSystems/usf-desi-cosmology}}
}
```

The analysis builds on the following foundational works:

- **Cobaya**: Torrado & Lewis (2021), *JCAP* 05, 057. [arXiv:2005.05290](https://arxiv.org/abs/2005.05290)
- **GetDist**: Lewis (2019). [arXiv:1910.13970](https://arxiv.org/abs/1910.13970)
- **DESI DR2 BAO**: DESI Collaboration (2025). [arXiv:2503.14738](https://arxiv.org/abs/2503.14738)
- **Pantheon+**: Scolnic et al. (2022), *ApJ* 938, 113. [arXiv:2112.03863](https://arxiv.org/abs/2112.03863)
- **SH0ES**: Riess et al. (2022), *ApJL* 934, L7. [arXiv:2112.04510](https://arxiv.org/abs/2112.04510)
- **Planck 2018**: Planck Collaboration (2020), *A&A* 641, A6. [arXiv:1807.06209](https://arxiv.org/abs/1807.06209)

---

## 11. License

This project is released under the [MIT License](LICENSE). You are free to use, modify, and distribute the code, provided that appropriate credit is given to the original authors.

---

*This repository ensures full reproducibility of the MCMC analysis in accordance with the open-science standards of the astronomical community. All parameter priors, convergence criteria, dataset versions, and software dependencies are pinned and documented here.*
