# Unified State Function — MCMC Cosmological Analysis

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)
[![Cobaya 3.6](https://img.shields.io/badge/cobaya-3.6.2-orange.svg)](https://cobaya.readthedocs.io/)
[![GetDist 1.7](https://img.shields.io/badge/getdist-1.7.7-green.svg)](https://getdist.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Public code, Cobaya configuration files and validation suite for the
> Bayesian analysis of the Unified State Function (USF) background
> parametrization against DESI DR2 BAO, Pantheon+ SNe Ia, Planck 2018
> distance priors and the SH0ES local H₀ measurement.**

> ⚠️ **Status (July 2026): pipeline corrected after referee review.**
> The first submission's numerical results were affected by implementation
> bugs in the sound-horizon sector (r_d inflated by +4.4%, which biased the
> ΛCDM reference to H₀ ≈ 64.2 km/s/Mpc). The pipeline has since been
> corrected, unified and validated against CAMB/astropy (see
> `PROGRESSO_CODIGO.md` and `pdfs/PLANO_CORRECAO_REFEREE1.md`). **All
> headline numbers from the earlier manuscript are superseded**. The chains
> and figures produced by the buggy pipeline have been **deleted** from the
> repository; everything under `chains/` and `results/` is now output of the
> corrected pipeline only.

---

## 1. The model

The USF framework is tested here as an **effective background
parametrization**: a redshift-dependent geometric correction added to the
flat ΛCDM Friedmann equation,

$$E^2(z) = \Omega_m (1+z)^3 + \Omega_r (1+z)^4 + \Omega_\Lambda + \Omega_\text{USF}(z),$$

$$\Omega_\text{USF}(z) = \frac{\alpha_q \, \sigma(z;\, z_\text{trans},\, \delta)\, (1+z)^3}{1 + (1+z)^2 / z_\text{trans}^2},$$

with a logistic suppression σ centred on z_trans (fixed width δ = 0.25) and
flatness enforced through Ω_Λ = 1 − Ω_m − Ω_r − Ω_USF(0), so E(0) = 1 by
construction. ΛCDM is recovered exactly at α_q = 0 (nested model).

Sampled parameters: **H₀, Ω_m (total matter), ω_b (BBN Gaussian prior),
α_q, z_trans** (the last two absent in the ΛCDM control). The physical
radiation density is fixed (ω_γ = 2.4728×10⁻⁵, N_eff = 3.044); Σm_ν =
0.06 eV — see `docs/NEUTRINOS.md`.

## 2. Pipeline validation (corrected, 2026-07)

| Check | Result |
| --- | --- |
| E(z), D_C, D_A vs astropy `FlatLambdaCDM` (ΛCDM limit) | < 0.05% over the full prior grid |
| r_drag vs CAMB (numerical integral + CAMB-calibrated z_drag) | **+0.006%** (was **+4.4%** before the fix) |
| D_M, D_A, H(z) vs CAMB (Σm_ν = 0.06 eV), z ≤ 5 | < 0.012% |
| DESI DR2 BAO likelihood at the published best fit | χ² = 10.3 / 13 points (was 227.5) |
| Planck 2018 distance priors at the Planck mean | χ² = 0.29 |

Run the full validation gate with:

```bash
.venv/bin/python -m pytest tests/ -q      # 44 tests
```

## 3. Datasets and likelihoods

| Likelihood | Data | Notes |
| --- | --- | --- |
| `theory_module.bao_desi_likelihood.BAODESILikelihood` | DESI **DR2** BAO (D_M/r_d, D_H/r_d, D_V/r_d; official mean + covariance) | r_d computed per sample from (ω_b, ω_m) with CAMB-calibrated z_drag; theory served via Cobaya provider |
| `sn.pantheonplus` (Cobaya built-in) | Pantheon+ (1701 SNe, STAT+SYS covariance) | Uses uncalibrated `m_b_corr`; **M marginalised analytically** → no H₀ information, no double counting (`docs/SN_CALIBRATION.md`) |
| `theory_module.cmb_distance_priors.CMBDistancePriors` | Planck 2018 (R, ℓ_A, ω_b) compressed likelihood (Chen, Huang & Wang 2019) | Theory→data mapping reproduces the compression convention exactly |
| `theory_module.sh0es_likelihood.SH0ESLikelihood` | SH0ES H₀ = 73.04 ± 1.04 (Riess et al. 2022) | Only in runs explicitly labelled `_sh0es` |
| BBN prior | ω_b = 0.02218 ± 0.00055 (Gaussian) | In the `params` block of every YAML |

## 4. Run grid (`cobaya_configs/`)

| Config | Data | Purpose |
| --- | --- | --- |
| `usf_bao_sn.yaml` / `lcdm_bao_sn.yaml` | BAO+BBN+SN | **Main runs** — no local H₀ information |
| `usf_bao.yaml`, `lcdm_bao.yaml`, `usf_sn.yaml`, `lcdm_sn.yaml` | single datasets | which dataset drives α_q |
| `usf_bao_sn_sh0es.yaml` / `lcdm_...` | + SH0ES prior | labelled variant |
| `*_cmb.yaml` (4) | + Planck distance priors | CMB geometric consistency |
| `*_polychord.yaml` (4) | nested sampling | Bayesian evidence (needs `cobaya-install polychord`) |

All configs carry fixed seeds; the USF prior on z_trans is U(0.5, 20)
(extended after the referee noted the earlier posterior was truncated at
the U(0.5, 5) boundary).

Statistical calibration of Δχ² under the ΛCDM null (z_trans profiled —
look-elsewhere effect):

```bash
.venv/bin/python analysis_scripts/run_mocks.py --n-mocks 300
```

## 5. Repository structure

```
usf-desi-cosmology/
├── cobaya_configs/                 # full run grid (see table above)
├── theory_module/
│   ├── background.py               # SINGLE SOURCE OF TRUTH: E(z), distances, r_d, z_drag/z*
│   ├── usf_background.py           # Cobaya theory class USF (provider)
│   ├── bao_desi_likelihood.py      # DESI DR2 BAO (via provider)
│   ├── cmb_distance_priors.py      # Planck 2018 compressed likelihood
│   ├── sh0es_likelihood.py         # Gaussian local H₀
│   └── data/zdrag_camb_grid.npz    # CAMB-calibrated z_drag/z* table
├── analysis_scripts/
│   ├── make_zdrag_grid.py          # regenerates the CAMB table
│   └── run_mocks.py                # null-hypothesis Δχ² calibration
├── plot_scripts/                   # all import theory_module.background
├── tests/                          # validation gate (astropy/CAMB/Cobaya)
├── docs/                           # NEUTRINOS.md, SN_CALIBRATION.md
├── usf_constants.py                # shared physical constants
├── run_pipeline.sh                 # tests → MCMC grid → figures
├── results/figures/<tag>/          # figures, one directory per run group
├── results/data/                   # numerical products (mocks, tables)
├── logs/<run>.log                  # one execution log per run (git-ignored)
├── PROGRESSO_CODIGO.md             # execution log of the correction plan
└── requirements.txt                # pinned dependencies (incl. camb, astropy, pytest)
```

## 6. Installation and data setup

```bash
git clone https://github.com/PantaleonSystems/usf-desi-cosmology.git
cd usf-desi-cosmology
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Pantheon+ catalogue + covariance (downloads into cobaya_packages/)
cobaya-install cobaya_configs/usf_bao_sn.yaml --packages-path ./cobaya_packages

# DESI DR2 BAO mean/covariance files
cobaya-install bao.desi_dr2 -p ./cobaya_packages

# (optional, for the *_polychord.yaml evidence runs)
cobaya-install polychord -p ./cobaya_packages
```

Required files after setup:
`cobaya_packages/data/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`
and
`cobaya_packages/data/sn_data/PantheonPlus/{Pantheon+SH0ES.dat, Pantheon+SH0ES_STAT+SYS.cov, config.dataset}`.

## 7. How to run

```bash
bash run_pipeline.sh                 # validation tests → main grid (BAO+SN) → figures
bash run_pipeline.sh --grid all      # full grid (main + datasets + SH0ES + CMB variants)
bash run_pipeline.sh --grid sh0es    # only the SH0ES-labelled variants
bash run_pipeline.sh --grid cmb      # only the Planck-distance-prior variants
bash run_pipeline.sh --plot-only     # regenerate figures from existing chains
bash run_pipeline.sh --clean         # wipe chains, logs and figures (keeps results/data)
bash run_pipeline.sh --clean-all     # idem, including results/data (mocks)
```

Every sampling session is preceded by the pytest validation gate (skippable
with `--skip-tests`, not recommended). Convergence: Gelman–Rubin
R−1 < 0.01 on all parameters; chains resume automatically (`--resume`).

**Testing a new model or parametrization.** Add
`cobaya_configs/<name>.yaml`, then run it without touching the standard
grids:

```bash
bash run_pipeline.sh --runs usf_newmodel_bao_sn,lcdm_bao_sn --tag newmodel --fresh
```

`--fresh` wipes that run's chain files before sampling instead of resuming.
Use it **whenever the YAML changes theory or likelihood blocks** — otherwise
Cobaya aborts with *"Old and new run information not compatible"* and the
stale chain silently stays in place. `--tag` keeps the figures in their own
`results/figures/<tag>/`, so runs never overwrite each other.

Output layout:

| Path | Content |
| --- | --- |
| `chains/<run>.*` | Cobaya samples and metadata |
| `logs/<run>.log` | one execution log per run (git-ignored) |
| `results/figures/<tag>/` | figures, one directory per run group |
| `results/data/` | numerical products (mocks, tables); survives `--clean` |

Individual runs:

```bash
cobaya-run cobaya_configs/usf_bao_sn.yaml
python plot_scripts/plot_triangle.py --usf-chain chains/usf_bao_sn \
    --lcdm-chain chains/lcdm_bao_sn --outdir results/figures/bao_sn
```

## 8. Authors & Contact

| Name | Affiliation | Contact |
| --- | --- | --- |
| **Efrain Marcelo Pulgar Pantaleon** | Graduate Program in Electrical Engineering and Computer Science, UFRN | [![GitHub](https://img.shields.io/badge/GitHub-efrainmpp1-181717?logo=github)](https://github.com/efrainmpp1) · [efrain.pulgar.110@ufrn.edu.br](mailto:efrain.pulgar.110@ufrn.edu.br) |
| **Efrain Pantaleon Matamoros** | School of Science and Technology, UFRN | [efrain.pantaleon@ufrn.br](mailto:efrain.pantaleon@ufrn.edu.br) |

For questions, reproducibility issues, or collaboration proposals, please
open a [GitHub Issue](../../issues) or contact the authors directly.

## 9. Citation

```bibtex
@misc{PantaleonSystems2026_USF,
  author = {Pulgar Pantaleon, E. M. and Pantaleon Matamoros, E.},
  title  = {{USF-desi-cosmology: MCMC pipeline for the Unified State
             Function background parametrization}},
  year   = {2026},
  note   = {\url{https://github.com/PantaleonSystems/usf-desi-cosmology}}
}
```

Foundational works: Cobaya ([arXiv:2005.05290](https://arxiv.org/abs/2005.05290)),
GetDist ([arXiv:1910.13970](https://arxiv.org/abs/1910.13970)),
DESI DR2 BAO ([arXiv:2503.14738](https://arxiv.org/abs/2503.14738)),
Pantheon+ ([arXiv:2112.03863](https://arxiv.org/abs/2112.03863)),
SH0ES ([arXiv:2112.04510](https://arxiv.org/abs/2112.04510)),
Planck 2018 ([arXiv:1807.06209](https://arxiv.org/abs/1807.06209)),
Planck distance priors ([arXiv:1808.05724](https://arxiv.org/abs/1808.05724)).

## 10. License

Released under the [MIT License](LICENSE).

---

_All priors, seeds, convergence criteria, dataset versions and software
dependencies are pinned and documented; the validation suite must pass
before any production sampling (`run_pipeline.sh` enforces this)._
