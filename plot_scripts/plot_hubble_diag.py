"""
Hubble diagram for the USF model vs. Pantheon+ data (with ΛCDM comparison).

Statistically consistent with the `sn.pantheonplus` likelihood used in the
MCMC (see docs/SN_CALIBRATION.md):

* uses the UNCALIBRATED corrected magnitudes `m_b_corr` (never `MU_SH0ES`,
  which carries the SH0ES calibration and made the original figure circular
  — Referee 1, Comment 3);
* applies the same z_HD > 0.01 cut as the likelihood (removes the Cepheid
  calibrators);
* for each model, fits the global magnitude offset M̂ (inverse-variance
  weighted mean of m_b_corr − μ_model), mimicking the analytic
  marginalisation over M performed by the likelihood.  The residual panel
  therefore shows only the *shape* information of d_L(z), which is what the
  data actually constrain — any constant offset is degenerate with M.

All background physics comes from theory_module.background (same functions
as the sampler).

Usage:
    Run from the repository root:
        python plot_scripts/plot_hubble_diag.py

Required data files:
    chains/feu_bao_sn_sh0es/          – USF MCMC chains
    chains/lcdm_bao_sn_sh0es/         – ΛCDM reference chains
    cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat

Output (saved to results/):
    Hubble_diagram_residuals.pdf / .png

Author: Pantaleon Systems
Date: 2026
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from getdist.mcsamples import loadMCSamples

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from theory_module import background as bg  # noqa: E402

os.makedirs('results', exist_ok=True)

Z_CUT = 0.01  # same mask as cobaya's sn.pantheonplus

parser = argparse.ArgumentParser()
parser.add_argument('--usf-chain', default='./chains/usf_bao_sn_sh0es')
parser.add_argument('--lcdm-chain', default='./chains/lcdm_bao_sn_sh0es')
args = parser.parse_args()

# ============================================================
# 1. Load MCMC chains and extract mean parameters
# ============================================================
feu = loadMCSamples(args.usf_chain)
lcdm = loadMCSamples(args.lcdm_chain)

feu_stats = feu.getMargeStats()
H0_feu = feu_stats.parWithName('H0').mean
Om_feu = feu_stats.parWithName('Omega_m').mean
aq_feu = feu_stats.parWithName('alpha_q').mean
zt_feu = feu_stats.parWithName('z_trans').mean

lcdm_stats = lcdm.getMargeStats()
H0_lcdm = lcdm_stats.parWithName('H0').mean
Om_lcdm = lcdm_stats.parWithName('Omega_m').mean

# ============================================================
# 2. Load Pantheon+ data (uncalibrated magnitudes, likelihood mask)
# ============================================================
pantheon_path = './cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat'
pdata = pd.read_csv(pantheon_path, comment='#', sep=r'\s+')
pdata = pdata[pdata['zHD'] > Z_CUT]

z_hd = pdata['zHD'].values
z_hel = pdata['zHEL'].values
mb_obs = pdata['m_b_corr'].values
mb_err = pdata['m_b_corr_err_DIAG'].values

print(f"Loaded Pantheon+ data: {len(z_hd)} supernovae with z_HD > {Z_CUT}")


# ============================================================
# 3. Model distance moduli and fitted magnitude offsets
# ============================================================
def mu_model(H0, Om, aq=0.0, zt=2.0):
    """Distance modulus μ(z) using d_L = (1+z_hel) · D_C(z_HD), as in the
    Cobaya SN likelihood convention."""
    dc = bg.comoving_distance(z_hd, H0, Om, aq, zt)
    dl = (1.0 + z_hel) * dc
    return 5.0 * np.log10(dl) + 25.0


def fit_offset(mb, mu, err):
    """Inverse-variance weighted global offset M̂ (mimics analytic
    marginalisation over the absolute magnitude M)."""
    w = 1.0 / err ** 2
    return np.sum(w * (mb - mu)) / np.sum(w)


mu_feu = mu_model(H0_feu, Om_feu, aq_feu, zt_feu)
mu_lcdm = mu_model(H0_lcdm, Om_lcdm)

M_feu = fit_offset(mb_obs, mu_feu, mb_err)
M_lcdm = fit_offset(mb_obs, mu_lcdm, mb_err)

print(f"Fitted offsets (m_b_corr − μ_model): "
      f"M̂_USF = {M_feu:.4f} mag, M̂_ΛCDM = {M_lcdm:.4f} mag")
print("(Any constant offset is degenerate with M — only residual *shape* "
      "is model-discriminating.)")

# Residuals after offset marginalisation
res_feu = mb_obs - mu_feu - M_feu
res_lcdm = mb_obs - mu_lcdm - M_lcdm

for name, res in [('USF', res_feu), ('ΛCDM', res_lcdm)]:
    chi2_diag = np.sum((res / mb_err) ** 2)
    print(f"{name}: diagonal chi²/N = {chi2_diag / len(res):.3f}, "
          f"weighted mean residual = "
          f"{np.average(res, weights=1 / mb_err ** 2) * 1e3:+.2f} mmag")

# ============================================================
# 4. Two-panel figure
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                               gridspec_kw={'height_ratios': [2, 1]},
                               sharex=True)

order = np.argsort(z_hd)

# --- Upper panel: Hubble diagram in offset-subtracted magnitudes ---
ax1.errorbar(z_hd, mb_obs - M_feu, yerr=mb_err, fmt='.', color='gray',
             alpha=0.35, label='Pantheon+ (m_b_corr − M̂_USF)')
ax1.plot(z_hd[order], mu_feu[order], 'r-', lw=1.4, alpha=0.9,
         label='USF best-fit')
ax1.plot(z_hd[order], mu_lcdm[order] + (M_lcdm - M_feu), 'b--', lw=1.2,
         alpha=0.9, label='ΛCDM best-fit (offset-aligned)')
ax1.set_ylabel('Distance modulus μ [mag]')
ax1.set_title('Hubble diagram – USF and ΛCDM vs Pantheon+ '
              '(M marginalised per model)')
ax1.set_xscale('log')
ax1.legend(markerscale=1.5, fontsize=9)
ax1.grid(True, alpha=0.2)

# --- Lower panel: shape-only residuals ---
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax2.errorbar(z_hd, res_lcdm, yerr=mb_err, fmt='.', color='blue', alpha=0.35,
             label='ΛCDM residuals (M̂ subtracted)')
ax2.errorbar(z_hd, res_feu, yerr=mb_err, fmt='.', color='red', alpha=0.5,
             label='USF residuals (M̂ subtracted)')
ax2.set_xlabel('Redshift z')
ax2.set_ylabel('Δm_b [mag]')
ax2.set_ylim(-1.0, 1.0)
ax2.legend(markerscale=1.5, fontsize=9)
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('results/Hubble_diagram_residuals.pdf')
plt.savefig('results/Hubble_diagram_residuals.png', dpi=300)

print("\nFigure saved to results/Hubble_diagram_residuals.pdf")
