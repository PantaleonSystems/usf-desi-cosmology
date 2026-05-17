"""
Hubble diagram for the USF model vs. Pantheon+SH0ES data (with ΛCDM comparison).

This script loads the MCMC chains of the USF and ΛCDM models and the
Pantheon+SH0ES supernova catalogue, then plots the distance moduli
predicted by the best-fit models against the observed values.  A residual
panel below the main diagram quantifies the improvement of the USF model.

Usage:
    Run from the repository root:
        python plot_scripts/plot_hubble_diag.py

Required data files:
    chains/feu_bao_sn_sh0es/          – USF MCMC chains
    chains/lcdm_bao_sn_sh0es/         – ΛCDM reference chains
    cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat

Output (saved to results/):
    Hubble_diagram_residuals.pdf
    Hubble_diagram_residuals.png

Author: Pantaleon Systems
Date: 2026
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from getdist.mcsamples import loadMCSamples

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from usf_constants import C_KM_S, OMEGA_GAMMA_H2, DELTA

os.makedirs('results', exist_ok=True)

# ============================================================
# 1. Load MCMC chains and extract best-fit parameters
# ============================================================
feu = loadMCSamples('./chains/feu_bao_sn_sh0es')
lcdm = loadMCSamples('./chains/lcdm_bao_sn_sh0es')

# USF parameters
feu_stats = feu.getMargeStats()
H0_feu = feu_stats.parWithName('H0').mean
Om_feu = feu_stats.parWithName('Omega_m').mean
aq_feu = feu_stats.parWithName('alpha_q').mean
zt_feu = feu_stats.parWithName('z_trans').mean

# ΛCDM parameters
lcdm_stats = lcdm.getMargeStats()
H0_lcdm = lcdm_stats.parWithName('H0').mean
Om_lcdm = lcdm_stats.parWithName('Omega_m').mean

# ============================================================
# 2. Distance functions (USF and ΛCDM)
# ============================================================

def H_feu(z, H0, Om, aq, zt):
    """USF Hubble parameter H(z) [km/s/Mpc], consistent with feu_background."""
    Omega_r = OMEGA_GAMMA_H2 / (H0 / 100.0) ** 2

    # USF correction at z=0 (required before setting Omega_L)
    supp0 = 1.0 / (1.0 + np.exp(-zt / DELTA))
    Oz0 = aq * supp0 / (1.0 + 1.0 / zt ** 2)

    # Flat-universe condition: Ω_m + Ω_r + Ω_Λ + Ω_USF(0) = 1
    Omega_L = 1.0 - Om - Omega_r - Oz0

    supp = 1.0 / (1.0 + np.exp(np.clip((z - zt) / DELTA, -500, 500)))
    Oz = aq * supp * (1.0 + z) ** 3 / (1.0 + (1.0 + z) ** 2 / zt ** 2)
    Esq = Om * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4 + Omega_L + Oz

    return H0 * np.sqrt(Esq)

def H_lcdm(z, H0, Om):
    """ΛCDM Hubble parameter H(z) [km/s/Mpc]."""
    Omega_r = OMEGA_GAMMA_H2 / (H0 / 100.0) ** 2
    Omega_L = 1.0 - Om - Omega_r
    Esq = Om * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4 + Omega_L
    return H0 * np.sqrt(Esq)

def dL(z, H_func, H0, Om, aq=None, zt=None, npts=200):
    """Luminosity distance d_L(z) [Mpc] for a given Hubble function."""
    z_int = np.linspace(0, z, npts)
    if aq is not None and zt is not None:
        H_int = H_func(z_int, H0, Om, aq, zt)
    else:
        H_int = H_func(z_int, H0, Om)
    dz = z_int[1] - z_int[0]
    integral = np.sum(1.0 / H_int) * dz
    return (1.0 + z) * C_KM_S * integral

# ============================================================
# 3. Load Pantheon+SH0ES data
# ============================================================
pantheon_path = './cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat'
pdata = pd.read_csv(pantheon_path, comment='#', sep=r'\s+')

z_sn   = pdata['zHD'].values
mu_obs = pdata['MU_SH0ES'].values
mu_err = pdata['MU_SH0ES_ERR_DIAG'].values

print(f"Loaded Pantheon+ data: {len(z_sn)} supernovae")

# ============================================================
# 4. Compute model predictions
# ============================================================
mu_feu  = 5 * np.log10([dL(z, H_feu, H0_feu, Om_feu, aq_feu, zt_feu) for z in z_sn]) + 25
mu_lcdm = 5 * np.log10([dL(z, H_lcdm, H0_lcdm, Om_lcdm) for z in z_sn]) + 25

# Residuals (data – model)
res_feu  = mu_obs - mu_feu
res_lcdm = mu_obs - mu_lcdm

# ============================================================
# 5. Create the two-panel figure
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                               gridspec_kw={'height_ratios': [2, 1]},
                               sharex=True)

# --- Upper panel: Hubble diagram ---
ax1.errorbar(z_sn, mu_obs, yerr=mu_err, fmt='.', color='gray', alpha=0.4,
             label='Pantheon+SH0ES')
ax1.plot(z_sn, mu_feu, 'r-', markersize=1.5, alpha=0.7, label='USF best‑fit')
ax1.plot(z_sn, mu_lcdm, 'b-', markersize=0.8, alpha=0.5, label='ΛCDM best‑fit')
ax1.set_ylabel('Distance modulus μ [mag]')
ax1.set_title('Hubble diagram – USF and ΛCDM vs Pantheon+SH0ES')
ax1.legend(markerscale=1.5, fontsize=9)
ax1.grid(True, alpha=0.2)

# --- Lower panel: Residuals ---
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax2.errorbar(z_sn, res_lcdm, yerr=mu_err, fmt='.', color='blue', alpha=0.4,
             label='ΛCDM residuals')
ax2.errorbar(z_sn, res_feu, yerr=mu_err, fmt='.', color='red', alpha=0.6,
             label='USF residuals')
ax2.set_xlabel('Redshift z')
ax2.set_ylabel('Δμ [mag]')
ax2.legend(markerscale=1.5, fontsize=9)
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('results/Hubble_diagram_residuals.pdf')
plt.savefig('results/Hubble_diagram_residuals.png', dpi=300)
plt.show()

print("\nFigure saved to results/Hubble_diagram_residuals.pdf")