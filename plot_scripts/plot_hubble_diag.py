"""
Hubble diagram for the USF model vs. Pantheon+SH0ES data.

This script loads the USF MCMC chain and the Pantheon+SH0ES supernova
catalogue, then plots the distance moduli predicted by the best-fit USF
model against the observed values.

Usage:
    Run from the repository root:
        python plot_scripts/plot_hubble_diag.py

Required data files:
    chains/feu_bao_sn/          – USF MCMC chains
    cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat

Output (saved to results/):
    Hubble_diagram_real.pdf   – Pantheon+ distance moduli vs USF best-fit

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
# 1. Load USF chain and extract best-fit parameters
# ============================================================
feu = loadMCSamples('./chains/feu_bao_sn')
feu_stats = feu.getMargeStats()

H0_feu = feu_stats.parWithName('H0').mean
Om_feu = feu_stats.parWithName('Omega_m').mean
aq_feu = feu_stats.parWithName('alpha_q').mean
zt_feu = feu_stats.parWithName('z_trans').mean

# ============================================================
# 2. USF distance functions (consistent with Cobaya theory module)
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
    # E²(0) = 1 by construction with the flat-universe Omega_L above.
    Esq = Om * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4 + Omega_L + Oz

    return H0 * np.sqrt(Esq)

def dL_feu(z, H0, Om, aq, zt, npts=200):
    """Luminosity distance d_L(z) [Mpc] for the USF model."""
    z_int = np.linspace(0, z, npts)
    H_int = H_feu(z_int, H0, Om, aq, zt)
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
# 4. Compute USF predictions and plot
# ============================================================
mu_feu = 5 * np.log10([dL_feu(z, H0_feu, Om_feu, aq_feu, zt_feu) for z in z_sn]) + 25

plt.figure(figsize=(8, 5))
plt.errorbar(z_sn, mu_obs, yerr=mu_err, fmt='.', color='gray', alpha=0.5,
             label='Pantheon+SH0ES')
plt.plot(z_sn, mu_feu, 'r.', markersize=1, alpha=0.8, label='USF best-fit')
plt.xlabel('Redshift z')
plt.ylabel('Distance modulus μ [mag]')
plt.title('Hubble diagram – USF vs Pantheon+SH0ES')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('results/Hubble_diagram_real.pdf')
plt.show()

print("\nFigure saved to results/Hubble_diagram_real.pdf")