"""
BAO observables and H0 comparison plot for the USF analysis.

This script loads the MCMC chains of the USF and ΛCDM models and generates:
    - a triangle plot comparing the H0 posteriors of both models,
    - plots of the BAO observables D_M/r_d and D_H/r_d vs. DESI data,
    - a parameter table printed to terminal.

Usage:
    Run from the repository root:
        python plot_scripts/plot_desi_bao.py

Required data files:
    chains/feu_bao_sn/          – USF MCMC chains
    chains/lcdm_bao_sn/         – ΛCDM reference chains
    cobaya_packages/data/bao_data/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt

Output (saved to results/):
    comparison_triangle.pdf       – H0 posterior comparison (FEU vs ΛCDM)
    BAO_observables_vs_DESI.pdf   – BAO observables D_M/r_d, D_H/r_d

Author: Pantaleon Systems
Date: 2026
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from getdist import plots
from getdist.mcsamples import loadMCSamples

# Ensure project root is on sys.path so usf_constants can be imported when
# the script is invoked as  python plot_scripts/plot_desi_bao.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from usf_constants import C_KM_S, OMEGA_B_H2, OMEGA_GAMMA_H2, Z_DRAG, DELTA

os.makedirs('results', exist_ok=True)

# ============================================================
# 1. Load MCMC chains
# ============================================================
feu = loadMCSamples('./chains/feu_bao_sn_sh0es')
lcdm = loadMCSamples('./chains/lcdm_bao_sn_sh0es')

feu.name = 'FEU'
lcdm.name = 'ΛCDM'

# ============================================================
# 2. Triangle plot: H0 comparison
# ============================================================
g = plots.get_subplot_plotter()
g.triangle_plot(
    [feu, lcdm],
    params=['H0'],
    filled=True,
    legend_labels=['FEU', 'ΛCDM'],
    title_fontsize=12
)
g.export('results/comparison_triangle.pdf')

# ============================================================
# 3. BAO observables D_M/r_d and D_H/r_d – real DESI data
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

def D_M_feu(z_val, H0, Om, aq, zt, npts=200):
    """Comoving distance D_M(z) [Mpc] for the USF model."""
    z_int = np.linspace(0, z_val, npts)
    H_int = H_feu(z_int, H0, Om, aq, zt)
    dz = z_int[1] - z_int[0]
    integral = np.sum(1.0 / H_int) * dz
    return C_KM_S * integral

def r_s_feu(H0, Om, aq, zt):
    """Sound horizon at the drag epoch (Mpc), consistent with bao_desi_likelihood."""
    from scipy import integrate

    def integrand(z):
        R = 3.0 * OMEGA_B_H2 / (4.0 * OMEGA_GAMMA_H2 * (1.0 + z))
        cs = C_KM_S / np.sqrt(3.0 * (1.0 + R))
        return cs / H_feu(z, H0, Om, aq, zt)

    result, _ = integrate.quad(integrand, Z_DRAG, np.inf, limit=200,
                               epsabs=1e-6, epsrel=1e-6)
    return result

# Best-fit parameters from USF chain
feu_stats = feu.getMargeStats()
H0_feu = feu_stats.parWithName('H0').mean
Om_feu = feu_stats.parWithName('Omega_m').mean
aq_feu = feu_stats.parWithName('alpha_q').mean
zt_feu = feu_stats.parWithName('z_trans').mean

# Compute rd consistently with the likelihood
rd_feu = r_s_feu(H0_feu, Om_feu, aq_feu, zt_feu)
print(f"Computed rd for best-fit USF: {rd_feu:.2f} Mpc")

# Load the DESI 2024 mean file (columns: z, value, quantity)
bao_file = './cobaya_packages/data/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt'
data = pd.read_csv(bao_file, comment='#', sep=r'\s+',
                   names=['z', 'value', 'quantity'])

dm = data[data['quantity'] == 'DM_over_rs']
dh = data[data['quantity'] == 'DH_over_rs']

# USF predictions for the BAO observables, using self-consistent rd
z_dm = dm['z'].values
obs_dm = dm['value'].values
pred_dm = np.array([D_M_feu(z, H0_feu, Om_feu, aq_feu, zt_feu) / rd_feu for z in z_dm])

z_dh = dh['z'].values
obs_dh = dh['value'].values
H_at_dh = H_feu(z_dh, H0_feu, Om_feu, aq_feu, zt_feu)
pred_dh = C_KM_S / (rd_feu * H_at_dh)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(z_dm, obs_dm, 'bo', label='DESI (D$_M$/r$_d$)')
ax1.plot(z_dm, pred_dm, 'r-', linewidth=2, label='USF best-fit')
ax1.set_xlabel('z')
ax1.set_ylabel('D$_M$/r$_d$')
ax1.set_title('Transverse comoving distance')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(z_dh, obs_dh, 'bo', label='DESI (D$_H$/r$_d$)')
ax2.plot(z_dh, pred_dh, 'r-', linewidth=2, label='USF best-fit')
ax2.set_xlabel('z')
ax2.set_ylabel('D$_H$/r$_d$')
ax2.set_title('Hubble distance')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.suptitle('BAO observables – USF vs DESI (real data)', fontsize=14)
plt.tight_layout()
plt.savefig('results/BAO_observables_vs_DESI.pdf')

# ============================================================
# 4. Parameter table (printed to terminal)
# ============================================================
print("\n" + "="*60)
print("PARAMETER TABLE – BAO (DESI) + SN (Pantheon+) + SH0ES")
print("="*60)
print(f"{'Parameter':<20} {'ΛCDM':<20} {'USF':<20}")
print("-"*60)

lcdm_stats = lcdm.getMargeStats()
lcdm_H0 = lcdm_stats.parWithName('H0')
lcdm_Om = lcdm_stats.parWithName('Omega_m')

feu_H0 = feu_stats.parWithName('H0')
feu_Om = feu_stats.parWithName('Omega_m')
feu_aq = feu_stats.parWithName('alpha_q')
feu_zt = feu_stats.parWithName('z_trans')

print(f"{'H0 [km/s/Mpc]':<20} {lcdm_H0.mean:.1f} ± {lcdm_H0.err:.1f}{'':>9} {feu_H0.mean:.1f} ± {feu_H0.err:.1f}")
print(f"{'Omega_m':<20} {lcdm_Om.mean:.4f} ± {lcdm_Om.err:.4f}{'':>5} {feu_Om.mean:.4f} ± {feu_Om.err:.4f}")
print(f"{'alpha_q':<20} {'0 (fixed)':<20} {feu_aq.mean:.4f} ± {feu_aq.err:.4f}")
print(f"{'z_trans':<20} {'— (fixed)':<20} {feu_zt.mean:.2f} ± {feu_zt.err:.2f}")

def read_best_fit_chi2(chain_prefix):
    """
    Returns the best-fit chi2 = 2 * min(minuslogpost) from a Cobaya chain.
    Assumes the chain file <chain_prefix>.1.txt with columns:
    weight  minuslogpost  H0  Omega_m  alpha_q  z_trans  chi2__SN ...
    """
    data = np.loadtxt(f"{chain_prefix}.1.txt")
    min_logpost = np.min(data[:, 1])  # column 1 = minuslogpost
    return 2.0 * min_logpost

chi2_feu  = read_best_fit_chi2("./chains/feu_bao_sn_sh0es")
chi2_lcdm = read_best_fit_chi2("./chains/lcdm_bao_sn_sh0es")

print("-"*60)
print(f"{'χ² (best‑fit)':<20} {chi2_lcdm:.2f}           {chi2_feu:.2f}")

plt.show()
print("\nAll figures saved to results/.")