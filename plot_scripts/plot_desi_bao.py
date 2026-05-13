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
g.export('results/comparison_triangle.png')

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

# --- ΛCDM Hubble and distances ---------------------------------
def H_lcdm(z, H0, Om):
    Omega_r = OMEGA_GAMMA_H2 / (H0 / 100.0) ** 2
    Omega_L = 1.0 - Om - Omega_r
    Esq = Om * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4 + Omega_L
    return H0 * np.sqrt(Esq)

def D_M_lcdm(z_val, H0, Om, npts=200):
    z_int = np.linspace(0, z_val, npts)
    H_int = H_lcdm(z_int, H0, Om)
    dz = z_int[1] - z_int[0]
    integral = np.sum(1.0 / H_int) * dz
    return C_KM_S * integral

def r_s_lcdm(H0, Om):
    from scipy import integrate
    def integrand(z):
        R = 3.0 * OMEGA_B_H2 / (4.0 * OMEGA_GAMMA_H2 * (1.0 + z))
        cs = C_KM_S / np.sqrt(3.0 * (1.0 + R))
        return cs / H_lcdm(z, H0, Om)
    result, _ = integrate.quad(integrand, Z_DRAG, np.inf, limit=200,
                               epsabs=1e-6, epsrel=1e-6)
    return result

# Best-fit parameters from USF chain
feu_stats = feu.getMargeStats()
H0_feu = feu_stats.parWithName('H0').mean
Om_feu = feu_stats.parWithName('Omega_m').mean
aq_feu = feu_stats.parWithName('alpha_q').mean
zt_feu = feu_stats.parWithName('z_trans').mean

lcdm_stats = lcdm.getMargeStats()
H0_lcdm = lcdm_stats.parWithName('H0').mean
Om_lcdm = lcdm_stats.parWithName('Omega_m').mean

# Compute rd consistently with the likelihood
# --- Sound horizons -------------------------------------------
rd_feu  = r_s_feu(H0_feu, Om_feu, aq_feu, zt_feu)
rd_lcdm = r_s_lcdm(H0_lcdm, Om_lcdm)
print(f"Computed rd for best‑fit USF: {rd_feu:.2f} Mpc")
print(f"Computed rd for best‑fit ΛCDM: {rd_lcdm:.2f} Mpc")

# Load the DESI 2024 mean file (columns: z, value, quantity)
# --- Load DESI data -------------------------------------------
bao_file = './cobaya_packages/data/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt'
data = pd.read_csv(bao_file, comment='#', sep=r'\s+',
                   names=['z', 'value', 'quantity'])
dm = data[data['quantity'] == 'DM_over_rs']
dh = data[data['quantity'] == 'DH_over_rs']

z_dm  = dm['z'].values
obs_dm = dm['value'].values
z_dh  = dh['z'].values
obs_dh = dh['value'].values

# --- Model predictions ----------------------------------------
# USF
pred_dm_feu = np.array([D_M_feu(z, H0_feu, Om_feu, aq_feu, zt_feu) / rd_feu for z in z_dm])
H_dh_feu   = H_feu(z_dh, H0_feu, Om_feu, aq_feu, zt_feu)
pred_dh_feu = C_KM_S / (rd_feu * H_dh_feu)

# ΛCDM
pred_dm_lcdm = np.array([D_M_lcdm(z, H0_lcdm, Om_lcdm) / rd_lcdm for z in z_dm])
H_dh_lcdm    = H_lcdm(z_dh, H0_lcdm, Om_lcdm)
pred_dh_lcdm = C_KM_S / (rd_lcdm * H_dh_lcdm)

# --- Residuals (percentage) -----------------------------------
res_dm_usf  = 100.0 * (obs_dm - pred_dm_feu) / obs_dm
res_dm_lcdm = 100.0 * (obs_dm - pred_dm_lcdm) / obs_dm
res_dh_usf  = 100.0 * (obs_dh - pred_dh_feu) / obs_dh
res_dh_lcdm = 100.0 * (obs_dh - pred_dh_lcdm) / obs_dh

# --- Plot (2 rows × 2 columns) --------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                         gridspec_kw={'height_ratios': [2, 1]},
                         sharex='col')
ax1, ax2 = axes[0]
ax1r, ax2r = axes[1]

# Upper‑left: D_M/r_d
ax1.plot(z_dm, obs_dm, 'ko', label='DESI')
ax1.plot(z_dm, pred_dm_feu, 'r-', linewidth=2, label='USF')
ax1.plot(z_dm, pred_dm_lcdm, 'b--', linewidth=2, label='ΛCDM')
ax1.set_ylabel('D$_M$/r$_d$')
ax1.set_title('Transverse comoving distance')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Lower‑left: D_M residuals
ax1r.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax1r.plot(z_dm, res_dm_lcdm, 'b--', linewidth=1.5, alpha=0.8, label='ΛCDM')
ax1r.plot(z_dm, res_dm_usf, 'r-', linewidth=1.5, alpha=0.8, label='USF')
ax1r.set_xlabel('z')
ax1r.set_ylabel('Δ (%)')
ax1r.legend(fontsize=8)
ax1r.grid(True, alpha=0.3)

# Upper‑right: D_H/r_d
ax2.plot(z_dh, obs_dh, 'ko', label='DESI')
ax2.plot(z_dh, pred_dh_feu, 'r-', linewidth=2, label='USF')
ax2.plot(z_dh, pred_dh_lcdm, 'b--', linewidth=2, label='ΛCDM')
ax2.set_ylabel('D$_H$/r$_d$')
ax2.set_title('Hubble distance')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Lower‑right: D_H residuals
ax2r.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax2r.plot(z_dh, res_dh_lcdm, 'b--', linewidth=1.5, alpha=0.8, label='ΛCDM')
ax2r.plot(z_dh, res_dh_usf, 'r-', linewidth=1.5, alpha=0.8, label='USF')
ax2r.set_xlabel('z')
ax2r.set_ylabel('Δ (%)')
ax2r.legend(fontsize=8)
ax2r.grid(True, alpha=0.3)

plt.suptitle('BAO observables – USF and ΛCDM vs DESI (real data)', fontsize=14)
plt.tight_layout()
plt.savefig('results/BAO_observables_vs_DESI.pdf')
plt.savefig('results/BAO_observables_vs_DESI.png', dpi=300)

# ============================================================
# 4. Parameter table (printed to terminal)
# ============================================================
print("\n" + "="*60)
print("PARAMETER TABLE – BAO (DESI) + SN (Pantheon+) + SH0ES")
print("="*60)
print(f"{'Parameter':<20} {'ΛCDM':<20} {'USF':<20}")
print("-"*60)

# Use the already existing lcdm_stats and feu_stats (do NOT redeclare them)
print(f"{'H0 [km/s/Mpc]':<20} {lcdm_stats.parWithName('H0').mean:.1f} ± {lcdm_stats.parWithName('H0').err:.1f}{'':>9} {feu_stats.parWithName('H0').mean:.1f} ± {feu_stats.parWithName('H0').err:.1f}")
print(f"{'Omega_m':<20} {lcdm_stats.parWithName('Omega_m').mean:.4f} ± {lcdm_stats.parWithName('Omega_m').err:.4f}{'':>5} {feu_stats.parWithName('Omega_m').mean:.4f} ± {feu_stats.parWithName('Omega_m').err:.4f}")
print(f"{'alpha_q':<20} {'0 (fixed)':<20} {feu_stats.parWithName('alpha_q').mean:.4f} ± {feu_stats.parWithName('alpha_q').err:.4f}")
print(f"{'z_trans':<20} {'— (fixed)':<20} {feu_stats.parWithName('z_trans').mean:.2f} ± {feu_stats.parWithName('z_trans').err:.2f}")

def read_best_fit_chi2(chain_prefix):
    data = np.loadtxt(f"{chain_prefix}.1.txt")
    min_logpost = np.min(data[:, 1])
    return 2.0 * min_logpost

chi2_feu  = read_best_fit_chi2("./chains/feu_bao_sn_sh0es")
chi2_lcdm = read_best_fit_chi2("./chains/lcdm_bao_sn_sh0es")

print("-"*60)
print(f"{'χ² (best‑fit)':<20} {chi2_lcdm:.2f}           {chi2_feu:.2f}")

plt.show()
print("\nAll figures saved to results/.")