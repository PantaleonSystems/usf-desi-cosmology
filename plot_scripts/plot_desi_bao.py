"""
BAO observables and H0 comparison plot for the USF analysis.

Loads the MCMC chains of the USF and ΛCDM models and generates:
    - a triangle plot comparing the H0 posteriors of both models,
    - plots of the BAO observables D_M/r_d and D_H/r_d vs. DESI DR2 data,
    - a parameter table printed to terminal.

All background physics (E(z), distances, r_d with CAMB-calibrated z_drag)
comes from theory_module.background — the exact functions used by the
sampler (Topic 9 of the correction plan; no local physics reimplementation).

Usage (from the repository root):
    python plot_scripts/plot_desi_bao.py \
        [--usf-chain chains/usf_bao_sn] \
        [--lcdm-chain chains/lcdm_bao_sn]

Output (saved to --outdir, default results/figures/):
    comparison_triangle.pdf       – H0 posterior comparison (USF vs ΛCDM)
    BAO_observables_vs_DESI.pdf   – BAO observables D_M/r_d, D_H/r_d

Author: Pantaleon Systems
Date: 2026
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from getdist import plots
from getdist.mcsamples import loadMCSamples

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from theory_module import background as bg  # noqa: E402
from usf_constants import C_KM_S, OMEGA_B_H2_BBN  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument('--usf-chain', default='./chains/usf_bao_sn')
parser.add_argument('--lcdm-chain', default='./chains/lcdm_bao_sn')
parser.add_argument('--outdir', default='results/figures',
                    help='destino das figuras (use um subdiretório por modelo/grade)')
args = parser.parse_args()

OUTDIR = args.outdir
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# 1. Load MCMC chains
# ============================================================
usf = loadMCSamples(args.usf_chain)
lcdm = loadMCSamples(args.lcdm_chain)

usf.name = 'USF'
lcdm.name = 'ΛCDM'


def mean_or(stats, name, default):
    par = stats.parWithName(name)
    return par.mean if par is not None else default


usf_stats = usf.getMargeStats()
H0_usf = mean_or(usf_stats, 'H0', 70.0)
Om_usf = mean_or(usf_stats, 'Omega_m', 0.31)
aq_usf = mean_or(usf_stats, 'alpha_q', 0.0)
zt_usf = mean_or(usf_stats, 'z_trans', 2.0)
wb_usf = mean_or(usf_stats, 'omega_b', OMEGA_B_H2_BBN)

lcdm_stats = lcdm.getMargeStats()
H0_lcdm = mean_or(lcdm_stats, 'H0', 68.0)
Om_lcdm = mean_or(lcdm_stats, 'Omega_m', 0.31)
wb_lcdm = mean_or(lcdm_stats, 'omega_b', OMEGA_B_H2_BBN)

# ============================================================
# 2. Triangle plot: H0 comparison
# ============================================================
g = plots.get_subplot_plotter()
g.triangle_plot([usf, lcdm], params=['H0'], filled=True,
                legend_labels=['USF', 'ΛCDM'], title_fontsize=12)
g.export(f'{OUTDIR}/comparison_triangle.pdf')
g.export(f'{OUTDIR}/comparison_triangle.png')

# ============================================================
# 3. BAO observables D_M/r_d and D_H/r_d – DESI DR2 data
# ============================================================
rd_usf = bg.r_drag(wb_usf, Om_usf * (H0_usf / 100.0) ** 2)
rd_lcdm = bg.r_drag(wb_lcdm, Om_lcdm * (H0_lcdm / 100.0) ** 2)
print(f"Computed rd for mean USF: {rd_usf:.2f} Mpc")
print(f"Computed rd for mean ΛCDM: {rd_lcdm:.2f} Mpc")

bao_file = ('./cobaya_packages/data/bao_data/desi_bao_dr2/'
            'desi_gaussian_bao_ALL_GCcomb_mean.txt')
data = pd.read_csv(bao_file, comment='#', sep=r'\s+',
                   names=['z', 'value', 'quantity'])
dm = data[data['quantity'] == 'DM_over_rs']
dh = data[data['quantity'] == 'DH_over_rs']

z_dm, obs_dm = dm['z'].values, dm['value'].values
z_dh, obs_dh = dh['z'].values, dh['value'].values

# --- Model predictions (shared background module) ---------------
pred_dm_usf = bg.comoving_distance(z_dm, H0_usf, Om_usf, aq_usf, zt_usf) / rd_usf
pred_dh_usf = C_KM_S / (rd_usf * bg.hubble(z_dh, H0_usf, Om_usf, aq_usf, zt_usf))

pred_dm_lcdm = bg.comoving_distance(z_dm, H0_lcdm, Om_lcdm) / rd_lcdm
pred_dh_lcdm = C_KM_S / (rd_lcdm * bg.hubble(z_dh, H0_lcdm, Om_lcdm))

# --- Residuals (percentage) -----------------------------------
res_dm_usf = 100.0 * (obs_dm - pred_dm_usf) / obs_dm
res_dm_lcdm = 100.0 * (obs_dm - pred_dm_lcdm) / obs_dm
res_dh_usf = 100.0 * (obs_dh - pred_dh_usf) / obs_dh
res_dh_lcdm = 100.0 * (obs_dh - pred_dh_lcdm) / obs_dh

# --- Plot (2 rows × 2 columns) --------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                         gridspec_kw={'height_ratios': [2, 1]},
                         sharex='col')
ax1, ax2 = axes[0]
ax1r, ax2r = axes[1]

ax1.plot(z_dm, obs_dm, 'ko', label='DESI DR2')
ax1.plot(z_dm, pred_dm_usf, 'r-', linewidth=2, label='USF')
ax1.plot(z_dm, pred_dm_lcdm, 'b--', linewidth=2, label='ΛCDM')
ax1.set_ylabel('D$_M$/r$_d$')
ax1.set_title('Transverse comoving distance')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax1r.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax1r.plot(z_dm, res_dm_lcdm, 'b--', linewidth=1.5, alpha=0.8, label='ΛCDM')
ax1r.plot(z_dm, res_dm_usf, 'r-', linewidth=1.5, alpha=0.8, label='USF')
ax1r.set_xlabel('z')
ax1r.set_ylabel('Δ (%)')
ax1r.legend(fontsize=8)
ax1r.grid(True, alpha=0.3)

ax2.plot(z_dh, obs_dh, 'ko', label='DESI DR2')
ax2.plot(z_dh, pred_dh_usf, 'r-', linewidth=2, label='USF')
ax2.plot(z_dh, pred_dh_lcdm, 'b--', linewidth=2, label='ΛCDM')
ax2.set_ylabel('D$_H$/r$_d$')
ax2.set_title('Hubble distance')
ax2.legend()
ax2.grid(True, alpha=0.3)

ax2r.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax2r.plot(z_dh, res_dh_lcdm, 'b--', linewidth=1.5, alpha=0.8, label='ΛCDM')
ax2r.plot(z_dh, res_dh_usf, 'r-', linewidth=1.5, alpha=0.8, label='USF')
ax2r.set_xlabel('z')
ax2r.set_ylabel('Δ (%)')
ax2r.legend(fontsize=8)
ax2r.grid(True, alpha=0.3)

plt.suptitle('BAO observables – USF and ΛCDM vs DESI DR2', fontsize=14)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/BAO_observables_vs_DESI.pdf')
plt.savefig(f'{OUTDIR}/BAO_observables_vs_DESI.png', dpi=300)

# ============================================================
# 4. Parameter table (printed to terminal)
# ============================================================
print("\n" + "=" * 60)
print("PARAMETER TABLE")
print("=" * 60)
print(f"{'Parameter':<20} {'ΛCDM':<20} {'USF':<20}")
print("-" * 60)
print(f"{'H0 [km/s/Mpc]':<20} "
      f"{H0_lcdm:.1f} ± {lcdm_stats.parWithName('H0').err:.1f}{'':>9} "
      f"{H0_usf:.1f} ± {usf_stats.parWithName('H0').err:.1f}")
print(f"{'Omega_m':<20} "
      f"{Om_lcdm:.4f} ± {lcdm_stats.parWithName('Omega_m').err:.4f}{'':>5} "
      f"{Om_usf:.4f} ± {usf_stats.parWithName('Omega_m').err:.4f}")
if usf_stats.parWithName('alpha_q') is not None:
    print(f"{'alpha_q':<20} {'0 (fixed)':<20} "
          f"{aq_usf:.4f} ± {usf_stats.parWithName('alpha_q').err:.4f}")
    print(f"{'z_trans':<20} {'— (fixed)':<20} "
          f"{zt_usf:.2f} ± {usf_stats.parWithName('z_trans').err:.2f}")


def read_best_fit_chi2(chain_prefix):
    data = np.loadtxt(f"{chain_prefix}.1.txt")
    min_logpost = np.min(data[:, 1])
    return 2.0 * min_logpost


chi2_usf = read_best_fit_chi2(args.usf_chain)
chi2_lcdm = read_best_fit_chi2(args.lcdm_chain)

print("-" * 60)
print(f"{'χ² (best-fit)':<20} {chi2_lcdm:.2f}           {chi2_usf:.2f}")
print(f"\nAll figures saved to {OUTDIR}/.")
