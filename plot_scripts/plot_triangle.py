"""
Triangle (corner) plot for the USF model parameter constraints.

Loads the USF and ΛCDM MCMC chains and generates:
    - an overlaid triangle (corner) plot of the free parameters,
    - a marginalised posterior histogram for H0,
    - printed parameter means and 1σ uncertainties.

Usage (from the repository root):
    python plot_scripts/plot_triangle.py \
        [--usf-chain chains/usf_bao_sn] \
        [--lcdm-chain chains/lcdm_bao_sn] \
        [--outdir results/figures]

Output (saved to --outdir, default results/figures/):
    Overlaid_triangle_USF_vs_LCDM.pdf/.png – corner plot
    USF_H0_posterior.pdf/.png              – H0 marginal posterior

Author: Pantaleon Systems
Date: 2026
"""

import argparse
import os

import matplotlib.pyplot as plt
from getdist import plots
from getdist.mcsamples import loadMCSamples

parser = argparse.ArgumentParser()
parser.add_argument('--usf-chain', default='./chains/usf_bao_sn')
parser.add_argument('--lcdm-chain', default='./chains/lcdm_bao_sn')
parser.add_argument('--outdir', default='results/figures',
                    help='destino das figuras (use um subdiretório por modelo/grade)')
args = parser.parse_args()

OUTDIR = args.outdir
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# 1. Load both chains (with burn-in)
# ============================================================
settings = {'ignore_rows': 0.3}   # discard first 30%

samples_usf = loadMCSamples(args.usf_chain, settings=settings)
samples_lcdm = loadMCSamples(args.lcdm_chain, settings=settings)

samples_usf.name = 'USF'
samples_lcdm.name = 'ΛCDM'

for s in (samples_usf, samples_lcdm):
    s.updateSettings({'fine_bins': 2048, 'fine_bins_2D': 1024})

# ============================================================
# 2. LaTeX parameter labels
# ============================================================
labels = {
    'H0': 'H_0',
    'Omega_m': '\\Omega_m',
    'omega_b': '\\omega_b',
    'alpha_q': '\\alpha_q',
    'z_trans': 'z_{\\rm trans}',
}

for name, label in labels.items():
    for s in (samples_usf, samples_lcdm):
        param = s.getParamNames().parWithName(name)
        if param is not None:
            param.label = label

# Plot only parameters actually present in the USF chain
plot_params = [p for p in ('H0', 'Omega_m', 'omega_b', 'alpha_q', 'z_trans')
               if samples_usf.getParamNames().parWithName(p) is not None]

# ============================================================
# 3. Overlaid triangle plot
# ============================================================
g = plots.get_subplot_plotter()
g.settings.legend_fontsize = 12

g.triangle_plot(
    [samples_lcdm, samples_usf],
    params=plot_params,
    filled=True,
    contour_colors=['#0072B2', '#D55E00'],    # blue ΛCDM, red USF
    legend_labels=['ΛCDM', 'USF'],
    markers={'H0': 73.04},                    # SH0ES reference line
    title_fontsize=14,
)

for ext in ('.pdf', '.png'):
    g.export(f'{OUTDIR}/Overlaid_triangle_USF_vs_LCDM{ext}',
             dpi=300 if ext == '.png' else None)

# ============================================================
# 4. Parameter statistics (printed to terminal)
# ============================================================
stats_usf = samples_usf.getMargeStats()
stats_lcdm = samples_lcdm.getMargeStats()

msg = ("USF  : H0 = {:.2f} ± {:.2f}, Ω_m = {:.3f} ± {:.3f}".format(
    stats_usf.parWithName('H0').mean, stats_usf.parWithName('H0').err,
    stats_usf.parWithName('Omega_m').mean, stats_usf.parWithName('Omega_m').err))
if stats_usf.parWithName('alpha_q') is not None:
    msg += ", α_q = {:.4f} ± {:.4f}, z_trans = {:.2f} ± {:.2f}".format(
        stats_usf.parWithName('alpha_q').mean, stats_usf.parWithName('alpha_q').err,
        stats_usf.parWithName('z_trans').mean, stats_usf.parWithName('z_trans').err)
print(msg)

print("ΛCDM : H0 = {:.2f} ± {:.2f}, Ω_m = {:.3f} ± {:.3f}".format(
    stats_lcdm.parWithName('H0').mean, stats_lcdm.parWithName('H0').err,
    stats_lcdm.parWithName('Omega_m').mean, stats_lcdm.parWithName('Omega_m').err))

# ============================================================
# 5. H0 posterior histogram (USF only)
# ============================================================
h0_usf = samples_usf.samples[:, samples_usf.index['H0']]
plt.figure(figsize=(6, 4))
plt.hist(h0_usf, bins=50, density=True, alpha=0.7, color='#D55E00')
plt.axvline(73.04, color='black', linestyle='--', linewidth=1.5, label='SH0ES')
plt.xlabel('H_0 [km/s/Mpc]')
plt.ylabel('Probability density')
plt.title('H_0 posterior – USF model')
plt.legend()
plt.tight_layout()
for ext in ('.pdf', '.png'):
    plt.savefig(f'{OUTDIR}/USF_H0_posterior{ext}', dpi=300 if ext == '.png' else None)

print(f"\nAll figures saved to {OUTDIR}/.")
