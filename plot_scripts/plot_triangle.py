"""
Triangle (corner) plot for the USF model parameter constraints.

This script loads the USF MCMC chain and generates:
    - a full triangle (corner) plot of the four free parameters,
    - a marginalised posterior histogram for H0,
    - printed parameter means and 1σ uncertainties.

Usage:
    Run from the repository root:
        python plot_scripts/plot_triangle.py

Required data files:
    chains/feu_bao_sn/   – USF MCMC chains

Output (saved to results/):
    Overlaid_triangle_USF_vs_LCDM.pdf – Corner plot (H0, Ω_m, α_q, z_trans)
    Overlaid_triangle_USF_vs_LCDM.png – H0 marginal posterior

Author: Pantaleon Systems
Date: 2026
"""

import os
import matplotlib.pyplot as plt
from getdist import plots
from getdist.mcsamples import loadMCSamples

os.makedirs('results', exist_ok=True)

# ============================================================
# 1. Load both chains (with burn‑in)
# ============================================================
settings = {'ignore_rows': 0.3}   # discard first 30%

samples_usf  = loadMCSamples('./chains/feu_bao_sn_sh0es', settings=settings)
samples_lcdm = loadMCSamples('./chains/lcdm_bao_sn_sh0es', settings=settings)

samples_usf.name  = 'USF'
samples_lcdm.name = 'ΛCDM'

# Increase bin resolution for smoother contours
for s in (samples_usf, samples_lcdm):
    s.updateSettings({'fine_bins': 2048, 'fine_bins_2D': 1024})

# ============================================================
# 2. Set LaTeX parameter labels (both chains share the same names)
# ============================================================
labels = {
    'H0':       'H_0',
    'Omega_m':  '\\Omega_m',
    'alpha_q':  '\\alpha_q',
    'z_trans':  'z_{\\rm trans}'
}

for name, label in labels.items():
    for s in (samples_usf, samples_lcdm):
        param = s.getParamNames().parWithName(name)
        if param is not None:
            param.label = label

# ============================================================
# 3. Overlaid triangle plot
# ============================================================
g = plots.get_subplot_plotter()
g.settings.legend_fontsize = 12

# Pass both chains; GetDist automatically plots only the common parameters
# for ΛCDM and the full set for USF.
g.triangle_plot(
    [samples_lcdm, samples_usf],
    params=['H0', 'Omega_m', 'alpha_q', 'z_trans'],
    filled=True,
    contour_colors=['#0072B2', '#D55E00'],    # blue for ΛCDM, red for USF
    legend_labels=['ΛCDM', 'USF'],
    markers={'H0': 73.04},                     # SH0ES reference line
    title_fontsize=14
)

# Export the overlaid figure
for ext in ('.pdf', '.png'):
    g.export(f'results/Overlaid_triangle_USF_vs_LCDM{ext}',
             dpi=300 if ext=='.png' else None)

# ============================================================
# 4. Parameter statistics (printed to terminal)
# ============================================================
stats_usf  = samples_usf.getMargeStats()
stats_lcdm = samples_lcdm.getMargeStats()

print("USF  : H0 = {:.2f} ± {:.2f}, Ω_m = {:.3f} ± {:.3f}, α_q = {:.4f} ± {:.4f}, z_trans = {:.2f} ± {:.2f}".format(
    stats_usf.parWithName('H0').mean, stats_usf.parWithName('H0').err,
    stats_usf.parWithName('Omega_m').mean, stats_usf.parWithName('Omega_m').err,
    stats_usf.parWithName('alpha_q').mean, stats_usf.parWithName('alpha_q').err,
    stats_usf.parWithName('z_trans').mean, stats_usf.parWithName('z_trans').err))

print("ΛCDM : H0 = {:.2f} ± {:.2f}, Ω_m = {:.3f} ± {:.3f}".format(
    stats_lcdm.parWithName('H0').mean, stats_lcdm.parWithName('H0').err,
    stats_lcdm.parWithName('Omega_m').mean, stats_lcdm.parWithName('Omega_m').err))

# ============================================================
# 5. H0 posterior histogram (USF only, kept for complement)
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
    plt.savefig(f'results/FEU_H0_posterior{ext}', dpi=300 if ext=='.png' else None)
plt.show()

print("\nAll figures saved to results/.")
