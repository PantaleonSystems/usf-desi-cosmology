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
    FEU_triangle_bao_sn_sh0es.pdf   – Corner plot (H0, Ω_m, α_q, z_trans)
    FEU_H0_posterior.pdf            – H0 marginal posterior

Author: Pantaleon Systems
Date: 2026
"""

import os
import matplotlib.pyplot as plt
from getdist import plots
from getdist.mcsamples import loadMCSamples

os.makedirs('results', exist_ok=True)

# ============================================================
# 1. Load USF chain
# ============================================================
feu_samples = loadMCSamples('./chains/feu_bao_sn_sh0es')
feu_samples.name = 'FEU'

# Increase bin resolution for smoother contours.
feu_samples.updateSettings({'fine_bins': 2048, 'fine_bins_2D': 1024})

# ============================================================
# 2. Set LaTeX parameter labels
# ============================================================
labels = {
    'H0':       'H_0',
    'Omega_m':  '\\Omega_m',
    'alpha_q':  '\\alpha_q',
    'z_trans':  'z_{\\rm trans}'
}

for name, label in labels.items():
    param = feu_samples.getParamNames().parWithName(name)
    if param is not None:
        param.label = label

# ============================================================
# 3. Triangle plot
# ============================================================
g = plots.get_subplot_plotter()
g.triangle_plot(
    feu_samples,
    params=['H0', 'Omega_m', 'alpha_q', 'z_trans'],
    filled=True,
    title_fontsize=14
)
g.export('results/FEU_triangle_bao_sn_sh0es.pdf')

# ============================================================
# 4. Parameter statistics (printed to terminal)
# ============================================================
stats = feu_samples.getMargeStats()

print(f"H0       = {stats.parWithName('H0').mean:.2f} ± {stats.parWithName('H0').err:.2f} km/s/Mpc")
print(f"Omega_m  = {stats.parWithName('Omega_m').mean:.3f} ± {stats.parWithName('Omega_m').err:.3f}")
print(f"alpha_q  = {stats.parWithName('alpha_q').mean:.4f} ± {stats.parWithName('alpha_q').err:.4f}")
print(f"z_trans  = {stats.parWithName('z_trans').mean:.2f} ± {stats.parWithName('z_trans').err:.2f}")

# ============================================================
# 5. H0 posterior histogram
# ============================================================
h0_samples = feu_samples.samples[:, feu_samples.index['H0']]
plt.figure(figsize=(6, 4))
plt.hist(h0_samples, bins=50, density=True, alpha=0.7, color='red')
plt.xlabel('H_0 [km/s/Mpc]')
plt.ylabel('Probability density')
plt.title('H_0 posterior – USF model')
plt.savefig('results/FEU_H0_posterior.pdf')
plt.show()
