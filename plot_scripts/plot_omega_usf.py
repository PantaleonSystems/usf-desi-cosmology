"""
Evolution of the effective geometric correction Ω_USF(z).

This script plots the best‑fit USF correction term as a function of
redshift, together with the ΛCDM asymptotic limit and the location
of the transition epoch.  The figure illustrates the three physical
regimes discussed in the paper:
    – primordial suppression (Ω_USF → 0 for z ≫ z_trans),
    – a negative dip near the transition (the “primordial brake”),
    – convergence to a small constant value at low redshift.

Usage:
    Run from the repository root:
        python plot_scripts/plot_omega_usf.py

Output (saved to results/):
    Omega_USF_evolution.pdf
    Omega_USF_evolution.png

Author: Pantaleon Systems
Date: 2026
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure project root is on sys.path so usf_constants can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from usf_constants import DELTA

os.makedirs('results', exist_ok=True)

# ============================================================
# Best‑fit parameters from the USF MCMC chain
# (mean marginalised values, as reported by plot_triangle.py)
# ============================================================
ALPHA_Q  = -0.1001
Z_TRANS  = 4.66
DELTA    = DELTA   # 0.25

# ============================================================
# USF correction function
# ============================================================
def omega_usf(z, alpha_q=ALPHA_Q, z_trans=Z_TRANS, delta=DELTA):
    """
    Effective geometric correction Ω_USF(z).

    :param z: redshift (float or array)
    :param alpha_q: amplitude of the correction
    :param z_trans: transition redshift
    :param delta: transition width
    :returns: Ω_USF(z)
    """
    term1 = (1.0 + z)**3 / (1.0 + (1.0 + z)**2 / z_trans**2)
    # clip to avoid overflow in the exponential
    term2 = 1.0 / (1.0 + np.exp(np.clip((z - z_trans) / delta, -500, 500)))
    return alpha_q * term1 * term2

# ============================================================
# Redshift grid and model curve
# ============================================================
z_vals = np.linspace(0, 10, 1000)
omega_vals = omega_usf(z_vals)

# ============================================================
# Plot
# ============================================================
plt.figure(figsize=(8, 5))

# USF curve
plt.plot(z_vals, omega_vals, color='red', linewidth=2, label='USF best‑fit')

# ΛCDM asymptotic limit (Ω = 0)
plt.axhline(0, color='gray', linestyle='--', linewidth=1.5,
            label=r'$\Lambda$CDM (asymptotic limit)')

# Transition epoch marker
plt.axvline(Z_TRANS, color='blue', linestyle=':', linewidth=1.5, alpha=0.6,
            label=rf'Transition epoch ($z_{{\rm trans}} = {Z_TRANS}$)')

# Shaded region to highlight the transition zone
z_fill = np.linspace(Z_TRANS - 0.5, Z_TRANS + 0.5, 200)
omega_fill = omega_usf(z_fill)
plt.fill_between(z_fill, omega_fill, alpha=0.1, color='blue')

# Labels and title
plt.xlabel(r'Redshift $z$', fontsize=14)
plt.ylabel(r'Effective geometric correction $\Omega_{\rm USF}(z)$', fontsize=14)
plt.title('Evolution of the USF geometric correction', fontsize=16)
plt.legend(loc='lower right', fontsize=12)

# Grid and fine‑tuning
plt.grid(alpha=0.3)
plt.tight_layout()

# Save figures
plt.savefig('results/Omega_USF_evolution.pdf', bbox_inches='tight')
plt.savefig('results/Omega_USF_evolution.png', dpi=300, bbox_inches='tight')
plt.show()

print("Figure saved to results/Omega_USF_evolution.pdf")