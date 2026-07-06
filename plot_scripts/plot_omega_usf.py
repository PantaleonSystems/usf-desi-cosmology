"""
Evolution of the effective geometric correction Ω_USF(z).

Plots the best-fit USF correction term as a function of redshift, together
with the ΛCDM asymptotic limit and the transition epoch.  The Ω_USF(z)
function is imported from theory_module.background — the exact expression
used by the sampler (Topic 9 of the correction plan).

Best-fit parameters are read from the USF chain when available; otherwise
the values passed on the command line are used.

Usage (from the repository root):
    python plot_scripts/plot_omega_usf.py [--chain chains/usf_bao_sn_sh0es]
        [--alpha-q -0.1] [--z-trans 4.66]

Output (saved to results/):
    Omega_USF_evolution.pdf / .png

Author: Pantaleon Systems
Date: 2026
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from theory_module import background as bg  # noqa: E402
from usf_constants import DELTA  # noqa: E402

os.makedirs('results', exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--chain', default='./chains/usf_bao_sn_sh0es')
parser.add_argument('--alpha-q', type=float, default=None)
parser.add_argument('--z-trans', type=float, default=None)
args = parser.parse_args()

alpha_q, z_trans = args.alpha_q, args.z_trans
if alpha_q is None or z_trans is None:
    try:
        from getdist.mcsamples import loadMCSamples
        stats = loadMCSamples(args.chain).getMargeStats()
        alpha_q = stats.parWithName('alpha_q').mean
        z_trans = stats.parWithName('z_trans').mean
        print(f"Chain means: alpha_q = {alpha_q:.4f}, z_trans = {z_trans:.2f}")
    except Exception as exc:  # chain absent: require explicit values
        raise SystemExit(
            f"Could not read chain '{args.chain}' ({exc}); "
            "pass --alpha-q and --z-trans explicitly."
        )

# ============================================================
# Model curve (shared background module)
# ============================================================
z_vals = np.linspace(0, 10, 1000)
omega_vals = bg.omega_usf(z_vals, alpha_q, z_trans, DELTA)

plt.figure(figsize=(8, 5))
plt.plot(z_vals, omega_vals, color='red', linewidth=2, label='USF best-fit')
plt.axhline(0, color='gray', linestyle='--', linewidth=1.5,
            label=r'$\Lambda$CDM (asymptotic limit)')
plt.axvline(z_trans, color='blue', linestyle=':', linewidth=1.5, alpha=0.6,
            label=rf'Transition epoch ($z_{{\rm trans}} = {z_trans:.2f}$)')

z_fill = np.linspace(z_trans - 0.5, z_trans + 0.5, 200)
plt.fill_between(z_fill, bg.omega_usf(z_fill, alpha_q, z_trans, DELTA),
                 alpha=0.1, color='blue')

plt.xlabel(r'Redshift $z$', fontsize=14)
plt.ylabel(r'Effective geometric correction $\Omega_{\rm USF}(z)$', fontsize=14)
plt.title('Evolution of the USF geometric correction', fontsize=16)
plt.legend(loc='lower right', fontsize=12)
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig('results/Omega_USF_evolution.pdf', bbox_inches='tight')
plt.savefig('results/Omega_USF_evolution.png', dpi=300, bbox_inches='tight')

print("Figure saved to results/Omega_USF_evolution.pdf")
