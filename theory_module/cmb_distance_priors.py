"""
Planck 2018 CMB distance priors (compressed likelihood) for the USF model.

Implements the acoustic scale ℓ_A, the shift parameter R and ω_b with the
central values, uncertainties and correlation matrix of

    Chen, Huang & Wang, JCAP 02 (2019) 028 (arXiv:1808.05724), Table I,
    ΛCDM row, Planck 2018 TT,TE,EE+lowE:

        R   = 1.7502  ± 0.0046
        ℓ_A = 301.471 ± 0.090
        ω_b = 0.02236 ± 0.00015
        corr(R,ℓ_A) = 0.46, corr(R,ω_b) = −0.66, corr(ℓ_A,ω_b) = −0.33

IMPORTANT — convention consistency (Referee 1, Comment 2): the theory→data
mapping below reproduces *exactly* the recipe used by Chen et al. to
compress the Planck chains (their Eqs. 1–6):

  * z* taken from the Planck/CAMB chains → here interpolated from our
    CAMB-calibrated table (theory_module/data/zdrag_camb_grid.npz);
  * E(z) with the *simplified* radiation density Ω_r = Ω_m/(1+z_eq),
    z_eq = 2.5×10⁴ ω_m (T_CMB/2.7)⁻⁴ (no massive-ν transition);
  * r_s with 3/(4ω_γ) = 31500 (T_CMB/2.7)⁻⁴, T_CMB = 2.7255 K.

This mapping is deliberately *not* the shared background module used for
BAO/SN: the data vector was compressed with this specific convention, and
using any other (e.g. the exact CAMB massive-ν background) would bias ℓ_A
by ≈ 2–3σ purely by definition.  Validated at the Planck 2018 mean:
prediction (R, ℓ_A) = (1.7503, 301.429) vs data (1.7502, 301.471) —
see tests/test_cmb_distance_priors.py.

The USF geometric term enters E(z) with the flat-universe normalisation
Ω_Λ = 1 − Ω_m − Ω_r − Ω_USF(0), so the distance to last scattering feels
the modification along the whole line of sight — precisely the effect the
referee asked to be constrained.

Usage (Cobaya YAML):
    likelihood:
        theory_module.cmb_distance_priors.CMBDistancePriors:
            python_path: .

Author: Pantaleon Systems
Date: 2026
"""

import numpy as np
from cobaya.likelihood import Likelihood
from scipy import integrate

from theory_module import background as bg
from usf_constants import C_KM_S

# --- Data vector and covariance (Chen, Huang & Wang 2019, Table I) ----------
_DATA = np.array([1.7502, 301.471, 0.02236])          # R, l_A, omega_b
_SIGMA = np.array([0.0046, 0.090, 0.00015])
_CORR = np.array([
    [1.00, 0.46, -0.66],
    [0.46, 1.00, -0.33],
    [-0.66, -0.33, 1.00],
])
_COV = _CORR * np.outer(_SIGMA, _SIGMA)

_T_FAC = (2.7255 / 2.7) ** -4  # (T_CMB/2.7 K)^-4


class CMBDistancePriors(Likelihood):
    """Gaussian likelihood on (R, ℓ_A, ω_b) from Planck 2018 TT,TE,EE+lowE."""

    # The CHW19 mapping needs the sampled parameters directly (its E(z)
    # convention differs, by construction, from the shared background).
    params = {
        "H0": None,
        "Omega_m": None,
        "omega_b": None,
        "Omega_TW": None,
        "n": None,
    }

    def initialize(self):
        self.cov_inv = np.linalg.inv(_COV)
        bg.z_star(0.02218, 0.14)  # warm the z* table

    # --- CHW19 background convention (their Eqs. 3, 5, 6) -----------------

    @staticmethod
    def _E_chw(z, Omega_m, omega_m, Omega_TW, n):
        z = np.asarray(z, dtype=np.float64)
        z_eq = 2.5e4 * omega_m * _T_FAC
        Omega_r = Omega_m / (1.0 + z_eq)
        Omega_de = 1.0 - Omega_m - Omega_r - Omega_TW
        Oz = Omega_TW * (1.0 + z) ** n
        return np.sqrt(Omega_r * (1.0 + z) ** 4 + Omega_m * (1.0 + z) ** 3
                       + Omega_de + Oz)

    def _predict(self, H0, Omega_m, omega_b, Omega_TW, n):
        omega_m = Omega_m * (H0 / 100.0) ** 2
        z_star = bg.z_star(omega_b, omega_m)

        # Comoving distance to z* (their Eq. 4, flat case)
        z_grid = np.linspace(0.0, z_star, 120000)
        DM = np.trapezoid(
            C_KM_S / (H0 * self._E_chw(z_grid, Omega_m, omega_m,
                                       Omega_TW, n)), z_grid)

        # Sound horizon at z* (their Eq. 3)
        fac = 31500.0 * omega_b * _T_FAC

        def integrand(z):
            return (C_KM_S
                    / (H0 * self._E_chw(z, Omega_m, omega_m, Omega_TW, n)
                       * np.sqrt(3.0 * (1.0 + fac / (1.0 + z)))))

        rs, _ = integrate.quad(integrand, z_star, np.inf,
                               limit=200, epsabs=1e-8, epsrel=1e-8)

        R = np.sqrt(omega_m) * 100.0 / C_KM_S * DM
        l_A = np.pi * DM / rs
        return R, l_A

    def logp(self, H0, Omega_m, omega_b, Omega_TW, n, **kwargs):
        R, l_A = self._predict(H0, Omega_m, omega_b, Omega_TW, n)
        delta = np.array([R, l_A, omega_b]) - _DATA
        return -0.5 * delta @ self.cov_inv @ delta
