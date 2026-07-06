"""
Single source of truth for the USF background physics.

Every consumer of the expansion history — the Cobaya theory class, the BAO
and CMB-distance likelihoods, the plotting scripts and the test suite —
imports the pure functions defined here, so the sampler and the
visualisation layer are guaranteed to evaluate identical physics
(Topics 2 and 9 of pdfs/PLANO_CORRECAO_REFEREE1.md).

Conventions
-----------
* ``Omega_m`` is the present-day **total** matter fraction (CDM + baryons +
  non-relativistic neutrinos), as constrained by late-time data.
* The radiation sector holds the **physical** density fixed:
  Omega_r(H0) = OMEGA_R_H2 / h², with ω_r = ω_γ (1 + 0.2271 N_eff).
* ``omega_b`` (= Ω_b h²) is a sampled parameter (BBN Gaussian prior); it
  only enters the early-time sound-horizon sector.
* The drag redshift is interpolated from a CAMB-calibrated table
  (theory_module/data/zdrag_camb_grid.npz).  The Eisenstein & Hu (1998)
  fitting formula is provided for reference/fallback but gives z_d ≈ 1021
  for Planck densities versus CAMB's ≈ 1060, which would bias r_d by
  ≈ +2.4%; it is therefore NOT used in production.
* r_drag is computed from the numerical integral of c_s/H using early-time
  physical densities (ω_cb + ω_r; neutrinos are relativistic for
  z > z_drag).  The Aubourg et al. (2015) fitting formula is provided as an
  independent cross-check (accurate to 0.021% versus CAMB).
"""

import os

import numpy as np
from scipy import integrate
from scipy.interpolate import RegularGridInterpolator

from usf_constants import (
    C_KM_S,
    DELTA,
    OMEGA_GAMMA_H2,
    OMEGA_NU_H2,
    OMEGA_R_H2,
)

_ZDRAG_TABLE = os.path.join(os.path.dirname(__file__), "data", "zdrag_camb_grid.npz")
_zdrag_interp = None
_zstar_interp = None


# ---------------------------------------------------------------------------
# Expansion history
# ---------------------------------------------------------------------------

def omega_usf(z, alpha_q, z_trans, delta=DELTA):
    """Effective geometric correction Ω_USF(z) (dimensionless density term)."""
    z = np.asarray(z, dtype=np.float64)
    arg = np.clip((z - z_trans) / delta, -500.0, 500.0)
    suppression = 1.0 / (1.0 + np.exp(arg))
    return (alpha_q * suppression * (1.0 + z) ** 3
            / (1.0 + (1.0 + z) ** 2 / z_trans ** 2))


def E2(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0, delta=DELTA):
    """Squared dimensionless Hubble rate E²(z) = (H(z)/H0)², with E(0) = 1.

    Flatness fixes Ω_Λ = 1 − Ω_m − Ω_r − Ω_USF(0): a nonzero Ω_USF(0) is
    compensated by Ω_Λ (not by H0 — see Topic 8 of the correction plan).
    """
    z = np.asarray(z, dtype=np.float64)
    Omega_r = OMEGA_R_H2 / (H0 / 100.0) ** 2
    Oz0 = omega_usf(0.0, alpha_q, z_trans, delta) if alpha_q != 0.0 else 0.0
    Omega_L = 1.0 - Omega_m - Omega_r - Oz0
    Oz = omega_usf(z, alpha_q, z_trans, delta) if alpha_q != 0.0 else 0.0
    return (Omega_m * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4
            + Omega_L + Oz)


def E(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0, delta=DELTA):
    """Dimensionless Hubble rate E(z) = H(z)/H0."""
    return np.sqrt(np.maximum(E2(z, H0, Omega_m, alpha_q, z_trans, delta), 1e-30))


def hubble(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0, delta=DELTA):
    """H(z) in km/s/Mpc."""
    return H0 * E(z, H0, Omega_m, alpha_q, z_trans, delta)


def comoving_distance(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0, delta=DELTA,
                      grid_size=4000):
    """Line-of-sight comoving distance D_C(z) [Mpc]; scalar or array z.

    Vectorised cumulative-trapezoid integration on a fine grid followed by
    interpolation (sub-0.01% accuracy for grid_size = 4000).
    """
    z_arr = np.atleast_1d(np.asarray(z, dtype=np.float64))
    z_max = float(z_arr.max())
    if z_max == 0.0:
        out = np.zeros_like(z_arr)
        return out[0] if np.isscalar(z) or np.ndim(z) == 0 else out

    z_grid = np.linspace(0.0, z_max, grid_size)
    integrand = C_KM_S / (H0 * E(z_grid, H0, Omega_m, alpha_q, z_trans, delta))
    dc_grid = integrate.cumulative_trapezoid(integrand, z_grid, initial=0.0)
    out = np.interp(z_arr, z_grid, dc_grid)
    return float(out[0]) if np.ndim(z) == 0 else out


def angular_diameter_distance(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0,
                              delta=DELTA):
    """Angular diameter distance D_A(z) = D_C(z)/(1+z) [Mpc]."""
    z_arr = np.asarray(z, dtype=np.float64)
    return comoving_distance(z, H0, Omega_m, alpha_q, z_trans, delta) / (1.0 + z_arr)


def luminosity_distance(z, H0, Omega_m, alpha_q=0.0, z_trans=2.0, delta=DELTA):
    """Luminosity distance D_L(z) = (1+z) D_C(z) [Mpc]."""
    z_arr = np.asarray(z, dtype=np.float64)
    return comoving_distance(z, H0, Omega_m, alpha_q, z_trans, delta) * (1.0 + z_arr)


# ---------------------------------------------------------------------------
# Early-universe scales (independent of H0 and of the USF late-time term)
# ---------------------------------------------------------------------------

def _load_tables():
    global _zdrag_interp, _zstar_interp
    if _zdrag_interp is None:
        tab = np.load(_ZDRAG_TABLE)
        _zdrag_interp = RegularGridInterpolator(
            (tab["omega_b"], tab["omega_m"]), tab["zdrag"],
            bounds_error=False, fill_value=None)
        _zstar_interp = RegularGridInterpolator(
            (tab["omega_b"], tab["omega_m"]), tab["zstar"],
            bounds_error=False, fill_value=None)


def z_drag(omega_b, omega_m):
    """Drag redshift z_d(ω_b, ω_m), interpolated from the CAMB table."""
    _load_tables()
    return float(_zdrag_interp([[omega_b, omega_m]])[0])


def z_star(omega_b, omega_m):
    """Recombination redshift z*(ω_b, ω_m), interpolated from the CAMB table."""
    _load_tables()
    return float(_zstar_interp([[omega_b, omega_m]])[0])


def z_drag_eh98(omega_b, omega_m):
    """Eisenstein & Hu (1998), Eq. 4 — reference only, NOT used in production.

    Gives z_d ≈ 1021 for Planck densities vs CAMB's ≈ 1060; kept as an
    importable fallback and to document the discrepancy quantitatively.
    """
    b1 = 0.313 * omega_m ** -0.419 * (1.0 + 0.607 * omega_m ** 0.674)
    b2 = 0.238 * omega_m ** 0.223
    return (1291.0 * omega_m ** 0.251 / (1.0 + 0.659 * omega_m ** 0.828)
            * (1.0 + b1 * omega_b ** b2))


def _E2_early(z, omega_cb, omega_l_h2):
    """h² E²(z) at early times, in physical-density form [units of h²].

    Radiation carries all N_eff neutrinos (relativistic for z > z_drag);
    the non-relativistic neutrino mass is *excluded* from the matter term
    (it belongs to Ω_m only at late times).  The Λ term is negligible at
    z > 1000 but kept for exactness.
    """
    return (omega_cb * (1.0 + z) ** 3 + OMEGA_R_H2 * (1.0 + z) ** 4
            + omega_l_h2)


def sound_horizon(z_end, omega_b, omega_m):
    """Comoving sound horizon r_s(z_end) [Mpc] = ∫_{z_end}^∞ c_s/H dz.

    Depends only on physical densities (ω_b, ω_m) — explicitly independent
    of H0 and of the USF late-time term, which vanishes for z ≫ z_trans.

    :param z_end: lower integration limit (z_drag for BAO, z_star for CMB).
    :param omega_b: physical baryon density Ω_b h².
    :param omega_m: physical *total* matter density Ω_m h².
    """
    omega_cb = omega_m - OMEGA_NU_H2  # neutrinos are radiation at z > z_end
    omega_l_h2 = 0.0  # negligible for z > 1000

    def integrand(z):
        R = 3.0 * omega_b / (4.0 * OMEGA_GAMMA_H2 * (1.0 + z))
        cs = C_KM_S / np.sqrt(3.0 * (1.0 + R))
        return cs / (100.0 * np.sqrt(_E2_early(z, omega_cb, omega_l_h2)))

    result, _ = integrate.quad(integrand, z_end, np.inf,
                               limit=200, epsabs=1e-8, epsrel=1e-8)
    return result


def r_drag(omega_b, omega_m):
    """Sound horizon at the drag epoch r_d [Mpc], CAMB-calibrated z_d."""
    return sound_horizon(z_drag(omega_b, omega_m), omega_b, omega_m)


def r_drag_aubourg15(omega_b, omega_m):
    """Aubourg et al. (2015), Eq. 16 — independent cross-check for r_d.

    Accurate to 0.021% against CAMB for standard N_eff and Σm_ν ≲ 0.6 eV.
    """
    omega_cb = omega_m - OMEGA_NU_H2
    return (55.154 * np.exp(-72.3 * (OMEGA_NU_H2 + 0.0006) ** 2)
            / (omega_cb ** 0.25351 * omega_b ** 0.12807))
