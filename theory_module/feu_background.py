"""
Unified State Function (USF) cosmological model – Cobaya theory module.

This module implements a phenomenological extension of ΛCDM in which
the large-scale cosmic expansion receives a small, effective geometric
correction governed by two additional parameters: α_q and z_trans.

The theory class provides the modified Hubble parameter H(z), the
comoving and angular diameter distances, and the sound horizon at the
drag epoch.  All quantities are computed internally and returned in the
format expected by Cobaya likelihoods.

Usage:
    Include this file in a package and refer to it in the Cobaya YAML
    configuration as:
        theory:
            theory_module.feu_background.FEU:

Author: Pantaleon Systems
Date: 2026
"""

import numpy as np
from scipy import integrate
from cobaya.theory import Theory
from usf_constants import C_KM_S, OMEGA_B_H2, OMEGA_GAMMA_H2, Z_DRAG, DELTA


class FEU(Theory):
    """
    USF effective geometric model for cosmic expansion.

    This theory introduces two extra parameters on top of the standard
    flat-ΛCDM background:
        - α_q  : amplitude of the geometric correction (dimensionless)
        - z_trans : redshift marking the onset of the correction

    The dimensionless Hubble ratio E(z) = H(z)/H0 is defined such that
    E(0) = 1 by construction, ensuring H(0) ≡ H0 regardless of the
    value of α_q.

    :ivar params: Dictionary of input parameters declared to Cobaya.
    :ivar _delta: Fixed transition width (see Sec. 2 of the paper).
    :ivar _c_km_s: Speed of light in km/s, used for distance calculations.
    """

    # Parameters that Cobaya will inject at runtime.
    params = {
        "H0": None,
        "Omega_m": None,
        "alpha_q": None,
        "z_trans": None,
    }

    _delta: float = DELTA
    _c_km_s: float = C_KM_S

    def initialize(self) -> None:
        """Perform any one-time setup (no heavy work needed here)."""
        pass

    # ------------------------------------------------------------------
    # Internal helper
    # ------------------------------------------------------------------

    def _E(self, z: np.ndarray) -> np.ndarray:
        """
        Dimensionless Hubble parameter E(z) = H(z) / H0.

        The function implements the USF-modified Friedmann equation,
        normalised to unity at z = 0.

        :param z: Redshift (float or array).
        :returns: E(z) values.
        """
        z = np.asarray(z, dtype=np.float64)
        Omega_r = OMEGA_GAMMA_H2 / (self.H0 / 100.0) ** 2

        # USF correction at z=0 (required before setting Omega_L)
        supp0 = 1.0 / (1.0 + np.exp(-self.z_trans / self._delta))
        Oz0 = self.alpha_q * supp0 / (1.0 + 1.0 / self.z_trans ** 2)

        # Flat-universe condition: Ω_m + Ω_r + Ω_Λ + Ω_USF(0) = 1
        Omega_L = 1.0 - self.Omega_m - Omega_r - Oz0

        arg = np.clip((z - self.z_trans) / self._delta, -500, 500)
        suppression = 1.0 / (1.0 + np.exp(arg))
        Oz = (self.alpha_q * suppression * (1.0 + z) ** 3 /
              (1.0 + (1.0 + z) ** 2 / self.z_trans ** 2))
        # With the flat-universe Omega_L, E²(0) = 1 by construction.
        Esq = (self.Omega_m * (1.0 + z) ** 3 + Omega_r * (1.0 + z) ** 4
               + Omega_L + Oz)

        return np.sqrt(np.maximum(Esq, 1e-30))

    # ------------------------------------------------------------------
    # Public interface expected by Cobaya likelihoods
    # ------------------------------------------------------------------

    def get_H(self, z: np.ndarray) -> np.ndarray:
        """
        Hubble parameter H(z) in km/s/Mpc.

        :param z: Redshift (float or array).
        :returns: H(z) array.
        """
        return self.H0 * self._E(z)

    def get_comoving_distance(self, z: float) -> float:
        """
        Line-of-sight comoving distance to redshift z.

        .. math:: D_C(z) = c \\int_0^z \\frac{dz'}{H(z')}

        :param z: Target redshift.
        :returns: Comoving distance in Mpc.
        """
        result, _ = integrate.quad(
            lambda zp: self._c_km_s / (self.H0 * self._E(zp)),
            0.0, z
        )
        return result

    def get_angular_diameter_distance(self, z):
        """
        Angular diameter distance to redshift z.

        .. math:: D_A(z) = \\frac{D_C(z)}{1+z}

        Accepts scalar or array z (Cobaya's sn.pantheonplus passes an array
        of ~1701 SN redshifts via the provider system).

        :param z: Target redshift(s) — float or array-like.
        :returns: Angular diameter distance(s) in Mpc, same shape as input.
        """
        z_arr = np.atleast_1d(np.asarray(z, dtype=np.float64))
        d_c = np.array([
            integrate.quad(
                lambda zp: self._c_km_s / (self.H0 * self._E(zp)),
                0.0, zi
            )[0]
            for zi in z_arr
        ])
        return d_c / (1.0 + z_arr)

    def get_r_s(self, z_drag: float = Z_DRAG) -> float:
        """
        Sound horizon at the drag epoch.
        """
        def integrand(z):
            R = 3.0 * OMEGA_B_H2 / (4.0 * OMEGA_GAMMA_H2 * (1.0 + z))
            cs = self._c_km_s / np.sqrt(3.0 * (1.0 + R))
            return cs / (self.H0 * self._E(z))

        result, _ = integrate.quad(integrand, z_drag, np.inf,
                                   limit=200, epsabs=1e-6, epsrel=1e-6)
        return result

    # ------------------------------------------------------------------
    # Cobaya internal methods
    # ------------------------------------------------------------------

    def calculate(self, state: dict, want_derived: bool = True,
                  **params_values_dict) -> None:
        """
        Store the current parameter values and populate the state dictionary.

        Cobaya passes sampled parameter values as kwargs each step.  We persist
        them as instance attributes so that provider methods like
        get_angular_diameter_distance (called after calculate returns) can use
        self.H0 / self.Omega_m / self.alpha_q / self.z_trans normally.

        :param state: Cobaya state dictionary (modified in-place).
        :param want_derived: Whether derived parameters are requested.
        :param params_values_dict: Current parameter values.
        """
        self.H0 = params_values_dict["H0"]
        self.Omega_m = params_values_dict["Omega_m"]
        self.alpha_q = params_values_dict["alpha_q"]
        self.z_trans = params_values_dict["z_trans"]
        state["Hubble"] = {"z": np.array([0.0]), "H": np.array([self.H0])}

    def get_requirements(self) -> dict:
        """
        Return an empty dictionary: this theory does not depend on any
        other component.

        :returns: Empty requirements dictionary.
        """
        return {}
