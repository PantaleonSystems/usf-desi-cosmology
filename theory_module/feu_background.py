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

    # Fixed physical constants and model settings.
    _delta: float = 0.25
    _c_km_s: float = 299792.458           # speed of light [km/s]

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
        Omega_L = 1.0 - self.Omega_m

        # Logistic suppression factor.
        arg = np.clip((z - self.z_trans) / self._delta, -500, 500)
        suppression = 1.0 / (1.0 + np.exp(arg))

        # Effective geometric correction term.
        Oz = (self.alpha_q * suppression * (1.0 + z)**3 /
              (1.0 + (1.0 + z)**2 / self.z_trans**2))
        Esq = self.Omega_m * (1.0 + z)**3 + Omega_L + Oz

        # Normalisation so that E(0) = 1.
        supp0 = 1.0 / (1.0 + np.exp(-self.z_trans / self._delta))
        Oz0 = self.alpha_q * supp0 / (1.0 + 1.0 / self.z_trans**2)
        Esq0 = self.Omega_m + Omega_L + Oz0

        return np.sqrt(np.maximum(Esq / Esq0, 1e-30))

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

    def get_angular_diameter_distance(self, z: float) -> float:
        """
        Angular diameter distance to redshift z.

        .. math:: D_A(z) = \\frac{D_C(z)}{1+z}

        :param z: Target redshift.
        :returns: Angular diameter distance in Mpc.
        """
        return self.get_comoving_distance(z) / (1.0 + z)

    def get_r_s(self, z_drag: float = 1059.94) -> float:
        """
        Sound horizon at the drag epoch.

        Uses the standard baryon-to-photon ratio and the adiabatic sound
        speed, integrated from z_drag to infinity.

        :param z_drag: Redshift of the drag epoch (default Planck 2018).
        :returns: Sound horizon r_d in Mpc.
        """
        # Fixed baryon and radiation parameters (Planck 2018).
        Omega_b = 0.0486
        rho_gamma_fac = 8.24e-5

        def integrand(z):
            R = 3.0 * Omega_b / (4.0 * rho_gamma_fac * (1.0 + z))
            cs = self._c_km_s / np.sqrt(3.0 * (1.0 + R))
            return cs / (self.H0 * self._E(z))

        result, _ = integrate.quad(
            integrand, z_drag, np.inf,
            limit=200, epsabs=1e-6, epsrel=1e-6
        )
        return result

    # ------------------------------------------------------------------
    # Cobaya internal methods
    # ------------------------------------------------------------------

    def calculate(self, state: dict, want_derived: bool = True,
                  **params_values_dict) -> None:
        """
        Provide the Hubble parameter at z=0 to satisfy Cobaya's requirements.

        This method is called by Cobaya to initialise the state dictionary
        with any derived quantities.  For the USF model we only need to
        ensure that the Hubble constant is correctly reported.

        :param state: Cobaya state dictionary (modified in-place).
        :param want_derived: Whether derived parameters are requested.
        :param params_values_dict: Current parameter values.
        """
        H0 = params_values_dict["H0"]
        state["Hubble"] = {"z": np.array([0.0]), "H": np.array([H0])}

    def get_requirements(self) -> dict:
        """
        Return an empty dictionary: this theory does not depend on any
        other component.

        :returns: Empty requirements dictionary.
        """
        return {}
