"""
DESI DR2 BAO likelihood for the Unified State Function (USF) model.

This module provides a custom Cobaya likelihood that evaluates the agreement
between the USF background expansion and the Baryon Acoustic Oscillation
measurements from the Dark Energy Spectroscopic Instrument (DESI) Data Release 2.

Usage:
    Place the directory `desi_bao_dr2/` (containing the mean and covariance
    files) inside the Cobaya packages data folder and set the `data_path`
    option accordingly in the YAML configuration.

Author: Pantaleon Systems
Date: 2026
"""

import os
import numpy as np
from scipy import integrate
from cobaya.likelihood import Likelihood


class BAODESILikelihood(Likelihood):
    """
    Likelihood for DESI DR2 BAO data within the USF effective geometric model.

    This likelihood computes the key BAO observables – transverse comoving
    distance D_M/r_d, Hubble distance D_H/r_d, and the angle-averaged
    distance D_V/r_d – using the USF expansion history, and compares them
    with the official DESI measurements.

    :ivar data_path: Path to the directory containing the DESI mean and
        covariance files.
    :ivar params: Dictionary of required parameter names (all mapped to
        the USF model).
    """

    # Path to the data directory; to be set in the YAML configuration.
    data_path: str = ""

    # Declaration of the input parameters expected by this likelihood.
    params = {
        "H0": None,
        "Omega_m": None,
        "alpha_q": None,
        "z_trans": None,
    }

    def initialize(self):
        """
        Called once at the start of the MCMC run.

        Reads the DESI mean vector and covariance matrix, stores the
        effective redshifts and the inverse covariance for fast evaluation
        of the log-likelihood.
        """
        mean_file = os.path.join(self.data_path,
                                 "desi_gaussian_bao_ALL_GCcomb_mean.txt")
        cov_file  = os.path.join(self.data_path,
                                 "desi_gaussian_bao_ALL_GCcomb_cov.txt")

        # Parse the mean file (columns: z, value, observable type).
        rows = []
        with open(mean_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                rows.append((float(parts[0]), float(parts[1]), parts[2]))

        self.z_eff     = np.array([r[0] for r in rows])
        self.obs_val   = np.array([r[1] for r in rows])
        self.obs_type  = [r[2] for r in rows]

        # Load covariance matrix and precompute its inverse.
        cov = np.loadtxt(cov_file)
        self.cov_inv = np.linalg.inv(cov)

        # Physical constants and fixed parameters.
        self._c_km_s = 299792.458          # speed of light in km/s
        self._delta  = 0.25                # transition width (Sec. 2)

    def _E(self, z, H0, Omega_m, alpha_q, z_trans):
        """
        Dimensionless Hubble parameter E(z) = H(z) / H0.

        :param z: Redshift (float or array).
        :param H0: Hubble constant (km/s/Mpc).
        :param Omega_m: Present-day matter density fraction.
        :param alpha_q: Amplitude of the geometric correction.
        :param z_trans: Transition redshift.
        :returns: E(z) array.
        """
        z = np.asarray(z, dtype=np.float64)
        Omega_L = 1.0 - Omega_m

        # Inside _E, after Omega_L = 1.0 - Omega_m
        Omega_r = 4.15e-5 / (H0 / 100.0)**2
        # Logistic suppression factor.
        suppression = 1.0 / (1.0 + np.exp(
            np.clip((z - z_trans) / self._delta, -500, 500)))
        Oz = (alpha_q * suppression * (1.0 + z)**3
              / (1.0 + (1.0 + z)**2 / z_trans**2))
        Esq = (Omega_m * (1.0 + z)**3 + Omega_r * (1.0 + z)**4 + Omega_L + Oz)
# also update Esq0 simi

        # Normalise so that E(z=0) = 1.
        supp0 = 1.0 / (1.0 + np.exp(-z_trans / self._delta))
        Oz0   = alpha_q * supp0 / (1.0 + 1.0 / z_trans**2)
        Esq0 = Omega_m + Omega_r + Omega_L + Oz0

        return np.sqrt(np.maximum(Esq / Esq0, 1e-30))

    def _comoving_distance(self, z, H0, Omega_m, alpha_q, z_trans):
        """
        Line-of-sight comoving distance to redshift z.

        .. math:: D_C(z) = c \\int_0^z \\frac{dz'}{H(z')}.

        :param z: Target redshift.
        :param H0: Hubble constant (km/s/Mpc).
        :param Omega_m: Matter density fraction.
        :param alpha_q: Geometric correction amplitude.
        :param z_trans: Transition redshift.
        :returns: Comoving distance in Mpc.
        """
        result, _ = integrate.quad(
            lambda zp: self._c_km_s / (H0 * self._E(zp, H0, Omega_m,
                                                     alpha_q, z_trans)),
            0.0, z
        )
        return result

    def _r_s(self, H0, Omega_m, alpha_q, z_trans, z_drag=1059.94):
        """
        Sound horizon at the drag epoch.

        The calculation assumes a standard baryon-to-photon ratio and uses
        the adiabatic sound speed.

        :param H0: Hubble constant (km/s/Mpc).
        :param Omega_m: Matter density fraction.
        :param alpha_q: Geometric correction amplitude.
        :param z_trans: Transition redshift.
        :param z_drag: Redshift of the drag epoch (default from Planck 2018).
        :returns: Sound horizon r_d in Mpc.
        """
        # Baryon and radiation densities (fixed to Planck 2018 best‑fit).
        Omega_b = 0.0486       # physical baryon density today (Omega_b h^2 / h^2)
        rho_gamma_fac = 8.24e-5  # factor converting Omega_b to the radiation epoch

        def integrand(z):
            R = 3.0 * Omega_b / (4.0 * rho_gamma_fac * (1.0 + z))
            cs = self._c_km_s / np.sqrt(3.0 * (1.0 + R))
            return cs / (H0 * self._E(z, H0, Omega_m, alpha_q, z_trans))

        result, _ = integrate.quad(integrand, z_drag, np.inf,
                                   limit=200, epsabs=1e-6, epsrel=1e-6)
        return result

    def logp(self, H0, Omega_m, alpha_q, z_trans, **kwargs):
        """
        Compute the log-likelihood of the DESI BAO data given the USF model.

        For each effective redshift, the corresponding BAO observable
        (D_M/r_d, D_H/r_d, or D_V/r_d) is computed and compared with the
        measured value using the inverse covariance matrix.

        :param H0: Hubble constant (km/s/Mpc).
        :param Omega_m: Present-day matter density fraction.
        :param alpha_q: Amplitude of the geometric correction.
        :param z_trans: Transition redshift.
        :returns: Log-likelihood value (float).
        """
        rs = self._r_s(H0, Omega_m, alpha_q, z_trans)

        pred = []
        for z, otype in zip(self.z_eff, self.obs_type):
            if otype == "DM_over_rs":
                DM = self._comoving_distance(z, H0, Omega_m, alpha_q, z_trans)
                pred.append(DM / rs)
            elif otype == "DH_over_rs":
                Hz = H0 * self._E(z, H0, Omega_m, alpha_q, z_trans)
                pred.append(self._c_km_s / (Hz * rs))
            elif otype == "DV_over_rs":
                DM = self._comoving_distance(z, H0, Omega_m, alpha_q, z_trans)
                Hz = H0 * self._E(z, H0, Omega_m, alpha_q, z_trans)
                DH = self._c_km_s / Hz
                DV = (z * DM**2 * DH) ** (1.0 / 3.0)
                pred.append(DV / rs)
            else:
                raise ValueError(f"Unknown BAO observable type: {otype}")

        pred  = np.array(pred)
        delta = pred - self.obs_val
        return -0.5 * delta @ self.cov_inv @ delta
