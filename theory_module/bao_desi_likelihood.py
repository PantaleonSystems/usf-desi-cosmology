"""
DESI DR2 BAO likelihood for the Unified State Function (USF) model.

Custom Cobaya likelihood comparing the USF background expansion with the
BAO measurements from DESI Data Release 2 (D_M/r_d, D_H/r_d, D_V/r_d).

All physics comes from :mod:`theory_module.background` — the same functions
used by the theory class, the plot scripts and the test suite.  The sound
horizon r_d is evaluated per sample from the physical densities
(ω_b sampled with a BBN prior; ω_m = Ω_m h²) with a CAMB-calibrated drag
redshift, so the BAO ruler is fully self-consistent with the sampled
cosmology (Referee 1, Comment 1).

Usage:
    Place the directory `desi_bao_dr2/` (mean and covariance files) inside
    the Cobaya packages data folder and set `data_path` in the YAML.

Author: Pantaleon Systems
Date: 2026
"""

import os

import numpy as np
from cobaya.likelihood import Likelihood

from theory_module import background as bg
from usf_constants import C_KM_S


class BAODESILikelihood(Likelihood):
    """Likelihood for DESI DR2 BAO data within the USF model.

    :ivar data_path: directory containing the DESI mean and covariance files.
    """

    data_path = ""

    # Declaration of the input parameters expected by this likelihood.
    params = {
        "H0": None,
        "Omega_m": None,
        "omega_b": None,
        "alpha_q": None,
        "z_trans": None,
    }

    def initialize(self):
        """Read the DESI mean vector and covariance matrix once."""
        mean_file = os.path.join(self.data_path,
                                 "desi_gaussian_bao_ALL_GCcomb_mean.txt")
        cov_file = os.path.join(self.data_path,
                                "desi_gaussian_bao_ALL_GCcomb_cov.txt")

        rows = []
        with open(mean_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                rows.append((float(parts[0]), float(parts[1]), parts[2]))

        self.z_eff = np.array([r[0] for r in rows])
        self.obs_val = np.array([r[1] for r in rows])
        self.obs_type = [r[2] for r in rows]

        cov = np.loadtxt(cov_file)
        self.cov_inv = np.linalg.inv(cov)

    def logp(self, H0, Omega_m, omega_b, alpha_q, z_trans, **kwargs):
        """Gaussian log-likelihood of the DESI BAO data given the USF model."""
        omega_m = Omega_m * (H0 / 100.0) ** 2
        rs = bg.r_drag(omega_b, omega_m)

        pred = []
        for z, otype in zip(self.z_eff, self.obs_type):
            if otype == "DM_over_rs":
                DM = bg.comoving_distance(z, H0, Omega_m, alpha_q, z_trans)
                pred.append(DM / rs)
            elif otype == "DH_over_rs":
                Hz = bg.hubble(z, H0, Omega_m, alpha_q, z_trans)
                pred.append(C_KM_S / (Hz * rs))
            elif otype == "DV_over_rs":
                DM = bg.comoving_distance(z, H0, Omega_m, alpha_q, z_trans)
                Hz = bg.hubble(z, H0, Omega_m, alpha_q, z_trans)
                DH = C_KM_S / Hz
                DV = (z * DM ** 2 * DH) ** (1.0 / 3.0)
                pred.append(DV / rs)
            else:
                raise ValueError(f"Unknown BAO observable type: {otype}")

        delta = np.array(pred) - self.obs_val
        return -0.5 * delta @ self.cov_inv @ delta
