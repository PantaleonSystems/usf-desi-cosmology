"""
DESI DR2 BAO likelihood for the Unified State Function (USF) model.

Custom Cobaya likelihood comparing the USF background expansion with the
BAO measurements from DESI Data Release 2 (D_M/r_d, D_H/r_d, D_V/r_d).

The expansion history and the sound horizon are obtained from the theory
component through Cobaya's requirements/provider system (H, comoving
distance and r_s are served by theory_module.usf_background.USF, which
delegates to theory_module.background — single source of truth).  r_d is
evaluated per sample from the physical densities (ω_b sampled with a BBN
prior; ω_m = Ω_m h²) with a CAMB-calibrated drag redshift, so the BAO ruler
is fully self-consistent with the sampled cosmology (Referee 1, Comment 1).

Usage:
    Place the directory `desi_bao_dr2/` (mean and covariance files) inside
    the Cobaya packages data folder and set `data_path` in the YAML.

Author: Pantaleon Systems
Date: 2026
"""

import os

import numpy as np
from cobaya.likelihood import Likelihood

from usf_constants import C_KM_S


class BAODESILikelihood(Likelihood):
    """Likelihood for DESI DR2 BAO data within the USF model.

    :ivar data_path: directory containing the DESI mean and covariance files.
    """

    data_path = ""

    def get_requirements(self):
        """Quantities served by the theory provider (USF get_* methods)."""
        return {"H": None, "comoving_distance": None, "r_s": None}

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

    def logp(self, **kwargs):
        """Gaussian log-likelihood of the DESI BAO data given the USF model.

        All theory quantities come from the provider (no local physics)."""
        rs = self.provider.get_r_s()

        pred = []
        for z, otype in zip(self.z_eff, self.obs_type):
            if otype == "DM_over_rs":
                DM = self.provider.get_comoving_distance(z)
                pred.append(DM / rs)
            elif otype == "DH_over_rs":
                Hz = self.provider.get_H(z)
                pred.append(C_KM_S / (Hz * rs))
            elif otype == "DV_over_rs":
                DM = self.provider.get_comoving_distance(z)
                Hz = self.provider.get_H(z)
                DH = C_KM_S / Hz
                DV = (z * DM ** 2 * DH) ** (1.0 / 3.0)
                pred.append(DV / rs)
            else:
                raise ValueError(f"Unknown BAO observable type: {otype}")

        delta = np.array(pred) - self.obs_val
        return -0.5 * delta @ self.cov_inv @ delta
