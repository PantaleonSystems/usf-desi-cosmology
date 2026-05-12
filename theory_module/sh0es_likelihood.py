"""
Gaussian likelihood for the local Hubble constant measurement (SH0ES).

This module provides a simple external likelihood for Cobaya that
encodes the 2023 SH0ES determination of H0 = 73.0 ± 1.04 km/s/Mpc.
It can be used as a prior-level constraint in any cosmological model.

Usage:
    Include this file in a package and add to the Cobaya YAML configuration:
        likelihood:
            theory_module.sh0es_likelihood.SH0ESLikelihood:

Author: Pantaleon Systems
Date: 2026
"""

import numpy as np
from cobaya.likelihood import Likelihood


class SH0ESLikelihood(Likelihood):
    """
    Gaussian likelihood representing the local measurement of H0.

    The adopted central value and uncertainty are taken from the
    SH0ES collaboration (Riess et al. 2023), based on Cepheid,
    TRGB, and Mira distance-ladder calibrations.

    :ivar H0_obs: Central value of the Hubble constant (km/s/Mpc).
    :ivar sigma_H0: 1σ uncertainty on H0 (km/s/Mpc).
    :ivar params: Dictionary declaring the required parameter(s).
    """

    # Central value and uncertainty (overridable in YAML).
    H0_obs:   float = 73.0
    sigma_H0: float = 1.04

    # The likelihood needs only the Hubble constant.
    params = {"H0": None}

    def initialize(self):
        """
        Perform any one-time setup (none required for a simple Gaussian).
        """
        pass

    def logp(self, H0: float, **kwargs) -> float:
        """
        Compute the log-likelihood for a given H0 value.

        Assumes a Gaussian form:
            log p = -0.5 * ((H0 - H0_obs) / sigma_H0)^2
        (normalisation constant omitted, as it cancels in MCMC).

        :param H0: Hubble constant value (km/s/Mpc).
        :returns: Log-likelihood value (float).
        """
        return -0.5 * ((H0 - self.H0_obs) / self.sigma_H0) ** 2
