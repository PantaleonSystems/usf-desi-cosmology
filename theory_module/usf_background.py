"""
Unified State Function (USF) cosmological model – Cobaya theory module.

Phenomenological extension of ΛCDM in which the large-scale cosmic expansion
receives a small, effective geometric correction governed by two additional
parameters: α_q and z_trans.

All physics is delegated to :mod:`theory_module.background` (single source of
truth shared with the likelihoods, the plot scripts and the test suite —
Topic 9 of pdfs/PLANO_CORRECAO_REFEREE1.md).  This class acts as a Cobaya
*provider*: its ``get_*`` methods are auto-discovered and serve both the
built-in SN likelihoods (angular_diameter_distance) and the custom BAO / CMB
distance-prior likelihoods (H, comoving_distance, r_s, sound scales).

Sampled parameters: H0, Omega_m (total matter), omega_b (= Ω_b h², with a BBN
Gaussian prior set in the YAML), alpha_q and z_trans.  The physical radiation
density is held fixed (ω_r = ω_γ(1 + 0.2271 N_eff), N_eff = 3.044); Σm_ν =
0.06 eV is counted inside Omega_m at late times (see docs/NEUTRINOS.md).

Usage (Cobaya YAML):
    theory:
        theory_module.usf_background.USF:

Author: Pantaleon Systems
Date: 2026
"""

import numpy as np
from cobaya.theory import Theory

from theory_module import background as bg
from usf_constants import DELTA


class USF(Theory):
    """USF effective geometric model for cosmic expansion.

    Extra parameters on top of flat ΛCDM:
        - alpha_q : amplitude of the geometric correction (dimensionless);
        - z_trans : redshift marking the onset of the correction.

    E(z) is normalised so that E(0) = 1 by construction (flatness fixes
    Ω_Λ = 1 − Ω_m − Ω_r − Ω_USF(0)); H(0) ≡ H0 for any value of alpha_q.
    """

    # Parameters that Cobaya will inject at runtime.
    params = {
        "H0": None,
        "Omega_m": None,
        "omega_b": None,
        "alpha_q": None,
        "z_trans": None,
    }

    _delta: float = DELTA

    def initialize(self) -> None:
        """One-time setup (warm the z_drag/z_star interpolation tables)."""
        bg.z_drag(0.02218, 0.14)

    # ------------------------------------------------------------------
    # Provider methods (auto-discovered by Cobaya's requirements system)
    # ------------------------------------------------------------------

    def get_H(self, z):
        """Hubble parameter H(z) in km/s/Mpc."""
        return bg.hubble(z, self.H0, self.Omega_m, self.alpha_q, self.z_trans)

    def get_comoving_distance(self, z):
        """Line-of-sight comoving distance D_C(z) [Mpc]."""
        return bg.comoving_distance(z, self.H0, self.Omega_m,
                                    self.alpha_q, self.z_trans)

    def get_angular_diameter_distance(self, z):
        """Angular diameter distance D_A(z) [Mpc]; scalar or array z."""
        return bg.angular_diameter_distance(z, self.H0, self.Omega_m,
                                            self.alpha_q, self.z_trans)

    def get_r_s(self):
        """Sound horizon at the drag epoch r_d [Mpc].

        Computed from the physical densities (ω_b sampled, ω_m = Ω_m h²)
        with the CAMB-calibrated drag redshift — independent of H0 and of
        the USF late-time term by construction.
        """
        return bg.r_drag(self.omega_b, self._omega_m_h2())

    def get_rs_zstar(self):
        """Sound horizon at recombination r_s(z*) [Mpc] (CMB priors)."""
        omega_m = self._omega_m_h2()
        return bg.sound_horizon(bg.z_star(self.omega_b, omega_m),
                                self.omega_b, omega_m)

    def get_z_star(self):
        """Recombination redshift z*(ω_b, ω_m) (CAMB-calibrated)."""
        return bg.z_star(self.omega_b, self._omega_m_h2())

    # ------------------------------------------------------------------
    # Cobaya internal methods
    # ------------------------------------------------------------------

    def _omega_m_h2(self):
        return self.Omega_m * (self.H0 / 100.0) ** 2

    def calculate(self, state: dict, want_derived: bool = True,
                  **params_values_dict) -> None:
        """Persist the sampled parameter values for the provider methods."""
        self.H0 = params_values_dict["H0"]
        self.Omega_m = params_values_dict["Omega_m"]
        self.omega_b = params_values_dict["omega_b"]
        self.alpha_q = params_values_dict["alpha_q"]
        self.z_trans = params_values_dict["z_trans"]
        state["Hubble"] = {"z": np.array([0.0]), "H": np.array([self.H0])}

    def get_requirements(self) -> dict:
        """This theory does not depend on any other component."""
        return {}
