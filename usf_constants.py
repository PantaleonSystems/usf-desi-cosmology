"""
Shared physical constants and fixed cosmological parameters for the USF analysis.

All modules (theory, likelihoods, plot scripts, tests) import from here so the
MCMC engine and the visualisation layer are guaranteed to use identical
numerical values.

Radiation-sector conventions (Referee 1 correction, Topic 2 of
pdfs/PLANO_CORRECAO_REFEREE1.md):

* ``OMEGA_GAMMA_H2`` is the **photon** density ω_γ = Ω_γ h² for
  T_CMB = 2.7255 K.  It enters the baryon-to-photon ratio
  R(z) = 3ω_b / (4ω_γ (1+z)) of the sound-speed in the r_d integral.
  (The previous revision wrongly used the *total* radiation density here,
  inflating r_d by +4.4% — see tests/test_sound_horizon_vs_camb.py.)
* ``OMEGA_R_H2`` is the **total** radiation density
  ω_r = ω_γ (1 + 0.2271 N_eff) with N_eff = 3.044, used in the (1+z)^4
  term of the Friedmann equation.
* Neutrino masses: Σm_ν = 0.06 eV is carried as ``OMEGA_NU_H2`` and counted
  as part of the (sampled) total matter density Ω_m at late times, while the
  early-time sound-horizon integral treats all neutrinos as relativistic
  (they are, for z > z_drag).  See docs/NEUTRINOS.md for the quantitative
  error budget of this approximation.

The baryon density is **not** fixed anymore: ω_b is a sampled parameter with
a BBN Gaussian prior (see cobaya_configs/*.yaml); ``OMEGA_B_H2_BBN`` /
``OMEGA_B_H2_BBN_SIGMA`` record the adopted BBN calibration
(Schöneberg 2024-type determination) and serve as reference/default values
for plotting scripts.
"""

# Speed of light [km/s]
C_KM_S: float = 299792.458

# --- Radiation sector -------------------------------------------------------
# Photon density for T_CMB = 2.7255 K (Fixsen 2009): Ω_γ h²
OMEGA_GAMMA_H2: float = 2.4728e-5
# Effective number of neutrino species (standard model, Froustey et al. 2020)
N_EFF: float = 3.044
# Total radiation density (photons + relativistic neutrinos): Ω_r h²
OMEGA_R_H2: float = OMEGA_GAMMA_H2 * (1.0 + 0.2271073 * N_EFF)

# --- Neutrino masses ---------------------------------------------------------
# Minimal normal-hierarchy sum of neutrino masses [eV]
SUM_MNU_EV: float = 0.06
# ω_ν = Σm_ν / 93.14 eV (non-relativistic today; part of Ω_m at late times)
OMEGA_NU_H2: float = SUM_MNU_EV / 93.14

# --- Baryons (BBN prior; ω_b is sampled in the MCMC) -------------------------
OMEGA_B_H2_BBN: float = 0.02218
OMEGA_B_H2_BBN_SIGMA: float = 0.00055

# Backwards-compatible alias used as a *default/reference* value only
# (e.g. quick plotting); the inference always samples omega_b.
OMEGA_B_H2: float = OMEGA_B_H2_BBN

# --- USF model ----------------------------------------------------------------
# USF logistic transition width (fixed by construction; see Sec. 2 of the paper)
DELTA: float = 0.25
