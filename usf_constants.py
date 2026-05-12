"""
Shared physical constants and fixed cosmological parameters for the USF analysis.

All modules (feu_background, bao_desi_likelihood, plot scripts) import from
here so the MCMC engine and the visualisation layer are guaranteed to use
identical numerical values.
"""

# Speed of light [km/s]
C_KM_S: float = 299792.458

# Physical baryon and photon densities (Planck 2018)
OMEGA_B_H2: float = 0.02237       # Ω_b h²
OMEGA_GAMMA_H2: float = 4.15e-5   # Ω_γ h²

# Sound-horizon drag redshift fixed to Planck 2018 best-fit for computational
# efficiency; the narrow priors used in our MCMC keep the true z_drag within
# ~0.5 of this value, preserving high accuracy (see Sec. 3 of the paper).
Z_DRAG: float = 1059.94

# USF logistic transition width (fixed by construction; see Sec. 2)
DELTA: float = 0.25
