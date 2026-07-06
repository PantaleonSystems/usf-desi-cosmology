"""
Topic 1 — External validation of the background in the ΛCDM limit (alpha_q = 0).

Compares E(z), the comoving distance and the angular diameter distance of the
shared background module (consumed by theory, likelihoods and plots) against
astropy.cosmology.FlatLambdaCDM over a grid of (H0, Omega_m).  The astropy
cosmology is built so that its *total* radiation density matches the
pipeline's Omega_r exactly (single effective relativistic fluid), isolating
the numerics from the radiation-content convention.

Tolerances: < 0.05 % on distances, < 1e-8 relative on E(z).
"""

import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM

from theory_module import background as bg
from usf_constants import OMEGA_R_H2
from tests.helpers import make_theory

GRID = [
    (60.0, 0.25),
    (67.36, 0.3153),
    (70.0, 0.30),
    (73.0, 0.35),
    (80.0, 0.40),
]

Z_TEST = np.array([0.01, 0.1, 0.295, 0.51, 0.706, 0.93, 1.317, 1.491, 2.33])


def matched_astropy(H0, Omega_m):
    """FlatLambdaCDM whose total radiation equals the pipeline's Omega_r."""
    h2 = (H0 / 100.0) ** 2
    omega_r_target = OMEGA_R_H2 / h2

    base = FlatLambdaCDM(H0=H0, Om0=Omega_m, Tcmb0=2.7255, Neff=0.0)
    # Choose Neff (massless) so that Ogamma0 + Onu0 == omega_r_target.
    neff = (omega_r_target / base.Ogamma0 - 1.0) / 0.22710731766
    cosmo = FlatLambdaCDM(H0=H0, Om0=Omega_m, Tcmb0=2.7255, Neff=neff)
    assert np.isclose(cosmo.Ogamma0 + cosmo.Onu0, omega_r_target, rtol=1e-10)
    return cosmo


@pytest.mark.parametrize("H0,Omega_m", GRID)
def test_E_of_z_matches_astropy(H0, Omega_m):
    cosmo = matched_astropy(H0, Omega_m)
    ours = bg.E(Z_TEST, H0, Omega_m, alpha_q=0.0)
    ref = cosmo.efunc(Z_TEST)
    np.testing.assert_allclose(ours, ref, rtol=1e-8)


@pytest.mark.parametrize("H0,Omega_m", GRID)
def test_comoving_distance_matches_astropy(H0, Omega_m):
    theo = make_theory(H0, Omega_m, alpha_q=0.0)
    cosmo = matched_astropy(H0, Omega_m)
    for z in Z_TEST:
        ours = theo.get_comoving_distance(float(z))
        ref = cosmo.comoving_distance(float(z)).value
        assert abs(ours / ref - 1.0) < 5e-4, (
            f"D_C({z}) off by {abs(ours / ref - 1) * 100:.4f}% "
            f"(ours={ours:.4f}, astropy={ref:.4f})"
        )


@pytest.mark.parametrize("H0,Omega_m", GRID)
def test_angular_diameter_distance_matches_astropy(H0, Omega_m):
    theo = make_theory(H0, Omega_m, alpha_q=0.0)
    cosmo = matched_astropy(H0, Omega_m)
    ours = np.asarray(theo.get_angular_diameter_distance(Z_TEST))
    ref = cosmo.angular_diameter_distance(Z_TEST).value
    np.testing.assert_allclose(ours, ref, rtol=5e-4)


def test_lcdm_limit_alpha_q_zero():
    """alpha_q = 0 must reduce the USF exactly to ΛCDM for any z_trans."""
    z = np.linspace(0.0, 10.0, 101)
    ref = bg.E(z, 70.0, 0.31, alpha_q=0.0, z_trans=2.0)
    for zt in (0.7, 4.66, 15.0):
        alt = bg.E(z, 70.0, 0.31, alpha_q=0.0, z_trans=zt)
        np.testing.assert_allclose(alt, ref, rtol=1e-14)


def test_E_normalised_at_zero():
    """Flatness bookkeeping: E(0) = 1 also with the USF term switched on."""
    for aq, zt in [(-0.1, 4.66), (0.05, 1.0), (-0.25, 12.0)]:
        e0 = float(bg.E(0.0, 72.0, 0.37, alpha_q=aq, z_trans=zt))
        assert abs(e0 - 1.0) < 1e-12
