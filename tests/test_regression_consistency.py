"""
Topic 9 — Regression: every consumer evaluates identical physics.

The theory class (Cobaya provider) must be a thin wrapper over
theory_module.background: bit-identical results for the same parameters.
This is what guarantees that the sampler, the likelihoods and the plotting
scripts (which import theory_module.background directly) can never drift
apart again.
"""

import numpy as np

from theory_module import background as bg
from tests.helpers import make_theory

PARAMS = dict(H0=72.1, Omega_m=0.373, Omega_TW=0.0998, n=1.5,
              omega_b=0.02218)
Z = np.array([0.0, 0.1, 0.51, 1.0, 2.33, 4.66, 8.0])


def test_theory_H_identical_to_background():
    theo = make_theory(**PARAMS)
    ours = theo.get_H(Z)
    ref = bg.hubble(Z, PARAMS["H0"], PARAMS["Omega_m"],
                    PARAMS["Omega_TW"], PARAMS["n"])
    np.testing.assert_array_equal(ours, ref)


def test_theory_distances_identical_to_background():
    theo = make_theory(**PARAMS)
    np.testing.assert_array_equal(
        np.asarray(theo.get_angular_diameter_distance(Z[1:])),
        bg.angular_diameter_distance(Z[1:], PARAMS["H0"], PARAMS["Omega_m"],
                                     PARAMS["Omega_TW"], PARAMS["n"]))
    assert theo.get_comoving_distance(1.0) == bg.comoving_distance(
        1.0, PARAMS["H0"], PARAMS["Omega_m"],
        PARAMS["Omega_TW"], PARAMS["n"])


def test_theory_rd_identical_to_background():
    theo = make_theory(**PARAMS)
    omega_m = PARAMS["Omega_m"] * (PARAMS["H0"] / 100.0) ** 2
    assert theo.get_r_s() == bg.r_drag(PARAMS["omega_b"], omega_m)


def test_omega_tw_normalisation():
    """The USF-26n density term equals Omega_TW at z=0."""
    got = float(bg.omega_tw(0.0, 0.0998, 1.5))
    assert abs(got - 0.0998) < 1e-14


def test_omega_tw_scales_with_redshift():
    """The USF-26n term follows its declared power-law scaling."""
    vals = bg.omega_tw(np.array([0.0, 1.0, 2.0]), 0.1, 2.0)
    np.testing.assert_allclose(vals, [0.1, 0.4, 0.9])
