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

PARAMS = dict(H0=72.1, Omega_m=0.373, alpha_q=-0.0998, z_trans=4.66,
              omega_b=0.02218)
Z = np.array([0.0, 0.1, 0.51, 1.0, 2.33, 4.66, 8.0])


def test_theory_H_identical_to_background():
    theo = make_theory(**PARAMS)
    ours = theo.get_H(Z)
    ref = bg.hubble(Z, PARAMS["H0"], PARAMS["Omega_m"],
                    PARAMS["alpha_q"], PARAMS["z_trans"])
    np.testing.assert_array_equal(ours, ref)


def test_theory_distances_identical_to_background():
    theo = make_theory(**PARAMS)
    np.testing.assert_array_equal(
        np.asarray(theo.get_angular_diameter_distance(Z[1:])),
        bg.angular_diameter_distance(Z[1:], PARAMS["H0"], PARAMS["Omega_m"],
                                     PARAMS["alpha_q"], PARAMS["z_trans"]))
    assert theo.get_comoving_distance(1.0) == bg.comoving_distance(
        1.0, PARAMS["H0"], PARAMS["Omega_m"],
        PARAMS["alpha_q"], PARAMS["z_trans"])


def test_theory_rd_identical_to_background():
    theo = make_theory(**PARAMS)
    omega_m = PARAMS["Omega_m"] * (PARAMS["H0"] / 100.0) ** 2
    assert theo.get_r_s() == bg.r_drag(PARAMS["omega_b"], omega_m)


def test_omega_usf_low_z_limit():
    """Ω_USF(0) → α_q/(1 + 1/z_trans²) — the low-z plateau quoted in the paper."""
    aq, zt = -0.0998, 4.66
    expected = aq / (1.0 + 1.0 / zt ** 2)
    got = float(bg.omega_usf(0.0, aq, zt))
    # supp0 factor deviates from 1 by exp(-z_trans/delta) ≈ 8e-9 — negligible
    assert abs(got / expected - 1.0) < 1e-6


def test_omega_usf_high_z_suppression():
    """Ω_USF must vanish for z ≫ z_trans (CMB-era safety)."""
    vals = bg.omega_usf(np.array([8.0, 10.0, 1100.0]), -0.1, 4.66)
    assert np.all(np.abs(vals) < 1e-4)
