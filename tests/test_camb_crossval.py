"""
Topic 9 — External cross-validation of the unified background module vs CAMB.

Pre-condition for any production run: distances and r_drag from
theory_module.background must match CAMB (full Boltzmann code with
Σm_ν = 0.06 eV) for a set of standard cosmologies covering the prior volume.

Documented tolerances (PLANO_CORRECAO_REFEREE1.md, Tópico 9):
    < 0.1 % on D_M, D_A, H(z) up to z = 5;
    < 0.2 % on r_drag.
Measured headroom (2026-07-06): ≤ 0.012 % and 0.006 % respectively —
see docs/NEUTRINOS.md.
"""

import numpy as np
import pytest

from theory_module import background as bg

MNU = 0.06
Z = np.array([0.1, 0.51, 1.0, 1.491, 2.33, 5.0])

# (ombh2, omch2, H0): Planck 2018, DESI DR2-like, and ±2σ-ish corners.
COSMOLOGIES = [
    (0.02237, 0.1200, 67.36),
    (0.02218, 0.1155, 68.5),
    (0.02160, 0.1300, 65.0),
    (0.02280, 0.1050, 74.0),
]

DIST_TOL = 1e-3
RD_TOL = 2e-3


@pytest.mark.parametrize("ombh2,omch2,h0", COSMOLOGIES)
def test_distances_and_rd_vs_camb(ombh2, omch2, h0):
    import camb

    pars = camb.set_params(ombh2=ombh2, omch2=omch2, H0=h0, mnu=MNU,
                           omk=0.0, As=2.1e-9, ns=0.9649, tau=0.0544)
    res = camb.get_background(pars)

    omega_m = ombh2 + omch2 + MNU / 93.14
    Om = omega_m / (h0 / 100.0) ** 2

    # Comoving distance
    dm_camb = res.comoving_radial_distance(Z)
    dm_ours = bg.comoving_distance(Z, h0, Om)
    np.testing.assert_allclose(dm_ours, dm_camb, rtol=DIST_TOL)

    # Angular diameter distance
    da_camb = res.angular_diameter_distance(Z)
    da_ours = bg.angular_diameter_distance(Z, h0, Om)
    np.testing.assert_allclose(da_ours, da_camb, rtol=DIST_TOL)

    # Hubble rate
    hz_camb = res.hubble_parameter(Z)
    hz_ours = bg.hubble(Z, h0, Om)
    np.testing.assert_allclose(hz_ours, hz_camb, rtol=DIST_TOL)

    # Sound horizon at the drag epoch
    rd_camb = res.get_derived_params()["rdrag"]
    rd_ours = bg.r_drag(ombh2, omega_m)
    assert abs(rd_ours / rd_camb - 1.0) < RD_TOL, (
        f"r_d off by {(rd_ours / rd_camb - 1) * 100:+.3f}%"
    )
