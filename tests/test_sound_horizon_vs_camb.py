"""
Topics 1–2 — Validation of the sound horizon r_d against CAMB.

Runs CAMB for the Planck 2018 baseline cosmology and compares its r_drag
with the pipeline's sound-horizon (numerical integral + CAMB-calibrated
z_drag).  Also cross-checks against the independent Aubourg et al. (2015)
fitting formula.

History: the Topic 1 audit (2026-07-06) measured a +4.401% deviation in the
original pipeline, caused by the total radiation density (4.15e-5) being
used as the photon density in the baryon-to-photon ratio R.  Fixed in
Topic 2 (see PROGRESSO_CODIGO.md); this test now guards the fix.
"""

import numpy as np
import pytest

from theory_module import background as bg
from tests.helpers import make_theory

# Planck 2018 TT,TE,EE+lowE+lensing baseline.
OMBH2 = 0.02237
OMCH2 = 0.1200
H0 = 67.36
MNU = 0.06  # eV, one massive eigenstate

RD_TOL = 0.003  # 0.3% on r_d vs CAMB
GRID = [
    (0.02218, 0.1200, 67.0),
    (0.02237, 0.1200, 67.36),
    (0.02218, 0.1432, 70.0),
    (0.02273, 0.1100, 68.0),
]


def camb_rdrag(ombh2, omch2, h0):
    import camb

    pars = camb.set_params(
        ombh2=ombh2, omch2=omch2, H0=h0, mnu=MNU, omk=0.0,
        As=2.1e-9, ns=0.9649, tau=0.0544,
    )
    return camb.get_background(pars).get_derived_params()


def test_rd_matches_camb_planck_baseline():
    derived = camb_rdrag(OMBH2, OMCH2, H0)
    rd_camb, zdrag_camb = derived["rdrag"], derived["zdrag"]

    omega_m = OMBH2 + OMCH2 + MNU / 93.14
    theo = make_theory(H0, omega_m / (H0 / 100.0) ** 2,
                       alpha_q=0.0, omega_b=OMBH2)

    rd_ours = theo.get_r_s()
    dev = rd_ours / rd_camb - 1.0
    print(
        f"\n[r_d] CAMB rdrag = {rd_camb:.3f} Mpc (zdrag={zdrag_camb:.2f}) | "
        f"pipeline r_d = {rd_ours:.3f} Mpc | deviation = {dev * 100:+.3f}%"
    )
    assert abs(dev) < RD_TOL


@pytest.mark.parametrize("ombh2,omch2,h0", GRID)
def test_rd_matches_camb_grid(ombh2, omch2, h0):
    derived = camb_rdrag(ombh2, omch2, h0)
    omega_m = ombh2 + omch2 + MNU / 93.14
    rd_ours = bg.r_drag(ombh2, omega_m)
    dev = rd_ours / derived["rdrag"] - 1.0
    assert abs(dev) < RD_TOL, (
        f"r_d off by {dev * 100:+.3f}% at ω_b={ombh2}, ω_c={omch2}"
    )


def test_rd_crosscheck_aubourg15():
    """Numerical integral vs Aubourg+15 fitting formula (independent check)."""
    for omega_b, omega_m in [(0.02218, 0.1400), (0.02237, 0.14301),
                             (0.02273, 0.1330)]:
        rd_int = bg.r_drag(omega_b, omega_m)
        rd_fit = bg.r_drag_aubourg15(omega_b, omega_m)
        assert abs(rd_int / rd_fit - 1.0) < 0.003, (
            f"integral {rd_int:.3f} vs Aubourg15 {rd_fit:.3f} at "
            f"ω_b={omega_b}, ω_m={omega_m}"
        )


def test_rd_independent_of_H0_and_usf():
    """r_d must depend only on physical densities, not on H0 or USF params."""
    rd_ref = bg.r_drag(0.02218, 0.1432)
    # r_drag takes only (omega_b, omega_m): recompute via theory at two H0
    t1 = make_theory(60.0, 0.1432 / 0.36, alpha_q=-0.1, z_trans=4.66)
    t2 = make_theory(80.0, 0.1432 / 0.64, alpha_q=0.0)
    assert abs(t1.get_r_s() / rd_ref - 1.0) < 1e-10
    assert abs(t2.get_r_s() / rd_ref - 1.0) < 1e-10


def test_eh98_discrepancy_documented():
    """EH98 z_drag ≈ 1021 vs CAMB ≈ 1060 — must NOT be used in production."""
    zd_eh98 = bg.z_drag_eh98(0.02237, 0.1432)
    zd_camb = bg.z_drag(0.02237, 0.1432)
    assert 1010 < zd_eh98 < 1035
    assert 1055 < zd_camb < 1065
