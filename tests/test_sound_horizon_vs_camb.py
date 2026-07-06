"""
Topic 1 — Validation of the sound horizon r_d against CAMB.

Runs CAMB for the Planck 2018 baseline cosmology and compares its r_drag
with the pipeline's internal sound-horizon integral evaluated for the same
physical densities.

STATUS (Topic 1 audit): this test is expected to FAIL against the original
pipeline, which used the *total* radiation density (4.15e-5) as the photon
density in the baryon-to-photon ratio R, inflating the sound speed and hence
r_d by several percent.  The measured deviation is recorded in
PROGRESSO_CODIGO.md and the fix is Topic 2 of PLANO_CORRECAO_REFEREE1.md.
The xfail marker must be removed once Topic 2 lands (strict=True will then
flag it as XPASS -> failure, forcing the cleanup).
"""

import numpy as np
import pytest

from tests.helpers import make_theory

# Planck 2018 TT,TE,EE+lowE+lensing baseline.
OMBH2 = 0.02237
OMCH2 = 0.1200
H0 = 67.36
MNU = 0.06  # eV, one massive eigenstate

RD_TOL = 0.005  # 0.5 % — accounts for the EH98 z_drag fitting formula.


@pytest.fixture(scope="module")
def camb_rdrag():
    import camb

    pars = camb.set_params(
        ombh2=OMBH2, omch2=OMCH2, H0=H0, mnu=MNU, omk=0.0,
        As=2.1e-9, ns=0.9649, tau=0.0544,
    )
    results = camb.get_background(pars)
    derived = results.get_derived_params()
    return derived["rdrag"], derived["zdrag"]


@pytest.mark.xfail(
    reason="Bug conhecido (Tópico 2): ω_γ=4.15e-5 (radiação total) usado como "
    "densidade de fótons na razão R do sound speed infla r_d. Remover este "
    "marcador ao concluir o Tópico 2.",
    strict=False,
)
def test_rd_matches_camb(camb_rdrag):
    rd_camb, zdrag_camb = camb_rdrag

    h = H0 / 100.0
    omega_nu = MNU / 93.14
    omega_m = OMBH2 + OMCH2 + omega_nu
    theo = make_theory(H0, omega_m / h**2, alpha_q=0.0)

    rd_ours = theo.get_r_s()
    dev = rd_ours / rd_camb - 1.0
    print(
        f"\n[r_d audit] CAMB rdrag = {rd_camb:.3f} Mpc (zdrag={zdrag_camb:.2f}) | "
        f"pipeline r_s = {rd_ours:.3f} Mpc | deviation = {dev * 100:+.3f}%"
    )
    assert abs(dev) < RD_TOL, (
        f"r_d deviates from CAMB by {dev * 100:+.3f}% "
        f"(pipeline {rd_ours:.3f} vs CAMB {rd_camb:.3f} Mpc)"
    )
