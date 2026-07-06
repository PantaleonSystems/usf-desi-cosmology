"""
Topic 1 — Isolated validation of the DESI DR2 BAO likelihood.

Evaluates the custom BAO likelihood at the published DESI DR2 flat-ΛCDM
best fit (Omega_m ≈ 0.2975; H0 ≈ 68.5 for BAO+BBN calibration).  A correct
implementation must yield a chi² compatible with the published fit quality
(chi² of order the number of data points, no systematic offset).

STATUS (Topic 1 audit): expected to FAIL with the original pipeline — the
inflated r_d (see test_sound_horizon_vs_camb.py) shifts every D/r_d
prediction coherently, producing a huge chi² at the published best fit.
The xfail marker must be removed once Topic 2 lands.
"""

import numpy as np
import pytest

from tests.helpers import make_bao_likelihood

# DESI DR2 (Abdul Karim et al. 2025) flat ΛCDM, BAO(+BBN):
DESI_OMEGA_M = 0.2975
DESI_H0 = 68.5


@pytest.fixture(scope="module")
def bao_like():
    return make_bao_likelihood()


def test_data_vector_loaded(bao_like):
    n = len(bao_like.z_eff)
    assert n >= 10, f"unexpectedly small DESI data vector ({n} points)"
    assert bao_like.cov_inv.shape == (n, n)
    assert set(bao_like.obs_type) <= {"DM_over_rs", "DH_over_rs", "DV_over_rs"}


@pytest.mark.xfail(
    reason="Bug conhecido (Tópico 2): r_d inflado desloca todas as predições "
    "D/r_d; chi² no best-fit publicado do DESI sai ordens de grandeza acima "
    "do aceitável. Remover este marcador ao concluir o Tópico 2.",
    strict=False,
)
def test_chi2_at_desi_bestfit(bao_like):
    logp = bao_like.logp(
        H0=DESI_H0, Omega_m=DESI_OMEGA_M, alpha_q=0.0, z_trans=2.0
    )
    chi2 = -2.0 * logp
    ndata = len(bao_like.z_eff)
    print(f"\n[BAO audit] chi² at DESI DR2 best fit = {chi2:.2f} for {ndata} points")
    assert chi2 < 2.5 * ndata, (
        f"chi²={chi2:.1f} at the published DESI best fit "
        f"({ndata} points) — BAO likelihood or r_d implementation inconsistent"
    )
