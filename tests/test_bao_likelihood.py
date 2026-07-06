"""
Topic 1 — Isolated validation of the DESI DR2 BAO likelihood.

Evaluates the custom BAO likelihood at the published DESI DR2 flat-ΛCDM
best fit (Omega_m ≈ 0.2975; H0 ≈ 68.5 for BAO+BBN calibration).  A correct
implementation must yield a chi² compatible with the published fit quality
(chi² of order the number of data points, no systematic offset).

History: the Topic 1 audit measured chi² = 227.5 for 13 points at the
published best fit, caused by the +4.4% inflated r_d (fixed in Topic 2).
This test now guards the fix.
"""

import numpy as np
import pytest

from tests.helpers import make_bao_likelihood, make_theory

# DESI DR2 (Abdul Karim et al. 2025) flat ΛCDM, BAO(+BBN):
DESI_OMEGA_M = 0.2975
DESI_H0 = 68.5


@pytest.fixture(scope="module")
def bao_like():
    theory = make_theory(DESI_H0, DESI_OMEGA_M, alpha_q=0.0,
                         z_trans=2.0, omega_b=0.02218)
    return make_bao_likelihood(theory)


def test_data_vector_loaded(bao_like):
    n = len(bao_like.z_eff)
    assert n >= 10, f"unexpectedly small DESI data vector ({n} points)"
    assert bao_like.cov_inv.shape == (n, n)
    assert set(bao_like.obs_type) <= {"DM_over_rs", "DH_over_rs", "DV_over_rs"}


def test_chi2_at_desi_bestfit(bao_like):
    logp = bao_like.logp()
    chi2 = -2.0 * logp
    ndata = len(bao_like.z_eff)
    print(f"\n[BAO audit] chi² at DESI DR2 best fit = {chi2:.2f} for {ndata} points")
    assert chi2 < 2.5 * ndata, (
        f"chi²={chi2:.1f} at the published DESI best fit "
        f"({ndata} points) — BAO likelihood or r_d implementation inconsistent"
    )
