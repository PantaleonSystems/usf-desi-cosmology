"""
Topic 5 — Validation of the Planck 2018 CMB distance-prior likelihood.

The theory→data mapping must reproduce the compression convention of
Chen, Huang & Wang (2019): evaluated at the Planck 2018 TT,TE,EE+lowE
posterior mean, the prediction must land on the published data vector
within a fraction of its 1σ.
"""

import numpy as np
import pytest

# Planck 2018 TT,TE,EE+lowE posterior means (base_plikHM_TTTEEE_lowl_lowE)
PLANCK_MEAN = dict(H0=67.27, Omega_m=0.3166, omega_b=0.02236,
                   alpha_q=0.0, z_trans=2.0)


@pytest.fixture(scope="module")
def like():
    from theory_module.cmb_distance_priors import CMBDistancePriors

    lk = object.__new__(CMBDistancePriors)
    lk.initialize()
    return lk


def test_prediction_reproduces_chw19_at_planck_mean(like):
    R, l_A = like._predict(**PLANCK_MEAN)
    # data: R = 1.7502 ± 0.0046, l_A = 301.471 ± 0.090
    assert abs(R - 1.7502) < 0.5 * 0.0046, f"R = {R:.4f}"
    assert abs(l_A - 301.471) < 0.6 * 0.090, f"l_A = {l_A:.3f}"


def test_chi2_small_at_planck_mean(like):
    chi2 = -2.0 * like.logp(**PLANCK_MEAN)
    print(f"\n[CMBdist] chi² at Planck 2018 mean = {chi2:.3f}")
    assert chi2 < 1.0


def test_lcdm_limit_independent_of_ztrans(like):
    p1 = like.logp(**{**PLANCK_MEAN, "z_trans": 0.7})
    p2 = like.logp(**{**PLANCK_MEAN, "z_trans": 4.9})
    assert abs(p1 - p2) < 1e-12


def test_usf_term_affects_distance_to_last_scattering(like):
    """Referee 1, Comment 2: the low-z modification must propagate to D(z*)."""
    chi2_lcdm = -2.0 * like.logp(**PLANCK_MEAN)
    chi2_usf = -2.0 * like.logp(**{**PLANCK_MEAN,
                                   "alpha_q": -0.1, "z_trans": 4.66})
    # A -0.1 geometric screening changes the line-of-sight integral: the
    # likelihood must register a substantially different chi².
    assert abs(chi2_usf - chi2_lcdm) > 1.0, (
        f"chi² insensitive to alpha_q: {chi2_lcdm:.3f} vs {chi2_usf:.3f}"
    )


def test_covariance_positive_definite():
    from theory_module.cmb_distance_priors import _COV

    eigvals = np.linalg.eigvalsh(_COV)
    assert np.all(eigvals > 0)
