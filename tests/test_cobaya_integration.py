"""
Topic 9 — End-to-end Cobaya integration test.

Builds the full Cobaya model from the production YAML (theory provider +
DESI BAO likelihood via requirements + Pantheon+ SN likelihood) and
evaluates the posterior at a reference point.  Validates the whole
provider wiring exactly as `cobaya-run` will exercise it.
"""

import os

import numpy as np
import pytest
import yaml

from tests.conftest import REPO_ROOT

REF_POINT = dict(H0=68.4, Omega_m=0.30, omega_b=0.02218,
                 alpha_q=0.0, z_trans=2.0)


@pytest.fixture(scope="module")
def model():
    from cobaya.model import get_model

    cfg_path = os.path.join(REPO_ROOT, "cobaya_configs", "usf_bao_sn.yaml")
    with open(cfg_path) as f:
        info = yaml.safe_load(f)
    info.pop("output", None)
    info.pop("sampler", None)
    info["packages_path"] = os.path.join(REPO_ROOT, "cobaya_packages")
    # python_path entries are relative to the repo root
    cwd = os.getcwd()
    os.chdir(REPO_ROOT)
    try:
        yield get_model(info)
    finally:
        os.chdir(cwd)


def test_posterior_finite_at_reference_point(model):
    post = model.logposterior(REF_POINT)
    assert np.isfinite(post.logpost)


def test_bao_chi2_reasonable_via_provider(model):
    post = model.logposterior(REF_POINT)
    loglikes = dict(zip([str(n) for n in model.likelihood], post.loglikes))
    bao_key = next(k for k in loglikes if "BAODESI" in k or "bao" in k.lower())
    chi2_bao = -2.0 * loglikes[bao_key]
    # 13 BAO points; the reference point is near the published best fit.
    assert chi2_bao < 40.0, f"BAO chi² through provider = {chi2_bao:.1f}"


def test_lcdm_limit_through_full_model(model):
    """alpha_q = 0 posterior must not depend on z_trans (nested model)."""
    p1 = model.logposterior({**REF_POINT, "z_trans": 0.7}).logpost
    p2 = model.logposterior({**REF_POINT, "z_trans": 4.9}).logpost
    assert abs(p1 - p2) < 1e-8
