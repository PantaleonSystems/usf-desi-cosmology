"""
Shared helpers for the validation suite.

Provides bare instances of the Cobaya components (theory / likelihood)
without going through Cobaya's full model machinery, so that the physics
methods can be unit-tested in isolation.
"""

import os

from tests.conftest import REPO_ROOT


def make_theory(H0, Omega_m, alpha_q=0.0, z_trans=2.0, omega_b=0.02218, **extra):
    """Return a theory instance with parameters set, bypassing Cobaya init."""
    from theory_module.usf_background import USF

    theo = object.__new__(USF)
    theo.H0 = H0
    theo.Omega_m = Omega_m
    theo.omega_b = omega_b
    theo.alpha_q = alpha_q
    theo.z_trans = z_trans
    for key, val in extra.items():
        setattr(theo, key, val)
    return theo


class StubProvider:
    """Minimal stand-in for Cobaya's provider, backed by a theory instance."""

    def __init__(self, theory):
        self._theory = theory

    def __getattr__(self, name):
        if name.startswith("get_"):
            return getattr(self._theory, name)
        raise AttributeError(name)


def make_bao_likelihood(theory=None):
    """Return an initialized BAO likelihood instance, bypassing Cobaya init.

    If a theory instance is given, a StubProvider is attached so that
    ``logp()`` can be called in isolation.
    """
    from theory_module.bao_desi_likelihood import BAODESILikelihood

    like = object.__new__(BAODESILikelihood)
    like.data_path = os.path.join(
        REPO_ROOT, "cobaya_packages", "data", "bao_data", "desi_bao_dr2"
    )
    like.initialize()
    if theory is not None:
        like.provider = StubProvider(theory)
    return like
