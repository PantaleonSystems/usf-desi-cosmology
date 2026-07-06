"""
Frequentist calibration of the USF vs ΛCDM Δχ² under the null hypothesis.

Because z_trans is not identifiable when α_q = 0, the naive |α_q|/σ cannot
be read as a Gaussian significance (Referee 1, Comment 8; look-elsewhere
effect).  This script builds the empirical distribution of

    Δχ² = χ²_min(ΛCDM) − χ²_min(USF)

under the ΛCDM null: it generates N mock realisations of the DESI DR2 BAO
vector and of the Pantheon+ magnitudes (full STAT+SYS covariance, M
marginalised analytically), fits both models to every mock — PROFILING over
a fixed z_trans grid, never fixing it — and stores the Δχ² sample.  The
observed Δχ² is then converted into a calibrated p-value / significance.

Usage (from the repository root; heavy — production numbers in the plan
call for N ≥ 300, ideally 500–1000):
    .venv/bin/python analysis_scripts/run_mocks.py --n-mocks 300 \
        --h0 68.4 --omega-m 0.30 --omega-b 0.02218 \
        [--with-sh0es] [--out results/mocks_null_dchi2.npz]

A quick smoke test: --n-mocks 2 --n-ztrans 4
"""

import argparse
import os
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from theory_module import background as bg  # noqa: E402
from usf_constants import (  # noqa: E402
    C_KM_S,
    OMEGA_B_H2_BBN,
    OMEGA_B_H2_BBN_SIGMA,
)

BAO_DIR = "cobaya_packages/data/bao_data/desi_bao_dr2"
SN_DAT = "cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES.dat"
SN_COV = "cobaya_packages/data/sn_data/PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov"
Z_CUT = 0.01
SH0ES_SIGMA = 1.04


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_bao():
    rows = []
    with open(os.path.join(BAO_DIR, "desi_gaussian_bao_ALL_GCcomb_mean.txt")) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            rows.append((float(p[0]), float(p[1]), p[2]))
    z = np.array([r[0] for r in rows])
    types = [r[2] for r in rows]
    cov = np.loadtxt(os.path.join(BAO_DIR, "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
    return z, types, cov


def load_sn():
    import pandas as pd

    df = pd.read_csv(SN_DAT, comment="#", sep=r"\s+")
    mask = (df["zHD"] > Z_CUT).values

    with open(SN_COV) as f:
        n = int(f.readline())
        cov = np.fromfile(f, sep=" ", count=n * n).reshape(n, n)
    cov = cov[np.ix_(mask, mask)]

    return (df["zHD"].values[mask], df["zHEL"].values[mask],
            df["m_b_corr"].values[mask], cov)


def marginalized_invcov(cov):
    """Project the constant (absolute-magnitude) mode out of C^-1 —
    identical to cobaya's sn.pantheonplus `_marginalize_abs_mag`."""
    invcov = np.linalg.inv(cov)
    deriv = np.ones((invcov.shape[0], 1))
    derivp = invcov @ deriv
    fisher = deriv.T @ derivp
    return invcov - derivp @ np.linalg.solve(fisher, derivp.T)


# ---------------------------------------------------------------------------
# Model predictions
# ---------------------------------------------------------------------------

def bao_prediction(z, types, H0, Om, wb, aq, zt):
    omega_m = Om * (H0 / 100.0) ** 2
    rs = bg.r_drag(wb, omega_m)
    pred = np.empty(len(z))
    for i, (zi, t) in enumerate(zip(z, types)):
        if t == "DM_over_rs":
            pred[i] = bg.comoving_distance(zi, H0, Om, aq, zt) / rs
        elif t == "DH_over_rs":
            pred[i] = C_KM_S / (bg.hubble(zi, H0, Om, aq, zt) * rs)
        elif t == "DV_over_rs":
            DM = bg.comoving_distance(zi, H0, Om, aq, zt)
            DH = C_KM_S / bg.hubble(zi, H0, Om, aq, zt)
            pred[i] = (zi * DM ** 2 * DH) ** (1.0 / 3.0) / rs
        else:
            raise ValueError(t)
    return pred


def sn_mu(z_hd, z_hel, H0, Om, aq, zt):
    dc = bg.comoving_distance(z_hd, H0, Om, aq, zt)
    return 5.0 * np.log10((1.0 + z_hel) * dc) + 25.0


class NullChi2:
    """χ² machinery shared by data and mocks."""

    def __init__(self, with_sh0es=False):
        self.z_bao, self.t_bao, self.cov_bao = load_bao()
        self.icov_bao = np.linalg.inv(self.cov_bao)
        self.chol_bao = np.linalg.cholesky(self.cov_bao)

        self.z_hd, self.z_hel, self.mb, cov_sn = load_sn()
        self.icov_sn = marginalized_invcov(cov_sn)
        self.chol_sn = np.linalg.cholesky(cov_sn)

        self.with_sh0es = with_sh0es

    def chi2(self, theta, bao_vec, sn_vec, sh0es_obs, aq, zt):
        H0, Om, wb = theta
        if not (55 < H0 < 85 and 0.15 < Om < 0.45 and 0.015 < wb < 0.03):
            return 1e10
        r_bao = bao_prediction(self.z_bao, self.t_bao, H0, Om, wb, aq, zt) - bao_vec
        chi2 = r_bao @ self.icov_bao @ r_bao
        r_sn = sn_vec - sn_mu(self.z_hd, self.z_hel, H0, Om, aq, zt)
        chi2 += r_sn @ self.icov_sn @ r_sn
        chi2 += ((wb - OMEGA_B_H2_BBN) / OMEGA_B_H2_BBN_SIGMA) ** 2
        if self.with_sh0es:
            chi2 += ((H0 - sh0es_obs) / SH0ES_SIGMA) ** 2
        return chi2

    def fit_lcdm(self, bao_vec, sn_vec, sh0es_obs, x0):
        res = optimize.minimize(
            self.chi2, x0, args=(bao_vec, sn_vec, sh0es_obs, 0.0, 2.0),
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-4, maxiter=2000))
        return res.fun, res.x

    def fit_usf_profiled(self, bao_vec, sn_vec, sh0es_obs, x0, ztrans_grid):
        """Profile over z_trans: minimise (H0, Ωm, ωb, α_q) at each grid node."""
        best = np.inf
        best_zt = None
        for zt in ztrans_grid:
            def chi2_usf(p):
                aq = p[3]
                if not -0.35 < aq < 0.35:
                    return 1e10
                return self.chi2(p[:3], bao_vec, sn_vec, sh0es_obs, aq, zt)

            res = optimize.minimize(
                chi2_usf, np.append(x0, 0.0), method="Nelder-Mead",
                options=dict(xatol=1e-4, fatol=1e-4, maxiter=3000))
            if res.fun < best:
                best, best_zt = res.fun, zt
        return best, best_zt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mocks", type=int, default=300)
    ap.add_argument("--h0", type=float, default=68.4)
    ap.add_argument("--omega-m", type=float, default=0.30)
    ap.add_argument("--omega-b", type=float, default=OMEGA_B_H2_BBN)
    ap.add_argument("--n-ztrans", type=int, default=12,
                    help="nodes of the z_trans profiling grid")
    ap.add_argument("--ztrans-max", type=float, default=20.0)
    ap.add_argument("--with-sh0es", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="results/mocks_null_dchi2.npz")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    machine = NullChi2(with_sh0es=args.with_sh0es)

    truth = (args.h0, args.omega_m, args.omega_b)
    bao_true = bao_prediction(machine.z_bao, machine.t_bao, *truth, 0.0, 2.0)
    mu_true = sn_mu(machine.z_hd, machine.z_hel, args.h0, args.omega_m, 0.0, 2.0)
    M_TRUE = -19.3  # arbitrary: marginalised away

    ztrans_grid = np.geomspace(0.5, args.ztrans_max, args.n_ztrans)
    x0 = np.array(truth)

    dchi2 = np.empty(args.n_mocks)
    zt_best = np.empty(args.n_mocks)
    for i in range(args.n_mocks):
        bao_mock = bao_true + machine.chol_bao @ rng.standard_normal(len(bao_true))
        sn_mock = (mu_true + M_TRUE
                   + machine.chol_sn @ rng.standard_normal(len(mu_true)))
        sh0es_obs = rng.normal(args.h0, SH0ES_SIGMA) if args.with_sh0es else None

        chi2_l, xfit = machine.fit_lcdm(bao_mock, sn_mock, sh0es_obs, x0)
        chi2_u, ztb = machine.fit_usf_profiled(bao_mock, sn_mock, sh0es_obs,
                                               xfit, ztrans_grid)
        dchi2[i] = chi2_l - chi2_u
        zt_best[i] = ztb
        print(f"mock {i + 1:4d}/{args.n_mocks}: "
              f"chi2_LCDM={chi2_l:9.2f} chi2_USF={chi2_u:9.2f} "
              f"dchi2={dchi2[i]:7.3f} (zt*={ztb:.2f})", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez(args.out, dchi2=dchi2, zt_best=zt_best,
             truth=np.array(truth), with_sh0es=args.with_sh0es,
             ztrans_grid=ztrans_grid, seed=args.seed)

    q = np.percentile(dchi2, [50, 90, 95, 99])
    print("\nNull Δχ² distribution:")
    print(f"  median = {q[0]:.2f} | 90% = {q[1]:.2f} | 95% = {q[2]:.2f} "
          f"| 99% = {q[3]:.2f}")
    print(f"Saved {args.out}.  Calibrated p-value of an observed Δχ²_obs: "
          "p = mean(dchi2 >= Δχ²_obs).")


if __name__ == "__main__":
    main()
