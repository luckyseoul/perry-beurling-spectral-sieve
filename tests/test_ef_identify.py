"""EF identification attack: shipped residual = modes + rem (multi-hypothesis)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.ef_identify import (  # noqa: E402
    HYPOTHESES,
    attack_one,
    hypothesis_residual,
    identify_ef,
    model_sanity_identify,
    multi_hypothesis_scan,
    summarize_attack,
)
from pbss.probes import explicit_formula_residual, primes_upto, sample_grid  # noqa: E402


def test_model_identity_peel_exact():
    """q = EF modes ⇒ identification remainder vanishes."""
    for T in (15.0, 25.0, 40.0):
        idn = model_sanity_identify(T=T, n_zeros=8, n_points=2048, degree=4)
        assert idn["model_identity_ok"] is True
        assert idn["frac_l2_remainder"] < 1e-10
        assert idn["triangle_holds"] is True


def test_model_identity_multi_N():
    u = sample_grid(2048)
    for N in (3, 6, 12):
        q, _, _ = explicit_formula_residual(u, T=30.0, n_zeros=N)
        idn = identify_ef(q, u, T=30.0, n_zeros=N, fit_scale=True, degree=3)
        assert idn["frac_l2_remainder"] < 1e-10


def test_hypotheses_build_finite_on_small_primes():
    primes = primes_upto(50_000)
    u = sample_grid(512)
    T = 10.0  # e^10 ~ 22026
    for h in HYPOTHESES:
        q, T_out, meta = hypothesis_residual(
            u, T=T, primes=primes, hypothesis=h, detrend="deg1"
        )
        assert T_out == T
        assert np.all(np.isfinite(q))
        assert q.shape == u.shape
        assert meta["kind"] in ("theta", "psi")


def test_attack_one_returns_metrics():
    primes = primes_upto(100_000)
    u = sample_grid(1024)
    row = attack_one(
        u,
        T=11.0,
        primes=primes,
        hypothesis="H_psi_sqrt",
        n_zeros=10,
        degree=4,
    )
    idn = row["identification"]
    assert 0.0 <= idn["frac_l2_remainder"] <= 2.0  # can exceed 1 if α overfits poorly
    assert idn["triangle_holds"] is True
    assert np.isfinite(idn["corr_q_modes"])
    assert idn["M5_bound_R_d_modes"] > 0.0


def test_multi_hypothesis_scan_shape():
    primes = primes_upto(80_000)
    rows = multi_hypothesis_scan(
        T_values=[9.0, 11.0],
        n_zeros_list=[5, 10],
        primes=primes,
        hypotheses=("H_theta_sqrt", "H_psi_sqrt"),
        n_points=512,
        degree=4,
    )
    # 2 T × 2 N × 2 H = 8
    assert len(rows) == 8
    summary = summarize_attack(rows)
    assert summary["best_hypothesis"] in ("H_theta_sqrt", "H_psi_sqrt")
    assert "hypothesis_means" in summary
    # on tiny x, we may or may not hit sharp block; just ensure structure
    assert summary["status"] in ("sharp_block", "partial_capture")


def test_m5_modes_decay_with_T_in_identification():
    """Mode block M5 majorant shrinks with T (proved-style)."""
    u = sample_grid(2048)
    q, _, _ = explicit_formula_residual(u, T=20.0, n_zeros=6)
    b20 = identify_ef(q, u, T=20.0, n_zeros=6)["M5_bound_R_d_modes"]
    q2, _, _ = explicit_formula_residual(u, T=40.0, n_zeros=6)
    b40 = identify_ef(q2, u, T=40.0, n_zeros=6)["M5_bound_R_d_modes"]
    assert b40 < b20
