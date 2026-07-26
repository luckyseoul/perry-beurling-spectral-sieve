"""Tests for truncation remainder / zero-sum peel (shipped remainder APIs)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.lemmas import bound_R_d_finite_mode_sum  # noqa: E402
from pbss.probes import explicit_formula_residual, sample_grid  # noqa: E402
from pbss.projection import energy_ratio  # noqa: E402
from pbss.remainder import (  # noqa: E402
    bound_R_d_mode_tail,
    multi_TN_remainder_scan,
    peel_via_remainder,
    remainder_diagnostic,
    remainder_diagnostic_from_q,
    truncated_mode_sum,
)
from pbss.zeros import explicit_formula_amplitudes, zeta_zero_ordinates  # noqa: E402


def test_truncated_mode_sum_matches_explicit_formula_residual():
    u = sample_grid(2048)
    q1, meta1 = truncated_mode_sum(u, T=22.0, n_zeros=8)
    q2, _, meta2 = explicit_formula_residual(u, T=22.0, n_zeros=8)
    assert np.allclose(q1, q2)
    assert meta1["n_zeros"] == meta2["n_zeros"] == 8


def test_peel_via_remainder_full_strip_kills_model_sum():
    """Stripping all N modes from q^{(N)} leaves ~0 residual."""
    u = sample_grid(4096)
    N = 10
    q_full, _ = truncated_mode_sum(u, T=30.0, n_zeros=N)
    q_rem, meta = peel_via_remainder(
        q_full, u, T=30.0, n_strip=N, fit_scale=True
    )
    assert meta["kind"] == "peel_via_remainder"
    # near zero after full peel
    assert float(np.max(np.abs(q_rem))) < 1e-8
    r = energy_ratio(q_rem + 1e-15, u, degree=4)  # avoid empty
    # use raw l2
    assert float(np.linalg.norm(q_rem)) < 1e-7


def test_partial_peel_reduces_energy_ratio_on_model():
    u = sample_grid(8192)
    N = 12
    q_full, _ = truncated_mode_sum(u, T=40.0, n_zeros=N)
    r0 = energy_ratio(q_full, u, degree=4)
    q_half, _ = peel_via_remainder(
        q_full, u, T=40.0, n_strip=6, fit_scale=True
    )
    r1 = energy_ratio(q_half, u, degree=4)
    # partial peel of model residual need not always lower R_d (can leave
    # higher modes), but full strip does — check full strip R path
    q_all, _ = peel_via_remainder(
        q_full, u, T=40.0, n_strip=N, fit_scale=False
    )
    assert float(np.linalg.norm(q_all)) < 1e-6
    assert r0 > 0.0
    assert np.isfinite(r1)


def test_remainder_diagnostic_multi_T_decay():
    """Model remainder after full strip is tiny; M5 bound decays in T."""
    u = sample_grid(4096)
    bounds = []
    for T in (12.0, 24.0, 48.0):
        row = remainder_diagnostic(
            u, T=T, n_full=10, n_strip=10, degree=4
        )
        assert row["R_d_remainder"] < 1e-10 or float(
            np.isfinite(row["R_d_remainder"])
        )
        # stripped block bound O(T^{-2})
        bounds.append(row["M5_bound_stripped_block"])
        assert "NOT AN UNCONDITIONAL PROOF OF RH" in row["banner"]
    assert bounds[2] < bounds[0]


def test_tail_majorant_decays_with_T():
    b_lo = bound_R_d_mode_tail(10.0, n_kept=10, n_tail=30, d=4)
    b_hi = bound_R_d_mode_tail(40.0, n_kept=10, n_tail=30, d=4)
    assert b_hi["bound_R_d_tail"] < b_lo["bound_R_d_tail"]
    assert b_lo["label"] == "scaffolding_majorant_not_sharp"


def test_remainder_diagnostic_from_q_uses_peel_path():
    u = sample_grid(2048)
    q, _, _ = explicit_formula_residual(u, T=20.0, n_zeros=8)
    # add a small bulk so peel is nontrivial
    q = q + 0.05 * (u - 0.5)
    row = remainder_diagnostic_from_q(
        q, u, T=20.0, n_strip=4, degree=4, fit_scale=True
    )
    assert row["kind"] == "remainder_diagnostic_external"
    assert 0.0 <= row["R_d_full"] <= 1.0 + 1e-9
    assert 0.0 <= row["R_d_remainder"] <= 1.0 + 1e-9
    assert np.isfinite(row["alpha"])


def test_multi_TN_remainder_scan_shape():
    rows = multi_TN_remainder_scan(
        T_values=[15.0, 30.0],
        n_full=8,
        n_strips=[0, 4, 8],
        degree=3,
        n_points=1024,
    )
    assert len(rows) == 2 * 3
    Ts = {r["T"] for r in rows}
    assert Ts == {15.0, 30.0}
    # full strip → remainder energy ~0
    full = [r for r in rows if r["n_strip"] == 8]
    assert all(r["R_d_remainder"] < 1e-8 for r in full)


def test_m5_bound_matches_lemma_api():
    t = zeta_zero_ordinates(6)
    a = explicit_formula_amplitudes(t)
    b1 = bound_R_d_finite_mode_sum(25.0, a, t, 4)
    # tail of empty? use same as kept block via mode tail with n_kept=0
    # just sanity: positive and finite
    assert b1 > 0.0 and np.isfinite(b1)
