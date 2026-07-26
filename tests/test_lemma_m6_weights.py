"""Lemma M6: admissible weights preserve model-mode O(T^{-2}) decay (shipped APIs)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.lemmas import (  # noqa: E402
    bound_R_d_weighted_finite_mode_sum,
    bound_R_d_weighted_sine_order,
    weighted_cl_R_d_order_T,
)
from pbss.probes import (  # noqa: E402
    explicit_formula_residual,
    probe_critical_line_mode,
    sample_grid,
)
from pbss.projection import energy_ratio  # noqa: E402
from pbss.weights import admissible_weight, apply_weight  # noqa: E402
from pbss.zeros import explicit_formula_amplitudes, zeta_zero_ordinates  # noqa: E402


def test_weighted_sine_majorant_decays_in_omega():
    b_lo = bound_R_d_weighted_sine_order(20.0, 4, w_linf=1.0)
    b_hi = bound_R_d_weighted_sine_order(80.0, 4, w_linf=1.0)
    assert b_hi < b_lo
    assert b_hi == pytest.approx(b_lo / 16.0, rel=1e-9)


def test_weighted_cl_empirical_below_majorant_and_decays():
    u = sample_grid(10000)
    t0 = 14.134725
    rows = []
    for T in (15.0, 30.0, 60.0):
        q = probe_critical_line_mode(u, T=T, t=t0)
        for name in ("hanning", "tukey"):
            w = admissible_weight(u, name=name, alpha=0.1)
            qw = apply_weight(q, w)
            r = energy_ratio(qw, u, degree=4)
            maj = bound_R_d_weighted_sine_order(
                t0 * T, 4, w_linf=1.0, wq_l2_floor=0.02
            )
            rows.append((T, name, r, maj))
            assert 0.0 <= r <= 1.0 + 1e-9
            assert r <= maj + 1e-6  # majorant above empirical
    # decay in T for hanning
    r15 = [r for T, n, r, m in rows if T == 15.0 and n == "hanning"][0]
    r60 = [r for T, n, r, m in rows if T == 60.0 and n == "hanning"][0]
    assert r60 < r15


def test_weighted_finite_ef_decays_and_majorant_scales():
    u = sample_grid(8192)
    t = zeta_zero_ordinates(8)
    a = explicit_formula_amplitudes(t)
    r_list = []
    for T in (12.0, 24.0, 48.0):
        q, _, _ = explicit_formula_residual(u, T=T, n_zeros=8)
        w = admissible_weight(u, name="tukey", alpha=0.12)
        r = energy_ratio(apply_weight(q, w), u, degree=4)
        maj = bound_R_d_weighted_finite_mode_sum(T, a, t, 4, w_linf=1.0)
        r_list.append(r)
        assert r <= maj + 0.05  # loose majorant room
    assert r_list[-1] < r_list[0]


def test_weighted_cl_order_T_stable_scaling():
    o1 = weighted_cl_R_d_order_T(20.0, d=4)
    o2 = weighted_cl_R_d_order_T(40.0, d=4)
    assert o1 == pytest.approx(o2, rel=1e-12)
