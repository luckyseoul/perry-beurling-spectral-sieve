"""
Tests that drive shipped projection against proved lemma predictions.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.lemmas import (
    bound_R_d_finite_mode_sum,
    continuous_R_d_orthogonal_defect,
    continuous_R_d_pure_mode,
    finite_mode_R_d_order_T,
    predicted_R_d_critical_scaling,
    synthetic_orthogonal_defect,
)
from pbss.probes import (
    finite_cl_superposition,
    probe_critical_line_mode,
    probe_low_degree,
    probe_off_critical_mode,
    sample_grid,
)
from pbss.projection import energy_ratio, project
from pbss.zeros import explicit_formula_amplitudes, zeta_zero_ordinates


def test_M1_pure_mode_in_subspace():
    u = sample_grid(3000)
    q = probe_low_degree(u, k=2)
    r = energy_ratio(q, u, degree=4)
    assert r == pytest.approx(continuous_R_d_pure_mode(2, 4), abs=0.03)


def test_M1_pure_mode_outside_subspace():
    u = sample_grid(4000)
    q = probe_low_degree(u, k=8)
    r = energy_ratio(q, u, degree=3)
    assert r == pytest.approx(continuous_R_d_pure_mode(8, 3), abs=0.05)


def test_M2_orthogonal_defect_formula():
    u = sample_grid(5000)
    eps = 0.4
    q = synthetic_orthogonal_defect(u, eps=eps, j=1, waves=100)
    r = energy_ratio(q, u, degree=4)
    assert r == pytest.approx(continuous_R_d_orthogonal_defect(eps), abs=0.05)


def test_M3_critical_line_R_d_decays_with_T():
    u = sample_grid(8000)
    r_small = energy_ratio(probe_critical_line_mode(u, T=5.0), u, degree=4)
    r_large = energy_ratio(probe_critical_line_mode(u, T=40.0), u, degree=4)
    assert r_large < r_small
    assert r_large < 0.05


def test_M3_matches_leading_R0_asymptotics():
    """Shipped R_0 for sin(t T u) ≈ closed-form predicted_R_d_critical_scaling."""
    u = sample_grid(10000)
    T = 25.0
    q = probe_critical_line_mode(u, T=T)
    r0 = energy_ratio(q, u, degree=0)
    pred = predicted_R_d_critical_scaling(T)
    assert r0 == pytest.approx(pred, rel=0.08, abs=0.01)


def test_M4_persistent_defect_independent_of_waves():
    u = sample_grid(4000)
    eps = 0.55
    r1 = energy_ratio(synthetic_orthogonal_defect(u, eps, waves=40), u, degree=3)
    r2 = energy_ratio(synthetic_orthogonal_defect(u, eps, waves=120), u, degree=3)
    assert r1 == pytest.approx(eps**2, abs=0.06)
    assert r2 == pytest.approx(eps**2, abs=0.06)


def test_off_critical_probe_finite():
    u = sample_grid(1024)
    q = probe_off_critical_mode(u, T=15.0, sigma=0.8)
    res = project(q, u, degree=4, T=15.0)
    assert np.isfinite(res.energy_ratio)
    assert 0.0 <= res.energy_ratio <= 1.0 + 1e-9


def test_M5_finite_mode_bound_is_O_T_minus_2():
    """Majorant scales as T^{-2}; T²·bound is T-independent (Lemma M5)."""
    t = zeta_zero_ordinates(5)
    a = explicit_formula_amplitudes(t)
    d = 4
    b10 = bound_R_d_finite_mode_sum(10.0, a, t, d)
    b40 = bound_R_d_finite_mode_sum(40.0, a, t, d)
    assert b40 < b10
    # (T_large/T_small)^2 · R_bound(T_large) ≈ R_bound(T_small)
    ratio = b10 / b40
    assert ratio == pytest.approx((40.0 / 10.0) ** 2, rel=1e-9)
    o10 = finite_mode_R_d_order_T(10.0, a, t, d)
    o40 = finite_mode_R_d_order_T(40.0, a, t, d)
    assert o10 == pytest.approx(o40, rel=1e-12)


def test_M5_finite_superposition_R_d_decays_with_T():
    """Shipped finite CL sum (N>1) has R_d decreasing in T — finite-mode A₀."""
    u = sample_grid(10000)
    t = zeta_zero_ordinates(4)
    a = explicit_formula_amplitudes(t)
    d = 4
    r_small = energy_ratio(
        finite_cl_superposition(u, T=8.0, amplitudes=a, ordinates=t, form="sin"),
        u,
        degree=d,
    )
    r_large = energy_ratio(
        finite_cl_superposition(u, T=48.0, amplitudes=a, ordinates=t, form="sin"),
        u,
        degree=d,
    )
    assert r_large < r_small
    assert r_large < 0.08
    # measured R stays under the crude M5 majorant at large T
    assert r_large <= bound_R_d_finite_mode_sum(48.0, a, t, d) + 1e-6


def test_M5_single_mode_recovers_M3_scale():
    """N=1 finite sum reduces to pure CL mode; R_d decays like M3."""
    u = sample_grid(10000)
    t = np.array([14.134725141734693])
    a = np.array([1.0])
    r5 = energy_ratio(
        finite_cl_superposition(u, T=5.0, amplitudes=a, ordinates=t), u, degree=4
    )
    r40 = energy_ratio(
        finite_cl_superposition(u, T=40.0, amplitudes=a, ordinates=t), u, degree=4
    )
    r_m3 = energy_ratio(probe_critical_line_mode(u, T=40.0), u, degree=4)
    assert r40 < r5
    assert r40 == pytest.approx(r_m3, rel=0.05, abs=1e-4)
