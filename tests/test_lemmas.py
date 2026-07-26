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
    continuous_R_d_orthogonal_defect,
    continuous_R_d_pure_mode,
    predicted_R_d_critical_scaling,
    synthetic_orthogonal_defect,
)
from pbss.probes import (
    probe_critical_line_mode,
    probe_low_degree,
    probe_off_critical_mode,
    sample_grid,
)
from pbss.projection import energy_ratio, project


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
