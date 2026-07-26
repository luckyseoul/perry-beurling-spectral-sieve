"""
Unit tests for the shipped PBSS projection path.

All assertions drive pbss.projection / pbss.basis / pbss.probes — no
re-implementation of the diagnostic inside the tests.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.basis import orthonormal_legendre_design, shifted_legendre_values
from pbss.probes import (
    normalize_l2,
    probe_critical_line_mode,
    probe_defective,
    probe_high_frequency,
    probe_low_degree,
    probe_prime_residual,
    sample_grid,
)
from pbss.projection import energy_ratio, project, projection_energy


def test_basis_orthonormal_quadrature():
    """Discrete trapezoid Gram matrix of φ_0..φ_4 ≈ I on a fine grid."""
    u = sample_grid(4000)
    Phi = orthonormal_legendre_design(4, u)
    du = u[1] - u[0]
    # trapezoid ≈ du * (½ ends + interior); use uniform weight du for dense grid
    G = (Phi.T * du) @ Phi
    # endpoints half-weight correction is O(1/n); check near-identity
    assert np.allclose(G, np.eye(5), atol=2e-2), G


def test_low_degree_probe_has_energy_ratio_near_one():
    """φ_2 lives in degree-≤4 space → R_4 ≈ 1."""
    u = sample_grid(2000)
    q = probe_low_degree(u, k=2)
    r = energy_ratio(q, u, degree=4)
    assert r == pytest.approx(1.0, abs=2e-2), r


def test_high_frequency_has_small_energy_ratio():
    """Rapid sinusoid is mostly orthogonal to degree ≤ 4."""
    u = sample_grid(4000)
    q = probe_high_frequency(u, waves=50)
    r = energy_ratio(q, u, degree=4)
    assert 0.0 <= r < 0.05, r


def test_defective_exceeds_high_frequency_energy_ratio():
    """Low-degree contamination must inflate R_d vs pure high-frequency."""
    u = sample_grid(4000)
    q_hf = probe_high_frequency(u, waves=50)
    q_def = probe_defective(u, waves=50, defect_degree=1, defect_weight=2.0)
    r_hf = energy_ratio(q_hf, u, degree=4)
    r_def = energy_ratio(q_def, u, degree=4)
    assert r_def > r_hf + 0.1, (r_hf, r_def)


def test_project_returns_finite_non_nan_on_tiny_synthetic():
    """Gating AC: real entry point on tiny q → finite non-NaN P."""
    u = sample_grid(64)
    q = probe_high_frequency(u, waves=8)
    res = project(q, u, degree=3, T=10.0)
    assert np.isfinite(res.energy)
    assert np.isfinite(res.energy_ratio)
    assert np.isfinite(res.scaled_strength)
    assert np.isfinite(res.P)
    assert res.P == res.scaled_strength
    assert res.coeffs.shape == (4,)
    assert res.n_points == 64


def test_unit_normalized_energy_equals_ratio():
    u = sample_grid(1000)
    q = normalize_l2(probe_defective(u), u)
    e = projection_energy(q, u, degree=3)
    r = energy_ratio(q, u, degree=3)
    assert e == pytest.approx(r, rel=1e-6)


def test_prime_residual_probe_runs_and_projects():
    """Prime-based RH-consistent probe: finite projection output."""
    u = sample_grid(1024)
    q, T = probe_prime_residual(u, x_max=1e4)
    assert T == pytest.approx(np.log(1e4))
    res = project(q, u, degree=4, T=T)
    assert np.isfinite(res.energy_ratio)
    assert 0.0 <= res.energy_ratio <= 1.0 + 1e-9
    assert np.isfinite(res.P)


def test_shifted_legendre_degree_zero_is_constant_one():
    u = sample_grid(50)
    phi0 = shifted_legendre_values(0, u)
    assert np.allclose(phi0, 1.0)


def test_critical_line_mode_below_defective():
    """RH-form oscillation has less low-degree energy than a defective probe."""
    u = sample_grid(4000)
    r_cl = energy_ratio(probe_critical_line_mode(u, T=20.0), u, degree=4)
    r_def = energy_ratio(probe_defective(u, waves=48, defect_weight=2.5), u, degree=4)
    assert r_cl < r_def
    assert r_cl < 0.15
