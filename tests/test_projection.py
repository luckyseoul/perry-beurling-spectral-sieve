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


# Regression checks use the discrete weighted inner product itself as oracle.
from pbss.projection import project_coefficients, scaled_projection_strength
from pbss.projection_backend import energy_ratio_auto, project_coefficients_auto


@pytest.mark.parametrize("u", [np.linspace(0.0, 1.0, 5), np.array([0.0, 0.01, 0.2, 0.4, 1.0])])
def test_coarse_polynomial_projection_is_exact_and_bounded(u):
    coefficients = np.array([0.3, -0.2, 0.7, 0.1, -0.05])
    design = orthonormal_legendre_design(4, u)
    q = design @ coefficients
    result = project(q, u, degree=4)
    assert result.energy_ratio == pytest.approx(1.0, abs=2e-14)
    assert result.energy_ratio <= 1.0
    assert result.coeffs == pytest.approx(coefficients, abs=2e-13)
    assert design @ result.coeffs == pytest.approx(q, abs=2e-13)


def test_weighted_fit_has_orthogonal_residual_and_pythagorean_energy():
    u = np.array([0.0, 0.03, 0.2, 0.45, 0.8, 1.0])
    w = np.array([0.2, 0.3, 1.7, 0.4, 0.8, 0.6])
    q = np.array([1.0, -0.3, 0.9, -1.4, 0.2, 1.8])
    design = orthonormal_legendre_design(2, u)
    expected, *_ = np.linalg.lstsq(np.sqrt(w)[:, None] * design, np.sqrt(w) * q, rcond=None)
    result = project(q, u, degree=2, weights=w)
    fitted = design @ result.coeffs
    residual = q - fitted
    assert result.coeffs == pytest.approx(expected, abs=2e-14)
    assert design.T @ (w * residual) == pytest.approx(np.zeros(3), abs=3e-14)
    assert result.energy == pytest.approx(np.dot(w, fitted**2), rel=2e-14)
    assert result.energy + np.dot(w, residual**2) == pytest.approx(result.l2_norm_sq, rel=2e-14)
    assert not np.isclose(result.energy, np.dot(result.coeffs, result.coeffs))
    assert result.energy_ratio == pytest.approx(result.energy / result.l2_norm_sq)
    assert energy_ratio(residual, u, degree=2, weights=w) < 1e-27


def test_gauss_weights_recover_continuous_legendre_coefficients():
    nodes, w = np.polynomial.legendre.leggauss(6)
    u, w = (nodes + 1.0) / 2.0, w / 2.0
    coefficients = np.array([0.5, -0.2, 0.3])
    q = orthonormal_legendre_design(2, u) @ coefficients
    result = project(q, u, degree=2, weights=w)
    assert result.coeffs == pytest.approx(coefficients, abs=2e-14)
    assert result.energy == pytest.approx(np.dot(coefficients, coefficients), rel=2e-14)


def test_zero_weights_ignore_samples_without_changing_fitted_polynomial():
    u = np.array([0.0, 0.2, 0.7, 1.0])
    q = np.array([1e300, 2.0, 3.0, 1e300])
    w = np.array([0.0, 2.0, 3.0, 0.0])
    result = project(q, u, degree=1, weights=w)
    fit = orthonormal_legendre_design(1, u) @ result.coeffs
    assert fit[1:3] == pytest.approx(q[1:3], abs=1e-14)
    assert result.l2_norm_sq == pytest.approx(35.0)
    assert result.energy_ratio == pytest.approx(1.0)


@pytest.mark.parametrize("q_scale", [1e-250, -1e-250, 1.0, 1e250, -1e250])
@pytest.mark.parametrize("w_scale", [1e-200, 1.0, 1e200])
def test_ratio_is_invariant_under_extreme_amplitude_and_weight_scaling(q_scale, w_scale):
    u = np.array([0.0, 0.02, 0.3, 0.8, 1.0])
    q = np.array([0.4, -0.5, 0.9, 0.1, -0.6])
    w = np.array([0.3, 0.7, 0.4, 0.9, 0.1])
    expected = energy_ratio(q, u, degree=2, weights=w)
    actual = energy_ratio(q * q_scale, u, degree=2, weights=w * w_scale)
    assert actual == pytest.approx(expected, rel=3e-14)


def test_nested_projection_energy_is_monotone():
    u = np.array([0.0, 0.02, 0.1, 0.3, 0.5, 0.8, 1.0])
    q = np.cos(13.0 * u) + u**2
    ratios = [energy_ratio(q, u, degree=d) for d in range(5)]
    assert np.all(np.diff(ratios) >= -2e-14)
    assert 0.0 <= ratios[0] <= ratios[-1] <= 1.0


@pytest.mark.parametrize("degree", [-1, 1.2, 2.0, True, np.nan, 10**30])
def test_invalid_projection_degree_is_rejected_before_design_allocation(degree):
    with pytest.raises(ValueError, match="degree"):
        project(np.ones(3), np.linspace(0.0, 1.0, 3), degree=degree)


@pytest.mark.parametrize("u", [
    [], [0.0], [0.0, 0.0, 1.0], [1.0, 0.5, 0.0], [-0.1, 0.5, 1.0],
    [0.0, 0.5, 1.1], [0.0, np.nan, 1.0], [0.0, np.inf, 1.0], [[0.0, 0.5, 1.0]],
])
def test_invalid_projection_grid_is_rejected(u):
    with pytest.raises(ValueError):
        energy_ratio(np.ones(np.asarray(u).shape), u, degree=0)


@pytest.mark.parametrize("q", [[1.0, np.nan, 2.0], [1.0, np.inf, 2.0], [[1.0, 2.0, 3.0]], [1.0, 2.0], [1.0, 2.0 + 1j, 3.0]])
def test_invalid_signal_is_rejected(q):
    with pytest.raises(ValueError, match="q"):
        energy_ratio(q, [0.0, 0.5, 1.0], degree=0)


@pytest.mark.parametrize("weights", [[1.0, -1.0, 1.0], [0.0, 0.0, 0.0], [1.0, np.nan, 1.0], [1.0, np.inf, 1.0], [1.0, 1.0], [[1.0, 1.0, 1.0]], [1.0, 1j, 1.0]])
def test_invalid_weights_are_rejected(weights):
    with pytest.raises(ValueError):
        project([1.0, 2.0, 3.0], [0.0, 0.5, 1.0], degree=0, weights=weights)


def test_insufficient_positive_weight_support_is_rejected():
    with pytest.raises(ValueError, match="positive-weight"):
        project([1.0, 2.0, 3.0], [0.0, 0.5, 1.0], degree=1, weights=[0.0, 1.0, 0.0])


def test_numerically_rank_deficient_grid_is_rejected():
    with pytest.raises(ValueError, match="rank deficient"):
        project(np.ones(8), np.linspace(0.0, 1e-14, 8), degree=3)


@pytest.mark.parametrize("T", [0.0, -1.0, np.nan, np.inf, -np.inf])
def test_invalid_T_is_rejected_consistently(T):
    q, u = np.ones(3), np.linspace(0.0, 1.0, 3)
    with pytest.raises(ValueError, match="T"):
        project(q, u, degree=1, T=T)
    with pytest.raises(ValueError, match="T"):
        scaled_projection_strength(q, u, degree=1, T=T)


def test_zero_signal_has_zero_coefficients_and_energy_but_no_ratio():
    q, u = np.zeros(4), np.linspace(0.0, 1.0, 4)
    assert project_coefficients(q, u, degree=2) == pytest.approx(np.zeros(3))
    assert projection_energy(q, u, degree=2) == 0.0
    with pytest.raises(ValueError, match="zero"):
        energy_ratio(q, u, degree=2)
    with pytest.raises(ValueError, match="zero"):
        project(q, u, degree=2)


def test_scaled_strength_avoids_intermediate_power_overflow():
    q = np.array([1.0, -1.0 + 1e-8, 1.0])
    u = np.array([0.0, 0.5, 1.0])
    ratio = energy_ratio(q, u, degree=0)
    strength = scaled_projection_strength(q, u, degree=0, T=1e160)
    expected = float(np.longdouble(1e160)**2 * ratio)
    assert np.isfinite(strength)
    assert strength == pytest.approx(expected, rel=3e-14)


@pytest.mark.parametrize("scale", [1e-250, 1e250])
def test_unrepresentable_absolute_energy_is_reported_without_breaking_ratio(scale):
    q, u = np.full(3, scale), np.linspace(0.0, 1.0, 3)
    assert energy_ratio(q, u, degree=0) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="float64"):
        project(q, u, degree=0)


@pytest.mark.parametrize("T", [1e-250, 1e250])
def test_unrepresentable_scaled_strength_is_reported(T):
    with pytest.raises(ValueError, match="scaled strength"):
        project(np.ones(3), np.linspace(0.0, 1.0, 3), degree=0, T=T)


def test_numpy_auto_backend_matches_custom_weighted_projection():
    u = np.array([0.0, 0.02, 0.3, 0.8, 1.0])
    q = np.sin(10.0 * u) + 0.2
    w = np.array([0.1, 0.2, 0.7, 0.3, 0.6])
    expected = project(q, u, degree=2, weights=w)
    coefficients, backend = project_coefficients_auto(q, u, 2, weights=w, prefer_gpu=False)
    ratio, ratio_backend = energy_ratio_auto(q, u, 2, weights=w, prefer_gpu=False)
    assert backend == ratio_backend == "numpy"
    assert coefficients == pytest.approx(expected.coeffs, abs=2e-14)
    assert ratio == pytest.approx(expected.energy_ratio, rel=2e-14)


def test_auto_backend_validates_before_attempting_gpu(monkeypatch):
    import pbss.projection_backend as backend

    def unexpected_gpu_probe():
        pytest.fail("invalid inputs must be rejected before touching GPU")

    monkeypatch.setattr(backend, "cupy_available", unexpected_gpu_probe)
    with pytest.raises(ValueError, match="q"):
        backend.energy_ratio_auto(np.full(8192, np.nan), np.linspace(0.0, 1.0, 8192), 2)


@pytest.mark.parametrize("k", [-1, 1.2, True])
def test_basis_rejects_invalid_degree(k):
    with pytest.raises(ValueError, match="degree"):
        shifted_legendre_values(k, np.array([0.0, 1.0]))


@pytest.mark.parametrize("u", [[0.0, np.nan], [0.0, np.inf], [-0.1, 1.0], [0.0, 1.1], [0.0, 1j]])
def test_basis_rejects_invalid_points(u):
    with pytest.raises(ValueError, match="u"):
        shifted_legendre_values(2, u)


def test_nonfinite_factorization_output_is_rejected(monkeypatch):
    original_qr = np.linalg.qr

    def broken_qr(*args, **kwargs):
        Q, R = original_qr(*args, **kwargs)
        Q[:] = np.nan
        return Q, R

    monkeypatch.setattr(np.linalg, "qr", broken_qr)
    with pytest.raises(ValueError, match="non-finite ratio"):
        energy_ratio(np.ones(3), np.linspace(0.0, 1.0, 3), degree=1)


def test_nonfinite_solver_output_is_rejected(monkeypatch):
    monkeypatch.setattr(np.linalg, "solve", lambda a, b: np.full_like(b, np.nan))
    with pytest.raises(ValueError, match="non-finite coefficients"):
        project(np.ones(3), np.linspace(0.0, 1.0, 3), degree=1)
