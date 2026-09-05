"""Independent checks of the continuous real-sine obstruction identities."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from pbss.off_critical_exact import (
    asymptotic_constants,
    detrended_ratio,
    off_critical_projection,
    partial_cancellation_ratio,
    phase_subsequence_T,
)


def _integrate(function):
    value, error = quad(function, 0.0, 1.0, epsabs=2e-13, epsrel=2e-12, limit=300)
    return value


@pytest.mark.parametrize("T", [0.001, 2.0, 20.0, 200.0])
def test_exact_identity_matches_independent_adaptive_integration(T):
    A, t, degree = 0.4, 2.0, 4
    result = off_critical_projection(T, A=A, t=t, degree=degree)

    def q(u):
        return np.exp(A * T * (u - 1.0)) * np.sin(t * T * u)

    coefficients = np.array([
        _integrate(lambda u: q(u) * np.sqrt(2 * k + 1) * eval_legendre(k, 2 * u - 1))
        for k in range(degree + 1)
    ])
    norm = _integrate(lambda u: q(u) ** 2)
    assert result.scaled_coeffs == pytest.approx(coefficients, abs=3e-13)
    assert result.scaled_norm_sq == pytest.approx(norm, rel=2e-11)
    assert result.energy_ratio == pytest.approx(np.dot(coefficients, coefficients) / norm, rel=2e-10)
    assert 0.0 <= result.energy_ratio <= 1.0


def test_large_exponential_envelope_needs_no_overflowing_samples():
    result = off_critical_projection(4000.0, A=0.4, t=2.0, degree=4)
    assert np.all(np.isfinite(result.scaled_coeffs))
    assert np.isfinite(result.scaled_norm_sq)
    assert 0.0 < result.energy_ratio < 0.01


@pytest.mark.parametrize("degree", [1, 4])
@pytest.mark.parametrize("kind,key,power", [
    ("maximal", "limsup_T_R", 1),
    ("leading_cancel", "leading_cancel_T3_R", 3),
    ("minimum_T3", "liminf_T3_R", 3),
])
def test_distinct_phase_subsequences_converge_to_proved_constants(degree, kind, key, power):
    A, t = 0.4, 2.0
    expected = asymptotic_constants(A=A, t=t, degree=degree)[key]
    errors = []
    for n in (200, 800, 3200):
        T = phase_subsequence_T(n, A=A, t=t, degree=degree, kind=kind)
        ratio = off_critical_projection(T, A=A, t=t, degree=degree).energy_ratio
        errors.append(abs(T**power * ratio / expected - 1.0))
    assert errors[-1] < 2e-3
    # The d=1 leading-cancel formula is already exact apart from exponentially
    # small endpoint terms; after those vanish, phase rounding need not improve.
    assert errors[-1] < max(errors[0], 1e-7)


def test_old_envelope_fails_as_a_lower_bound_even_without_cancellation():
    T = phase_subsequence_T(800, A=0.4, t=2.0, degree=4)
    result = off_critical_projection(T, A=0.4, t=2.0, degree=4)
    assert result.energy_ratio < 1e-4 * result.phase_envelope_ratio


def test_degree_zero_phase_cancellation_has_exact_exponential_formula():
    A, t = 0.4, 2.0
    T = phase_subsequence_T(5, A=A, t=t)
    expected = 4 * A / ((A * A + t * t) * T * np.expm1(2 * A * T))
    actual = off_critical_projection(T, A=A, t=t, degree=0).energy_ratio
    assert actual == pytest.approx(expected, rel=2e-10)


@pytest.mark.parametrize("alpha", [0.0, 0.3, 0.8, 1.0])
def test_partial_cancellation_accounts_for_changed_total_norm(alpha):
    T, A, t, degree = 10.0, 0.4, 2.0, 3
    result = off_critical_projection(T, A=A, t=t, degree=degree)

    def remainder(u):
        low = sum(c * np.sqrt(2 * k + 1) * eval_legendre(k, 2 * u - 1)
                  for k, c in enumerate(result.scaled_coeffs))
        return np.exp(A * T * (u - 1)) * np.sin(t * T * u) - alpha * low

    norm_after = _integrate(lambda u: remainder(u)**2)
    direct = (1 - alpha)**2 * result.scaled_projection_energy / norm_after
    assert partial_cancellation_ratio(result.energy_ratio, alpha) == pytest.approx(direct, rel=2e-11, abs=1e-14)


def test_partial_cancellation_zero_signal_has_no_ratio():
    with pytest.raises(ValueError, match="zero norm"):
        partial_cancellation_ratio(1.0, 1.0)


def test_detrending_uses_nested_orthogonal_projection_identity():
    args = dict(T=20.0, A=0.4, t=2.0)
    full = off_critical_projection(**args, degree=4).energy_ratio
    removed = off_critical_projection(**args, degree=1).energy_ratio
    expected = (full - removed) / (1 - removed)
    assert detrended_ratio(**args, degree=4, remove_degree=1) == pytest.approx(expected, rel=1e-12)
    assert detrended_ratio(**args, degree=2, remove_degree=4) == 0.0


@pytest.mark.parametrize("parameter,value", [("T", 0), ("T", np.inf), ("A", -0.1), ("A", np.nan), ("t", 0), ("degree", -1), ("degree", 1.2)])
def test_invalid_model_parameters_are_rejected(parameter, value):
    args = dict(T=10.0, A=0.4, t=2.0, degree=4)
    args[parameter] = value
    with pytest.raises(ValueError):
        off_critical_projection(**args)
