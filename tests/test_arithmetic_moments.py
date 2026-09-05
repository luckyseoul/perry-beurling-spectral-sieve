"""Continuous arithmetic checks with independent adaptive integration."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from pbss.arithmetic_moments import arithmetic_moments
from pbss.probes import primes_upto


def _adaptive_reference(T, degree):
    # Direct integration of the actual step function, split at every prime.
    primes = primes_upto(int(math.ceil(math.exp(T))))
    logs = np.log(primes.astype(float))
    logs = logs[logs < T]
    edges = np.concatenate(([0.0], logs, [T]))
    theta = np.concatenate(([0.0], np.cumsum(logs)))
    coefficients = np.zeros(degree + 1)
    norm = 0.0
    for a, b, value in zip(edges[:-1], edges[1:], theta):
        q = lambda y: value * np.exp(-y / 2) - np.exp(y / 2)
        norm += quad(lambda y: q(y)**2, a, b, epsabs=2e-13, epsrel=2e-13)[0] / T
        for k in range(degree + 1):
            coefficients[k] += quad(
                lambda y: q(y) * np.sqrt(2 * k + 1) * eval_legendre(k, 2 * y / T - 1),
                a, b, epsabs=2e-13, epsrel=2e-13,
            )[0] / T
    return coefficients, norm


@pytest.mark.parametrize("x", [2.0, 10.0, 100.0, 1000.0])
def test_prime_interval_moments_and_norm_match_direct_integration(x):
    T = math.log(x)
    result = arithmetic_moments(T, degree=4)
    coefficients, norm = _adaptive_reference(T, 4)
    assert result.coeffs == pytest.approx(coefficients, abs=3e-12)
    assert result.l2_norm_sq == pytest.approx(norm, rel=3e-12)
    assert result.energy_ratio == pytest.approx(np.dot(coefficients, coefficients) / norm, rel=3e-12)
    expected_affine_norm = norm - np.dot(coefficients[:2], coefficients[:2])
    assert result.affine_detrended_l2_norm_sq == pytest.approx(expected_affine_norm, rel=2e-8, abs=3e-12)


def test_window_below_first_prime_is_not_clamped_to_two():
    T = np.longdouble("0.25")
    result = arithmetic_moments(T, degree=2)
    expected_c0 = -2 * np.expm1(T / 2) / T
    expected_norm = np.expm1(T) / T
    assert result.metadata["n_primes"] == 0
    assert result.coeffs[0] == pytest.approx(expected_c0, rel=1e-15)
    assert result.l2_norm_sq == pytest.approx(expected_norm, rel=1e-15)


def test_degree_zero_matches_global_finite_prime_sum():
    T = math.log(1000)
    primes = primes_upto(1000)
    prime_sum = math.fsum(math.log(int(p)) * (1 / math.sqrt(int(p)) - math.exp(-T / 2)) for p in primes)
    expected = 2 / T * (prime_sum - math.expm1(T / 2))
    result = arithmetic_moments(T, degree=0)
    assert result.coeffs[0] == pytest.approx(expected, abs=3e-13)
    assert result.affine_detrended_energy == 0


def test_continuous_affine_detrending_preserves_higher_absolute_moments():
    result = arithmetic_moments(math.log(1000), degree=6)
    assert np.array_equal(result.affine_detrended_coeffs[2:], result.coeffs[2:])
    assert np.array_equal(result.affine_detrended_coeffs[:2], np.zeros(2))
    assert result.affine_detrended_l2_norm_sq == pytest.approx(
        result.l2_norm_sq - np.dot(result.coeffs[:2], result.coeffs[:2]), rel=1e-14)
    assert result.affine_detrended_energy_ratio == pytest.approx(
        np.dot(result.coeffs[2:], result.coeffs[2:]) / result.affine_detrended_l2_norm_sq, rel=1e-14)


def test_integrals_are_continuous_across_a_prime_jump():
    boundary = np.log(np.longdouble(29))
    below = np.nextafter(boundary, np.longdouble(-np.inf))
    above = np.nextafter(boundary, np.longdouble(np.inf))
    left, right = arithmetic_moments(below), arithmetic_moments(above)
    assert right.metadata["n_primes"] == left.metadata["n_primes"] + 1
    assert left.coeffs == pytest.approx(right.coeffs, abs=2e-16)
    assert left.l2_norm_sq == pytest.approx(right.l2_norm_sq, abs=2e-16)


def test_supplied_complete_larger_prefix_matches_generated_result():
    T = math.log(20)
    generated = arithmetic_moments(T)
    supplied = arithmetic_moments(T, primes=primes_upto(100))
    assert np.array_equal(generated.coeffs, supplied.coeffs)
    assert supplied.metadata["prime_prefix_status"] == "supplied_prefix_verified_against_complete_finite_sieve"


@pytest.mark.parametrize("primes", [
    [2, 3, 5, 7, 13, 17, 19],
    [2, 3, 5, 7, 9, 11, 13, 17, 19],
    [2, 3, 5, 5, 7, 11, 13, 17, 19],
    [3, 2, 5, 7, 11, 13, 17, 19],
    [1, 2, 3, 5, 7, 11, 13, 17, 19],
    [[2, 3, 5]],
    [2.0, 3.0, 5.0],
])
def test_incomplete_or_malformed_prime_prefix_is_rejected(primes):
    with pytest.raises(ValueError, match="prime"):
        arithmetic_moments(math.log(20), primes=primes)


@pytest.mark.parametrize("T", [0, -1, np.inf, np.nan, 1j, np.complex128(1 + 2j),
                              np.array(1 + 0j), [1.0], np.array([1.0, 2.0]), "bad", 1000])
def test_invalid_or_out_of_budget_window_is_rejected(T):
    with pytest.raises(ValueError):
        arithmetic_moments(T)


@pytest.mark.parametrize("degree", [-1, 1.5, True, 10**30])
def test_invalid_or_unsupported_degree_is_rejected_before_sieve(degree):
    with pytest.raises(ValueError, match="degree"):
        arithmetic_moments(1.0, degree=degree)


def test_explicit_prime_prefix_budget_is_enforced():
    with pytest.raises(ValueError, match="max_prime_limit"):
        arithmetic_moments(math.log(1000), max_prime_limit=100)


def test_unresolved_affine_norm_is_exposed_without_losing_absolute_coefficients():
    result = arithmetic_moments(1e-6, degree=2)
    assert np.all(np.isfinite(result.coeffs))
    assert result.l2_norm_sq > 0
    assert result.affine_detrended_l2_norm_sq is None
    assert result.affine_detrended_energy_ratio is None
    assert result.metadata["affine_detrend_norm_status"] == "unresolved_by_working_precision_subtraction"


def test_small_window_keeps_absolute_moment_when_affine_ratio_is_unresolved():
    T = np.longdouble("0.001")
    result = arithmetic_moments(T, degree=2)
    # Rodrigues gives c2=-sqrt(5)*T^2/8 * integral exp(Tu/2)*u^2*(1-u)^2 du.
    scale = np.sqrt(np.longdouble(5)) * T**2 / 240
    assert -scale * np.exp(T / 2) <= result.coeffs[2] <= -scale
    assert result.l2_norm_sq > 0
    assert np.array_equal(result.affine_detrended_coeffs[2:], result.coeffs[2:])
    if result.affine_detrended_l2_norm_sq is None:
        assert result.affine_detrended_energy_ratio is None
        assert result.metadata["affine_detrend_norm_status"] == "unresolved_by_working_precision_subtraction"
    else:
        assert 0 <= result.affine_detrended_energy_ratio <= 1


def test_json_report_preserves_absolute_values_and_honest_precision_status():
    result = arithmetic_moments(math.log(100), degree=4)
    report = json.loads(json.dumps(result.to_dict(), allow_nan=False))
    reconstructed = np.asarray(report["coefficients_decimal"], dtype=np.longdouble)
    assert np.array_equal(reconstructed, result.coeffs)
    assert report["coefficient_signs"] == np.sign(result.coeffs).astype(int).tolist()
    assert np.asarray(report["log_abs_over_T_decimal"], dtype=np.longdouble) == pytest.approx(
        np.log(np.abs(result.coeffs)) / result.T, rel=1e-15)
    assert report["metadata"]["certified_interval_result"] is False
    assert report["metadata"]["B_ABS_growth_bound_proved"] is False
    assert report["metadata"]["RH_proved"] is False
    assert report["metadata"]["working_precision_bits"] == np.finfo(np.longdouble).nmant + 1
    assert np.all(result.kernel_series_tail_estimates < result.coefficient_roundoff_indicators)
