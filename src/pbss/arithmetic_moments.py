"""Absolute continuous theta moments from complete finite prime intervals.

This evaluates analytic integral identities, not a sampled or interpolated
residual. All absolute coefficients and norms are retained in longdouble.
Floating-point values, series-tail estimates and roundoff indicators are not
interval certificates and cannot prove an all-window RH-equivalent bound.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import factorial
from numbers import Integral
from typing import Any

import numpy as np

from .basis import _validate_degree
from .probes import primes_upto

_LD = np.longdouble
_SERIES_ORDER = 40


def _decimal(value):
    if value is None or not np.isfinite(value):
        return None
    return np.format_float_scientific(_LD(value), unique=True, trim="-")


@dataclass(frozen=True)
class ArithmeticMoments:
    T: np.longdouble
    degree: int
    coeffs: np.ndarray
    affine_detrended_coeffs: np.ndarray
    l2_norm_sq: np.longdouble
    affine_detrended_l2_norm_sq: np.longdouble | None
    energy: np.longdouble
    affine_detrended_energy: np.longdouble
    energy_ratio: float
    affine_detrended_energy_ratio: float | None
    coefficient_roundoff_indicators: np.ndarray
    coefficient_cancellation_indicators: np.ndarray
    kernel_series_tail_estimates: np.ndarray
    norm_roundoff_indicator: np.longdouble
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Strict-JSON-compatible report retaining absolute values as decimals."""
        logs = [np.log(abs(c)) if c != 0 else None for c in self.coeffs]
        return {
            "T_decimal": _decimal(self.T),
            "degree": self.degree,
            "coefficients_decimal": [_decimal(c) for c in self.coeffs],
            "coefficient_signs": [int(np.sign(c)) for c in self.coeffs],
            "log_abs_coefficients_decimal": [_decimal(v) for v in logs],
            "log_abs_over_T_decimal": [_decimal(v / self.T) if v is not None else None for v in logs],
            "affine_detrended_coefficients_decimal": [_decimal(c) for c in self.affine_detrended_coeffs],
            "l2_norm_sq_decimal": _decimal(self.l2_norm_sq),
            "affine_detrended_l2_norm_sq_decimal": _decimal(self.affine_detrended_l2_norm_sq),
            "projection_energy_decimal": _decimal(self.energy),
            "affine_detrended_projection_energy_decimal": _decimal(self.affine_detrended_energy),
            "energy_ratio": self.energy_ratio,
            "affine_detrended_energy_ratio": self.affine_detrended_energy_ratio,
            "coefficient_roundoff_indicators_decimal": [_decimal(v) for v in self.coefficient_roundoff_indicators],
            "coefficient_cancellation_indicators_decimal": [_decimal(v) for v in self.coefficient_cancellation_indicators],
            "kernel_series_tail_estimates_decimal": [_decimal(v) for v in self.kernel_series_tail_estimates],
            "norm_roundoff_indicator_decimal": _decimal(self.norm_roundoff_indicator),
            "metadata": self.metadata,
        }


def _verified_prefix(T, primes, max_prime_limit):
    if isinstance(max_prime_limit, (bool, np.bool_)) or not isinstance(max_prime_limit, Integral) or max_prime_limit < 2:
        raise ValueError("max_prime_limit must be an integer at least 2")
    max_prime_limit = int(max_prime_limit)
    try:
        if np.iscomplexobj(T) or np.asarray(T).ndim != 0:
            raise ValueError("T must be a finite positive real scalar")
        T = _LD(T)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("T must be finite and positive") from exc
    if not np.isfinite(T) or T <= 0:
        raise ValueError("T must be finite and positive")
    # A tiny tolerance allows a binary64 log(integer cap) as the window input.
    if T > np.log(_LD(max_prime_limit)) + _LD("1e-12"):
        raise ValueError("exp(T) exceeds max_prime_limit; increase the explicit finite-prefix limit")
    endpoint = np.exp(T)
    limit = min(max_prime_limit, int(np.ceil(endpoint)))
    expected = primes_upto(limit)
    expected = expected[np.log(expected.astype(_LD)) <= T]
    if primes is None:
        prefix = expected
        provenance = "generated_by_complete_finite_sieve"
    else:
        supplied = np.asarray(primes)
        if supplied.ndim != 1 or (supplied.size and supplied.dtype.kind not in "iu"):
            raise ValueError("primes must be a one-dimensional integer array")
        if np.any(supplied < 2) or np.any(supplied[1:] <= supplied[:-1]):
            raise ValueError("primes must be strictly increasing integers at least 2")
        prefix = supplied[supplied <= limit].astype(np.int64)
        prefix = prefix[np.log(prefix.astype(_LD)) <= T]
        if not np.array_equal(prefix, expected):
            raise ValueError("supplied primes are not the complete prime prefix for this window")
        provenance = "supplied_prefix_verified_against_complete_finite_sieve"
    logs = np.log(prefix.astype(_LD))
    use = logs <= T
    return T, prefix[use], logs[use], limit, endpoint, provenance


def _exponential_kernels(h: np.ndarray, degree: int):
    """Integrals on local v in [0,1] of v^j exp(-hv/2) and 2v^j sinh(hv/2).

    The series are exact identities; truncation after order N has the explicit
    tail estimate documented in EXACT_ARITHMETIC_MOMENTS.md. Complete ordinary
    prime intervals have h<=log(2), making N=40 ample in working precision.
    """
    j = np.arange(degree + 1, dtype=_LD)[None, :]
    f = np.broadcast_to(1 / (j + 1), (h.size, degree + 1)).copy()
    g = np.zeros_like(f)
    power = np.ones(h.size, dtype=_LD)
    x = h / 2
    for n in range(1, _SERIES_ORDER + 1):
        power *= x / n
        term = power[:, None] / (j + n + 1)
        f += (-1 if n % 2 else 1) * term
        if n % 2:
            g += 2 * term
    tail = (np.exp(x) * x ** (_SERIES_ORDER + 1)
            / _LD(factorial(_SERIES_ORDER + 1)))[:, None] / (j + _SERIES_ORDER + 2)
    return f, g, tail


def _positive_norm_intervals(x: np.ndarray, delta: np.ndarray, h: np.ndarray):
    """Exact interval norm identity evaluated without subtracting its large terms."""
    a = -np.expm1(-h)
    b = np.zeros_like(h)
    c = np.zeros_like(h)
    term = np.ones_like(h)
    for n in range(1, _SERIES_ORDER + 1):
        term *= h / n
        if n >= 2:
            b += (-1 if n % 2 else 1) * term
        if n >= 3 and n % 2:
            c += 2 * term
    variance = c - b * b / a
    if np.any(variance < 0):
        raise ValueError("interval norm is numerically unresolved")
    return x * (a * (delta - b / a) ** 2 + variance)


def arithmetic_moments(T: float, degree: int = 4, *, primes=None,
                       max_prime_limit: int = 1_000_000) -> ArithmeticMoments:
    """Integrate c_0,...,c_d and the true squared norm on exact prime intervals.

    Q(y)=exp(-y/2)*(theta(exp(y))-exp(y)). No clamping below x=2,
    smoothing, sampled detrending, or quadrature-grid interpolation is used.
    Supplied integer prime data are checked against a complete finite sieve.
    ``max_prime_limit`` makes the intended finite-prefix scope explicit.

    Continuous affine detrending removes coefficients 0 and 1 exactly; higher
    moments are unchanged. Its norm is the orthogonal Pythagorean difference.
    If that difference is below the working roundoff indicator, its norm and
    ratio are returned as None rather than as a spurious nonnegative number.
    """
    degree = _validate_degree(degree)
    # This evaluator is for a fixed low degree. Prevent accidental allocations
    # from an unbounded CLI/library integer and ill-conditioned high-degree jets.
    if degree > 32:
        raise ValueError("degree exceeds the supported low-degree range 0..32")
    T, prefix, logs, limit, endpoint, provenance = _verified_prefix(T, primes, max_prime_limit)
    internal_degree = max(degree, 1)
    left = np.concatenate((np.zeros(1, dtype=_LD), logs))
    right = np.concatenate((logs, np.array([T], dtype=_LD)))
    h = right - left
    x = np.concatenate((np.ones(1, dtype=_LD), prefix.astype(_LD)))
    theta = np.concatenate((np.zeros(1, dtype=_LD), np.cumsum(logs, dtype=_LD)))
    nonempty = h > 0
    left, h, x, theta = left[nonempty], h[nonempty], x[nonempty], theta[nonempty]
    if np.any(h > np.log(_LD(2)) + 32 * np.finfo(_LD).eps):
        raise ValueError("prime-interval width exceeds the complete-prefix bound")
    delta = (theta - x) / x
    f, g, kernel_tail = _exponential_kernels(h, internal_degree)
    integrals = delta[:, None] * f - g
    prefactor = np.sqrt(x) * h / T
    coordinate = 2 * left / T - 1
    step = 2 * h / T

    coefficients, absolute_sums, tail_estimates = [], [], []
    previous = np.zeros((h.size, 1), dtype=_LD)
    polynomial = np.ones((h.size, 1), dtype=_LD)
    for k in range(internal_degree + 1):
        basis_scale = np.sqrt(_LD(2 * k + 1))
        terms = (basis_scale * polynomial * integrals[:, :k + 1]) * prefactor[:, None]
        coefficients.append(np.sum(terms, dtype=_LD))
        absolute_sums.append(np.sum(np.abs(terms), dtype=_LD))
        tails = np.abs(polynomial) * kernel_tail[:, :k + 1] * (np.abs(delta) + 2)[:, None]
        tail_estimates.append(basis_scale * np.sum(prefactor[:, None] * tails, dtype=_LD))
        if k < internal_degree:
            # Legendre recurrence with local argument coordinate + step*v.
            nxt = np.zeros((h.size, k + 2), dtype=_LD)
            nxt[:, :k + 1] += (2 * k + 1) * coordinate[:, None] * polynomial
            nxt[:, 1:k + 2] += (2 * k + 1) * step[:, None] * polynomial
            nxt[:, :previous.shape[1]] -= k * previous
            previous, polynomial = polynomial, nxt / (k + 1)

    all_coeffs = np.asarray(coefficients, dtype=_LD)
    absolute_sums = np.asarray(absolute_sums, dtype=_LD)
    tails = np.asarray(tail_estimates, dtype=_LD)
    norm = np.sum(_positive_norm_intervals(x, delta, h), dtype=_LD) / T
    eps = np.finfo(_LD).eps
    theta_error_indicator = eps * (prefix.size + 2) * np.sum(logs, dtype=_LD)
    coefficient_error = (eps * (64 + 8 * internal_degree) * absolute_sums
                         + 2 * np.sqrt(2 * np.arange(internal_degree + 1, dtype=_LD) + 1) / T * theta_error_indicator)
    e_minus_norm = -np.expm1(-T) / T
    norm_error = (128 * eps * norm
                  + 2 * theta_error_indicator * np.sqrt(norm * e_minus_norm)
                  + theta_error_indicator**2 * e_minus_norm)
    if not np.all(np.isfinite(all_coeffs)) or not np.isfinite(norm) or norm <= 0:
        raise ValueError("absolute moments or norm are unresolved in working precision")
    raw = all_coeffs[:degree + 1].copy()
    detrended = raw.copy()
    detrended[:min(2, degree + 1)] = 0
    energy = np.dot(raw, raw)
    affine_energy = np.dot(detrended, detrended)
    affine_norm = norm - all_coeffs[0]**2 - all_coeffs[1]**2
    subtraction_indicator = norm_error + np.sum(2 * np.abs(all_coeffs[:2]) * coefficient_error[:2])
    if energy > norm * (1 + max(_LD("1e-10"), 128 * eps)):
        raise ValueError("projection energy exceeds the norm; use a lower degree or higher precision")
    if affine_norm <= subtraction_indicator:
        affine_norm = None
        affine_ratio = None
        affine_status = "unresolved_by_working_precision_subtraction"
    else:
        affine_ratio = float(affine_energy / affine_norm)
        if affine_ratio > 1 + max(1e-10, 128 * float(eps)):
            affine_norm = None
            affine_ratio = None
            affine_status = "unresolved_by_working_precision_subtraction"
        else:
            affine_ratio = min(1.0, affine_ratio)
            affine_status = "continuous_orthogonal_P1_difference"
    cancellation = np.divide(absolute_sums[:degree + 1], np.abs(raw),
                             out=np.full(degree + 1, np.inf, dtype=_LD), where=raw != 0)
    metadata = {
        "object": "exact_continuous_theta_sqrt_residual",
        "evaluation": "prime_interval_analytic_identities_in_longdouble",
        "certified_interval_result": False,
        "error_metadata_role": "roundoff_indicators_and_series_tail_estimates_not_certified_bounds",
        "prime_prefix_status": provenance,
        "prime_limit": limit,
        "n_primes": int(prefix.size),
        "n_positive_intervals": int(h.size),
        "endpoint_x_decimal": _decimal(endpoint),
        "max_interval_log_width_decimal": _decimal(np.max(h)),
        "working_precision_bits": int(np.finfo(_LD).nmant + 1),
        "working_dtype": str(np.dtype(_LD)),
        "exponential_series_order": _SERIES_ORDER,
        "norm_series_tail_estimate_decimal": _decimal(np.sum(
            x * (2 * np.abs(delta) + 2) * np.exp(h) * h ** (_SERIES_ORDER + 1)
            / _LD(factorial(_SERIES_ORDER + 1)), dtype=_LD) / T),
        "moments_below_roundoff_indicator": [int(k) for k in range(degree + 1)
                                               if abs(raw[k]) <= coefficient_error[k]],
        "affine_detrend_norm_status": affine_status,
        "affine_norm_subtraction_indicator_decimal": _decimal(subtraction_indicator),
        "RH_proved": False,
        "B_ABS_growth_bound_proved": False,
    }
    return ArithmeticMoments(T, degree, raw, detrended, norm, affine_norm,
                             energy, affine_energy, min(1.0, float(energy / norm)),
                             affine_ratio, coefficient_error[:degree + 1], cancellation,
                             tails[:degree + 1], norm_error, metadata)
