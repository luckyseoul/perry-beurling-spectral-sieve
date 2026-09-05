"""Continuous, phase-aware checks for one exponentially weighted real sine.

For fixed A,t>0, q_T(u)=exp(A*T*u)*sin(t*T*u) and the continuous Legendre
projector, the exact identities below imply R_d=O(1/T), uniformly in phase.
The old phase-independent expression is the asymptotic upper envelope:

    limsup T*R_d = 4*A*(d+1)**2/(A*A+t*t),    liminf T*R_d = 0.

These are model statements, not an arithmetic counterexample or a solution of
B-RES or RH. "Exact" refers to analytic identities; this module evaluates them
in float64, without interval certification. The optional experiment checks an
independent finite antiderivative formula with mpmath at higher precision.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import factorial

import numpy as np
from scipy.special import hyp1f1, ive

from .basis import _validate_degree


@dataclass(frozen=True)
class OffCriticalProjection:
    T: float
    A: float
    t: float
    degree: int
    scaled_coeffs: np.ndarray  # exp(-A*T) * continuous Legendre coefficients
    scaled_norm_sq: float  # exp(-2*A*T) * ||q_T||_2**2
    scaled_projection_energy: float
    energy_ratio: float
    phase_leading_ratio: float  # additive error O(T^-2), not relative at zeros
    phase_envelope_ratio: float  # the limsup constant divided by T
    evaluation: str = "analytic_identities_evaluated_in_float64"


def _parameters(T: float, A: float, t: float, degree: int):
    degree = _validate_degree(degree)
    try:
        T, A, t = float(T), float(A), float(t)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("T, A and t must be finite and positive") from exc
    if not all(np.isfinite(v) and v > 0.0 for v in (T, A, t)):
        raise ValueError("T, A and t must be finite and positive")
    a, b = A * T, t * T
    if not np.isfinite(a) or not np.isfinite(b):
        raise ValueError("A*T and t*T must be finite")
    return T, A, t, degree, a, b


def _phase_terms(A: float, t: float, theta: float):
    radius = float(np.hypot(A, t))
    s, c = np.sin(theta), np.cos(theta)
    numerator = (A / radius) * s - (t / radius) * c
    B = numerator / radius
    # Equivalent to 1/A-Re(exp(2j*theta)/(A+it)), without cancellation
    # between two positive terms when t/A is small.
    H = (numerator * numerator + s * s) / A
    return B, H, radius


def _scaled_norm(a: float, b: float) -> float:
    if max(a, b) < 0.5:
        # sin²(bu) series avoids cancellation in the two-exponential formula
        # for small T. Integrals of u^(2j)*exp(2au) use Kummer's identity.
        total = 0.0
        for j in range(1, 32):
            term = ((-1.0) ** (j + 1) * (2.0 * b) ** (2 * j)
                    / (2.0 * factorial(2 * j) * (2 * j + 1))
                    * float(hyp1f1(2 * j + 1, 2 * j + 2, 2.0 * a)))
            total += term
            if abs(term) <= np.finfo(float).eps * abs(total):
                break
        norm = float(np.exp(-2.0 * a) * total)
    else:
        _, H, radius = _phase_terms(a, b, b)
        H0 = (b / radius) ** 2 / a
        norm = float((H - np.exp(-2.0 * a) * H0) / 4.0)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("norm is unresolved in float64; use the high-precision checker")
    return norm


def off_critical_projection(T: float, *, A: float, t: float, degree: int) -> OffCriticalProjection:
    """Evaluate the continuous projection using exponentially scaled identities.

    With z=T(A+it), the complex Legendre moment is
      sqrt(2k+1)*exp(z/2)*sqrt(pi/z)*I_(k+1/2)(z/2).
    ``ive`` removes exp(Re(z)/2), so no exp(A*T) is formed. Taking the
    imaginary part yields the real sine coefficients. This is independent of
    the finite endpoint-derivative formula used by the mpmath checker.
    """
    T, A, t, degree, a, b = _parameters(T, A, t, degree)
    z = complex(a, b)
    k = np.arange(degree + 1, dtype=float)
    moments = (np.sqrt(2.0 * k + 1.0) * np.sqrt(np.pi / z)
               * np.exp(0.5j * b) * ive(k + 0.5, 0.5 * z))
    coeffs = np.asarray(moments.imag, dtype=float)
    norm = _scaled_norm(a, b)
    energy = float(np.dot(coeffs, coeffs))
    ratio = energy / norm
    if not np.all(np.isfinite(coeffs)) or not np.isfinite(ratio) or ratio > 1.0 + 1e-10:
        raise ValueError("projection is unresolved in float64; use the high-precision checker")
    B, H, radius = _phase_terms(A, t, b)
    leading = 4.0 * (degree + 1) ** 2 * B * B / (T * H)
    envelope = 4.0 * (A / radius) * (degree + 1) ** 2 / (radius * T)
    return OffCriticalProjection(T, A, t, degree, coeffs, norm, energy,
                                 min(1.0, ratio), leading, envelope)


def asymptotic_constants(*, A: float, t: float, degree: int) -> dict[str, float]:
    """Sharp phase subsequence constants for the continuous real-sine model.

    For d>=1, liminf(T^3*R_d) is positive. For d=0, exact zero-coefficient
    phases exist and the liminf is zero. ``leading_cancel_T3`` describes
    T=(arg(A+it)+n*pi)/t; it is four times the optimized liminf for d>=1.
    """
    _, A, t, degree, _, _ = _parameters(1.0, A, t, degree)
    radius = float(np.hypot(A, t))
    factor = float(degree * (degree + 1) * (degree + 2))
    minimum = (A / radius) * (factor / radius) ** 2 / (3.0 * radius)
    return {
        "limsup_T_R": 4.0 * (A / radius) * (degree + 1) ** 2 / radius,
        "liminf_T_R": 0.0,
        "leading_cancel_T3_R": 4.0 * minimum,
        "liminf_T3_R": minimum,
    }


def phase_subsequence_T(n: int, *, A: float, t: float, degree: int = 0,
                        kind: str = "leading_cancel") -> float:
    """Return a window on a stated phase subsequence (n is a positive integer)."""
    n = _validate_degree(n)
    if n == 0:
        raise ValueError("n must be a positive integer")
    _, A, t, degree, _, _ = _parameters(1.0, A, t, degree)
    if kind == "maximal":
        return float(n * np.pi / t)
    theta = float(np.arctan2(t, A) + n * np.pi)
    if kind == "leading_cancel":
        return theta / t
    if kind == "minimum_T3":
        radius = float(np.hypot(A, t))
        mu = degree * (degree + 2) / 2.0
        discriminant = theta * theta - 4.0 * (t / radius) ** 2 * mu
        if discriminant <= 0.0:
            raise ValueError("n is too small for the minimum-T3 phase sequence")
        return (theta + np.sqrt(discriminant)) / (2.0 * t)
    raise ValueError("kind must be maximal, leading_cancel or minimum_T3")


def partial_cancellation_ratio(ratio: float, alpha: float) -> float:
    """R(q-alpha*Pq)=(1-alpha)^2*R/[1-(2alpha-alpha^2)*R]."""
    if not np.isfinite(ratio) or not 0.0 <= ratio <= 1.0:
        raise ValueError("ratio must be finite and in [0,1]")
    if not np.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be finite and in [0,1]")
    numerator = (1.0 - alpha) ** 2 * ratio
    denominator = (1.0 - ratio) + numerator
    if denominator == 0.0:
        raise ValueError("full cancellation leaves zero norm")
    return numerator / denominator


def detrended_ratio(T: float, *, A: float, t: float, degree: int, remove_degree: int) -> float:
    """Continuous orthogonal detrending: R_d((I-P_m)q).

    This identity assumes the same continuous inner product for detrending and
    projection; it does not assert that an arbitrary sampled regression does.
    """
    degree, remove_degree = _validate_degree(degree), _validate_degree(remove_degree)
    result = off_critical_projection(T, A=A, t=t, degree=max(degree, remove_degree))
    removed = float(np.dot(result.scaled_coeffs[:remove_degree + 1], result.scaled_coeffs[:remove_degree + 1]))
    remaining_norm = result.scaled_norm_sq - removed
    if remaining_norm <= 0.0:
        raise ValueError("detrended norm is unresolved in float64")
    retained = result.scaled_coeffs[remove_degree + 1:degree + 1]
    return float(np.dot(retained, retained) / remaining_norm)
