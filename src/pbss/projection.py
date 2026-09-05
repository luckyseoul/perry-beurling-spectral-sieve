"""Weighted orthogonal projection onto shifted Legendre polynomials.

The continuous basis is orthonormal for du on [0, 1], but its sampled columns
are generally *not* orthonormal for a finite quadrature. We therefore fit the
polynomial by weighted QR, using the same weights for projection and norms::

    a = argmin_a sum_i w_i (q_i - (Phi @ a)_i)**2
    E_d = ||Phi @ a||_w**2 = a.T @ (Phi.T @ W @ Phi) @ a
    R_d = E_d / ||q||_w**2, in [0, 1]
    S_d = T**(2*(d+1)) * R_d

``coeffs`` are the expansion coefficients a, not raw quadrature moments and
not coordinates in the QR basis. Consequently E_d is not generally sum(a**2).
At quadrature rules integrating the basis products exactly, these agree with
the continuous orthonormal formulas. Trapezoid weights are the default;
nonnegative custom weights need not sum to one.

S_d is an archive-inspired heuristic scale, not a proved RH criterion. The
legacy P(q) normalization and classification threshold are unrecovered.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .basis import _validate_degree, orthonormal_legendre_design


@dataclass(frozen=True)
class ProjectionResult:
    """A discrete weighted projection; ``Phi @ coeffs`` reconstructs its fit."""

    degree: int
    coeffs: np.ndarray  # Legendre expansion coefficients a_0 .. a_d
    energy: float  # E_d = ||Phi @ coeffs||_w**2, not sum(coeffs**2)
    l2_norm_sq: float  # ||q||_w**2 via the same quadrature
    energy_ratio: float  # R_d = E_d / ||q||_w**2
    scaled_strength: float  # S_d = T**(2*(d+1)) * R_d
    T: float
    n_points: int

    @property
    def P(self) -> float:
        """Working heuristic strength S_d; not the recovered legacy P(q)."""
        return self.scaled_strength


def _real_vector(value: np.ndarray, name: str) -> np.ndarray:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real")
    try:
        out = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite one-dimensional array") from exc
    if out.ndim != 1 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite one-dimensional array")
    return out


def _validate_grid(u: np.ndarray) -> np.ndarray:
    u = _real_vector(u, "u")
    if u.size < 2:
        raise ValueError("need at least 2 sample points")
    if np.any((u < 0.0) | (u > 1.0)) or np.any(np.diff(u) <= 0.0):
        raise ValueError("u must be strictly increasing in [0, 1]")
    return u


def _trapezoid_weights(u: np.ndarray) -> np.ndarray:
    """Composite trapezoid weights on a strictly increasing grid in [0,1]."""
    u = _validate_grid(u)
    du = np.diff(u)
    w = np.empty(u.size, dtype=np.float64)
    w[0], w[-1] = 0.5 * du[0], 0.5 * du[-1]
    w[1:-1] = 0.5 * (du[:-1] + du[1:])
    return w


def _prepare_projection_inputs(q, u, degree, weights=None):
    degree = _validate_degree(degree)
    u = _validate_grid(u)
    q = _real_vector(q, "q")
    if q.shape != u.shape:
        raise ValueError("q and u must have the same shape")
    w = _trapezoid_weights(u) if weights is None else _real_vector(weights, "weights")
    if w.shape != u.shape:
        raise ValueError("weights must match u")
    if np.any(w < 0.0):
        raise ValueError("weights must be nonnegative")
    if np.count_nonzero(w > 0.0) < degree + 1:
        raise ValueError("need at least degree + 1 positive-weight sample points")
    return q, u, degree, w


@dataclass(frozen=True)
class _ProjectionFit:
    coeffs_unit: np.ndarray
    ratio: float
    signal_norm: np.longdouble  # ||sqrt(w / weight_scale) * q||_2
    weight_scale: float

    def coefficients(self) -> np.ndarray:
        coeffs = self.coeffs_unit.astype(np.longdouble) * self.signal_norm
        if np.any(np.abs(coeffs) > np.finfo(np.float64).max):
            raise ValueError("projection coefficients exceed float64 range; rescale q")
        return coeffs.astype(np.float64)

    def norm_sq(self) -> float:
        return _finite_energy(self.signal_norm**2 * self.weight_scale, "L2 norm squared")

    def energy(self) -> float:
        return _finite_energy(self.signal_norm**2 * self.weight_scale * self.ratio, "projection energy")


def _finite_energy(value: np.longdouble, name: str) -> float:
    if not np.isfinite(value) or value > np.finfo(np.float64).max:
        raise ValueError(f"{name} exceeds float64 range; rescale q or weights")
    result = float(value)
    if value > 0.0 and result == 0.0:
        raise ValueError(f"{name} underflows float64; rescale q or weights")
    return result


def _fit_projection(q, u, degree, weights, *, xp=np) -> _ProjectionFit:
    """QR fit of validated inputs on NumPy or CuPy, with a unit-norm target.

    Scaling the target and weights before factorization keeps ratios valid even
    when the unscaled squared norm would overflow or underflow float64.
    """
    positive = weights > 0.0
    w = weights[positive]
    values = q[positive]
    weight_scale = float(np.max(w))
    # sqrt first preserves a wider range than sqrt(w / max(w)).
    sqrt_w = np.sqrt(w) / np.sqrt(weight_scale)
    design = orthonormal_legendre_design(degree, u[positive])
    A = xp.asarray(design) * xp.asarray(sqrt_w)[:, None]
    try:
        Q, R = xp.linalg.qr(A, mode="reduced")
        singular = xp.linalg.svd(R, compute_uv=False)
    except np.linalg.LinAlgError as exc:
        raise ValueError("weighted polynomial factorization failed") from exc
    if not bool(xp.all(xp.isfinite(singular))):
        raise ValueError("weighted polynomial factorization produced non-finite values")
    cutoff = np.finfo(np.float64).eps * max(A.shape) * float(singular[0])
    if float(singular[-1]) <= cutoff:
        raise ValueError("weighted polynomial design is numerically rank deficient; reduce degree or improve grid/weights")

    q_scale = float(np.max(np.abs(values)))
    if q_scale == 0.0:
        return _ProjectionFit(np.zeros(degree + 1), 0.0, np.longdouble(0.0), weight_scale)
    target = xp.asarray(sqrt_w) * xp.asarray(values / q_scale)
    target_scale = float(xp.max(xp.abs(target)))
    if target_scale == 0.0:
        raise ValueError("weighted signal is numerically unresolved; rescale weights")
    target = target / target_scale
    target_norm = float(xp.linalg.norm(target))
    coordinates = Q.T @ (target / target_norm)
    ratio = float(xp.dot(coordinates, coordinates))
    if not np.isfinite(ratio):
        raise ValueError("weighted projection produced a non-finite ratio")
    # QR is orthonormal; this only removes floating-point overshoot at 1.
    ratio = min(1.0, max(0.0, ratio))
    try:
        coefficients = xp.linalg.solve(R, coordinates)
    except np.linalg.LinAlgError as exc:
        raise ValueError("weighted polynomial solve failed") from exc
    if not bool(xp.all(xp.isfinite(coefficients))):
        raise ValueError("weighted polynomial solve produced non-finite coefficients")
    if xp is not np:
        coefficients = xp.asnumpy(coefficients)
    signal_norm = np.longdouble(q_scale) * target_scale * target_norm
    return _ProjectionFit(coefficients, ratio, signal_norm, weight_scale)


def _validate_T(T: float) -> float:
    try:
        T = float(T)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("T must be finite and positive") from exc
    if not np.isfinite(T) or T <= 0.0:
        raise ValueError("T must be finite and positive")
    return T


def _scaled_ratio(ratio: float, degree: int, T: float) -> float:
    if ratio == 0.0:
        return 0.0
    # Combine logs before exponentiation: T**exponent may overflow even when S
    # itself is representable after multiplying a small R.
    log_strength = np.longdouble(2 * (degree + 1)) * np.log(np.longdouble(T)) + np.log(np.longdouble(ratio))
    if log_strength > np.log(np.longdouble(np.finfo(np.float64).max)):
        raise ValueError("scaled strength exceeds float64 range; reduce T or degree")
    with np.errstate(under="ignore"):
        result = float(np.exp(log_strength))
    if not np.isfinite(result):
        raise ValueError("scaled strength exceeds float64 range; reduce T or degree")
    if result == 0.0:
        raise ValueError("scaled strength underflows float64; change T or degree")
    return result


def project_coefficients(q: np.ndarray, u: np.ndarray, degree: int, weights: Optional[np.ndarray] = None) -> np.ndarray:
    """Weighted least-squares Legendre coefficients; reconstruct with Phi @ a.

    These are expansion coefficients, not the raw moments Phi.T @ (w*q).
    The grid is finite, one-dimensional, strictly increasing and in [0,1].
    Weights are finite and nonnegative, with at least degree+1 positive nodes.
    """
    inputs = _prepare_projection_inputs(q, u, degree, weights)
    return _fit_projection(*inputs).coefficients()


def projection_energy(q: np.ndarray, u: np.ndarray, degree: int, weights: Optional[np.ndarray] = None) -> float:
    """Weighted squared norm of the fitted polynomial (zero for a zero signal)."""
    inputs = _prepare_projection_inputs(q, u, degree, weights)
    return _fit_projection(*inputs).energy()


def energy_ratio(q: np.ndarray, u: np.ndarray, degree: int, weights: Optional[np.ndarray] = None) -> float:
    """Scale-invariant R_d in [0,1]; undefined for a signal of zero weighted norm."""
    inputs = _prepare_projection_inputs(q, u, degree, weights)
    fit = _fit_projection(*inputs)
    if fit.signal_norm == 0.0:
        raise ValueError("||q||^2 is zero; cannot form energy ratio")
    return fit.ratio


def scaled_projection_strength(q: np.ndarray, u: np.ndarray, degree: int, T: float, weights: Optional[np.ndarray] = None) -> float:
    """S_d = T**(2*(d+1))*R_d, the archive-inspired heuristic scale."""
    T = _validate_T(T)
    ratio = energy_ratio(q, u, degree, weights)
    return _scaled_ratio(ratio, _validate_degree(degree), T)


def project(q: np.ndarray, u: np.ndarray, degree: int, T: float = 1.0, weights: Optional[np.ndarray] = None) -> ProjectionResult:
    """Return weighted projection coefficients, E_d, R_d and heuristic S_d.

    Ratios remain stable under large/small amplitude changes. If the requested
    absolute energies or S_d cannot be represented as finite float64 values,
    raise ValueError; ``energy_ratio`` can still be used independently.
    """
    T = _validate_T(T)
    inputs = _prepare_projection_inputs(q, u, degree, weights)
    q, u, degree, _ = inputs
    fit = _fit_projection(*inputs)
    if fit.signal_norm == 0.0:
        raise ValueError("||q||^2 is zero; cannot form energy ratio")
    return ProjectionResult(
        degree=degree,
        coeffs=fit.coefficients(),
        energy=fit.energy(),
        l2_norm_sq=fit.norm_sq(),
        energy_ratio=fit.ratio,
        scaled_strength=_scaled_ratio(fit.ratio, degree, T),
        T=T,
        n_points=int(u.size),
    )
