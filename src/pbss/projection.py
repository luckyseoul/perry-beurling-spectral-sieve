"""
Projection of a density perturbation q onto low-degree shifted Legendre space.

Formulas (unit log-window [0,1], measure du)
--------------------------------------------
Let {φ_k}_{k=0}^d be the orthonormal shifted Legendre basis (see basis.py).
For a (discrete) sample of q at quadrature nodes u_i with weights w_i
approximating ∫_0^1 · du:

  c_k  = ∑_i w_i q(u_i) φ_k(u_i)     ≈ ⟨q, φ_k⟩
  E_d  = ∑_{k=0}^d |c_k|²            projection energy  ||P_d q||²
  R_d  = E_d / ||q||²                 L² energy ratio ∈ [0,1]
  S_d  = T^{2(d+1)} R_d               scaled strength (Thm A heuristic scale)

Under the conditional decay R_d(q_T) = O(T^{-2(d+1)}) suggested by the
archive notes (RH ⇒ projection energy dies), S_d stays O(1). Defective
systems inject low-degree mass and inflate R_d and S_d.

P(q) naming
-----------
The archive quotes P(q)≈3.92 for zeta and a classification threshold ≈29.5.
Those figures used lost high-precision scripts and an incompletely documented
normalization. This module therefore exposes the **honest** primitives E_d,
R_d, S_d and does **not** hard-code legacy numerics. Callers may treat S_d
(or E_d after unit-normalizing q) as the working "projection strength."
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .basis import orthonormal_legendre_design


@dataclass(frozen=True)
class ProjectionResult:
    """Output of a finite Legendre projection of q."""

    degree: int
    coeffs: np.ndarray  # c_0 .. c_d
    energy: float  # E_d = sum |c_k|²
    l2_norm_sq: float  # ||q||² via the same quadrature
    energy_ratio: float  # R_d = E_d / ||q||²
    scaled_strength: float  # S_d = T^{2(d+1)} R_d
    T: float  # window parameter used for scaling
    n_points: int

    @property
    def P(self) -> float:
        """
        Working projection-strength scalar for this reconstruction.

        Defined as the scaled L² energy ratio S_d = T^{2(d+1)} R_d,
        which is O(1) under the archive's RH decay heuristic and larger
        when low-degree mass is present.
        """
        return self.scaled_strength


def _trapezoid_weights(u: np.ndarray) -> np.ndarray:
    """Composite trapezoid weights for nonuniform or uniform samples on [0,1]."""
    u = np.asarray(u, dtype=np.float64)
    n = u.size
    if n < 2:
        raise ValueError("need at least 2 sample points")
    w = np.zeros(n, dtype=np.float64)
    # general nonuniform trapezoid
    du = np.diff(u)
    w[0] = 0.5 * du[0]
    w[-1] = 0.5 * du[-1]
    if n > 2:
        w[1:-1] = 0.5 * (du[:-1] + du[1:])
    return w


def project_coefficients(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Discrete inner products c_k = ⟨q, φ_k⟩ via quadrature.

    Parameters
    ----------
    q : samples of the density perturbation at abscissae u
    u : sample points in [0, 1] (log-window coordinate)
    degree : max polynomial degree d (uses φ_0..φ_d)
    weights : optional quadrature weights; default trapezoid on u
    """
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if q.shape != u.shape:
        raise ValueError("q and u must have the same shape")
    if weights is None:
        weights = _trapezoid_weights(u)
    else:
        weights = np.asarray(weights, dtype=np.float64).ravel()
        if weights.shape != u.shape:
            raise ValueError("weights must match u")

    Phi = orthonormal_legendre_design(degree, u)  # (n, d+1)
    # c = Phi^T (w ⊙ q)
    return Phi.T @ (weights * q)


def projection_energy(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    weights: Optional[np.ndarray] = None,
) -> float:
    """E_d = ∑_{k=0}^d |c_k|²."""
    c = project_coefficients(q, u, degree, weights=weights)
    return float(np.dot(c, c))


def energy_ratio(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    weights: Optional[np.ndarray] = None,
) -> float:
    """R_d = E_d / ||q||² with the same quadrature for ||q||²."""
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if weights is None:
        weights = _trapezoid_weights(u)
    else:
        weights = np.asarray(weights, dtype=np.float64).ravel()
    l2 = float(np.sum(weights * q * q))
    if l2 <= 0.0:
        raise ValueError("||q||^2 is zero; cannot form energy ratio")
    e = projection_energy(q, u, degree, weights=weights)
    return e / l2


def scaled_projection_strength(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    T: float,
    weights: Optional[np.ndarray] = None,
) -> float:
    """
    S_d = T^{2(d+1)} R_d.

    Scaling matches the archive note that under RH the projection energy
    decays like O(T^{-2(d+1)}); S_d is then the natural O(1) statistic.
    """
    if T <= 0:
        raise ValueError("T must be positive")
    r = energy_ratio(q, u, degree, weights=weights)
    return float((T ** (2 * (degree + 1))) * r)


def project(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    T: float = 1.0,
    weights: Optional[np.ndarray] = None,
) -> ProjectionResult:
    """Full projection package: coefficients, E_d, R_d, S_d."""
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if weights is None:
        weights = _trapezoid_weights(u)
    else:
        weights = np.asarray(weights, dtype=np.float64).ravel()

    c = project_coefficients(q, u, degree, weights=weights)
    e = float(np.dot(c, c))
    l2 = float(np.sum(weights * q * q))
    if l2 <= 0.0:
        raise ValueError("||q||^2 is zero")
    r = e / l2
    s = float((T ** (2 * (degree + 1))) * r)
    return ProjectionResult(
        degree=degree,
        coeffs=c,
        energy=e,
        l2_norm_sq=l2,
        energy_ratio=r,
        scaled_strength=s,
        T=float(T),
        n_points=int(u.size),
    )
