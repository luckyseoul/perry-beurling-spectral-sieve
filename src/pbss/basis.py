"""
Orthonormal shifted Legendre basis on the unit log-window [0, 1].

Standard Legendre polynomials L_k are orthogonal on [-1, 1] with weight 1
and ||L_k||_{L2[-1,1]}^2 = 2/(2k+1).

Map u ∈ [0,1] → t = 2u-1 ∈ [-1,1]. Then
  φ_k(u) = √(2k+1) · L_k(2u-1)
are orthonormal on [0,1] with measure du:
  ∫_0^1 φ_j(u) φ_k(u) du = δ_{jk}.

This is the finite polynomial space used for the P(q) / energy-ratio diagnostic.
"""
from __future__ import annotations

import numpy as np
from scipy.special import eval_legendre
from numbers import Integral


def _validate_degree(degree: int) -> int:
    if isinstance(degree, (bool, np.bool_)) or not isinstance(degree, Integral) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    return int(degree)


def _basis_points(u: np.ndarray) -> np.ndarray:
    if np.iscomplexobj(u):
        raise ValueError("u must be real")
    u = np.asarray(u, dtype=np.float64)
    if not np.all(np.isfinite(u)) or np.any((u < 0.0) | (u > 1.0)):
        raise ValueError("u must contain finite points in [0, 1]")
    return u


def shifted_legendre_values(k: int, u: np.ndarray) -> np.ndarray:
    """
    Evaluate φ_k(u) = √(2k+1) L_k(2u-1) at points u ∈ [0,1].
    """
    k = _validate_degree(k)
    u = _basis_points(u)
    t = 2.0 * u - 1.0
    return np.sqrt(2.0 * k + 1.0) * eval_legendre(k, t)


def orthonormal_legendre_design(degree: int, u: np.ndarray) -> np.ndarray:
    """
    Design matrix Φ of shape (n_points, degree+1) with columns φ_0,...,φ_d
    evaluated at sample abscissae u.
    """
    degree = _validate_degree(degree)
    u = _basis_points(u)
    if u.ndim != 1 or u.size == 0:
        raise ValueError("u must be a nonempty one-dimensional array")
    cols = [shifted_legendre_values(k, u) for k in range(degree + 1)]
    return np.column_stack(cols)
