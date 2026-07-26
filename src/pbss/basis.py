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


def shifted_legendre_values(k: int, u: np.ndarray) -> np.ndarray:
    """
    Evaluate φ_k(u) = √(2k+1) L_k(2u-1) at points u ∈ [0,1].
    """
    u = np.asarray(u, dtype=np.float64)
    t = 2.0 * u - 1.0
    return np.sqrt(2.0 * k + 1.0) * eval_legendre(k, t)


def orthonormal_legendre_design(degree: int, u: np.ndarray) -> np.ndarray:
    """
    Design matrix Φ of shape (n_points, degree+1) with columns φ_0,...,φ_d
    evaluated at sample abscissae u.
    """
    if degree < 0:
        raise ValueError("degree must be >= 0")
    u = np.asarray(u, dtype=np.float64)
    cols = [shifted_legendre_values(k, u) for k in range(degree + 1)]
    return np.column_stack(cols)
