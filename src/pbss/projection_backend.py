"""
Optional CuPy-accelerated projection with NumPy CPU fallback.

Large grids: build design matrix and inner products on GPU when available.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .basis import orthonormal_legendre_design
from .projection import _trapezoid_weights

_CUPY = None
_CUPY_ERR = None


def cupy_available() -> bool:
    global _CUPY, _CUPY_ERR
    if _CUPY is False:
        return False
    if _CUPY is not None:
        return True
    try:
        import cupy as cp

        # touch device
        cp.cuda.Device(0).use()
        _ = cp.zeros(1)
        _CUPY = cp
        return True
    except Exception as e:
        _CUPY = False
        _CUPY_ERR = str(e)
        return False


def project_coefficients_auto(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    weights: Optional[np.ndarray] = None,
    prefer_gpu: bool = True,
) -> Tuple[np.ndarray, str]:
    """
    Returns (coeffs, backend_name) with backend in {'numpy','cupy'}.
    """
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if weights is None:
        weights = _trapezoid_weights(u)
    else:
        weights = np.asarray(weights, dtype=np.float64).ravel()

    if prefer_gpu and q.size >= 8192:
        try:
            # Re-check in this process (forked workers must re-init CUDA)
            import cupy as cp

            Phi = orthonormal_legendre_design(degree, u)
            Phi_g = cp.asarray(Phi)
            wq = cp.asarray(weights * q)
            c = Phi_g.T @ wq
            return cp.asnumpy(c), "cupy"
        except Exception:
            pass  # fall through to NumPy

    Phi = orthonormal_legendre_design(degree, u)
    return Phi.T @ (weights * q), "numpy"


def energy_ratio_auto(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    prefer_gpu: bool = True,
) -> Tuple[float, str]:
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    w = _trapezoid_weights(u)
    c, backend = project_coefficients_auto(q, u, degree, weights=w, prefer_gpu=prefer_gpu)
    e = float(np.dot(c, c))
    l2 = float(np.sum(w * q * q))
    if l2 <= 0:
        raise ValueError("zero L2 norm")
    return e / l2, backend
