"""Optional CuPy weighted QR projection, with an identical NumPy fallback."""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .projection import _fit_projection, _prepare_projection_inputs

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

        cp.cuda.Device(0).use()
        _ = cp.zeros(1)
        _CUPY = cp
        return True
    except Exception as exc:
        _CUPY = False
        _CUPY_ERR = str(exc)
        return False


def _fit_auto(q, u, degree, weights, prefer_gpu):
    # Validate before GPU fallback so malformed input has the same error on both
    # paths and never gets silently reinterpreted by a backend.
    inputs = _prepare_projection_inputs(q, u, degree, weights)
    if prefer_gpu and inputs[0].size >= 8192 and cupy_available():
        try:
            return _fit_projection(*inputs, xp=_CUPY), "cupy"
        except ValueError:
            raise
        except Exception:
            # CUDA memory/runtime failures may fall back; the reported backend
            # identifies the engine that actually completed the factorization.
            pass
    return _fit_projection(*inputs), "numpy"


def project_coefficients_auto(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    weights: Optional[np.ndarray] = None,
    prefer_gpu: bool = True,
) -> Tuple[np.ndarray, str]:
    """Return weighted least-squares Legendre coefficients and backend name."""
    fit, backend = _fit_auto(q, u, degree, weights, prefer_gpu)
    return fit.coefficients(), backend


def energy_ratio_auto(
    q: np.ndarray,
    u: np.ndarray,
    degree: int,
    prefer_gpu: bool = True,
    *,
    weights: Optional[np.ndarray] = None,
) -> Tuple[float, str]:
    """Return the true weighted projection ratio in [0,1] and backend name."""
    fit, backend = _fit_auto(q, u, degree, weights, prefer_gpu)
    if fit.signal_norm == 0.0:
        raise ValueError("||q||^2 is zero; cannot form energy ratio")
    return fit.ratio, backend
