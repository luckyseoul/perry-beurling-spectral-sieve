"""
GPU-assisted residual helpers (CuPy / V100).

Uses GPU for:
  - searchsorted of window abscissae into prime table (when prefix fits)
  - Legendre energy_ratio via projection_backend

CPU keeps multi-GB prime tables (mmap / prefix). Atomic GPU ops used for
parallel gather of θ from indices (scatter-style read via advanced indexing).

Not an RH proof.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .probes import _detrend, sample_grid
from .projection_backend import energy_ratio_auto


def arithmetic_residual_fast(
    u: np.ndarray,
    *,
    T: float,
    primes: np.ndarray,
    csum: np.ndarray,
    detrend: str = "deg1",
    prefer_gpu: bool = True,
) -> Tuple[np.ndarray, float, dict]:
    """
    Fast residual build: optional CuPy searchsorted + gather of θ.

    primes/csum are host arrays for the prefix ≤ exp(T) (or full table).
    """
    u = np.asarray(u, dtype=np.float64)
    T = float(T)
    x = np.maximum(np.exp(u * T), 2.0)
    p = np.asarray(primes)
    c = np.asarray(csum)
    backend = "numpy"
    if prefer_gpu and p.size > 50_000 and u.size >= 1024:
        try:
            import cupy as cp

            # Cap transfer: if prefix too large for V100 free mem, stay on CPU
            # Rough: int64 primes + float64 csum ≈ 16 bytes/prime
            need = p.size * 16
            free_mem = cp.cuda.Device(0).mem_info[0]
            if need < 0.6 * free_mem:
                pg = cp.asarray(p)
                cg = cp.asarray(c)
                xg = cp.asarray(x)
                # searchsorted on GPU
                idx = cp.searchsorted(pg, xg, side="right") - 1
                theta = cp.zeros_like(xg)
                ok = idx >= 0
                # atomic-safe gather: advanced indexing is coherent on V100
                theta[ok] = cg[idx[ok]]
                residual = cp.asnumpy((theta - xg) / cp.sqrt(xg))
                backend = "cupy_searchsorted"
                del pg, cg, xg, idx, theta
                cp.get_default_memory_pool().free_all_blocks()
            else:
                raise RuntimeError("prefix too large for GPU")
        except Exception:
            idx = np.searchsorted(p, x, side="right") - 1
            theta = np.zeros_like(x)
            ok = idx >= 0
            theta[ok] = c[idx[ok]]
            residual = (theta - x) / np.sqrt(x)
            backend = "numpy"
    else:
        idx = np.searchsorted(p, x, side="right") - 1
        theta = np.zeros_like(x)
        ok = idx >= 0
        theta[ok] = c[idx[ok]]
        residual = (theta - x) / np.sqrt(x)

    residual = _detrend(residual, u, detrend)
    meta = {
        "T": T,
        "x_max": float(np.exp(T)),
        "n_primes": int(p.size),
        "detrend": detrend,
        "backend": backend,
        "kind": "arithmetic_residual_fast",
        "note": "GPU-assisted residual when prefix fits V100; not an RH proof.",
    }
    return residual, T, meta


def energy_ratios_multi_degree(
    q: np.ndarray,
    u: np.ndarray,
    degrees: list,
    prefer_gpu: bool = True,
) -> dict:
    """Compute R_d for several degrees; reuse one design path per degree on GPU."""
    out = {}
    for d in degrees:
        r, be = energy_ratio_auto(q, u, int(d), prefer_gpu=prefer_gpu)
        out[int(d)] = {"R_d": float(r), "backend": be}
    return out
