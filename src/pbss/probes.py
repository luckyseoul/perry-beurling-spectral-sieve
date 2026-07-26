"""
Synthetic and prime-based density-perturbation probes on the unit log-window.

Coordinate: u ∈ [0,1] stands for the normalized log-abscissa
  u = log(x) / log(X_max)  with  x ∈ [1, X_max],  T := log(X_max)
(so T is the logarithmic window length used for scaled strength S_d).

Probes
------
- low_degree: pure φ_m (should put nearly all energy in degree ≤ m)
- high_frequency: rapid oscillation (RH-like: mostly orthogonal to low poly space)
- defective: high-frequency + large low-degree contamination (non-RH control)
- prime_residual: Chebyshev-style residual from real primes (RH-consistent probe)
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .basis import shifted_legendre_values


def sample_grid(n: int = 2048) -> np.ndarray:
    """Uniform grid on [0,1] with n ≥ 2 points (endpoints included)."""
    if n < 2:
        raise ValueError("n >= 2 required")
    return np.linspace(0.0, 1.0, n, dtype=np.float64)


def normalize_l2(q: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Return q / ||q||_L2 so energy ratios are pure subspace fractions."""
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    du = np.diff(u)
    w = np.zeros_like(q)
    w[0] = 0.5 * du[0]
    w[-1] = 0.5 * du[-1]
    if q.size > 2:
        w[1:-1] = 0.5 * (du[:-1] + du[1:])
    nrm = np.sqrt(float(np.sum(w * q * q)))
    if nrm <= 0:
        raise ValueError("zero L2 norm")
    return q / nrm


def probe_low_degree(u: np.ndarray, k: int = 2) -> np.ndarray:
    """Exact basis mode φ_k — projection of degree ≥ k captures all energy."""
    return shifted_legendre_values(k, u)


def probe_high_frequency(u: np.ndarray, waves: int = 40, phase: float = 0.0) -> np.ndarray:
    """
    Rapid sinusoid on [0,1]: mostly orthogonal to low-degree polynomials
    (RH-like / high-frequency spectral signature).
    """
    u = np.asarray(u, dtype=np.float64)
    return np.sin(2.0 * np.pi * waves * u + phase)


def probe_critical_line_mode(
    u: np.ndarray,
    T: float = 20.0,
    t: float = 14.134725,  # Im(ρ) of the first zeta zero
    phase: float = 0.0,
) -> np.ndarray:
    """
    Model residual mode cos/sin(t log x) on the log window — the shape
    associated with a critical-line zero at 1/2 + it (RH-consistent form).

    With x = exp(u T),  t log x = t T u.
    """
    u = np.asarray(u, dtype=np.float64)
    return np.sin(t * T * u + phase)


def probe_off_critical_mode(
    u: np.ndarray,
    T: float = 20.0,
    t: float = 14.134725,
    sigma: float = 0.75,
    phase: float = 0.0,
) -> np.ndarray:
    """
    Model contribution of a zero at σ + it with σ ≠ 1/2:

      x^{σ - 1/2} sin(t log x)  on  x = exp(u T)

    becomes  exp(α u) sin(ω u) with α = T(σ - 1/2), ω = t T.

    This is the standard off-line envelope in the log-window coordinate.
    """
    u = np.asarray(u, dtype=np.float64)
    alpha = T * (sigma - 0.5)
    omega = t * T
    return np.exp(alpha * u) * np.sin(omega * u + phase)


def probe_persistent_defect(
    u: np.ndarray,
    eps: float = 0.5,
    j: int = 0,
    waves: int = 80,
) -> np.ndarray:
    """
    Lemma M2 family: fixed low-degree mass ε independent of frequency.
    Uses lemmas.synthetic_orthogonal_defect.
    """
    from .lemmas import synthetic_orthogonal_defect

    return synthetic_orthogonal_defect(u, eps=eps, j=j, waves=waves)


def probe_defective(
    u: np.ndarray,
    waves: int = 40,
    defect_degree: int = 1,
    defect_weight: float = 2.0,
) -> np.ndarray:
    """
    High-frequency oscillation contaminated by a strong low-degree mode.
    Models a spectral defect (extra low-frequency density mass) — non-RH control.
    """
    u = np.asarray(u, dtype=np.float64)
    return probe_high_frequency(u, waves=waves) + defect_weight * shifted_legendre_values(
        defect_degree, u
    )


def primes_upto(n: int, segment_size: int = 10_000_000) -> np.ndarray:
    """
    Sieve of Eratosthenes; returns primes ≤ n as int64 array.

    For n ≥ 50_000_000 uses a segmented sieve so peak RAM stays ~O(√n + segment)
    rather than a full bool array of length n (important for x_max ~ 1e8–1e9).
    """
    n = int(n)
    if n < 2:
        return np.array([], dtype=np.int64)

    # Small-n: classic dense sieve
    if n < 50_000_000:
        sieve = np.ones(n + 1, dtype=bool)
        sieve[0:2] = False
        r = int(n**0.5)
        for p in range(2, r + 1):
            if sieve[p]:
                sieve[p * p :: p] = False
        return np.nonzero(sieve)[0].astype(np.int64)

    # Segmented sieve
    r = int(n**0.5) + 1
    base = primes_upto(r, segment_size=segment_size)  # recursive small path
    out: list[np.ndarray] = [base]
    seg = int(segment_size)
    low = r + 1
    while low <= n:
        high = min(low + seg - 1, n)
        mark = np.ones(high - low + 1, dtype=bool)
        for p in base:
            p = int(p)
            # first multiple of p at or above low
            start = ((low + p - 1) // p) * p
            if start < p * p:
                start = p * p
            if start > high:
                continue
            mark[start - low : high - low + 1 : p] = False
        out.append((np.nonzero(mark)[0] + low).astype(np.int64))
        low = high + 1
    return np.concatenate(out)


# backward-compatible alias
_primes_upto = primes_upto


def _detrend(residual: np.ndarray, u: np.ndarray, mode: str) -> np.ndarray:
    """
    Remove slow bulk from arithmetic residual without projecting out all of V_d.

    - none:   raw residual
    - deg0:   subtract mean
    - deg1:   subtract least-squares a + b u  (default for multi-T)
    """
    mode = mode.lower()
    r = np.asarray(residual, dtype=np.float64).copy()
    u = np.asarray(u, dtype=np.float64)
    if mode in ("none", "raw"):
        return r
    if mode in ("deg0", "mean"):
        return r - float(np.mean(r))
    if mode in ("deg1", "linear"):
        # least squares fit a + b u
        A = np.column_stack([np.ones_like(u), u])
        coef, _, _, _ = np.linalg.lstsq(A, r, rcond=None)
        return r - (A @ coef)
    raise ValueError(f"unknown detrend mode: {mode}")


def arithmetic_residual(
    u: np.ndarray,
    *,
    T: Optional[float] = None,
    x_max: Optional[float] = None,
    primes: Optional[np.ndarray] = None,
    detrend: str = "deg1",
    smooth: int = 1,
) -> Tuple[np.ndarray, float, dict]:
    """
    Arithmetic Chebyshev residual on the unit log-window (shipped entry point).

    Window: x = exp(u * T) for u in [0,1], with T = log(x_max).
      θ(x) = ∑_{p ≤ x} log p
      raw(u) = (θ(x) - x) / √x
      q(u)   = detrend(raw)   [optional light moving-average smooth first]

    Parameters
    ----------
    u : sample grid on [0,1]
    T, x_max : specify one; T = log(x_max)
    primes : optional precomputed primes (must cover ≤ x_max); speeds multi-T
    detrend : 'none' | 'deg0' | 'deg1' (default deg1 removes slow linear bulk)
    smooth : moving-average width in samples (1 = off)

    Returns
    -------
    q : residual samples
    T : window length
    meta : dict with x_max, n_primes, detrend, smooth
    """
    u = np.asarray(u, dtype=np.float64)
    if T is None and x_max is None:
        raise ValueError("provide T or x_max")
    if T is None:
        x_max = float(x_max)
        if x_max < 3:
            raise ValueError("x_max too small")
        T = float(np.log(x_max))
    else:
        T = float(T)
        if T <= 0:
            raise ValueError("T must be positive")
        x_max = float(np.exp(T))

    x_max_i = int(np.floor(x_max))
    if primes is None:
        primes = primes_upto(x_max_i)
    else:
        primes = np.asarray(primes, dtype=np.int64)
        # keep only primes in range
        primes = primes[primes <= x_max_i]

    x = np.exp(u * T)
    x = np.maximum(x, 2.0)

    if primes.size == 0:
        raise ValueError("no primes in range")
    logp = np.log(primes.astype(np.float64))
    csum = np.cumsum(logp)
    idx = np.searchsorted(primes, x, side="right") - 1
    theta = np.zeros_like(x)
    ok = idx >= 0
    theta[ok] = csum[idx[ok]]

    residual = (theta - x) / np.sqrt(x)
    if smooth > 1:
        kernel = np.ones(smooth, dtype=np.float64) / float(smooth)
        residual = np.convolve(residual, kernel, mode="same")
    residual = _detrend(residual, u, detrend)

    meta = {
        "x_max": x_max,
        "T": T,
        "n_primes": int(primes.size),
        "detrend": detrend,
        "smooth": int(smooth),
    }
    return residual, T, meta


def probe_prime_residual(
    u: np.ndarray,
    x_max: float = 1e5,
    smooth: int = 3,
    detrend: str = "deg1",
    primes: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, float]:
    """
    Backward-compatible wrapper around ``arithmetic_residual``.

    Returns (q, T) only. Prefer ``arithmetic_residual`` for new code.
    """
    q, T, _ = arithmetic_residual(
        u, x_max=x_max, primes=primes, detrend=detrend, smooth=smooth
    )
    return q, T


def default_T_from_xmax(x_max: float) -> float:
    return float(np.log(x_max))


def xmax_from_T(T: float) -> float:
    return float(np.exp(T))