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


def _primes_upto(n: int) -> np.ndarray:
    """Simple sieve of Eratosthenes; returns primes ≤ n."""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def probe_prime_residual(
    u: np.ndarray,
    x_max: float = 1e5,
    smooth: int = 3,
) -> Tuple[np.ndarray, float]:
    """
    Build a Chebyshev-function residual probe from real primes.

    On the multiplicative scale x = exp(u * log x_max):
      θ(x) = ∑_{p≤x} log p
      q(u) ∝ (θ(x) - x) / x^{1/2}     (normalized later)

    Under RH, θ(x)-x = O(x^{1/2} log² x); the residual is oscillatory and
    high-frequency on the log window — an RH-consistent probe.

    Returns
    -------
    q : residual samples on u
    T : logarithmic window length log(x_max)
    """
    u = np.asarray(u, dtype=np.float64)
    if x_max < 3:
        raise ValueError("x_max too small")
    T = float(np.log(x_max))
    x = np.exp(u * T)
    # avoid x<2
    x = np.maximum(x, 2.0)

    primes = _primes_upto(int(x_max))
    logp = np.log(primes.astype(np.float64))
    # θ at each x via searchsorted cumulative sum
    csum = np.cumsum(logp)
    idx = np.searchsorted(primes, x, side="right") - 1
    theta = np.where(idx >= 0, csum[np.clip(idx, 0, len(csum) - 1)], 0.0)
    theta = np.where(idx < 0, 0.0, theta)

    residual = (theta - x) / np.sqrt(x)
    # light moving average to reduce pure step noise (optional)
    if smooth > 1:
        kernel = np.ones(smooth) / smooth
        residual = np.convolve(residual, kernel, mode="same")
    # remove mean (degree-0 bulk) so the diagnostic sees oscillatory content;
    # raw staircase bias otherwise floods low-degree Legendre modes at modest x_max
    residual = residual - float(np.mean(residual))
    return residual, T


def default_T_from_xmax(x_max: float) -> float:
    return float(np.log(x_max))
