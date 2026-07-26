"""
Beurling / generalized prime systems for PBSS battery scoring.

Two constructive families (finite, offline, no network):

1. **RH-like / PNT-friendly (integers):** ordinary primes — the classical system.
   Use arithmetic residual path (caller supplies primes).

2. **Defective / non-RH-like (gaps):** generalized primes with large regular gaps
   so θ_B(x) - x has a strong low-frequency / low-degree component on the
   log-window (persistent bulk).

3. **Smooth integers (highly composite-ish density):** primes thinned by keeping
   only every k-th ordinary prime — changes density and injects structure.

These are **diagnostic controls**, not a classification of all Beurling systems.
Not a proof of RH.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .probes import _detrend, sample_grid


def ordinary_primes_system(primes: np.ndarray) -> Dict:
    """Metadata wrapper: classical integers' primes."""
    p = np.asarray(primes)
    return {
        "name": "ordinary_primes",
        "kind": "rh_like",
        "description": "Classical rational primes (PNT / RH ambient system)",
        "n_primes": int(p.size),
        "p_max": float(p[-1]) if p.size else 0.0,
    }


def gapped_beurling_primes(
    x_max: float,
    gap: float = 3.0,
    p0: float = 2.0,
) -> np.ndarray:
    """
    Generalized primes: arithmetic progression-like sequence
      p_n = p0 + n * gap
    truncated at x_max. Strong regular spacing → defective residual structure.
    """
    x_max = float(x_max)
    gap = float(gap)
    if gap <= 0 or x_max < p0:
        return np.array([], dtype=np.float64)
    # number of terms
    n = int(np.floor((x_max - p0) / gap)) + 1
    p = p0 + gap * np.arange(n, dtype=np.float64)
    return p[p <= x_max]


def thinned_ordinary_primes(primes: np.ndarray, keep_every: int = 3) -> np.ndarray:
    """Keep every k-th ordinary prime (density defect)."""
    p = np.asarray(primes)
    k = int(keep_every)
    if k < 1:
        raise ValueError("keep_every >= 1")
    return p[::k].astype(np.float64, copy=False)


def beurling_theta_residual(
    u: np.ndarray,
    primes_b: np.ndarray,
    *,
    T: Optional[float] = None,
    x_max: Optional[float] = None,
    detrend: str = "deg1",
) -> Tuple[np.ndarray, float, dict]:
    """
    Chebyshev-style residual for a generalized prime sequence B:

      θ_B(x) = ∑_{p∈B, p≤x} log p
      q(u)   = detrend((θ_B(x) - x) / √x),  x = exp(u T)

    For non-integer primes, log p is still used (Beurling convention).
    """
    u = np.asarray(u, dtype=np.float64)
    p = np.asarray(primes_b, dtype=np.float64)
    p = p[np.isfinite(p) & (p >= 2.0)]
    p = np.sort(p)
    if T is None and x_max is None:
        raise ValueError("provide T or x_max")
    if T is None:
        x_max = float(x_max)
        T = float(np.log(x_max))
    else:
        T = float(T)
        x_max = float(np.exp(T))

    if p.size == 0:
        raise ValueError("empty Beurling prime set")

    # restrict to window
    hi = int(np.searchsorted(p, x_max, side="right"))
    p = p[:hi]
    if p.size == 0:
        raise ValueError("no Beurling primes ≤ x_max")

    csum = np.cumsum(np.log(p))
    x = np.maximum(np.exp(u * T), 2.0)
    idx = np.searchsorted(p, x, side="right") - 1
    theta = np.zeros_like(x)
    ok = idx >= 0
    theta[ok] = csum[idx[ok]]
    residual = (theta - x) / np.sqrt(x)
    residual = _detrend(residual, u, detrend)
    meta = {
        "T": T,
        "x_max": x_max,
        "n_primes": int(p.size),
        "detrend": detrend,
        "kind": "beurling_theta_residual",
        "note": "Generalized-prime residual; not an RH proof.",
    }
    return residual, T, meta


def default_battery_specs() -> List[dict]:
    """
    Fixed battery definitions for multi-T scoring.

    - ordinary: classical primes (rh_like)
    - gapped: regular gaps (defective)
    - thinned: every 3rd prime (defective density)
    """
    return [
        {
            "name": "ordinary_primes",
            "kind": "rh_like",
            "builder": "ordinary",
        },
        {
            "name": "gapped_gap3",
            "kind": "defective",
            "builder": "gapped",
            "gap": 3.0,
            "p0": 2.0,
        },
        {
            "name": "thinned_every3",
            "kind": "defective",
            "builder": "thinned",
            "keep_every": 3,
        },
    ]


def marathon_battery_specs(n_systems: int = 100) -> List[dict]:
    """
    Expand to ≥ n_systems Beurling constructions for overnight / open-plateau battery.

    Includes one ordinary (rh_like) plus many gapped / thinned / shifted-gap
    defective systems. Deterministic, offline, no network. Scales to tens of
    thousands of systems for deep multi-T stress.
    """
    n_systems = max(3, int(n_systems))
    specs: List[dict] = [
        {"name": "ordinary_primes", "kind": "rh_like", "builder": "ordinary"}
    ]
    # Gapped family: dense gap grid × many p0 seeds
    p0_seeds = (
        2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0,
        31.0, 37.0, 41.0, 43.0, 47.0, 53.0, 59.0, 61.0, 67.0, 71.0,
    )
    g = 1.25
    while len(specs) < n_systems and g < 5000.0:
        for p0 in p0_seeds:
            if len(specs) >= n_systems:
                break
            specs.append(
                {
                    "name": f"gapped_g{g:g}_p{p0:g}",
                    "kind": "defective",
                    "builder": "gapped",
                    "gap": float(g),
                    "p0": float(p0),
                }
            )
        if g < 5:
            g += 0.05
        elif g < 20:
            g += 0.1
        elif g < 100:
            g += 0.25
        elif g < 500:
            g += 0.5
        else:
            g += 1.0
    # Thinned ordinary (every k-th prime)
    k = 2
    while len(specs) < n_systems and k <= 100_000:
        specs.append(
            {
                "name": f"thinned_every{k}",
                "kind": "defective",
                "builder": "thinned",
                "keep_every": int(k),
            }
        )
        k += 1
    # Extra rh_like variants if still short (should be rare)
    extra = 0
    while len(specs) < n_systems and extra < 10:
        specs.append(
            {
                "name": f"ordinary_primes_v{extra}",
                "kind": "rh_like",
                "builder": "ordinary",
            }
        )
        extra += 1
    if len(specs) < n_systems:
        raise ValueError(f"could only build {len(specs)} systems < {n_systems}")
    return specs[:n_systems]


def build_system_primes(
    spec: dict,
    ordinary_primes: np.ndarray,
    x_max: float,
) -> np.ndarray:
    """Materialize generalized primes for one battery spec."""
    b = spec.get("builder", "ordinary")
    if b == "ordinary":
        p = np.asarray(ordinary_primes)
        hi = int(np.searchsorted(p, x_max, side="right"))
        return p[:hi].astype(np.float64, copy=False)
    if b == "gapped":
        return gapped_beurling_primes(
            x_max, gap=float(spec.get("gap", 3.0)), p0=float(spec.get("p0", 2.0))
        )
    if b == "thinned":
        p = thinned_ordinary_primes(
            ordinary_primes, keep_every=int(spec.get("keep_every", 3))
        )
        hi = int(np.searchsorted(p, x_max, side="right"))
        return p[:hi]
    raise ValueError(f"unknown builder {b}")
