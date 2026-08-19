"""
Gamma-shaped measure sensitivity for off-line vs on-line residual discrimination.

Reproduces the project-record sensitivity experiment:
  T_w(q) = ∫_0^1 w(s) q(s)^2 ds
with flat w≡1 vs Gamma bump w(s)=s^{k-1} e^{-σ s} (normalized), using Fisher's d′
between synthetic on-line and off-line ensembles.

Historical claim (~53% gain on a noisy Beurling-like ensemble; idealized
seed=20260522, k=4, σ=6 → flat d′≈2.37, Gamma d′≈5.53) is re-run here.

**Not RH. Engineering sensitivity result only.**
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

from .probes import (
    probe_critical_line_mode,
    probe_defective,
    probe_high_frequency,
    probe_off_critical_mode,
    sample_grid,
)

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

# Project-record idealized reference (docs/PROJECT_RECORD § survivor note)
HISTORICAL_IDEALIZED = {
    "seed": 20260522,
    "k": 4,
    "sigma": 6.0,
    "flat_dprime": 2.37,
    "gamma_dprime": 5.53,
    "source": "docs/PROJECT_RECORD_2025-11_to_2026-06.md (historical note)",
}


def gamma_weight(
    u: np.ndarray,
    *,
    k: float = 4.0,
    sigma: float = 6.0,
    normalize: bool = True,
) -> np.ndarray:
    """
    Gamma-shaped weight w(u) = u^{k-1} exp(-σ u) on [0,1].

    Peaks at u⋆ = (k-1)/σ when that lies in (0,1) (e.g. k=4, σ=6 → u⋆=0.5).
    """
    u = np.asarray(u, dtype=np.float64).ravel()
    k = float(k)
    sigma = float(sigma)
    if k < 1.0:
        raise ValueError("k must be >= 1")
    if sigma <= 0.0:
        raise ValueError("sigma must be > 0")
    # avoid 0^{k-1} issues for k<1; for k>=1, u=0 is fine (0 if k>1)
    w = np.power(np.maximum(u, 0.0), k - 1.0) * np.exp(-sigma * u)
    w = np.asarray(w, dtype=np.float64)
    if normalize:
        du = float(u[1] - u[0]) if u.size > 1 else 1.0
        mass = float(np.sum(w) * du)
        if mass <= 0.0:
            raise ValueError("gamma weight has zero mass")
        w = w / mass
    return w


def weighted_energy(q: np.ndarray, u: np.ndarray, w: np.ndarray) -> float:
    """T_w(q) = ∫ w(u) q(u)^2 du (trapezoid)."""
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    w = np.asarray(w, dtype=np.float64).ravel()
    if q.shape != u.shape or w.shape != u.shape:
        raise ValueError("q, u, w must match shape")
    du = float(u[1] - u[0]) if u.size > 1 else 1.0
    # trapezoid with endpoint half-weights
    tw = np.ones_like(u)
    if u.size >= 2:
        tw[0] = tw[-1] = 0.5
    return float(np.sum(tw * w * q * q) * du)


def fisher_dprime(a: np.ndarray, b: np.ndarray) -> float:
    """Fisher discriminability d′ = |μa−μb| / sqrt((σa²+σb²)/2)."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    ma, mb = float(np.mean(a)), float(np.mean(b))
    va, vb = float(np.var(a, ddof=1)), float(np.var(b, ddof=1))
    den = np.sqrt(0.5 * (va + vb))
    if den <= 1e-30:
        return 0.0 if abs(ma - mb) < 1e-30 else float("inf")
    return abs(ma - mb) / den


def _make_ensemble(
    u: np.ndarray,
    *,
    kind: str,
    n: int,
    noise: float,
    rng: np.random.Generator,
    T: float = 20.0,
) -> np.ndarray:
    """
    Return (n, len(u)) residual matrix.

    kind:
      online  — high-frequency / CL-like (RH-like signature)
      offline — off-critical envelope or low-degree defect (non-RH-like)
    """
    n = int(n)
    out = np.zeros((n, u.size), dtype=np.float64)
    for i in range(n):
        if kind == "online":
            # mix CL mode + mild HF
            q = probe_critical_line_mode(u, T=T, t=14.134725 + 0.01 * i)
            q = q + 0.15 * probe_high_frequency(u, waves=40 + (i % 7))
        elif kind == "offline":
            # off-critical envelope OR defective low-degree mass
            if i % 2 == 0:
                q = probe_off_critical_mode(u, T=T, sigma=0.75 + 0.05 * (i % 3), t=14.13)
            else:
                q = probe_defective(u, waves=40, defect_degree=1, defect_weight=1.8)
        else:
            raise ValueError("kind must be online or offline")
        q = q + noise * rng.normal(size=u.size)
        # unit L2 normalize for fair Tw comparison
        nrm = np.sqrt(float(np.mean(q * q)) + 1e-30)
        out[i] = q / nrm
    return out


def sensitivity_experiment(
    *,
    n_per_class: int = 200,
    n_points: int = 1024,
    noise: float = 0.35,
    k: float = 4.0,
    sigma: float = 6.0,
    seed: int = 20260522,
    T: float = 20.0,
) -> Dict[str, Any]:
    """
    Run flat vs Gamma d′ on online/offline ensembles.

    Returns d′ values and relative gain (gamma - flat) / flat.
    """
    rng = np.random.default_rng(int(seed))
    u = sample_grid(int(n_points))
    w_flat = np.ones_like(u)
    w_gamma = gamma_weight(u, k=k, sigma=sigma, normalize=True)

    online = _make_ensemble(
        u, kind="online", n=n_per_class, noise=noise, rng=rng, T=T
    )
    offline = _make_ensemble(
        u, kind="offline", n=n_per_class, noise=noise, rng=rng, T=T
    )

    def scores(w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        s_on = np.array([weighted_energy(online[i], u, w) for i in range(n_per_class)])
        s_off = np.array([weighted_energy(offline[i], u, w) for i in range(n_per_class)])
        return s_on, s_off

    on_f, off_f = scores(w_flat)
    on_g, off_g = scores(w_gamma)
    d_flat = fisher_dprime(on_f, off_f)
    d_gamma = fisher_dprime(on_g, off_g)
    gain = (d_gamma - d_flat) / d_flat if d_flat > 1e-12 else float("nan")

    return {
        "banner": BANNER,
        "rh_claimed": False,
        "seed": int(seed),
        "k": float(k),
        "sigma": float(sigma),
        "peak_u_star": (k - 1.0) / sigma,
        "n_per_class": int(n_per_class),
        "noise": float(noise),
        "n_points": int(n_points),
        "flat_dprime": float(d_flat),
        "gamma_dprime": float(d_gamma),
        "relative_gain": float(gain),
        "relative_gain_percent": float(100.0 * gain) if gain == gain else float("nan"),
        "historical_idealized": HISTORICAL_IDEALIZED,
        "note": (
            "Gain is ensemble-dependent. Project record: ~53% on a noisy Beurling "
            "ensemble; idealized seed=20260522 often shows larger gain."
        ),
    }


def confirm_sensitivity_claim(
    *,
    min_gain: float = 0.53,
    noisy_noise: float = 0.55,
    clean_noise: float = 0.08,
    seed: int = 20260522,
) -> Dict[str, Any]:
    """
    Confirm the project-record sensitivity claim in two regimes:

    1) **noisy** Beurling-like ensemble — ask whether relative gain ≥ 53%.
    2) **clean/idealized** — report d′ and compare directionally to historical 2.37/5.53.

    Returns a structured verdict without hard-coding the experimental d′ values.
    """
    noisy = sensitivity_experiment(
        n_per_class=300,
        noise=noisy_noise,
        seed=seed,
        k=4.0,
        sigma=6.0,
    )
    clean = sensitivity_experiment(
        n_per_class=300,
        noise=clean_noise,
        seed=seed,
        k=4.0,
        sigma=6.0,
    )
    noisy_ok = bool(noisy["relative_gain"] >= min_gain - 1e-9)
    clean_improves = bool(clean["gamma_dprime"] > clean["flat_dprime"])
    # Idealized historical: gamma much larger than flat (more than double in record)
    clean_strong = bool(clean["gamma_dprime"] >= 1.5 * clean["flat_dprime"])

    return {
        "banner": BANNER,
        "rh_claimed": False,
        "claim": (
            f"Gamma weight improves offline/online discriminability by at least "
            f"{100*min_gain:.0f}% relative gain on a noisy ensemble; clean regime "
            f"improves further (historical idealized ~2.37→5.53)."
        ),
        "min_gain_required": float(min_gain),
        "noisy": noisy,
        "clean": clean,
        "noisy_meets_53pct": noisy_ok,
        "clean_improves": clean_improves,
        "clean_strong_gain": clean_strong,
        "confirmed": bool(noisy_ok and clean_improves),
        "verdict": (
            "CONFIRMED"
            if (noisy_ok and clean_improves)
            else (
                "PARTIAL"
                if clean_improves
                else "NOT_CONFIRMED"
            )
        ),
    }
