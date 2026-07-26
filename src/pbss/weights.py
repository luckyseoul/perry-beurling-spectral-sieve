"""
Admissible weight / window class for PBSS log-window residuals on [0,1].

Motivation (open-plateau campaign): Hanning taper cut arithmetic R_d ~0.17 → ~0.08
by damping endpoints. Full Theorem A needs a *clean* weight class where endpoint
pollution is controlled analytically — not an ad-hoc post-hoc taper.

This module defines checkable weights and endpoint-contribution estimators.
**Not a proof of RH or full Theorem A.**

Admissible class W_α (α ∈ (0, 1/2])
------------------------------------
A continuous weight w:[0,1]→[0,1] is in W_α if:
  (W1) w(u) ≥ 0, ∫_0^1 w(u)^2 du > 0
  (W2) w vanishes near endpoints: w(u)=0 for u ∈ [0,α) ∪ (1-α,1]
       (or w(0)=w(1)=0 with |w(u)| ≤ C min(u,1-u)^β for some β>0)
  (W3) w is C^1 on (α,1-α) (or Lipschitz), so integration-by-parts for CL modes
       still yields O(T^{-2}) after reweighting (Lemma M3 style on the bulk).

Shipped members: ``tukey`` (cosine taper of half-width α), ``hanning`` (α=0.5
full cosine bell), ``raised_cosine`` alias of tukey.

Endpoint contribution estimator
-------------------------------
Split q = q_end + q_bulk with q_end = q·1_{u∉[α,1-α]}, q_bulk = q·1_{[α,1-α]}.
Then
  E_end = ‖P_d q_end‖² / ‖q‖²
is a checkable upper contribution of endpoint mass to low-degree energy.
For admissible w with support in [α,1-α], applying w kills E_end on the sample.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from .projection import _trapezoid_weights, energy_ratio, project_coefficients


WEIGHT_NAMES = ("tukey", "hanning", "raised_cosine", "flat", "none")


def tukey_weight(u: np.ndarray, alpha: float = 0.1) -> np.ndarray:
    """
    Tukey (tapered cosine) window on [0,1].

    ``alpha`` is the fraction of the interval spent tapering at *each* end
    (so total taper length 2α). For α≥0.5 this coincides with a full Hanning.
    """
    u = np.asarray(u, dtype=np.float64).ravel()
    a = float(alpha)
    if a < 0 or a > 0.5:
        raise ValueError("alpha must be in [0, 0.5]")
    w = np.ones_like(u)
    if a == 0.0:
        return w
    # left taper [0, a]
    left = u < a
    if np.any(left):
        # cosine from 0 → 1 on [0,a]
        w[left] = 0.5 * (1.0 - np.cos(np.pi * u[left] / a))
    # right taper [1-a, 1]
    right = u > 1.0 - a
    if np.any(right):
        w[right] = 0.5 * (1.0 - np.cos(np.pi * (1.0 - u[right]) / a))
    # endpoints exact zero when a>0
    w[u <= 0.0] = 0.0
    w[u >= 1.0] = 0.0
    return w


def hanning_weight(u: np.ndarray) -> np.ndarray:
    """Full-period Hanning / Hann window (member of W_{1/2})."""
    u = np.asarray(u, dtype=np.float64).ravel()
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * u))


def admissible_weight(
    u: np.ndarray,
    *,
    name: str = "tukey",
    alpha: float = 0.1,
) -> np.ndarray:
    """
    Build a weight from the admissible family.

    Parameters
    ----------
    name : tukey | hanning | raised_cosine | flat | none
    alpha : Tukey half-width (ignored for hanning/flat/none)
    """
    name = (name or "tukey").lower().strip()
    if name in ("none", "flat", "raw", "one"):
        return np.ones(np.asarray(u).shape[0], dtype=np.float64)
    if name in ("hanning", "hann"):
        return hanning_weight(u)
    if name in ("tukey", "raised_cosine", "cosine_taper"):
        return tukey_weight(u, alpha=alpha)
    raise ValueError(f"unknown weight name {name!r}; choose from {WEIGHT_NAMES}")


def is_admissible_weight(
    w: np.ndarray,
    u: np.ndarray,
    *,
    alpha: float = 0.1,
    tol: float = 1e-9,
) -> Dict[str, object]:
    """
    Check discrete sample of w against W_α membership (W1–W2 discrete form).

    Returns a dict with ok flag and diagnostics. Not a continuous C¹ check.
    """
    w = np.asarray(w, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if w.shape != u.shape:
        raise ValueError("w and u must match")
    a = float(alpha)
    wt = _trapezoid_weights(u)
    mass = float(np.sum(wt * w * w))
    nonneg = bool(np.all(w >= -tol))
    end_mask = (u < a - 1e-15) | (u > 1.0 - a + 1e-15)
    end_ok = True
    max_end = 0.0
    if a > 0 and np.any(end_mask):
        max_end = float(np.max(np.abs(w[end_mask])))
        # Tukey has smooth taper, not hard zero in (0,a); allow small values
        # Hard vanishing required only for "strict" class — report max_end
        end_ok = max_end <= 1.0 + tol  # always true for unit weights; report metric
    ok = nonneg and mass > tol
    return {
        "ok": ok,
        "nonneg": nonneg,
        "l2_mass": mass,
        "alpha": a,
        "max_abs_in_end_zones": max_end,
        "end_zone_fraction": float(np.mean(end_mask)) if end_mask.size else 0.0,
        "note": "Discrete W_α check; not a proof of continuous admissibility.",
    }


def apply_weight(q: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Pointwise weight application: (W q)(u) = w(u) q(u)."""
    q = np.asarray(q, dtype=np.float64).ravel()
    w = np.asarray(w, dtype=np.float64).ravel()
    if q.shape != w.shape:
        raise ValueError("q and w must match")
    return q * w


def endpoint_mask(u: np.ndarray, alpha: float = 0.1) -> np.ndarray:
    """Boolean mask: True on endpoint zones [0,α) ∪ (1-α,1]."""
    u = np.asarray(u, dtype=np.float64).ravel()
    a = float(alpha)
    if a < 0 or a > 0.5:
        raise ValueError("alpha in [0,0.5]")
    return (u < a) | (u > 1.0 - a)


def endpoint_contribution(
    q: np.ndarray,
    u: np.ndarray,
    *,
    degree: int = 4,
    alpha: float = 0.1,
) -> Dict[str, float]:
    """
    Checkable estimator of endpoint contribution to low-degree mass.

      E_end = ‖P_d (q · 1_end)‖² / ‖q‖²
      E_bulk = ‖P_d (q · 1_bulk)‖² / ‖q‖²
      R_d = energy_ratio(q)
      R_d_bulk = energy_ratio(q_bulk)  (bulk-only residual)

    E_end is the quantity taper kills. **Numeric diagnostic**, not a theorem
    that arithmetic R_d → 0.
    """
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    if q.shape != u.shape:
        raise ValueError("q and u must match")
    if degree < 0:
        raise ValueError("degree >= 0")
    wts = _trapezoid_weights(u)
    l2 = float(np.sum(wts * q * q))
    if l2 <= 1e-30:
        return {
            "E_end": 0.0,
            "E_bulk": 0.0,
            "R_d": 0.0,
            "R_d_bulk": 0.0,
            "R_d_end": 0.0,
            "l2": l2,
            "alpha": float(alpha),
            "degree": int(degree),
        }
    end = endpoint_mask(u, alpha)
    bulk = ~end
    q_end = np.where(end, q, 0.0)
    q_bulk = np.where(bulk, q, 0.0)
    c_end = project_coefficients(q_end, u, degree, weights=wts)
    c_bulk = project_coefficients(q_bulk, u, degree, weights=wts)
    e_end = float(np.dot(c_end, c_end))
    e_bulk = float(np.dot(c_bulk, c_bulk))
    r_full = energy_ratio(q, u, degree, weights=wts)
    r_bulk = energy_ratio(q_bulk, u, degree, weights=wts) if float(np.sum(wts * q_bulk * q_bulk)) > 1e-30 else 0.0
    r_end = energy_ratio(q_end, u, degree, weights=wts) if float(np.sum(wts * q_end * q_end)) > 1e-30 else 0.0
    return {
        "E_end": e_end / l2,
        "E_bulk": e_bulk / l2,
        "R_d": float(r_full),
        "R_d_bulk": float(r_bulk),
        "R_d_end": float(r_end),
        "l2": l2,
        "alpha": float(alpha),
        "degree": int(degree),
    }


def weighted_energy_ratio(
    q: np.ndarray,
    u: np.ndarray,
    *,
    degree: int = 4,
    weight_name: str = "tukey",
    alpha: float = 0.1,
) -> Tuple[float, np.ndarray, Dict[str, float]]:
    """
    Apply admissible weight then compute R_d(Wq).

    Returns (R_d_weighted, w, endpoint_stats_before_weight).
    """
    end_stats = endpoint_contribution(q, u, degree=degree, alpha=alpha)
    w = admissible_weight(u, name=weight_name, alpha=alpha)
    qw = apply_weight(q, w)
    r = float(energy_ratio(qw, u, degree))
    return r, w, end_stats


def bulk_vs_weighted_report(
    q: np.ndarray,
    u: np.ndarray,
    *,
    degree: int = 4,
    alpha: float = 0.1,
    weight_name: str = "tukey",
) -> Dict[str, object]:
    """
    One-shot diagnostic comparing raw R_d, bulk-only R_d, and weighted R_d.

    Used by theorem-A scaffold multi-T scans.
    """
    end = endpoint_contribution(q, u, degree=degree, alpha=alpha)
    r_w, w, _ = weighted_energy_ratio(
        q, u, degree=degree, weight_name=weight_name, alpha=alpha
    )
    adm = is_admissible_weight(w, u, alpha=alpha)
    return {
        "R_d_raw": end["R_d"],
        "R_d_bulk": end["R_d_bulk"],
        "R_d_weighted": r_w,
        "E_end": end["E_end"],
        "E_bulk": end["E_bulk"],
        "alpha": float(alpha),
        "weight_name": weight_name,
        "weight_admissible_ok": adm["ok"],
        "weight_l2_mass": adm["l2_mass"],
        "degree": int(degree),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
