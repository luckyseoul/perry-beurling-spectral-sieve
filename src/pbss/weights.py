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
  (W3) w is fixed and absolutely continuous on [0,1], with finite ‖w'‖₁,
       so integration-by-parts for CL modes yields O(T^{-2}) after reweighting
       when the weighted norm has a positive lower bound (Lemma M6).

Shipped members: ``tukey`` (cosine taper of half-width α), ``hanning`` (α=0.5
full cosine bell), ``raised_cosine`` alias of tukey.

Endpoint contribution estimator
-------------------------------
Split q = q_end + q_bulk with q_end = q·1_{u∉[α,1-α]}, q_bulk = q·1_{[α,1-α]}.
Then
  E_end = ‖P_d q_end‖² / ‖q‖²
is the endpoint-only projected energy. Projection mixes the disjoint supports:
R_d(q) = E_end + E_bulk + 2 Re⟨P_d q_end,P_d q_bulk⟩ / ‖q‖².
The cross term can have either sign, so E_end is not an additive attribution.
For admissible w with support in [α,1-α], applying w kills E_end on the sample.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from .basis import orthonormal_legendre_design
from .projection import _finite_energy, _trapezoid_weights, _validate_grid, energy_ratio, project_coefficients


WEIGHT_NAMES = ("tukey", "hanning", "raised_cosine", "flat", "none")


def _alpha(alpha: float) -> float:
    a = float(alpha)
    if not np.isfinite(a) or not 0.0 <= a <= 0.5:
        raise ValueError("alpha must be finite and in [0, 0.5]")
    return a


def tukey_weight(u: np.ndarray, alpha: float = 0.1) -> np.ndarray:
    """
    Tukey (tapered cosine) window on [0,1].

    ``alpha`` is the fraction of the interval spent tapering at *each* end
    (so total taper length 2α). For α≥0.5 this coincides with a full Hanning.
    """
    u = _validate_grid(u)
    a = _alpha(alpha)
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
    u = _validate_grid(u)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * u))


def admissible_weight(
    u: np.ndarray,
    *,
    name: str = "tukey",
    alpha: float = 0.1,
) -> np.ndarray:
    """
    Build an endpoint taper, or the flat/none comparison baseline.

    Parameters
    ----------
    name : tukey | hanning | raised_cosine | flat | none
    alpha : Tukey half-width (ignored for hanning/flat/none)
    """
    name = (name or "tukey").lower().strip()
    if name in ("none", "flat", "raw", "one"):
        return np.ones(_validate_grid(u).size, dtype=np.float64)
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
    Check necessary sampled conditions for endpoint-vanishing W_α weights.

    ``ok`` requires values in [0,1], positive mass, and sampled endpoints
    at 0 and 1 with zero weight. ``strict_support_ok`` additionally checks
    the hard support condition; Hann/Tukey generally satisfy only the taper
    alternative. Neither flag proves regularity between sample points.
    """
    w = np.asarray(w, dtype=np.float64)
    u = _validate_grid(u)
    if w.ndim != 1 or w.shape != u.shape:
        raise ValueError("w and u must be matching 1D arrays")
    a = _alpha(alpha)
    tol = float(tol)
    if not np.isfinite(tol) or tol < 0.0:
        raise ValueError("tol must be finite and nonnegative")
    wt = _trapezoid_weights(u)
    finite = bool(np.all(np.isfinite(w)))
    mass = float(np.sum(wt * w * w)) if finite else 0.0
    nonneg = bool(np.all(w >= -tol))
    bounded = bool(np.all(w <= 1.0 + tol))
    end_mask = (u < a - 1e-15) | (u > 1.0 - a + 1e-15)
    max_end = float(np.max(np.abs(w[end_mask]))) if finite and np.any(end_mask) else None
    endpoints_sampled = bool(abs(u[0]) <= tol and abs(u[-1] - 1.0) <= tol)
    endpoints_zero = bool(finite and endpoints_sampled and abs(w[0]) <= tol and abs(w[-1]) <= tol)
    ok = finite and nonneg and bounded and mass > 0.0 and endpoints_zero
    strict_support_ok = bool(ok and (max_end is None or max_end <= tol))
    return {
        "ok": ok,
        "finite": finite,
        "nonneg": nonneg,
        "bounded_by_one": bounded,
        "endpoints_sampled": endpoints_sampled,
        "endpoints_zero": endpoints_zero,
        "strict_support_ok": strict_support_ok,
        "continuous_regularity_verified": False,
        "l2_mass": mass,
        "alpha": a,
        "max_abs_in_end_zones": max_end,
        "end_zone_fraction": float(np.mean(end_mask)) if end_mask.size else 0.0,
        "note": "Necessary sampled conditions only; not a proof of continuous admissibility or M6 regularity.",
    }


def apply_weight(q: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Pointwise weight application: (W q)(u) = w(u) q(u)."""
    q = np.asarray(q)
    w = np.asarray(w, dtype=np.float64)
    if q.ndim != 1 or w.ndim != 1 or q.shape != w.shape:
        raise ValueError("q and w must be matching 1D arrays")
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(w)):
        raise ValueError("q and w must be finite")
    return q * w


def endpoint_mask(u: np.ndarray, alpha: float = 0.1) -> np.ndarray:
    """Boolean mask: True on endpoint zones [0,α) ∪ (1-α,1]."""
    u = _validate_grid(u)
    a = _alpha(alpha)
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
      E_cross = 2 Re⟨P_d q_end,P_d q_bulk⟩ / ‖q‖²
      R_d = E_end + E_bulk + E_cross
      R_d_bulk = energy_ratio(q_bulk)  (bulk-only residual)

    E_cross can be negative. A hard support mask removes q_end; a smooth
    taper attenuates it and can also change the bulk. These energies alone
    do not predict the change in the normalized ratio under tapering.
    A zero residual has no energy ratio. If the requested absolute l2
    cannot be represented in float64, rescale q before requesting this report.
    """
    q = np.asarray(q)
    u = _validate_grid(u)
    if q.ndim != 1 or q.shape != u.shape or not np.all(np.isfinite(q)):
        raise ValueError("q must be a finite 1D array matching u")
    end = endpoint_mask(u, alpha)
    wts = _trapezoid_weights(u)
    # Validate the projection space even for a zero residual.
    project_coefficients(q, u, degree, weights=wts)
    scale = float(np.max(np.abs(q)))
    if scale == 0.0:
        raise ValueError("||q||^2 is zero; cannot form energy ratio")
    # Ratios use a scaled residual to preserve amplitude invariance.
    q_scaled = q / scale
    l2_scaled = float(np.sum(wts * np.abs(q_scaled) ** 2))
    l2 = _finite_energy(np.longdouble(scale)**2 * np.longdouble(l2_scaled), "L2 norm squared")
    bulk = ~end
    q_end = np.where(end, q_scaled, 0.0)
    q_bulk = np.where(bulk, q_scaled, 0.0)
    c_end = project_coefficients(q_end, u, degree, weights=wts)
    c_bulk = project_coefficients(q_bulk, u, degree, weights=wts)
    Phi = orthonormal_legendre_design(degree, u)
    p_end, p_bulk = Phi @ c_end, Phi @ c_bulk
    e_end = float(np.sum(wts * np.abs(p_end) ** 2)) / l2_scaled
    e_bulk = float(np.sum(wts * np.abs(p_bulk) ** 2)) / l2_scaled
    e_cross = float(2.0 * np.real(np.sum(wts * np.conj(p_end) * p_bulk))) / l2_scaled
    r_full = energy_ratio(q_scaled, u, degree, weights=wts)
    r_bulk = energy_ratio(q_bulk, u, degree, weights=wts) if np.any(q_bulk) else 0.0
    r_end = energy_ratio(q_end, u, degree, weights=wts) if np.any(q_end) else 0.0
    return {
        "E_end": e_end,
        "E_bulk": e_bulk,
        "E_cross": e_cross,
        "projection_triangle_upper": (np.sqrt(e_end) + np.sqrt(e_bulk)) ** 2,
        "decomposition_error": abs(float(r_full) - e_end - e_bulk - e_cross),
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
        "E_cross": end["E_cross"],
        "projection_triangle_upper": end["projection_triangle_upper"],
        "decomposition_error": end["decomposition_error"],
        "alpha": float(alpha),
        "weight_name": weight_name,
        "weight_admissible_ok": adm["ok"],
        "weight_l2_mass": adm["l2_mass"],
        "degree": int(degree),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
