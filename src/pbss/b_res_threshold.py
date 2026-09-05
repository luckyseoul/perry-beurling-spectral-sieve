"""
Rank-5: B-RES as a threshold / obstruction problem (not full RH).

Records a sufficient total-projection hypothesis and finite sampled controls.
The old component-only hypothesis omitted projected remainder cancellation.
Even an intact continuous off-critical mode has R_d = O(1/T), with deeper
phase-dependent troughs. For the exact continuous arithmetic residual, the
positive-floor version of B-RES is equivalent to RH; it is not established here.

**Does not solve B-RES for arithmetic ζ.** Does not claim RH.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from .probes import (
    probe_critical_line_mode,
    sample_grid,
)
from .projection import energy_ratio, project

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"
B_RES_ID = "B-RES"


def _scaled_mode(u: np.ndarray, *, T: float, sigma: float, t: float) -> np.ndarray:
    """A common exponential scale leaves the ratio unchanged and avoids overflow."""
    T, sigma, t = float(T), float(sigma), float(t)
    if not all(np.isfinite(x) for x in (T, sigma, t)) or T <= 0 or t <= 0:
        raise ValueError("T and t must be positive; T, sigma and t must be finite")
    exponent, frequency = T * (sigma - 0.5), T * t
    if not np.isfinite(exponent) or not np.isfinite(frequency):
        raise ValueError("model exponent and frequency must be finite")
    return np.exp(exponent * u - max(exponent, 0.0)) * np.sin(frequency * u)


def off_critical_rd_lower_model(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = 14.134725,
    degree: int = 4,
    n_points: int = 2048,
) -> Dict[str, float]:
    """Finite sampled ratio; the legacy name does not assert a lower bound."""
    u = sample_grid(int(n_points))
    q = _scaled_mode(u, T=T, sigma=sigma, t=t)
    r = float(energy_ratio(q, u, int(degree)))
    r_cl = float(
        energy_ratio(probe_critical_line_mode(u, T=float(T), t=float(t)), u, int(degree))
    )
    return {
        "T": float(T),
        "sigma": float(sigma),
        "R_d_off": r,
        "R_d_cl": r_cl,
        "ratio": r / max(r_cl, 1e-30),
    }


def cancelled_off_critical_rd(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = 14.134725,
    degree: int = 4,
    n_points: int = 2048,
) -> Dict[str, float]:
    """
    Model counterexample *below* the threshold hypothesis:

    Take off-critical mode and subtract its own projection onto V_d (remove all
    low-degree mass by force). Then R_d = 0 by construction even though σ≠1/2
    was present before cancellation.

    Shows: without a non-cancellation / injection hypothesis, off-critical origin
    alone does not force liminf R_d > 0.
    """
    u = sample_grid(int(n_points))
    q = _scaled_mode(u, T=T, sigma=sigma, t=t)
    # subtract P_d q
    pr = project(q, u, degree=int(degree), T=1.0)
    # reconstruct low part from coefficients
    from .basis import shifted_legendre_values

    low = np.zeros_like(u)
    for k, c in enumerate(pr.coeffs):
        low = low + float(c) * shifted_legendre_values(k, u)
    q_orth = q - low
    r = float(energy_ratio(q_orth, u, int(degree)))
    return {
        "T": float(T),
        "sigma": float(sigma),
        "R_d_after_Vd_kill": r,
        "threshold_hypothesis_holds": False,
        "note": "Forced V_d removal ⇒ R_d≈0; model of cancellation below B-RES threshold",
    }


def threshold_hypothesis_statement() -> Dict[str, str]:
    """Sufficient projection margin, including the missing remainder term."""
    return {
        "id": "H_star_injection",
        "statement": (
            "For q_arith=q_off+q_rest, a sufficient H* is "
            "liminf_T (||P_d q_off||-||P_d q_rest||)/||q_arith|| ≥ ε > 0. "
            "A lower bound on ||P_d q_off|| alone is insufficient."
        ),
        "implies": (
            "The reverse triangle inequality gives liminf R_d(q_arith) ≥ ε^2. "
            "The strengthened component margin is sufficient, not necessary."
        ),
        "counterexample_below": (
            "An intact continuous off-mode already has R_d=O(1/T). Also "
            "q_rest=-P_d q_off+h_perp cancels its projected contribution exactly."
        ),
        "status_for_zeta": (
            "For exact continuous degree-one detrending and fixed d>=2, "
            "liminf R_d=0 unconditionally. The assertion off-critical zero => "
            "liminf R_d>0 is therefore RH-equivalent and remains open."
        ),
        "scope": "Exact continuous theta residual, du projection; no sampled transfer asserted.",
    }


def b_res_threshold_report(
    *,
    T: float = 20.0,
    sigma: float = 0.9,
    degree: int = 4,
) -> Dict[str, Any]:
    """Rank-5 entry: threshold package + model obstruction + cancellation counterexample."""
    off = off_critical_rd_lower_model(T, sigma=sigma, degree=degree)
    cancel = cancelled_off_critical_rd(T, sigma=sigma, degree=degree)
    # sanity: pure off has larger R_d than cancelled
    pure_above_cancel = off["R_d_off"] > cancel["R_d_after_Vd_kill"] + 1e-6
    return {
        "banner": BANNER,
        "rh_claimed": False,
        "b_res_solved": False,
        "normalized_b_res_floor": "rh_equivalent_open",
        "normalized_recurrence": "proved_unconditional_continuous",
        "model_asymptotic": "vanishing_with_phase_dependent_troughs",
        "rank": 5,
        "title": "B-RES threshold / obstruction package",
        "B_RES_id": B_RES_ID,
        "threshold_hypothesis": threshold_hypothesis_statement(),
        "model_off_critical": off,
        "model_cancellation_counterexample": cancel,
        "pure_above_cancel": bool(pure_above_cancel),
        "conclusion": (
            "Finite sampled controls do not prove a projection floor. Exact "
            "continuous off-modes vanish; the arithmetic liminf is also zero "
            "unconditionally. B-ABS proves an absolute-moment RH equivalence; "
            "its unconditional growth bound and the normalized converse remain open."
        ),
    }
