"""
Rank-5: B-RES as a threshold / obstruction problem (not full RH).

Formalizes the *weakest model hypothesis* that makes an off-critical mode force
nonvanishing R_d, and a model counterexample when that hypothesis fails
(cancellation / projection orthogonalization).

**Does not solve B-RES for arithmetic ζ.** Does not claim RH.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from .probes import (
    probe_critical_line_mode,
    probe_off_critical_mode,
    sample_grid,
)
from .projection import energy_ratio, project

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"
B_RES_ID = "B-RES"


def off_critical_rd_lower_model(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = 14.134725,
    degree: int = 4,
    n_points: int = 2048,
) -> Dict[str, float]:
    """R_d of pure off-critical model mode (envelope intact)."""
    u = sample_grid(int(n_points))
    q = probe_off_critical_mode(u, T=float(T), sigma=float(sigma), t=float(t))
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
    q = probe_off_critical_mode(u, T=float(T), sigma=float(sigma), t=float(t))
    # subtract P_d q
    pr = project(q, u, degree=int(degree), T=float(T))
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
    """
    Weakest *model* hypothesis H* sufficient for nonvanishing R_d from an
    off-critical contribution — the shape of what B-RES must prove for ζ.
    """
    return {
        "id": "H_star_injection",
        "statement": (
            "H*: After all EF main terms, secondary terms, and admissible weights, "
            "an off-critical zero ρ=σ+it contributes a residual component q_off whose "
            "correlation with V_d stays bounded below: "
            "liminf_T ||P_d q_off|| / ||q_arith|| ≥ ε(σ,t,d) > 0 "
            "(no total cancellation into V_d^⊥)."
        ),
        "implies": (
            "Under H*, liminf R_d(q_arith) ≥ ε^2 > 0 whenever such a zero exists, "
            "hence rapid R_d→0 forces no off-critical zeros (model Full B)."
        ),
        "counterexample_below": (
            "If H* fails (perfect V_d-orthogonalization of the off-critical piece), "
            "R_d may vanish — see cancelled_off_critical_rd."
        ),
        "status_for_zeta": (
            "B-RES is exactly the claim that H* holds for the true arithmetic residual "
            "of ζ. Open / RH-hard. Not proved here."
        ),
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
        "rank": 5,
        "title": "B-RES threshold / obstruction package",
        "B_RES_id": B_RES_ID,
        "threshold_hypothesis": threshold_hypothesis_statement(),
        "model_off_critical": off,
        "model_cancellation_counterexample": cancel,
        "pure_above_cancel": bool(pure_above_cancel),
        "conclusion": (
            "Full B reduces to B-RES = H* for arithmetic ζ. Model modes satisfy a "
            "positive R_d lower direction; forced V_d kill shows H* is necessary. "
            "Do not attack full RH; only prove H* or a weaker unconditional lemma."
        ),
    }
