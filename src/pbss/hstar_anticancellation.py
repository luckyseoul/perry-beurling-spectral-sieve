"""
Historical L* proxy, corrected after the real-sine phase obstruction.

The proposed positive coercivity interpretation is retracted. For every fixed
A,t>0 and degree d, the intact continuous real off-mode already has R_d→0.
The old phase-independent expression is a limsup envelope, not a lower bound.

This module retains finite-grid projection, cancellation, and injected-bump
checks as numerical diagnostics. It evaluates no Guinand-Weil quadratic form
and proves neither arithmetic B-RES nor RH. Rank-4 proportion remains STOP.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .b_res_threshold import _scaled_mode, cancelled_off_critical_rd, off_critical_rd_lower_model
from .basis import shifted_legendre_values
from .probes import (
    arithmetic_residual,
    arithmetic_zero_peel,
    probe_critical_line_mode,
    sample_grid,
)
from .projection import _trapezoid_weights, energy_ratio, project_coefficients
from .off_critical_exact import off_critical_projection, partial_cancellation_ratio

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"
LEMMA_ID = "L_star_local_anticancellation"
CANDIDATE_STATUS = "retracted_coercivity_proxy"

# Default campaign grid (overridden by callers / CLI)
DEFAULT_T = (8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0)
DEFAULT_SIGMA = (0.60, 0.70, 0.80, 0.90)
DEFAULT_DEGREE = (2, 4, 6)

DEFAULT_T_FIRST_ZERO = 14.134725


def _scaled_off_mode(u, *, T, sigma, t):
    """An amplitude-equivalent off-mode without forming exp(T*(sigma-1/2))."""
    sigma = float(sigma)
    if not np.isfinite(sigma) or sigma <= 0.5:
        raise ValueError("sigma must be finite and greater than 1/2")
    return _scaled_mode(np.asarray(u), T=T, sigma=sigma, t=t)


def analytic_off_critical_Rd_leading(
    T: float,
    *,
    sigma: float,
    degree: int,
    t: float = DEFAULT_T_FIRST_ZERO,
) -> Dict[str, Any]:
    """Phase-aware continuous asymptotic and exact-identity evaluation.

    R_d_leading now includes the real endpoint phase. The old expression is
    retained as R_d_limsup_envelope. eps_star_analytic is a legacy half-envelope
    reference only: eps_star_analytic_is_lower_bound is explicitly false.
    See OFF_CRITICAL_PHASE_OBSTRUCTION.md for the analytic model proof.
    """
    exact = off_critical_projection(T, A=float(sigma) - 0.5, t=t, degree=degree)
    leading = exact.phase_leading_ratio
    envelope = exact.phase_envelope_ratio
    return {
        "T": exact.T,
        "sigma": float(sigma),
        "degree": float(exact.degree),
        "t": exact.t,
        "a": exact.T * exact.A,
        "b": exact.T * exact.t,
        "R_d_leading": leading,
        "R_d_phase_leading": leading,
        "R_d_limsup_envelope": envelope,
        "R_d_continuous_exact": exact.energy_ratio,
        "C_d_leading": float(np.sqrt(leading)),
        "eps_star_analytic": 0.5 * float(np.sqrt(envelope)),
        "eps_star_analytic_is_lower_bound": False,
        "coercivity_status": CANDIDATE_STATUS,
        "continuous_Rd_limit": 0.0,
        "analytic_evaluation": exact.evaluation,
    }


def analytic_vs_numeric_row(
    T: float,
    *,
    sigma: float = 0.9,
    degree: int = 4,
    t: float = DEFAULT_T_FIRST_ZERO,
    n_points: int = 4096,
) -> Dict[str, Any]:
    """Compare sampled QR projection with the exact continuous model identity.

    Half-envelope flags describe this finite sample only; they certify no
    coercivity inequality. The leading asymptotic has additive O(T^-2) error.
    """
    lead = analytic_off_critical_Rd_leading(T, sigma=sigma, degree=degree, t=t)
    u = sample_grid(int(n_points))
    q = _scaled_off_mode(u, T=T, sigma=sigma, t=t)
    r = float(energy_ratio(q, u, degree))
    c = float(np.sqrt(r))
    exact = lead["R_d_continuous_exact"]
    relative_discretization_error = abs(r - exact) / max(exact, 1e-30)
    above_reference = float(c >= lead["eps_star_analytic"])
    return {
        **lead,
        "R_d_numeric": r,
        "C_d_numeric": c,
        "rel_err_Rd": relative_discretization_error,
        "relative_discretization_error": relative_discretization_error,
        "relative_phase_leading_error": abs(lead["R_d_phase_leading"] - exact) / max(exact, 1e-30),
        "numeric_above_half_leading": above_reference,
        "numeric_above_legacy_half_envelope": above_reference,
        "half_envelope_role": "finite_sample_reference_not_a_lower_bound",
    }


def arithmetic_injection_bridge_row(
    x_max: float,
    *,
    amp: float = 0.1,
    sigma: float = 0.9,
    degree: int = 4,
    n_strip: int = 30,
    n_points: int = 2048,
    t: float = DEFAULT_T_FIRST_ZERO,
) -> Dict[str, float]:
    """
    Calibration on an arithmetic background:
      1) build a peeled arithmetic residual and remove V_d;
      2) add an explicitly chosen off-mode with norm amp*||q_orth||;
      3) verify C_d(mix)=C_d(off)*amp*||q_orth||/||mix||.

    Every norm uses the projector's trapezoid inner product. The equality is
    linearity plus P_d q_orth=0. No arithmetic zero forces the injected bump;
    this is not evidence of arithmetic non-cancellation or a Weil theorem.
    """
    u = sample_grid(int(n_points))
    x_max = float(x_max)
    if n_strip <= 0:
        q, T, meta = arithmetic_residual(u, x_max=x_max, detrend="deg1")
    else:
        q, T, meta = arithmetic_zero_peel(
            u, x_max=x_max, n_strip=int(n_strip), fit_scale=True, detrend="deg1"
        )
    q_orth = partial_cancel(q, u, degree=int(degree), T=float(T), alpha=1.0)
    c_orth = coercivity_Cd(q_orth, u, degree=int(degree), T=float(T))
    if not np.isfinite(amp) or amp <= 0:
        raise ValueError("amp must be finite and positive")
    weights = _trapezoid_weights(u)
    norm = lambda v: float(np.sqrt(np.dot(weights, v * v)))
    n_orth = norm(q_orth)
    qoff = _scaled_off_mode(u, T=T, sigma=sigma, t=t)
    qoff = qoff / max(norm(qoff), 1e-30) * (float(amp) * n_orth)
    c_pure = coercivity_Cd(qoff, u, degree=int(degree), T=float(T))
    qmix = q_orth + qoff
    c_mix = coercivity_Cd(qmix, u, degree=int(degree), T=float(T))
    pred = c_pure * (float(amp) * n_orth) / max(norm(qmix), 1e-30)
    lead = analytic_off_critical_Rd_leading(float(T), sigma=float(sigma), degree=int(degree), t=t)
    return {
        "x_max": x_max,
        "T": float(T),
        "n_strip": float(n_strip),
        "amp": float(amp),
        "sigma": float(sigma),
        "degree": float(degree),
        "n_primes": float(meta.get("n_primes", 0)),
        "C_d_arith": float(coercivity_Cd(q, u, degree=int(degree), T=float(T))),
        "C_d_orth": float(c_orth),
        "C_d_pure_bump": float(c_pure),
        "C_d_mix": float(c_mix),
        "C_d_pred": float(pred),
        "pred_ratio": float(c_mix / max(pred, 1e-30)),
        "eps_star_analytic": float(lead["eps_star_analytic"]),
        "eps_star_analytic_is_lower_bound": False,
        "injected_model_component": True,
        "arithmetic_off_zero_forced": False,
        "bridge_role": "projection_identity_calibration_only",
        "bridge_ok": float(1.0 if (c_orth < 1e-4 and 0.85 <= c_mix / max(pred, 1e-30) <= 1.15) else 0.0),
    }


def _bridge_job(args: Tuple[float, float, float, int, int]) -> Dict[str, float]:
    x_max, amp, sigma, degree, n_strip = args
    return arithmetic_injection_bridge_row(
        x_max, amp=amp, sigma=sigma, degree=degree, n_strip=n_strip
    )


def arithmetic_bridge_campaign(
    *,
    x_max_values: Sequence[float] = (1e5, 1e6, 3e6, 1e7),
    amps: Sequence[float] = (0.05, 0.1, 0.2),
    sigma_values: Sequence[float] = (0.7, 0.9),
    degrees: Sequence[int] = (4,),
    n_strip: int = 30,
    max_workers: Optional[int] = None,
) -> List[Dict[str, float]]:
    jobs = [
        (float(x), float(a), float(s), int(d), int(n_strip))
        for x in x_max_values
        for a in amps
        for s in sigma_values
        for d in degrees
    ]
    if max_workers is None:
        max_workers = min(len(jobs), max(1, (os.cpu_count() or 4) - 2))
    if len(jobs) <= 2 or max_workers <= 1:
        return [_bridge_job(j) for j in jobs]
    rows: List[Dict[str, float]] = []
    with ProcessPoolExecutor(max_workers=int(max_workers)) as ex:
        futs = [ex.submit(_bridge_job, j) for j in jobs]
        for fut in as_completed(futs):
            rows.append(fut.result())
    rows.sort(key=lambda r: (r["x_max"], r["sigma"], r["amp"]))
    return rows


def _analytic_job(args: Tuple[float, float, int, int]) -> Dict[str, float]:
    T, s, d, n_points = args
    return analytic_vs_numeric_row(T, sigma=s, degree=d, n_points=n_points)


def analytic_envelope_campaign(
    *,
    T_values: Sequence[float] = (16.0, 24.0, 32.0, 48.0, 64.0, 96.0, 128.0),
    sigma_values: Sequence[float] = DEFAULT_SIGMA,
    degrees: Sequence[int] = DEFAULT_DEGREE,
    n_points: int = 4096,
    max_workers: Optional[int] = None,
) -> List[Dict[str, float]]:
    jobs = [
        (float(T), float(s), int(d), int(n_points))
        for T in T_values
        for s in sigma_values
        for d in degrees
    ]
    if max_workers is None:
        max_workers = min(len(jobs), max(1, (os.cpu_count() or 4) - 2))
    if max_workers <= 1 or len(jobs) <= 2:
        return [_analytic_job(j) for j in jobs]
    rows = []
    with ProcessPoolExecutor(max_workers=int(max_workers)) as ex:
        futs = [ex.submit(_analytic_job, j) for j in jobs]
        for fut in as_completed(futs):
            rows.append(fut.result())
    rows.sort(key=lambda r: (r["sigma"], r["degree"], r["T"]))
    return rows



def lemma_statement() -> Dict[str, str]:
    """Historical L* label with its unsupported coercivity interpretation retracted."""
    return {
        "id": LEMMA_ID,
        "name": "L*",
        "doc": "docs/HSTAR_ANTICANCELLATION.md",
        "status": CANDIDATE_STATUS,
        "hypotheses": (
            "Historical finite-polynomial proxy C_d=||P_d q||/||q||. "
            "No Guinand-Weil quadratic form is evaluated or identified with this norm."
        ),
        "conclusion": (
            "For fixed A,t>0 and d, the intact real off-mode has R_d=O(1/T)→0. "
            "The old phase-independent expression is the limsup envelope for T R_d; "
            "phase troughs are order T^-3 for d>=1. Positive limiting pure-mode "
            "coercivity is retracted. Orthogonal cancellation identities remain valid."
        ),
        "not_a_proportion": "No zero-proportion bound; Rank-4 proportion remains STOP.",
        "gaps_to_full_H_star": (
            "Identify and control the total arithmetic residual and denominator; "
            "prove non-cancellation for actual arithmetic contributions. An injected "
            "bump calibrates a projection identity and does not establish either premise."
        ),
    }


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-30 or nb < 1e-30:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def off_critical_signature(
    q: np.ndarray,
    u: np.ndarray,
    *,
    T: float,
    sigma: float,
    t: float,
) -> float:
    """Shape correlation with off-mode minus |corr| with CL mode (same T,t)."""
    q = np.asarray(q, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    q_off = _scaled_off_mode(u, T=T, sigma=sigma, t=t)
    q_cl = probe_critical_line_mode(u, T=float(T), t=float(t))
    return _corr(q, q_off) - abs(_corr(q, q_cl))


def coercivity_Cd(q: np.ndarray, u: np.ndarray, *, degree: int, T: float) -> float:
    """W = C_d = ||P_d q|| / ||q|| = sqrt(R_d) on the shipped projector."""
    r = float(energy_ratio(q, u, int(degree)))
    return float(np.sqrt(max(r, 0.0)))


def _vd_component(q: np.ndarray, u: np.ndarray, *, degree: int, T: float) -> np.ndarray:
    coefficients = project_coefficients(q, u, degree)
    low = np.zeros_like(u)
    for k, c in enumerate(coefficients):
        low = low + float(c) * shifted_legendre_values(k, u)
    return low


def partial_cancel(
    q: np.ndarray,
    u: np.ndarray,
    *,
    degree: int,
    T: float,
    alpha: float,
) -> np.ndarray:
    """Remove fraction α∈[0,1] of the V_d component: q ↦ q - α P_d q."""
    alpha = float(alpha)
    if not np.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be finite and in [0,1]")
    low = _vd_component(q, u, degree=degree, T=T)
    return np.asarray(q) - alpha * low


def weil_local_proxy_row(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = DEFAULT_T_FIRST_ZERO,
    degree: int = 4,
    n_points: int = 2048,
    eps_C: Optional[float] = None,
) -> Dict[str, float]:
    """Historical name: finite sampled separation, with no Weil-form assertion."""
    u = sample_grid(int(n_points))
    q_pure = _scaled_off_mode(u, T=T, sigma=sigma, t=t)
    q_cancel = partial_cancel(q_pure, u, degree=degree, T=T, alpha=1.0)
    q_cl = probe_critical_line_mode(u, T=float(T), t=float(t))

    s_pure = off_critical_signature(q_pure, u, T=T, sigma=sigma, t=t)
    s_cancel = off_critical_signature(q_cancel, u, T=T, sigma=sigma, t=t)
    c_pure = coercivity_Cd(q_pure, u, degree=degree, T=T)
    c_cancel = coercivity_Cd(q_cancel, u, degree=degree, T=T)
    c_cl = coercivity_Cd(q_cl, u, degree=degree, T=T)
    r_pure = float(c_pure**2)
    r_cancel = float(c_cancel**2)
    gap = r_pure - r_cancel

    # Finite-row comparisons only; none is an asymptotic coercivity statement.
    beats_cl = c_pure > c_cl * 1.05 if c_cl > 1e-12 else c_pure > 1e-4
    cancel_killed = c_cancel < max(1e-4, 0.02 * c_pure)
    shape_off = s_pure > 0.02
    if eps_C is None:
        # Per-row soft floor: half the observed pure coercivity (always holds for pure;
        # cancelled must fall below it).
        eps_use = 0.5 * c_pure
    else:
        eps_use = float(eps_C)
        if not np.isfinite(eps_use) or eps_use < 0:
            raise ValueError("eps_C must be finite and nonnegative")
    above_reference = c_pure >= eps_use * 0.99
    below_reference_after_cancel = c_cancel < eps_use
    sample_separation = bool(
        shape_off and beats_cl and cancel_killed and above_reference and below_reference_after_cancel and gap > 1e-12
    )
    return {
        "T": float(T),
        "sigma": float(sigma),
        "degree": float(degree),
        "t": float(t),
        "S_off_pure": float(s_pure),
        "S_off_cancelled": float(s_cancel),
        "C_d_pure": float(c_pure),
        "C_d_cancelled": float(c_cancel),
        "C_d_cl": float(c_cl),
        "R_d_pure": float(r_pure),
        "R_d_cancelled": float(r_cancel),
        "R_d_ratio_off_over_cl": float(r_pure / max(c_cl**2, 1e-30)),
        "proxy_gap_Rd": float(gap),
        "eps_C_used": float(eps_use),
        "constraint_signal": 0.0,  # retired coercivity claim; no theorem established
        "sample_separation_signal": float(sample_separation),
        "reference_floor_source": "half_observed_pure" if eps_C is None else "external_numeric_reference",
        "weil_form_evaluated": False,
        "coercivity_status": CANDIDATE_STATUS,
    }


def _row_job(args: Tuple[float, float, int, float, int]) -> Dict[str, float]:
    T, sigma, degree, t, n_points = args
    return weil_local_proxy_row(T, sigma=sigma, degree=degree, t=t, n_points=n_points)


def partial_cancel_curve(
    T: float,
    *,
    sigma: float = 0.9,
    degree: int = 4,
    t: float = DEFAULT_T_FIRST_ZERO,
    n_points: int = 2048,
    alphas: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
) -> List[Dict[str, float]]:
    """W(α) = C_d(q - α P_d q) should decrease in α for pure off-modes."""
    u = sample_grid(int(n_points))
    q = _scaled_off_mode(u, T=T, sigma=sigma, t=t)
    base_ratio = energy_ratio(q, u, degree)
    rows = []
    for a in alphas:
        predicted = partial_cancellation_ratio(base_ratio, float(a))
        qa = partial_cancel(q, u, degree=degree, T=T, alpha=float(a))
        c = coercivity_Cd(qa, u, degree=degree, T=T)
        rows.append(
            {
                "T": float(T),
                "sigma": float(sigma),
                "degree": float(degree),
                "alpha": float(a),
                "C_d": float(c),
                "R_d": float(c**2),
                "R_d_predicted": predicted,
                "C_d_predicted": float(np.sqrt(predicted)),
                "identity_absolute_error": abs(float(c**2) - predicted),
            }
        )
    return rows


def multi_grid_campaign(
    *,
    T_values: Sequence[float] = DEFAULT_T,
    sigma_values: Sequence[float] = DEFAULT_SIGMA,
    degrees: Sequence[int] = DEFAULT_DEGREE,
    t: float = DEFAULT_T_FIRST_ZERO,
    n_points: int = 2048,
    max_workers: Optional[int] = None,
) -> List[Dict[str, float]]:
    """
    ProcessPool sweep over (T, σ, d). Independent units — use full_workers.
    """
    jobs = [
        (float(T), float(sig), int(d), float(t), int(n_points))
        for T in T_values
        for sig in sigma_values
        for d in degrees
    ]
    if max_workers is None:
        max_workers = min(len(jobs), max(1, (os.cpu_count() or 4) - 2))
    if len(jobs) == 1 or max_workers <= 1:
        return [_row_job(j) for j in jobs]
    rows: List[Dict[str, float]] = []
    with ProcessPoolExecutor(max_workers=int(max_workers)) as ex:
        futs = {ex.submit(_row_job, j): j for j in jobs}
        for fut in as_completed(futs):
            rows.append(fut.result())
    rows.sort(key=lambda r: (r["sigma"], r["degree"], r["T"]))
    return rows


def multi_T_proxy_campaign(
    T_values: Sequence[float] = (8.0, 12.0, 16.0, 24.0, 32.0),
    *,
    sigma: float = 0.9,
    degree: int = 4,
    n_points: int = 2048,
) -> List[Dict[str, float]]:
    """Backward-compatible single-(σ,d) multi-T slice."""
    return [
        weil_local_proxy_row(float(T), sigma=sigma, degree=degree, n_points=n_points)
        for T in T_values
    ]


def empirical_eps_star(rows: Sequence[Dict[str, float]]) -> Dict[str, Any]:
    """Finite-sample minima, grouped by (σ,d); not an asymptotic lower bound."""
    from collections import defaultdict

    groups: Dict[Tuple[float, int], List[float]] = defaultdict(list)
    for r in rows:
        groups[(float(r["sigma"]), int(r["degree"]))].append(float(r["C_d_pure"]))
    envelope = {}
    for (sig, d), vals in sorted(groups.items()):
        envelope[f"sigma={sig:g}|d={d}"] = {
            "min_C_d_pure": float(min(vals)),
            "median_C_d_pure": float(np.median(vals)),
            "n": len(vals),
            "eps_star_half_min": float(0.5 * min(vals)),
        }
    global_min = float(min(r["C_d_pure"] for r in rows)) if rows else 0.0
    return {
        "global_min_C_d_pure": global_min,
        "eps_star_global_half_min": 0.5 * global_min,
        "by_sigma_degree": envelope,
        "scope": "observed_finite_grid_only",
        "is_uniform_lower_bound": False,
    }


def write_campaign_artifacts(
    out_dir: Path | str,
    *,
    T_values: Sequence[float] = DEFAULT_T,
    sigma_values: Sequence[float] = DEFAULT_SIGMA,
    degrees: Sequence[int] = DEFAULT_DEGREE,
    max_workers: Optional[int] = None,
    make_plot: bool = True,
) -> Dict[str, str]:
    """Write multi_T.json, grid.json, partial_cancel.json, SUMMARY, optional PNG."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    grid_rows = multi_grid_campaign(
        T_values=T_values,
        sigma_values=sigma_values,
        degrees=degrees,
        max_workers=max_workers,
    )
    eps = empirical_eps_star(grid_rows)
    # Compare observed samples with a finite-grid reference, not a theorem.
    eps_floor = float(eps["eps_star_global_half_min"])
    n_ok = 0
    for r in grid_rows:
        # cancelled must sit below empirical floor; pure above it
        ok = (r["C_d_pure"] >= eps_floor) and (r["C_d_cancelled"] < eps_floor) and (
            r["sample_separation_signal"] >= 0.5
        )
        r["eps_star_floor"] = eps_floor
        r["meets_eps_star"] = 0.0  # retired uniform-coercivity interpretation
        r["meets_sample_reference"] = float(ok)
        n_ok += int(ok)

    # Backward-compatible multi_T slice at σ=0.9,d=4
    slice_rows = [r for r in grid_rows if abs(r["sigma"] - 0.9) < 1e-9 and int(r["degree"]) == 4]
    if not slice_rows:
        slice_rows = multi_T_proxy_campaign(T_values=T_values, sigma=0.9, degree=4)

    partial = []
    for T in (16.0, 32.0, 64.0):
        partial.extend(partial_cancel_curve(T, sigma=0.9, degree=4))

    grid_payload = {
        "banner": BANNER,
        "rh_claimed": False,
        "lemma_id": LEMMA_ID,
        "n_rows": len(grid_rows),
        "T_values": list(map(float, T_values)),
        "sigma_values": list(map(float, sigma_values)),
        "degrees": list(map(int, degrees)),
        "eps_star": eps,
        "candidate_status": CANDIDATE_STATUS,
        "n_meets_eps_star": 0,
        "fraction_meets_eps_star": 0.0,
        "n_meets_sample_reference": n_ok,
        "fraction_meets_sample_reference": n_ok / max(len(grid_rows), 1),
        "rows": grid_rows,
        "note": "Finite-grid model diagnostics; coercivity interpretation retracted; not arithmetic B-RES or RH.",
    }
    multi_payload = {
        "banner": BANNER,
        "rh_claimed": False,
        "lemma_id": LEMMA_ID,
        "sigma": 0.9,
        "degree": 4,
        "rows": slice_rows,
        "candidate_status": CANDIDATE_STATUS,
        "all_constraint_signals": False,
        "all_sample_separation_signals": bool(slice_rows) and all(r["sample_separation_signal"] >= 0.5 for r in slice_rows),
        "eps_star_floor": eps_floor,
        "note": (
            "Finite σ=0.9,d=4 sample of projection and cancellation. "
            "No Weil form is evaluated and no uniform visibility floor follows."
        ),
    }
    partial_payload = {
        "banner": BANNER,
        "rh_claimed": False,
        "rows": partial,
        "note": "Orthogonal-projector identity: R_alpha=(1-alpha)^2 R/[1-(2alpha-alpha^2)R].",
    }

    paths = {}
    paths["grid"] = str(out / "grid.json")
    Path(paths["grid"]).write_text(json.dumps(grid_payload, indent=2))
    paths["multi_T"] = str(out / "multi_T.json")
    Path(paths["multi_T"]).write_text(json.dumps(multi_payload, indent=2))
    paths["partial"] = str(out / "partial_cancel.json")
    Path(paths["partial"]).write_text(json.dumps(partial_payload, indent=2))

    analytic_rows = analytic_envelope_campaign(max_workers=max_workers)
    n_above = sum(1 for r in analytic_rows if r["numeric_above_half_leading"] >= 0.5)
    analytic_payload = {
        "banner": BANNER,
        "rh_claimed": False,
        "candidate_status": CANDIDATE_STATUS,
        "formula": "R_d=4(d+1)^2 B(theta)^2/[T H(theta)]+O(T^-2), theta=tT",
        "limsup_formula": "limsup T R_d=4 A(d+1)^2/(A^2+t^2), A=sigma-1/2",
        "half_envelope_is_lower_bound": False,
        "median_rel_err_role": "sampled_vs_exact_continuous_projection",
        "n_rows": len(analytic_rows),
        "n_numeric_above_half_leading": n_above,
        "fraction_above": n_above / max(len(analytic_rows), 1),
        "median_rel_err": float(np.median([r["rel_err_Rd"] for r in analytic_rows])),
        "rows": analytic_rows,
    }
    paths["analytic"] = str(out / "analytic_envelope.json")
    Path(paths["analytic"]).write_text(json.dumps(analytic_payload, indent=2))

    bridge_rows = arithmetic_bridge_campaign(max_workers=max_workers)
    n_bridge_ok = sum(1 for r in bridge_rows if r["bridge_ok"] >= 0.5)
    bridge_payload = {
        "banner": BANNER,
        "rh_claimed": False,
        "b_res_solved": False,
        "n_rows": len(bridge_rows),
        "n_bridge_ok": n_bridge_ok,
        "fraction_ok": n_bridge_ok / max(len(bridge_rows), 1),
        "mean_pred_ratio": float(np.mean([r["pred_ratio"] for r in bridge_rows])),
        "rows": bridge_rows,
        "note": (
            "After V_d-kill, an inserted model bump verifies projection linearity "
            "with matching weighted norms. It is not forced by an arithmetic zero."
        ),
    }
    paths["bridge"] = str(out / "arithmetic_bridge.json")
    Path(paths["bridge"]).write_text(json.dumps(bridge_payload, indent=2))

    summary = out / "SUMMARY.txt"
    lines = [
        BANNER,
        f"lemma={LEMMA_ID} status={CANDIDATE_STATUS}",
        f"grid_rows={len(grid_rows)} meets_sample_reference={n_ok}/{len(grid_rows)} "
        f"({100*n_ok/max(len(grid_rows),1):.1f}%)",
        f"finite_grid_half_min_reference={eps_floor:.6e} (not a uniform lower bound)",
        f"slice_sigma0.9_d4 n={len(slice_rows)} all_constraint="
        f"{multi_payload['all_constraint_signals']}",
    ]
    for r in slice_rows:
        lines.append(
            f"T={r['T']:.1f}  R_pure={r['R_d_pure']:.4e}  R_cancel={r['R_d_cancelled']:.4e}  "
            f"C_pure={r['C_d_pure']:.4e}  sample_signal={int(r['sample_separation_signal'])}"
        )
    # Monotone check on partial curves
    if "analytic" in paths:
        lines.append(
            f"analytic_envelope n={analytic_payload['n_rows']} above_half={n_above}/{analytic_payload['n_rows']} "
            f"median_rel_err={analytic_payload['median_rel_err']:.4f}"
        )
    if "bridge" in paths:
        lines.append(
            f"arith_bridge n={bridge_payload['n_rows']} ok={n_bridge_ok}/{bridge_payload['n_rows']} "
            f"mean_pred_ratio={bridge_payload['mean_pred_ratio']:.4f}"
        )
    for T in (16.0, 32.0, 64.0):
        curve = [r for r in partial if abs(r["T"] - T) < 1e-9]
        cds = [r["C_d"] for r in curve]
        mono = all(cds[i] >= cds[i + 1] - 1e-9 for i in range(len(cds) - 1))
        lines.append(f"partial_cancel T={T:.0f} monotone_decreasing={mono} C_d={cds}")
    summary.write_text("\n".join(lines) + "\n")
    paths["summary"] = str(summary)

    if make_plot:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            ax = axes[0]
            for sig in sigma_values:
                series = [
                    r
                    for r in grid_rows
                    if abs(r["sigma"] - sig) < 1e-9 and int(r["degree"]) == 4
                ]
                series.sort(key=lambda r: r["T"])
                ax.plot(
                    [r["T"] for r in series],
                    [r["C_d_pure"] for r in series],
                    "o-",
                    label=f"σ={sig:g} pure",
                )
                ax.plot(
                    [r["T"] for r in series],
                    [r["C_d_cancelled"] for r in series],
                    "x--",
                    label=f"σ={sig:g} cancel",
                )
            ax.axhline(eps_floor, color="k", ls=":", label="finite-grid half-min reference")
            ax.set_xlabel("T")
            ax.set_ylabel("C_d")
            ax.set_title("Polynomial C_d: pure vs V_d-kill (d=4)")
            ax.legend(fontsize=7, ncol=2)
            ax.grid(True, alpha=0.3)

            ax2 = axes[1]
            for T in (16.0, 32.0, 64.0):
                curve = [r for r in partial if abs(r["T"] - T) < 1e-9]
                curve.sort(key=lambda r: r["alpha"])
                ax2.plot(
                    [r["alpha"] for r in curve],
                    [r["C_d"] for r in curve],
                    "o-",
                    label=f"T={T:.0f}",
                )
            ax2.set_xlabel("cancel fraction α")
            ax2.set_ylabel("C_d(q-α P_d q)")
            ax2.set_title("Partial V_d kill (σ=0.9, d=4)")
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            fig.suptitle("Retracted coercivity proxy: model diagnostics")
            fig.tight_layout()
            plot_path = out / "hstar_proxy.png"
            fig.savefig(plot_path, dpi=120)
            plt.close(fig)
            paths["plot"] = str(plot_path)
        except Exception as exc:  # pragma: no cover
            paths["plot_error"] = str(exc)

    return paths


def write_multi_T_artifacts(
    out_dir: Path | str,
    *,
    T_values: Sequence[float] = (8.0, 12.0, 16.0, 24.0, 32.0),
    sigma: float = 0.9,
    degree: int = 4,
) -> Path:
    """Compat wrapper: full campaign with at least the requested T slice."""
    paths = write_campaign_artifacts(
        out_dir,
        T_values=tuple(sorted(set(list(T_values) + list(DEFAULT_T)))),
        sigma_values=DEFAULT_SIGMA if sigma in DEFAULT_SIGMA else DEFAULT_SIGMA + (sigma,),
        degrees=DEFAULT_DEGREE if degree in DEFAULT_DEGREE else DEFAULT_DEGREE + (degree,),
        max_workers=max(1, (os.cpu_count() or 4) - 2),
        make_plot=True,
    )
    return Path(paths["multi_T"])


def anticancellation_candidates() -> List[Dict[str, str]]:
    """Separate axis from Rank-4 incremental_inequality_candidates()."""
    return [
        {
            "id": LEMMA_ID,
            "idea": (
                "Historical local coercivity proxy retracted: intact real off-modes "
                "already have C_d→0; projection-cancellation identities remain valid."
            ),
            "status": CANDIDATE_STATUS,
            "axis": "anti_cancellation_toward_H_star",
            "doc": "docs/HSTAR_ANTICANCELLATION.md",
        }
    ]


def hstar_anticancellation_report(
    *,
    T_values: Sequence[float] = DEFAULT_T,
    sigma_values: Sequence[float] = DEFAULT_SIGMA,
    degrees: Sequence[int] = DEFAULT_DEGREE,
    sigma: float = 0.9,
    degree: int = 4,
    write_results: bool = False,
    results_dir: Path | str | None = None,
    max_workers: Optional[int] = None,
    heavy: bool = False,
) -> Dict[str, Any]:
    """
    Shipped report.

    Default: full (T,σ,d) grid via ProcessPool when heavy=True or write_results=True;
    otherwise a fast σ/d slice for CLI/tests.
    """
    if len(T_values) == 0:
        raise ValueError("T_values must be nonempty")
    lemma = lemma_statement()
    if heavy or write_results:
        workers = max_workers
        if workers is None:
            workers = max(1, (os.cpu_count() or 4) - 2)
        grid_rows = multi_grid_campaign(
            T_values=T_values,
            sigma_values=sigma_values,
            degrees=degrees,
            max_workers=workers,
        )
        eps = empirical_eps_star(grid_rows)
        slice_rows = [
            r
            for r in grid_rows
            if abs(r["sigma"] - float(sigma)) < 1e-9 and int(r["degree"]) == int(degree)
        ]
        if not slice_rows:
            slice_rows = multi_T_proxy_campaign(T_values=T_values, sigma=sigma, degree=degree)
        partial = partial_cancel_curve(32.0, sigma=sigma, degree=degree)
        analytic_rows = analytic_envelope_campaign(max_workers=workers)
        bridge_rows = arithmetic_bridge_campaign(max_workers=workers)
        artifacts = None
        if write_results:
            root = Path(results_dir) if results_dir else Path("results/hstar_corrected")
            artifacts = write_campaign_artifacts(
                root,
                T_values=T_values,
                sigma_values=sigma_values,
                degrees=degrees,
                max_workers=workers,
                make_plot=True,
            )
        all_ok = bool(slice_rows) and all(r["sample_separation_signal"] >= 0.5 for r in slice_rows)
        frac = float(np.mean([r["sample_separation_signal"] for r in grid_rows])) if grid_rows else 0.0
    else:
        slice_rows = multi_T_proxy_campaign(T_values=T_values[:5], sigma=sigma, degree=degree)
        grid_rows = slice_rows
        eps = empirical_eps_star(slice_rows)
        partial = partial_cancel_curve(32.0, sigma=sigma, degree=degree)
        analytic_rows = []
        bridge_rows = []
        artifacts = None
        all_ok = bool(slice_rows) and all(r["sample_separation_signal"] >= 0.5 for r in slice_rows)
        frac = float(np.mean([r["sample_separation_signal"] for r in slice_rows])) if slice_rows else 0.0

    mid_T = float(T_values[len(T_values) // 2])
    off = off_critical_rd_lower_model(mid_T, sigma=sigma, degree=degree)
    cancel = cancelled_off_critical_rd(mid_T, sigma=sigma, degree=degree)

    # Partial-cancel monotone at T=32
    cds = [r["C_d"] for r in partial]
    monotone = all(cds[i] >= cds[i + 1] - 1e-9 for i in range(len(cds) - 1))

    return {
        "banner": BANNER,
        "rh_claimed": False,
        "b_res_solved": False,
        "full_H_star_proved": False,
        "axis": "anti_cancellation_toward_H_star",
        "title": "L* proxy retraction and phase-aware model diagnostics",
        "weil_form_evaluated": False,
        "positive_limiting_model_coercivity": False,
        "proved_model_obstruction": "docs/OFF_CRITICAL_PHASE_OBSTRUCTION.md",
        "lemma": lemma,
        "candidate_status": CANDIDATE_STATUS,
        "candidates": anticancellation_candidates(),
        "multi_T_proxy": slice_rows,
        "grid_n_rows": len(grid_rows),
        "grid_constraint_fraction": 0.0,
        "grid_sample_separation_fraction": frac,
        "eps_star": eps,
        "partial_cancel_T32": partial,
        "partial_cancel_monotone": bool(monotone),
        "analytic_envelope": (
            {
                "n": len(analytic_rows),
                "half_envelope_is_lower_bound": False,
                "median_rel_err_role": "sampled_vs_exact_continuous_projection",
                "fraction_above_half_leading": float(np.mean([r["numeric_above_half_leading"] for r in analytic_rows])),
                "median_rel_err": float(np.median([r["rel_err_Rd"] for r in analytic_rows])),
                "sample": analytic_rows[:3],
            }
            if heavy or write_results
            else None
        ),
        "arithmetic_bridge": (
            {
                "n": len(bridge_rows),
                "role": "projection_identity_calibration_only",
                "arithmetic_off_zero_forced": False,
                "fraction_ok": float(np.mean([r["bridge_ok"] for r in bridge_rows])),
                "mean_pred_ratio": float(np.mean([r["pred_ratio"] for r in bridge_rows])),
                "sample": bridge_rows[:3],
            }
            if heavy or write_results
            else None
        ),
        "all_constraint_signals": False,
        "all_sample_separation_signals": bool(all_ok),
        "rank5_crosscheck": {
            "model_off_critical": off,
            "model_cancellation_counterexample": cancel,
            "pure_above_cancel": bool(off["R_d_off"] > cancel["R_d_after_Vd_kill"] + 1e-6),
        },
        "rank4_proportion_axis": "STOP (unchanged; see zero_proportion_feasibility)",
        "artifacts": artifacts,
        "assimilated": [],
        "retained_model_identities": [
            "exact orthogonal-projector partial-cancellation identity",
            "phase-aware continuous coefficients, norm, and sharp subsequence limits",
            "injected-bump calibration using the same weighted inner product",
        ],
        "rejected": [
            "phase-independent real-sine asymptotic or positive limiting C_d floor",
            "identifying polynomial projection energy with a Weil quadratic form",
            "deriving arithmetic non-cancellation from an injected model bump",
            "BGST/Bombieri proportion / pair-correlation as R_d engine",
            "Anthropic 41.6%→67.2% reproduction",
            "full Guinand–Weil positivity cascade as RH proof",
        ],
        "conclusion": (
            "The L* coercivity proxy is retracted. Intact continuous real off-modes "
            "have R_d→0 at fixed degree; the old formula is a phase limsup envelope. "
            "Finite-grid separation and injected bumps prove no arithmetic visibility "
            "floor. Orthogonal projection identities remain valid; B-RES and RH stay open."
        ),
    }
