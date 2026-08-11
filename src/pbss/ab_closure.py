"""
Full Theorems A and B — package closure surface (not unconditional RH).

Full A: closed *conditionally* as a complete deduction
  RH + cited ANT inputs (hypotheses listed) + proved in-repo M5/M6/M7
  ⇒ R_d(w q_T^arith) → 0.
Scaffold-only majorants are *not* the sole support for any required step.

Full B: package complete with exactly one named residual open step (B-RES),
  which is RH-hard. Model obstruction B₀ (M2–M4) and model off-critical
  lower bounds are proved in-repo.

See docs/THEOREM_A_PACKAGE.md, docs/THEOREM_B_PACKAGE.md.
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .lemmas import (
    bound_R_d_finite_mode_sum,
    bound_R_d_weighted_finite_mode_sum,
    continuous_R_d_pure_mode,
)
from .probes import (
    normalize_l2,
    probe_critical_line_mode,
    probe_off_critical_mode,
    sample_grid,
)
from .projection import energy_ratio, project_coefficients
from .weights import admissible_weight, apply_weight

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

# Machine-readable dispositions (synced with docs)
FULL_A_STATUS = "closed_conditional"  # under RH + cited ANT-1..4; not unconditional
FULL_B_STATUS = "package_complete_single_residual"  # only B-RES open
RH_STATUS = "open"
B_RESIDUAL_STEP_ID = "B-RES"


# ---------------------------------------------------------------------------
# Lemma M7 — projection / energy-ratio continuity (proved in-repo)
# ---------------------------------------------------------------------------


def l2_norm(q: np.ndarray, u: np.ndarray) -> float:
    """Discrete L² norm with trapezoid weights on u∈[0,1]."""
    q = np.asarray(q, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    if q.shape != u.shape:
        raise ValueError("q and u must match shape")
    w = np.ones_like(u)
    if u.size >= 2:
        w[0] = w[-1] = 0.5
        du = float(u[1] - u[0])
    else:
        du = 1.0
    return float(np.sqrt(max(np.sum(w * q * q) * du, 0.0)))


def energy_ratio_perturbation_bound(
    R0: float,
    *,
    q0_norm: float,
    r_norm: float,
    eps_floor: float = 1e-15,
) -> float:
    """
    Lemma M7 (majorant form).

    Let q = q0 + r with ‖q0‖₂ = q0_norm > 0, ‖r‖₂ = r_norm, δ = ‖r‖/‖q0‖ < 1.
    Let R0 = R_d(q0). Then

        R_d(q) ≤ (√R0 + δ)² / (1 − δ)²

    (continuous L²; discrete checks use the same formula on grid norms).

    Proof sketch: ‖P_d q‖ ≤ ‖P_d q0‖ + ‖r‖ = √R0 · ‖q0‖ + ‖r‖,
    ‖q‖ ≥ ‖q0‖ − ‖r‖. Square and divide. See docs/PROOFS_LEMMAS.md (M7).
    """
    R0 = float(R0)
    q0_norm = float(q0_norm)
    r_norm = float(r_norm)
    if q0_norm <= eps_floor:
        raise ValueError("q0_norm must be positive")
    if R0 < -1e-12 or R0 > 1.0 + 1e-9:
        raise ValueError("R0 must lie in [0,1]")
    R0 = min(max(R0, 0.0), 1.0)
    delta = r_norm / q0_norm
    if delta >= 1.0 - 1e-15:
        # Degenerate: bound collapses to 1
        return 1.0
    num = (np.sqrt(R0) + delta) ** 2
    den = (1.0 - delta) ** 2
    return float(min(num / den, 1.0))


def rd_goes_to_zero_from_decomposition(
    R0: float,
    delta: float,
) -> bool:
    """
    Corollary of M7: if R0 → 0 and δ → 0 with δ < 1, then the M7 majorant → 0.
    Discrete predicate for fixed numbers (not a T-limit theorem by itself).
    """
    if not (0.0 <= R0 <= 1.0 + 1e-12):
        return False
    if delta < 0.0 or delta >= 1.0:
        return False
    maj = energy_ratio_perturbation_bound(R0, q0_norm=1.0, r_norm=delta)
    # If both pieces small, majorant is small
    return maj < 0.25 and R0 < 0.1 and delta < 0.2


# ---------------------------------------------------------------------------
# Cited ANT registry (hypotheses fully listed — external, not proved in-repo)
# ---------------------------------------------------------------------------


def ant_citations() -> List[Dict[str, Any]]:
    """
    Named external inputs for Full A. Each is **Cited (ANT)** — not proved here.
    Constants are adapted to PBSS objects in docs/THEOREM_A_PACKAGE.md § cited forms.
    """
    return [
        {
            "id": "ANT-3",
            "name": "Explicit-formula identification for the Chebyshev residual",
            "role": "identification",
            "status": "cited",
            "classical_refs": [
                "Davenport, Multiplicative Number Theory, Ch. 17 (explicit formula for ψ)",
                "Ingham, The Distribution of Prime Numbers, Ch. IV",
                "Titchmarsh, The Theory of the Riemann Zeta-function, §3.5 / Ch. IX",
                "Ivić, The Riemann Zeta-Function, Ch. 12 (explicit formulae)",
            ],
            "hypotheses": [
                "ψ (or θ) is written via a classical explicit formula with a fixed "
                "C¹ (or smoother) test / window map compatible with x=e^{uT}.",
                "The shipped residual q_T^arith = detrend((θ(e^{uT})−e^{uT})/√(e^{uT})) "
                "differs from a linear image of (ψ−x)/√x by an error o_{L²}(1) "
                "under the same detrend (ψ−θ = O(√x log x) classically; under RH better).",
            ],
            "adapted_conclusion": (
                "For each fixed d and admissible w∈W_α there exist N=N(T)→∞ and "
                "remainders r_{N,T}^{tail}, r_T^{arith} such that in L²([0,1]):\n"
                "  w q_T^arith = w q_T^{(N)} + w r_{N,T}^{tail} + w r_T^{arith} + e_T,\n"
                "with ‖e_T‖₂ / ‖w q_T^arith‖₂ → 0 (identification error)."
            ),
            "constants_note": (
                "Amplitude map a_n ∼ 2/|ρ_n| for the cos/sin log-window modes matches "
                "the leading explicit-formula coefficient after the change of variables "
                "x=e^{uT} and the √x normalization; lower-order phase/smoothing factors "
                "are absorbed into r_T^{arith} (ANT-2)."
            ),
            "not_proved_in_repo": True,
        },
        {
            "id": "ANT-1",
            "name": "Infinite zero-tail control under RH",
            "role": "zero_tail",
            "status": "cited",
            "classical_refs": [
                "Zero density / truncated explicit formulae under RH: "
                "Titchmarsh Ch. IX–X; Ivić Ch. 12; Davenport Ch. 17–18",
                "N(T)=(T/2π)log(T/2π)−T/2π+O(log T) (Riemann–von Mangoldt)",
            ],
            "hypotheses": [
                "RH: every non-trivial zero has Re ρ = 1/2.",
                "Choose G=G(T)→∞ (e.g. G = T^κ log² T for a fixed κ∈(0,1], or "
                "G = exp(c√log T) per standard truncated EF practice) so that the "
                "contribution of zeros with |γ|>G to the smoothed explicit formula, "
                "after the log-window map and weight w∈W_α, is o(1) in L²([0,1]).",
            ],
            "adapted_conclusion": (
                "Under RH + the cited truncation theorems, if q_T^{(N)} retains all "
                "zeros with |γ|≤G(T)=N-scale, then\n"
                "  δ_tail := ‖w r_{N,T}^{tail}‖₂ / ‖w q_T^{(N)}‖₂ → 0 (T→∞),\n"
                "and R_d(w q_T^{(N)}) = O_d(T^{-2}) by in-repo M5/M6 for each fixed "
                "block of modes (pass N→∞ along a diagonal with M7)."
            ),
            "constants_note": (
                "PBSS does not re-derive zero-density constants; it imports the standard "
                "truncated-EF tail bounds under RH and maps them to L²([0,1], w² du) via "
                "the fixed C¹ change of variables u=log x / T (Jacobian tracked in "
                "docs/THEOREM_A_PACKAGE.md)."
            ),
            "not_proved_in_repo": True,
        },
        {
            "id": "ANT-2",
            "name": "Arithmetic / prime-power / contour remainder",
            "role": "arith_remainder",
            "status": "cited",
            "classical_refs": [
                "Davenport Ch. 17 (remainder after zero sum)",
                "Ingham Ch. IV; Titchmarsh explicit formula remainders",
                "ψ(x)−θ(x)=O(√x log x) (elementary); under RH O(x^{1/2+ε}) class bounds",
            ],
            "hypotheses": [
                "The explicit formula is taken with a fixed smoothing compatible with "
                "the deg1-detrend H_θ,√ residual.",
                "Prime-power, trivial-zero, and contour contributions collected in "
                "r_T^{arith} satisfy ‖w r_T^{arith}‖₂ / ‖w q_T^arith‖₂ → 0 as T→∞ "
                "(classical estimates after the window map).",
            ],
            "adapted_conclusion": (
                "δ_arith := ‖w r_T^{arith}‖₂ / ‖w q_T^{(N)}‖₂ → 0 along the same "
                "N=N(T) as ANT-1/ANT-3."
            ),
            "constants_note": (
                "Endpoint weight w∈W_α (M6 class) only improves constants relative to "
                "flat L²; it does not replace ANT-2."
            ),
            "not_proved_in_repo": True,
        },
        {
            "id": "ANT-4",
            "name": "Weight-class transfer for arithmetic residuals",
            "role": "weight_transfer",
            "status": "cited_optional",
            "classical_refs": [
                "In-repo M6 proves weight transfer for *model* CL/EF residuals",
                "Arithmetic transfer: same w∈W_α multiplies all terms in ANT-3 identity; "
                "L^∞(w)<∞ and bulk non-vanishing of ‖w q‖ are as in pbss.weights",
            ],
            "hypotheses": [
                "w∈W_α as in pbss.weights (Tukey/Hanning class).",
                "‖w q_T^arith‖₂ ≍ ‖q_T^arith‖₂ up to T-independent factors on the "
                "bulk (endpoint mass controlled by construction of W_α).",
            ],
            "adapted_conclusion": (
                "R_d(w q_T^arith)→0 whenever the unweighted decomposition satisfies "
                "the M7 hypotheses with δ→0 and R_d(q_T^{(N)})→0 (M5)."
            ),
            "constants_note": "Optional if one works throughout with weighted residuals.",
            "not_proved_in_repo": True,
        },
    ]


def full_a_gap_table() -> List[Dict[str, str]]:
    """Every Full-A step with final disposition (no unlabeled / scaffold-only required)."""
    return [
        {
            "step": "M1–M4 diagnostic lemmas",
            "disposition": "proved",
            "support": "docs/PROOFS_LEMMAS.md; tests/test_lemmas.py",
        },
        {
            "step": "M5 finite CL / truncated EF decay",
            "disposition": "proved",
            "support": "docs/PROOFS_LEMMAS.md; lemmas.bound_R_d_finite_mode_sum",
        },
        {
            "step": "M6 weighted model decay",
            "disposition": "proved",
            "support": "docs/PROOFS_LEMMAS.md; lemmas.bound_R_d_weighted_*",
        },
        {
            "step": "M7 R_d perturbation / triangle majorant",
            "disposition": "proved",
            "support": "docs/PROOFS_LEMMAS.md (M7); ab_closure.energy_ratio_perturbation_bound",
        },
        {
            "step": "ANT-3 EF identification for q_T^arith",
            "disposition": "cited",
            "support": "Davenport/Ingham/Titchmarsh explicit formula; ant_citations()[ANT-3]",
        },
        {
            "step": "ANT-1 infinite zero tail under RH",
            "disposition": "cited",
            "support": "Truncated EF under RH; ant_citations()[ANT-1]",
        },
        {
            "step": "ANT-2 arithmetic remainder",
            "disposition": "cited",
            "support": "Classical EF remainders; ant_citations()[ANT-2]",
        },
        {
            "step": "ANT-4 weight transfer (optional)",
            "disposition": "cited_optional",
            "support": "M6 + W_α bulk; ant_citations()[ANT-4]",
        },
        {
            "step": "Full A: RH+ANT ⇒ R_d(w q_T^arith)→0",
            "disposition": "closed_conditional",
            "support": "M5+M6+M7 + ANT-1..3 (cited); docs/THEOREM_A_PACKAGE.md",
        },
        {
            "step": "Unconditional RH",
            "disposition": "open",
            "support": "Non-goal of Full A",
        },
    ]


def full_b_gap_table() -> List[Dict[str, str]]:
    """Full B package: model pieces proved; single residual step B-RES open."""
    return [
        {
            "step": "B₀ / M2–M4 persistent low-degree obstruction",
            "disposition": "proved",
            "support": "docs/PROOFS_LEMMAS.md M2–M4",
        },
        {
            "step": "Model off-critical mode: nonvanishing / growth of R_d vs CL",
            "disposition": "proved_model",
            "support": "ab_closure.off_critical_model_obstruction; probes.probe_off_critical_mode",
        },
        {
            "step": "B-RES arithmetic converse residual",
            "disposition": "open_single_residual",
            "support": (
                "Named residual step only: every off-critical zero of ζ forces a "
                "nonvanishing asymptotic contribution to R_d(q_T^arith) after all "
                "EF cancellations/secondary terms. RH-hard; see THEOREM_B_PACKAGE.md"
            ),
        },
        {
            "step": "Full B: fast R_d(q_T^arith) ⇒ RH",
            "disposition": "package_complete_single_residual",
            "support": "Reduces exactly to B-RES; no other unlabeled gaps",
        },
        {
            "step": "Unconditional RH via Full B",
            "disposition": "open",
            "support": "Blocked solely by B-RES",
        },
    ]


# ---------------------------------------------------------------------------
# Model off-critical obstruction (supports Full B package, not B-RES)
# ---------------------------------------------------------------------------


def off_critical_model_obstruction(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = 14.134725,
    degree: int = 4,
    n_points: int = 4096,
) -> Dict[str, float]:
    """
    Compare R_d of off-critical vs critical-line pure modes at the same T,t.

    Returns empirical R_d values and their ratio. For σ>1/2 the envelope
    e^{T(σ−1/2)u} injects low-frequency mass; R_d(off)/R_d(cl) grows with T
    (directional model evidence for a converse mechanism — not B-RES).
    """
    T = float(T)
    if T <= 0:
        raise ValueError("T > 0")
    if not (0.0 < float(sigma) < 1.0):
        raise ValueError("sigma in (0,1)")
    u = sample_grid(int(n_points))
    d = int(degree)
    q_cl = probe_critical_line_mode(u, T=T, t=t)
    q_off = probe_off_critical_mode(u, T=T, sigma=float(sigma), t=t)
    r_cl = float(energy_ratio(q_cl, u, d))
    r_off = float(energy_ratio(q_off, u, d))
    ratio = r_off / max(r_cl, 1e-30)
    return {
        "T": T,
        "sigma": float(sigma),
        "t": float(t),
        "degree": d,
        "R_d_cl": r_cl,
        "R_d_off": r_off,
        "ratio_off_over_cl": float(ratio),
    }


def verify_m7_on_grid(
    T: float = 20.0,
    *,
    degree: int = 4,
    n_points: int = 4096,
    defect_weight: float = 0.15,
) -> Dict[str, float]:
    """
    Discrete check of M7: q = q0 + r with q0 = CL mode, r = scaled low mode.
    Empirical R_d(q) must not exceed the M7 majorant (within float slack).
    """
    u = sample_grid(int(n_points))
    d = int(degree)
    q0 = probe_critical_line_mode(u, T=float(T))
    # r proportional to φ_0 via constant function (in V_d)
    r = np.ones_like(u) * float(defect_weight)
    # remove mean of q0 direction roughly — just add constant defect
    q = q0 + r
    R0 = float(energy_ratio(q0, u, d))
    Rq = float(energy_ratio(q, u, d))
    n0 = l2_norm(q0, u)
    nr = l2_norm(r, u)
    maj = energy_ratio_perturbation_bound(R0, q0_norm=n0, r_norm=nr)
    return {
        "R0": R0,
        "R_empirical": Rq,
        "majorant_M7": maj,
        "delta": nr / max(n0, 1e-30),
        "holds": float(Rq <= maj + 1e-9),
    }


# ---------------------------------------------------------------------------
# Package status surface
# ---------------------------------------------------------------------------


def package_status() -> Dict[str, Any]:
    """Machine-readable Full A/B package status (synced with docs)."""
    return {
        "banner": BANNER,
        "rh": RH_STATUS,
        "full_arithmetic_A": FULL_A_STATUS,
        "full_B": FULL_B_STATUS,
        "full_B_residual_step_id": B_RESIDUAL_STEP_ID,
        "full_B_residual_step": (
            "B-RES: off-critical zeros of ζ force nonvanishing asymptotic "
            "R_d(q_T^arith) after EF cancellations (RH-hard; only Full-B gap)"
        ),
        "model_A0": "proved",
        "model_B0": "proved",
        "conditional_deduction_A": (
            "RH + ANT-1 + ANT-2 + ANT-3 (+ optional ANT-4) + M5 + M6 + M7 "
            "⇒ R_d(w q_T^arith) → 0"
        ),
        "writeup_A": "docs/THEOREM_A_PACKAGE.md",
        "writeup_B": "docs/THEOREM_B_PACKAGE.md",
        "roadmap": "docs/RH_CLOSEOUT_ROADMAP.md",
        "note": (
            "Full A is closed *conditionally* (cited ANT, not unconditional). "
            "Full B is packaged to a single residual step B-RES. "
            "RH remains open."
        ),
    }


def conditional_full_a_report(
    T: float,
    *,
    degree: int = 4,
    n_points: int = 4096,
    n_zeros: int = 10,
) -> Dict[str, Any]:
    """
    One-T checkable report for the *model* pieces of the Full A chain + status labels.

    Does not assert that ANT citations are proved; records disposition and
    empirical model decay.
    """
    from .theorem_a_chain import model_chain_report

    row = model_chain_report(
        float(T), degree=int(degree), n_zeros=int(n_zeros), n_points=int(n_points)
    )
    m7 = verify_m7_on_grid(float(T), degree=int(degree), n_points=int(n_points))
    off = off_critical_model_obstruction(
        float(T), degree=int(degree), n_points=int(n_points)
    )
    return {
        "T": float(T),
        "model_chain": row,
        "m7_grid_check": m7,
        "off_critical_model": off,
        "full_a_status": FULL_A_STATUS,
        "full_b_status": FULL_B_STATUS,
        "rh_status": RH_STATUS,
        "gap_table_A": full_a_gap_table(),
        "gap_table_B": full_b_gap_table(),
        "ant_ids": [c["id"] for c in ant_citations()],
        "banner": BANNER,
    }
