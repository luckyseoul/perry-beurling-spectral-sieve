"""Proof-status registry and sampled model checks; not an unconditional RH proof.

Continuous A: RH implies R_d=O(T^-2) for the exact theta residual with the
continuous orthogonal affine detrend, fixed d>=2, and no additional weight.
B-ABS: one absolute Legendre moment has subexponential growth iff RH.
Signed moments: every off-line zero forces both signs of exponential excursions;
either one-sided subexponential bound is RH-equivalent and unproved.
Normalized recurrence: liminf R_d=0 unconditionally for that continuous detrend.
An attained rightmost off-line zero implies limsup T*R_d>0. The case of an
unattained rightmost real-part supremum remains open.
The positive-floor version of B-RES is consequently RH-equivalent and open.

The sampled functions below check model calculations, not these limit theorems.
See CONTINUOUS_THEOREM_A.md, ABSOLUTE_MOMENT_CONVERSE.md,
SIGNED_MOMENT_OSCILLATION.md, NORMALIZED_RECURRENCE_THEOREM.md, and
ATTAINED_BOUNDARY_CONVERSE.md for the actual proofs and hypotheses.
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
    sample_grid,
)
from .projection import energy_ratio, project_coefficients
from .weights import admissible_weight, apply_weight

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

# Machine-readable dispositions (synced with docs)
FULL_A_STATUS = "proved_conditional_continuous"
FULL_B_STATUS = "open_normalized_converse"
RH_STATUS = "open"
B_RESIDUAL_STEP_ID = "B-RES"

CONTINUOUS_SCOPE = {
    "residual": "Q(y)=exp(-y/2)*(theta(exp(y))-exp(y))",
    "window": "q_T(u)=(I-P_1)Q(Tu), u in [0,1]",
    "measure": "continuous Lebesgue du",
    "detrend": "continuous orthogonal affine projection P_1",
    "degree": "fixed integer d>=2",
    "weight": "none beyond Lebesgue measure",
    "sampled_transfer_proved": False,
    "weighted_transfer_proved": False,
}


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
    """Exact external inputs and their proved continuous adaptations.

    This is a source registry, not an executable verification of the citations.
    The former generic ANT adaptations are superseded by the cited equations
    and the deductions in docs/CONTINUOUS_THEOREM_A.md.
    """
    source = "https://www.cs.uleth.ca/~akbary/Akbary-Ng-Shahabi.pdf"
    return [
        {
            "id": "ANT-3",
            "name": "Exact normalized psi explicit formula",
            "role": "identification",
            "status": "cited_exact_statement",
            "classical_refs": [
                "Akbary-Ng-Shahabi, arXiv:1306.1657v1, equation (1.5), p. 2",
            ],
            "source_url": source,
            "hypotheses": ["RH", "y>=log(2), cutoff G>1, exact psi residual"],
            "adapted_conclusion": (
                "F_psi(y)=Re sum_{0<gamma<=G}[-2 m_gamma exp(i gamma y)/(1/2+i gamma)] "
                "+ O(exp(y/2)*(y+log G)^2/G + y exp(-y/2))."
            ),
            "constants_note": "The negative sign, phases and multiplicities are retained.",
            "not_proved_in_repo": True,
            "adaptation_proved_in_repo": True,
            "proof": "docs/CONTINUOUS_THEOREM_A.md, section 1",
        },
        {
            "id": "ANT-1",
            "name": "Uniform control of the growing zero-block moments",
            "role": "moment_uniformity",
            "status": "proved_continuous_adaptation",
            "classical_refs": [
                "Akbary-Ng-Shahabi equation (1.5); classical N(H)=O(H log(H+2))",
            ],
            "source_url": source,
            "hypotheses": [
                "RH", "G=exp(T)", "fixed polynomial test g on [0,1]",
            ],
            "adapted_conclusion": (
                "Integral_0^1 F_psi(Tu)g(u)du=O_g(1/T), uniformly in the retained "
                "zero block because sum m_gamma/(gamma*abs(1/2+i gamma)) converges."
            ),
            "constants_note": (
                "Uses direct moment bounds; does not assert a relative L2 tail "
                "estimate or use fixed-N M5 constants on an arbitrary N(T)."
            ),
            "not_proved_in_repo": False,
            "adaptation_proved_in_repo": True,
            "proof": "docs/CONTINUOUS_THEOREM_A.md, section 1",
        },
        {
            "id": "ANT-2",
            "name": "Explicit prime-square correction and integrable secondary error",
            "role": "arith_remainder",
            "status": "proved_continuous_adaptation",
            "classical_refs": [
                "Exact identity psi(x)-theta(x)=sum_{j>=2} theta(x^(1/j)); RH error estimate",
            ],
            "hypotheses": ["RH", "exact theta and psi", "continuous affine detrend"],
            "adapted_conclusion": (
                "Q_theta(y)=F_psi(y)-1+h(y), h in L1 intersect L2. "
                "The detrend removes -1; h(Tu) contributes O(1/T) moments and O(1/T) squared norm."
            ),
            "constants_note": "An O(sqrt(x) log x) bound alone would not prove this conclusion.",
            "not_proved_in_repo": False,
            "adaptation_proved_in_repo": True,
            "proof": "docs/CONTINUOUS_THEOREM_A.md, section 2",
        },
        {
            "id": "ANT-MS",
            "name": "Positive mean-square limit for the exact psi residual",
            "role": "denominator",
            "status": "cited_exact_statement",
            "classical_refs": [
                "Akbary-Ng-Shahabi, arXiv:1306.1657v1, Corollary 1.15(i), p. 11, L=zeta",
            ],
            "source_url": source,
            "hypotheses": ["RH", "exact psi residual; no simplicity or linear independence assumption"],
            "adapted_conclusion": (
                "(1/T) Integral_0^T F_psi(y)^2 dy -> V="
                "2 sum_{gamma>0} m_gamma^2/(1/4+gamma^2), 0<V<infinity. "
                "The continuous affine-detrended theta residual has the same norm limit."
            ),
            "constants_note": "This supplies the positive denominator; it is not a guessed model norm floor.",
            "not_proved_in_repo": True,
            "adaptation_proved_in_repo": True,
            "proof": "docs/CONTINUOUS_THEOREM_A.md, section 3",
        },
        {
            "id": "ANT-4",
            "name": "Optional transfer to a weighted arithmetic residual",
            "role": "weight_transfer",
            "status": "not_established_for_this_extension",
            "classical_refs": ["M6 covers fixed absolutely continuous weights on model modes"],
            "hypotheses": [
                "A separate weighted numerator estimate and positive weighted denominator are required.",
            ],
            "adapted_conclusion": (
                "Not included in the continuous arithmetic theorem. Bounded weights "
                "alone do not ensure comparable weighted and unweighted norms."
            ),
            "constants_note": "Optional model tapers remain diagnostic tools.",
            "not_proved_in_repo": True,
            "adaptation_proved_in_repo": False,
            "proof": "docs/CONTINUOUS_THEOREM_A.md, scope",
        },
    ]


def full_a_gap_table() -> List[Dict[str, Any]]:
    """Separate the proved continuous theorem from optional and sampled transfers."""
    return [
        {
            "step": "M1-M7 diagnostic model lemmas",
            "disposition": "proved_model",
            "scope": "continuous projection lemmas with their stated hypotheses",
            "required_for_continuous_A": False,
            "support": "docs/PROOFS_LEMMAS.md",
        },
        {
            "step": "ANT-3 exact psi explicit formula",
            "disposition": "cited_exact_statement",
            "scope": "continuous, under RH",
            "required_for_continuous_A": True,
            "support": "Akbary-Ng-Shahabi equation (1.5)",
        },
        {
            "step": "ANT-1 uniform polynomial moments and ANT-2 prime-square correction",
            "disposition": "proved_continuous_adaptation",
            "scope": "continuous, under RH",
            "required_for_continuous_A": True,
            "support": "docs/CONTINUOUS_THEOREM_A.md sections 1-2",
        },
        {
            "step": "ANT-MS positive mean-square denominator",
            "disposition": "cited_exact_statement",
            "scope": "continuous, under RH",
            "required_for_continuous_A": True,
            "support": "Akbary-Ng-Shahabi Corollary 1.15(i); continuous adaptation section 3",
        },
        {
            "step": "Full A: RH implies R_d((I-P_1)Q(Tu))=O_d(T^-2)",
            "disposition": FULL_A_STATUS,
            "scope": "exact continuous theta residual, fixed d>=2, no added weight",
            "required_for_continuous_A": True,
            "support": "docs/CONTINUOUS_THEOREM_A.md",
        },
        {
            "step": "ANT-4 optional weighted arithmetic transfer",
            "disposition": "open_transfer",
            "scope": "optional extension beyond the proved theorem",
            "required_for_continuous_A": False,
            "support": "Requires weighted numerator and denominator estimates",
        },
        {
            "step": "Sampled grid and regression transfer",
            "disposition": "open_transfer",
            "scope": "sampled diagnostic, not the continuous theorem",
            "required_for_continuous_A": False,
            "support": "Requires quantitative quadrature and detrend errors as T grows",
        },
        {
            "step": "Unconditional RH",
            "disposition": "open",
            "scope": "arithmetic zeta",
            "required_for_continuous_A": False,
            "support": "The proved arithmetic decay theorem assumes RH",
        },
    ]


def full_b_gap_table() -> List[Dict[str, str]]:
    """Record distinct normalized and absolute converse statements."""
    return [
        {
            "step": "B0 / M2-M4 persistent orthogonal low-degree defect",
            "disposition": "proved_model",
            "support": "docs/PROOFS_LEMMAS.md M2-M4; this is not an arithmetic converse",
        },
        {
            "step": "Pure off-critical mode",
            "disposition": "phase_dependent_model",
            "support": "docs/OFF_CRITICAL_PHASE_OBSTRUCTION.md; empirical ratios are not a uniform positive floor",
        },
        {
            "step": "B-ABS absolute Legendre moment converse",
            "disposition": "proved_rh_equivalence",
            "support": "docs/ABSOLUTE_MOMENT_CONVERSE.md; equivalent absolute growth bound is unproved unconditionally",
        },
        {
            "step": "Signed absolute-moment oscillation from each off-line zero",
            "disposition": "proved_arithmetic",
            "support": "docs/SIGNED_MOMENT_OSCILLATION.md; either one-sided subexponential bound is RH-equivalent and open",
        },
        {
            "step": "Normalized converse with an attained rightmost off-line zero",
            "disposition": "proved_conditional",
            "support": "docs/ATTAINED_BOUNDARY_CONVERSE.md; limsup T*R_d>0 without a gap or finite dominant block",
        },
        {
            "step": "Normalized converse with an unattained rightmost real-part supremum",
            "disposition": "open",
            "support": "docs/ATTAINED_BOUNDARY_CONVERSE.md section 4; the energy little-o estimate does not close this case",
        },
        {
            "step": "Normalized recurrence: liminf R_d=0",
            "disposition": "proved_unconditional_continuous",
            "support": "docs/NORMALIZED_RECURRENCE_THEOREM.md; no full limit or rate",
        },
        {
            "step": "B-RES positive-floor formulation",
            "disposition": "rh_equivalent_open",
            "support": "By normalized recurrence, the off-zero implies positive-floor assertion is equivalent to RH",
        },
        {
            "step": "Full B: normalized fast decay implies RH",
            "disposition": FULL_B_STATUS,
            "support": "The attained-boundary case is proved; the unattained-boundary case remains open",
        },
        {
            "step": "Normalized converse for general positive discrete Beurling primes",
            "disposition": "refuted_by_discrete_counterexample",
            "support": "docs/BEURLING_NORMALIZED_COUNTEREXAMPLE.md; does not settle ordinary zeta",
        },
        {
            "step": "Unconditional RH",
            "disposition": "open",
            "support": "An RH-equivalent condition is not an unconditional proof of that condition",
        },
    ]


# ---------------------------------------------------------------------------
# Sampled off-critical comparison (not an arithmetic converse)
# ---------------------------------------------------------------------------


def off_critical_model_obstruction(
    T: float,
    *,
    sigma: float = 0.9,
    t: float = 14.134725,
    degree: int = 4,
    n_points: int = 4096,
) -> Dict[str, Any]:
    """
    Compare R_d of off-critical vs critical-line pure modes at the same T,t.

    Returns sampled R_d values and their ratio at this T. The real sinusoid
    has a phase-dependent leading term, so this does not establish monotone
    growth or a uniform positive lower envelope as T varies.
    """
    from .b_res_threshold import _scaled_mode

    T = float(T)
    if not np.isfinite(T) or T <= 0:
        raise ValueError("T must be positive and finite")
    if not (0.0 < float(sigma) < 1.0):
        raise ValueError("sigma in (0,1)")
    u = sample_grid(int(n_points))
    d = int(degree)
    q_off = _scaled_mode(u, T=T, sigma=float(sigma), t=t)
    q_cl = probe_critical_line_mode(u, T=T, t=t)
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
        "evidence_kind": "sampled_model_comparison",
        "uniform_lower_bound_proved": False,
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
    """Theorem registry; calling this function does not verify a proof or an RH hypothesis."""
    return {
        "banner": BANNER,
        "rh": RH_STATUS,
        "full_arithmetic_A": FULL_A_STATUS,
        "full_B": FULL_B_STATUS,
        "full_B_residual_step_id": B_RESIDUAL_STEP_ID,
        "full_B_residual_step": (
            "The positive-floor B-RES formulation is RH-equivalent and open; "
            "continuous normalized recurrence gives liminf R_d=0 unconditionally."
        ),
        "continuous_scope": dict(CONTINUOUS_SCOPE),
        "absolute_moment_converse": "proved_rh_equivalence",
        "absolute_moment_growth_bound": "open_unconditional",
        "signed_moment_oscillation": "proved_arithmetic",
        "one_sided_moment_bound": "rh_equivalent_open",
        "normalized_attained_boundary_converse": "proved_conditional",
        "normalized_unattained_boundary_case": "open",
        "beurling_normalized_converse": "refuted_by_discrete_counterexample",
        "normalized_recurrence": "proved_unconditional_continuous",
        "normalized_b_res_floor": "rh_equivalent_open",
        "normalized_full_limit": "open_unconditional",
        "normalized_rate": "proved_only_under_rh",
        "rh_claimed": False,
        "evidence_kind": "theorem_registry",
        "model_A0": "proved",
        "model_B0": "proved",
        "conditional_deduction_A": (
            "RH + exact explicit formula + positive mean-square theorem "
            "imply R_d((I-P_1)Q(Tu))=O_d(T^-2), fixed d>=2, continuous unweighted residual."
        ),
        "writeup_A": "docs/CONTINUOUS_THEOREM_A.md",
        "writeup_B": "docs/THEOREM_B_PACKAGE.md",
        "writeup_absolute_moments": "docs/ABSOLUTE_MOMENT_CONVERSE.md",
        "writeup_signed_moments": "docs/SIGNED_MOMENT_OSCILLATION.md",
        "writeup_attained_boundary": "docs/ATTAINED_BOUNDARY_CONVERSE.md",
        "writeup_normalized_recurrence": "docs/NORMALIZED_RECURRENCE_THEOREM.md",
        "roadmap": "docs/RH_CLOSEOUT_ROADMAP.md",
        "note": (
            "Continuous A is proved under RH. B-ABS is an equivalence and "
            "either one-sided moment bound remains unproved. Normalized recurrence "
            "proves only liminf zero. The normalized converse is proved when the "
            "rightmost off-line boundary is attained; its unattained case and RH remain open."
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

    Finite-grid model checks accompany the theorem registry. They neither
    verify the cited analytic inputs nor prove a continuous limit theorem.
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
        "evidence_kind": "sampled_model_report",
        "continuous_scope": dict(CONTINUOUS_SCOPE),
        "continuous_limit_verified_by_this_run": False,
        "rh_verified_by_this_run": False,
        "gap_table_A": full_a_gap_table(),
        "gap_table_B": full_b_gap_table(),
        "ant_ids": [c["id"] for c in ant_citations()],
        "banner": BANNER,
    }
