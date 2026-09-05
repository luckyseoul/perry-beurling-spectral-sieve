"""Registry of analytic interfaces for the exact continuous Theorem A.

The registry records the proof and source dispositions; this function does
not execute a proof checker or certify sampled asymptotic convergence.
Historical generic ANT matches are replaced by explicit formula/moment/norm
inputs and separate unresolved transfer obligations. RH remains open.
"""
from __future__ import annotations

from typing import Any, Dict, List

from .ab_closure import FULL_A_STATUS, ant_citations, full_a_gap_table, package_status

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

PBSS_OBJECTS = {
    "residual": "q_T=(I-P_1)[exp(-Tu/2)*(theta(exp(Tu))-exp(Tu))]",
    "window": "u in [0,1], continuous Lebesgue du",
    "detrend": "continuous orthogonal affine projection P_1",
    "degree": "fixed integer d>=2",
    "metric": "R_d=||P_d q||_2^2/||q||_2^2, shifted Legendre V_d",
    "weight": "no added weight in the proved arithmetic theorem",
    "mode_sum": "Re sum[-2 m_gamma exp(i gamma y)/(1/2+i gamma)], with correct sign and phases",
    "sampled_residual": "probes.arithmetic_residual returns finite samples with sample least-squares detrend",
}


def interface_checklist() -> List[Dict[str, Any]]:
    """Registered proof obligations, including unresolved optional transfers."""
    return [
        {
            "id": "ANT-3",
            "pbss_target": "exact normalized psi formula before the explicit theta correction",
            "required_matches": ["RH", "ANS equation (1.5)", "negative coefficient sign and full complex phases"],
            "uniformity_notes": ["Choose one cutoff G=exp(T) on log(2)<=y<=T"],
            "status": "cited_exact_statement",
            "required_for_continuous_A": True,
            "action": "use_exact_source_statement",
        },
        {
            "id": "ANT-1",
            "pbss_target": "O(1/T) moments of the full psi residual",
            "required_matches": [
                "sum m_gamma/(gamma*abs(rho)) converges",
                "integrated explicit-formula error is O(1/T)",
            ],
            "uniformity_notes": ["The zero-block constant is independent of G; no arbitrary-N M5 inference"],
            "status": "proved_continuous_adaptation",
            "required_for_continuous_A": True,
            "action": "use_continuous_theorem_section_1",
        },
        {
            "id": "ANT-2",
            "pbss_target": "Q_theta=F_psi-1+h with h in L1 intersect L2",
            "required_matches": ["prime-square constant isolated", "continuous affine detrend removes -1"],
            "uniformity_notes": ["h(Tu) has O(1/T) moments and O(1/T) squared norm"],
            "status": "proved_continuous_adaptation",
            "required_for_continuous_A": True,
            "action": "use_continuous_theorem_section_2",
        },
        {
            "id": "ANT-MS",
            "pbss_target": "positive limiting norm of continuous q_T",
            "required_matches": ["ANS Corollary 1.15(i), L=zeta", "RH", "retain zero multiplicities"],
            "uniformity_notes": ["No simplicity or linear independence assumption; V>0"],
            "status": "cited_exact_statement",
            "required_for_continuous_A": True,
            "action": "use_mean_square_theorem_and_section_3",
        },
        {
            "id": "ANT-4",
            "pbss_target": "optional weighted arithmetic residual",
            "required_matches": ["weighted numerator bound", "positive weighted denominator"],
            "uniformity_notes": ["M6 proves model transfer only; a bounded taper does not prove norm comparability"],
            "status": "gap_named",
            "required_for_continuous_A": False,
            "action": "separate_optional_extension",
        },
        {
            "id": "SAMPLED",
            "pbss_target": PBSS_OBJECTS["sampled_residual"],
            "required_matches": ["grid refinement as T grows", "quadrature and detrend errors on the claimed scale"],
            "uniformity_notes": ["Fixed samples do not inherit the continuous limit theorem"],
            "status": "gap_named",
            "required_for_continuous_A": False,
            "action": "separate_sampled_transfer",
        },
        {
            "id": "M7",
            "pbss_target": "projection-ratio perturbation bound",
            "required_matches": ["R(q0+r)<=(sqrt(R(q0))+delta)^2/(1-delta)^2 when delta<1"],
            "uniformity_notes": ["A finite numerical check does not prove that delta tends to zero"],
            "status": "proved_model",
            "required_for_continuous_A": False,
            "action": "retain_model_lemma",
        },
    ]


def ant_interface_audit() -> Dict[str, Any]:
    """Summarize registered interfaces without upgrading models or open transfers."""
    cites = {c["id"]: c for c in ant_citations()}
    rows = interface_checklist()
    for row in rows:
        if row["id"] in cites:
            source = cites[row["id"]]
            row["citation_name"] = source["name"]
            row["classical_refs"] = source["classical_refs"]
            row["adapted_conclusion"] = source["adapted_conclusion"]
            row["source_url"] = source.get("source_url")
            row["proof"] = source["proof"]

    known = {"cited_exact_statement", "proved_continuous_adaptation", "proved_model", "gap_named"}
    closed = known - {"gap_named"}
    unlabeled = [r for r in rows if r["status"] not in known]
    gaps = [r for r in rows if r["status"] == "gap_named"]
    required_gaps = [r for r in rows if r["required_for_continuous_A"] and r["status"] not in closed]
    pkg = package_status()
    continuous_complete = not required_gaps and pkg["full_arithmetic_A"] == FULL_A_STATUS
    return {
        "banner": BANNER,
        "rh_claimed": False,
        "rank": 3,
        "title": "continuous theorem and sampled-interface audit",
        "evidence_kind": "theorem_and_source_registry",
        "analytic_inputs_verified_by_this_call": False,
        "sampled_transfer_proved": False,
        "weighted_transfer_proved": False,
        "pbss_objects": dict(PBSS_OBJECTS),
        "checklist": rows,
        "unlabeled_count": len(unlabeled),
        "gap_named_count": len(gaps),
        "required_gap_count": len(required_gaps),
        "full_a_status": pkg["full_arithmetic_A"],
        "freeze_continuous_theorem_packaging": continuous_complete,
        # Legacy broad freeze flag must not hide even a named unresolved transfer.
        "freeze_full_a_packaging": continuous_complete and not gaps and not unlabeled,
        "recommendation": (
            "Use the exact continuous theorem and its primary sources. Keep optional weighted "
            "and sampled transfers separate; finite model checks do not establish a T-limit."
        ),
        "gap_table_A": full_a_gap_table(),
    }
