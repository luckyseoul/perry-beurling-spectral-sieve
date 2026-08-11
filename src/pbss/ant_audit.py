"""
Rank-3: bounded ANT interface audit for Full A citations.

Maps each cited ANT-1..3 input to the exact PBSS objects it must control.
Does **not** re-prove classical theorems; records a checkable interface checklist.

**Not RH.** Full A remains closed_conditional under these interfaces.
"""
from __future__ import annotations

from typing import Any, Dict, List

from .ab_closure import ant_citations, full_a_gap_table, package_status

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

# PBSS residual / weight / metric contracts (must match docs/THEOREM_A_PACKAGE.md)
PBSS_OBJECTS = {
    "residual": "q_T^arith = detrend_deg1((θ(e^{uT})-e^{uT})/√(e^{uT})) on u∈[0,1]",
    "window": "x = e^{uT}, Lebesgue du on [0,1]",
    "metric": "R_d = ||P_d q||_2^2 / ||q||_2^2, shifted Legendre V_d",
    "weight": "optional w ∈ W_α (Tukey/Hanning class, pbss.weights)",
    "mode_sum": "q_T^{(N)} from pbss.probes.explicit_formula_residual / finite CL sum",
}


def interface_checklist() -> List[Dict[str, Any]]:
    """
    Fixed audit rows: each ANT id → required PBSS match items + status.

    Status values:
      - matched: citation hypotheses cover the PBSS object as adapted in package
      - gap_named: missing uniformity / constant — recorded, not unlabeled
      - not_required: optional ANT-4
    """
    return [
        {
            "id": "ANT-3",
            "pbss_target": PBSS_OBJECTS["residual"] + " ↔ " + PBSS_OBJECTS["mode_sum"],
            "required_matches": [
                "classical EF for ψ or θ with C1 window map compatible with x=e^{uT}",
                "√x normalization + deg1 detrend absorbed into remainder or main terms",
                "identification error e_T with ||e_T||/||w q|| → 0",
            ],
            "uniformity_notes": [
                "Need T→∞ with N=N(T); constants may depend on d and w class",
            ],
            "status": "matched",
            "action": "stop_reproof",
        },
        {
            "id": "ANT-1",
            "pbss_target": "tail r_{N,T}^{tail} in L2([0,1], w^2 du) under RH",
            "required_matches": [
                "RH: Re ρ = 1/2 for all non-trivial zeros",
                "height cutoff G(T)→∞ so |γ|>G contributes o(1) after window map",
                "compatible with M5/M6 on the retained finite mode block",
            ],
            "uniformity_notes": [
                "Scaffold bound_infinite_zero_tail_scaffold is diagnostic only",
            ],
            "status": "matched",
            "action": "stop_reproof",
        },
        {
            "id": "ANT-2",
            "pbss_target": "arithmetic remainder r_T^{arith} (prime powers, contours, trivial zeros)",
            "required_matches": [
                "classical EF remainder after zero sum",
                "ψ−θ elementary / RH-scale bounds mapped to log-window L2",
            ],
            "uniformity_notes": [
                "Endpoint weight w∈W_α improves constants; does not replace ANT-2",
            ],
            "status": "matched",
            "action": "stop_reproof",
        },
        {
            "id": "ANT-4",
            "pbss_target": PBSS_OBJECTS["weight"],
            "required_matches": [
                "w∈W_α multiplies all terms in ANT-3 identity",
                "bulk ||w q|| ≍ ||q|| up to T-independent factors",
            ],
            "uniformity_notes": ["Optional if working unweighted"],
            "status": "not_required",
            "action": "optional",
        },
        {
            "id": "M7",
            "pbss_target": PBSS_OBJECTS["metric"] + " perturbation continuity",
            "required_matches": [
                "in-repo proved: R_d(q0+r) ≤ (√R0+δ)^2/(1-δ)^2",
            ],
            "uniformity_notes": [],
            "status": "matched",
            "action": "proved_in_repo",
        },
    ]


def ant_interface_audit() -> Dict[str, Any]:
    """
    Rank-3 entry: full audit report. Freeze Full A packaging after this.

    No unlabeled gaps: every ANT row has status matched|not_required|gap_named.
    """
    cites = {c["id"]: c for c in ant_citations()}
    rows = interface_checklist()
    # attach citation titles
    for row in rows:
        cid = row["id"]
        if cid in cites:
            row["citation_name"] = cites[cid]["name"]
            row["classical_refs"] = cites[cid]["classical_refs"]
            row["adapted_conclusion"] = cites[cid]["adapted_conclusion"]

    unlabeled = [
        r for r in rows if r["status"] not in ("matched", "not_required", "gap_named")
    ]
    gaps = [r for r in rows if r["status"] == "gap_named"]
    pkg = package_status()

    return {
        "banner": BANNER,
        "rh_claimed": False,
        "rank": 3,
        "title": "bounded ANT interface audit",
        "pbss_objects": PBSS_OBJECTS,
        "checklist": rows,
        "unlabeled_count": len(unlabeled),
        "gap_named_count": len(gaps),
        "full_a_status": pkg.get("full_arithmetic_A"),
        "freeze_full_a_packaging": len(unlabeled) == 0 and pkg.get("full_arithmetic_A")
        == "closed_conditional",
        "recommendation": (
            "Stop Full A packaging work. Do not re-prove classical EF theory unless "
            "a new quantitative theorem needs explicit constants."
        ),
        "gap_table_A": full_a_gap_table(),
    }
