"""
Rank-4: narrow feasibility check for Anthropic-style zero-proportion methods.

Compares theorem *classes* only. Ships **no** new Weil/BGST inequality.
Stop condition: if no incremental inequality is proposed, record STOP.

**Not RH. Not a claim that 41.6%→67.2% is reproduced here.**
"""
from __future__ import annotations

from typing import Any, Dict, List

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"


def theorem_class_comparison() -> List[Dict[str, str]]:
    """Structural comparison: Anthropic proportion bound vs PBSS residual diagnostic."""
    return [
        {
            "axis": "object",
            "anthropic": "lower bound on fraction of zeros with Re=1/2",
            "pbss": "R_d of density residual on log-window / Beurling systems",
        },
        {
            "axis": "tools",
            "anthropic": "Weil quadratic forms; BGST; Bombieri; rank inequalities",
            "pbss": "shifted Legendre projection; M1–M7; EF peel; Beurling battery",
        },
        {
            "axis": "unconditional progress type",
            "anthropic": "partial zero-count proportion (41.6%→67.2% claimed upstream)",
            "pbss": "model decay + conditional Full A + limit theorems (Jensen-blindness)",
        },
        {
            "axis": "implies RH?",
            "anthropic": "no (explicit)",
            "pbss": "no (explicit)",
        },
        {
            "axis": "overlap with B-RES",
            "anthropic": "none direct — more on-line zeros ≠ residual injection theorem",
            "pbss": "B-RES is residual converse; different class",
        },
    ]


def incremental_inequality_candidates() -> List[Dict[str, str]]:
    """
    Only candidates that would justify *continuing* this rank.
    None are implemented — feasibility says STOP unless one is adopted later.
    """
    return [
        {
            "idea": "PBSS R_d bound that implies a positive on-line zero density",
            "status": "not_formulated",
            "why_hard": "R_d is a residual energy ratio, not a zero counter",
        },
        {
            "idea": "Import BGST pair correlation without RH into a residual weight",
            "status": "not_formulated",
            "why_hard": "no natural map from quadratic-form ranks to V_d energy",
        },
    ]


def zero_proportion_feasibility_report() -> Dict[str, Any]:
    """Rank-4 entry: feasibility report with explicit STOP."""
    cands = incremental_inequality_candidates()
    any_ready = any(c["status"] == "ready_to_implement" for c in cands)
    return {
        "banner": BANNER,
        "rh_claimed": False,
        "rank": 4,
        "title": "Anthropic-style zero-proportion feasibility",
        "upstream_ref": "docs/related/anthropic-riemann-zeta/",
        "upstream_blog": "https://www.anthropic.com/research/riemann-zeta",
        "comparison": theorem_class_comparison(),
        "incremental_candidates": cands,
        "decision": "STOP" if not any_ready else "CONTINUE",
        "rationale": (
            "No incremental inequality tying Weil/BGST rank methods to PBSS R_d is "
            "formulated. Reproducing 41.6%→67.2% would be a separate research program, "
            "not a continuation of residual diagnostics. Park this rank."
        ),
        "do_not": [
            "rebuild Spectral Sieve operators",
            "claim Anthropic methods close B-RES",
            "run large zero-count campaigns without a new theorem statement",
        ],
    }
