"""
Rank-2: arithmetic plateau analysis via secondary explicit-formula-style terms.

Drives the **shipped** path ``ef_identify.multi_N_enrich_scan`` / ``identify_ef``.
Predeclared predictions are checked as live inequalities.

**Not RH. Not Full A.** Finite-T mechanism study for the soft R_d / Ed plateau.
"""
from __future__ import annotations

from typing import Any, Dict, List, Sequence

import numpy as np

from .ef_identify import multi_N_enrich_scan
from .probes import primes_upto

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

PREDICTIONS = {
    "P1_smooth_beats_zeros": (
        "zeros_smooth mean Ed_r/||q||^2 is strictly below zeros-only "
        "(secondary exp columns reduce low-degree mass)."
    ),
    "P2_smooth_not_full_kill": (
        "zeros_smooth leaves positive Ed (not a free V_d oracle)."
    ),
    "P3_Vd_is_oracle": (
        "zeros_Vd drives Ed near 0 (remaining mass after deg1 is in V_d)."
    ),
    "P4_flat_in_N_zeros": (
        "zeros-only Ed does not fall by more than 10% from min-N to max-N "
        "(more CL modes alone do not eat the plateau mass)."
    ),
}


def enrich_scan_small(
    *,
    T: float = 12.0,
    n_zeros_list: Sequence[int] = (5, 10, 20),
    degree: int = 4,
    x_max: float = 5e4,
    n_points: int = 1024,
    enrichments: Sequence[str] = ("zeros", "zeros_smooth", "zeros_Vd"),
) -> List[dict]:
    """Live multi-N enrichment scan via shipped multi_N_enrich_scan."""
    primes = primes_upto(float(x_max))
    return multi_N_enrich_scan(
        T=float(T),
        n_zeros_list=tuple(int(n) for n in n_zeros_list),
        primes=primes,
        enrichments=tuple(enrichments),
        n_points=int(n_points),
        degree=int(degree),
        detrend="deg1",
    )


def evaluate_predictions(rows: List[dict]) -> Dict[str, Any]:
    """Check predeclared predictions on live scan rows (Ed_r_over_l2q)."""
    by_en: Dict[str, List[float]] = {}
    by_en_N: Dict[str, Dict[int, float]] = {}
    for r in rows:
        en = str(r["m_enrich"])
        ed = float(r["Ed_r_over_l2q"])
        by_en.setdefault(en, []).append(ed)
        by_en_N.setdefault(en, {})[int(r["n_zeros"])] = ed

    means = {en: float(np.mean(v)) for en, v in by_en.items()}
    results: Dict[str, Any] = {"means": means, "predictions": {}}

    if "zeros" in means and "zeros_smooth" in means:
        results["predictions"]["P1_smooth_beats_zeros"] = {
            "pass": bool(means["zeros_smooth"] < means["zeros"] - 1e-6),
            "zeros_mean": means["zeros"],
            "smooth_mean": means["zeros_smooth"],
            "statement": PREDICTIONS["P1_smooth_beats_zeros"],
        }
    if "zeros_smooth" in means:
        results["predictions"]["P2_smooth_not_full_kill"] = {
            "pass": bool(means["zeros_smooth"] > 1e-4),
            "smooth_mean": means["zeros_smooth"],
            "statement": PREDICTIONS["P2_smooth_not_full_kill"],
        }
    if "zeros_Vd" in means:
        results["predictions"]["P3_Vd_is_oracle"] = {
            "pass": bool(means["zeros_Vd"] < 1e-3),
            "Vd_mean": means["zeros_Vd"],
            "statement": PREDICTIONS["P3_Vd_is_oracle"],
        }
    if "zeros" in by_en_N and len(by_en_N["zeros"]) >= 2:
        ns = sorted(by_en_N["zeros"])
        eds = [by_en_N["zeros"][n] for n in ns]
        results["predictions"]["P4_flat_in_N_zeros"] = {
            "pass": bool(not (eds[-1] < eds[0] * 0.9)),
            "N_list": ns,
            "Ed_list": eds,
            "statement": PREDICTIONS["P4_flat_in_N_zeros"],
        }

    preds = results["predictions"]
    results["all_tested_pass"] = bool(preds) and all(v["pass"] for v in preds.values())
    return results


def plateau_secondary_report(
    *,
    T: float = 14.0,
    x_max: float = 2e6,
    n_points: int = 1024,
    degree: int = 4,
) -> Dict[str, Any]:
    """Rank-2 entry point: secondary-term plateau analysis report.

    Default x_max=2e6 is large enough for P1 (smooth < zeros) on deg1 residual;
    tiny tables (x≲5e4) can invert P1 and should not be used as the default check.
    """
    rows = enrich_scan_small(
        T=T, x_max=x_max, n_points=n_points, degree=degree
    )
    pred = evaluate_predictions(rows)
    compact = [
        {
            "N": r["n_zeros"],
            "enrich": r["m_enrich"],
            "Ed_r_over_l2q": r["Ed_r_over_l2q"],
            "frac_l2": r["frac_l2"],
        }
        for r in rows
    ]
    return {
        "banner": BANNER,
        "rh_claimed": False,
        "rank": 2,
        "title": "arithmetic plateau via secondary EF-style terms",
        "T": float(T),
        "x_max": float(x_max),
        "degree": int(degree),
        "rows": compact,
        "prediction_eval": pred,
        "mechanism_summary": (
            "zeros-only leaves roughly flat low-degree Ed in N; zeros_smooth cuts Ed "
            "without free V_d; zeros_Vd kills Ed (oracle). Plateau mass is largely "
            "V_d-shaped secondary structure, not missing low critical-line zeros."
        ),
        "stop_conditions": (
            "Stop claiming a robust secondary mechanism if P1–P3 fail under "
            "T/x_max/degree changes; do not increase x alone without a transfer theorem."
        ),
    }
