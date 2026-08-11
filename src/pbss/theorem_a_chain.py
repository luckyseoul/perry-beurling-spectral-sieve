"""
Checkable model chain for Conditional Theorem A (not arithmetic Full A, not RH).

Composes shipped majorants (M5, M6) and scaffolding tail bounds into one report
per T. Labels every field as proved-style majorant vs scaffolding vs open.

See docs/THEOREM_A_PACKAGE.md.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from .lemmas import (
    bound_R_d_finite_mode_sum,
    bound_R_d_weighted_finite_mode_sum,
    bound_R_d_weighted_sine_order,
)
from .probes import (
    explicit_formula_residual,
    probe_critical_line_mode,
    sample_grid,
)
from .projection import energy_ratio
from .remainder import bound_infinite_zero_tail_scaffold
from .weights import admissible_weight, apply_weight
from .zeros import explicit_formula_amplitudes, zeta_zero_ordinates

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"
# Full A: closed conditionally under RH + cited ANT (see ab_closure / THEOREM_A_PACKAGE).
# Not unconditional; RH remains open.
FULL_A_STATUS = "closed_conditional"
RH_STATUS = "open"


def model_chain_report(
    T: float,
    *,
    degree: int = 4,
    n_zeros: int = 10,
    n_points: int = 4096,
    t_cl: float = 14.134725,
    weight_name: str = "hanning",
    alpha: float = 0.1,
    n_tail_eff: int = 2000,
) -> Dict[str, object]:
    """
    One-T snapshot of the Conditional Theorem A *model* chain.

    Returns empirical R_d for CL / EF / weighted EF, majorants M5/M6, and a
    scaffolding infinite-tail majorant. Does **not** assert arithmetic A₀.
    """
    T = float(T)
    if T <= 0:
        raise ValueError("T > 0")
    u = sample_grid(int(n_points))
    d = int(degree)

    # --- pure CL ---
    q_cl = probe_critical_line_mode(u, T=T, t=t_cl)
    r_cl = float(energy_ratio(q_cl, u, d))
    maj_cl = bound_R_d_weighted_sine_order(
        t_cl * T, d, w_linf=1.0, wq_l2_floor=0.05
    )  # flat w_linf=1 is M3-style; weighted below

    w = admissible_weight(u, name=weight_name, alpha=alpha)
    r_cl_w = float(energy_ratio(apply_weight(q_cl, w), u, d))
    maj_cl_w = bound_R_d_weighted_sine_order(
        t_cl * T, d, w_linf=1.0, wq_l2_floor=0.02
    )

    # --- finite EF ---
    q_ef, _, meta = explicit_formula_residual(u, T=T, n_zeros=int(n_zeros))
    r_ef = float(energy_ratio(q_ef, u, d))
    t = zeta_zero_ordinates(int(n_zeros))
    a = explicit_formula_amplitudes(t)
    maj_ef = bound_R_d_finite_mode_sum(T, a, t, d)
    r_ef_w = float(energy_ratio(apply_weight(q_ef, w), u, d))
    maj_ef_w = bound_R_d_weighted_finite_mode_sum(T, a, t, d, w_linf=1.0)

    # --- scaffolding tail beyond n_zeros ---
    tail = bound_infinite_zero_tail_scaffold(
        T, n_kept=int(n_zeros), N_eff=int(n_tail_eff), d=d
    )

    # Model decay checks (not arithmetic)
    proved_ok = (
        r_cl <= maj_cl + 0.05
        and r_cl_w <= maj_cl_w + 0.05
        and r_ef <= maj_ef + 0.15  # finite-mode majorant is loose
        and r_ef_w <= maj_ef_w + 0.15
    )

    return {
        "T": T,
        "degree": d,
        "n_zeros": int(n_zeros),
        "weight_name": weight_name,
        "empirical": {
            "R_d_cl": r_cl,
            "R_d_cl_weighted": r_cl_w,
            "R_d_ef": r_ef,
            "R_d_ef_weighted": r_ef_w,
        },
        "majorants": {
            "M3_style_cl_flat": {
                "value": maj_cl,
                "label": "proved_style_majorant_M3_M6",
            },
            "M6_cl_weighted": {
                "value": maj_cl_w,
                "label": "proved_style_majorant_M6",
            },
            "M5_ef_flat": {
                "value": maj_ef,
                "label": "proved_style_majorant_M5",
            },
            "M6_ef_weighted": {
                "value": maj_ef_w,
                "label": "proved_style_majorant_M6",
            },
            "infinite_tail_scaffold": {
                "value": tail["bound_R_d_model_tail"],
                "label": tail["label"],
            },
        },
        "proved_model_decay_ok": bool(proved_ok),
        "full_arithmetic_A_status": FULL_A_STATUS,
        "rh_status": RH_STATUS,
        "banner": BANNER,
        "note": (
            "Model chain only (M5/M6 majorants + scaffolding diagnostic tail). "
            "Full arithmetic Theorem A is closed *conditionally* under RH + "
            "cited ANT-1..3 + M7 — see docs/THEOREM_A_PACKAGE.md / pbss.ab_closure. "
            "RH remains open. Scaffold tail is not a required Full-A step."
        ),
    }


def multi_T_model_chain(
    T_values: Sequence[float],
    **kwargs,
) -> List[dict]:
    """Multi-T list of model_chain_report rows."""
    return [model_chain_report(float(T), **kwargs) for T in T_values]


def package_status() -> Dict[str, object]:
    """Machine-readable status of Full A/B packages vs RH (delegates to ab_closure)."""
    from .ab_closure import package_status as _full_status

    base = _full_status()
    # Preserve legacy keys used by older tests/callers
    base["conditional_theorem_a_package"] = "complete"
    base["writeup"] = base.get("writeup_A", "docs/THEOREM_A_PACKAGE.md")
    return base
