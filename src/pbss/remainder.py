"""
Truncated explicit-formula / zero-sum remainder scaffolding for Theorem A.

Under RH the arithmetic residual is formally a sum over critical-line modes.
Lemma M5 proves R_d → 0 for any *finite* truncation. Full Theorem A needs
control of the **tail** (zeros beyond N) and arithmetic remainder terms.

This module ships:
  - truncated mode sum q^{(N)} (via probes.explicit_formula_residual)
  - peel-through-remainder path (q - α q^{(N)})
  - explicit majorant for the *mode tail* ∑_{n>N} under amplitude model a_n = 2/|ρ_n|
  - multi-(T,N) remainder size diagnostics

**Not a proof of RH or full arithmetic Theorem A.** Tail majorants use the
built-in finite zero table and a crude a_n ∼ 2/t_n model for hypothetical
further zeros — labeled as scaffolding bounds, not sharp AN T estimates.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .lemmas import bound_R_d_finite_mode_sum
from .probes import (
    explicit_formula_residual,
    peel_residual,
    sample_grid,
)
from .projection import energy_ratio, _trapezoid_weights
from .zeros import (
    ZETA_ZERO_ORDINATES_50,
    explicit_formula_amplitudes,
    zeta_zero_ordinates,
)


def truncated_mode_sum(
    u: np.ndarray,
    *,
    T: float,
    n_zeros: int,
    form: str = "cos",
    ordinates: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, dict]:
    """
    First-N explicit-formula mode sum q_T^{(N)} on the log window.

    Thin wrapper around shipped ``explicit_formula_residual`` (real entry path).
    """
    q, T_out, meta = explicit_formula_residual(
        u,
        T=T,
        n_zeros=int(n_zeros),
        ordinates=ordinates,
        form=form,
        bulk="none",
        bulk_scale=0.0,
    )
    meta = {
        **meta,
        "kind": "truncated_mode_sum",
        "note": "Finite CL truncation; not full arithmetic residual; not an RH proof.",
    }
    return q, meta


def peel_via_remainder(
    q_full: np.ndarray,
    u: np.ndarray,
    *,
    T: float,
    n_strip: int,
    form: str = "cos",
    mode_scale: float = 1.0,
    fit_scale: bool = False,
    ordinates: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, dict]:
    """
    Zero-sum peel driven through the truncation path:

      q_rem = q_full - α · q_T^{(N)}

    Uses shipped ``peel_residual`` / explicit-formula builder (not reimplemented).
    If ``fit_scale``, α is L² least-squares against the mode sum.
    """
    u = np.asarray(u, dtype=np.float64)
    q_full = np.asarray(q_full, dtype=np.float64)
    if n_strip < 0:
        raise ValueError("n_strip >= 0")
    if n_strip == 0:
        return q_full.copy(), {
            "n_strip": 0,
            "alpha": 0.0,
            "fit_scale": False,
            "kind": "peel_via_remainder",
            "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        }

    if fit_scale:
        q_modes, meta_m = truncated_mode_sum(
            u, T=T, n_zeros=n_strip, form=form, ordinates=ordinates
        )
        w = _trapezoid_weights(u)
        num = float(np.sum(w * q_full * q_modes))
        den = float(np.sum(w * q_modes * q_modes))
        alpha = num / den if den > 1e-30 else 0.0
        q_rem = q_full - alpha * q_modes
        return q_rem, {
            "n_strip": int(n_strip),
            "alpha": float(alpha),
            "fit_scale": True,
            "mode_meta": meta_m,
            "kind": "peel_via_remainder",
            "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        }

    q_rem, meta = peel_residual(
        q_full,
        u,
        T,
        n_strip,
        form=form,
        mode_scale=mode_scale,
        ordinates=ordinates,
    )
    meta = {
        **meta,
        "alpha": float(mode_scale),
        "fit_scale": False,
        "kind": "peel_via_remainder",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    return q_rem, meta


def tail_amplitude_majorant(
    n_start: int,
    n_end: int,
    *,
    t_min: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Scaffolding amplitudes for zeros with indices in [n_start, n_end)
    (1-based index of the zero).

    Uses the built-in table when available; beyond the table uses a crude
    model t_n ≈ t_table_last * (n / n_table) and a_n = 2/t_n.

    Returns (amplitudes, ordinates) for the tail segment.
    """
    n_start = int(n_start)
    n_end = int(n_end)
    if n_start < 1 or n_end <= n_start:
        raise ValueError("need 1 <= n_start < n_end")
    table = np.asarray(ZETA_ZERO_ORDINATES_50, dtype=np.float64)
    n_tab = int(table.size)
    ordinates: List[float] = []
    for n in range(n_start, n_end):
        if n <= n_tab:
            ordinates.append(float(table[n - 1]))
        else:
            # crude linear density extrapolation from last tabulated zero
            t_last = float(table[-1])
            ordinates.append(t_last * (n / n_tab))
    t = np.asarray(ordinates, dtype=np.float64)
    if t_min is not None:
        t = np.maximum(t, float(t_min))
    a = 2.0 / t  # |ρ|≈t model
    return a, t


def bound_R_d_mode_tail(
    T: float,
    *,
    n_kept: int,
    n_tail: int = 50,
    d: int = 4,
    l2_floor: Optional[float] = None,
) -> Dict[str, float]:
    """
    Explicit M5-style majorant for a *hypothetical* tail of ``n_tail`` modes
    after the first ``n_kept`` zeros.

    bound = bound_R_d_finite_mode_sum(T, a_tail, t_tail, d)

    This bounds only the **model tail** under a_n=2/t_n, not the full arithmetic
    explicit-formula remainder (prime powers, contour integrals, …).

    Labels: scaffolding majorant — not sharp, not RH.
    """
    T = float(T)
    if T <= 0:
        raise ValueError("T > 0")
    n_kept = int(n_kept)
    n_tail = int(n_tail)
    if n_kept < 0 or n_tail < 1:
        raise ValueError("n_kept >= 0 and n_tail >= 1")
    a, t = tail_amplitude_majorant(n_kept + 1, n_kept + 1 + n_tail)
    b = bound_R_d_finite_mode_sum(T, a, t, d, l2_floor=l2_floor)
    return {
        "T": T,
        "n_kept": n_kept,
        "n_tail": n_tail,
        "degree": int(d),
        "bound_R_d_tail": float(b),
        "sum_abs_a_over_tT": float(np.sum(np.abs(a) / (t * T))),
        "label": "scaffolding_majorant_not_sharp",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }


def remainder_diagnostic(
    u: np.ndarray,
    *,
    T: float,
    n_full: int,
    n_strip: int,
    degree: int = 4,
    form: str = "cos",
    fit_scale: bool = True,
    n_tail_bound: int = 40,
) -> Dict[str, object]:
    """
    Multi-purpose remainder diagnostic at one (T, N_strip):

    1. Build full truncated sum q^{(n_full)}
    2. Peel first n_strip modes → q_rem
    3. Report R_d(q_full), R_d(q_rem), M5 bound on kept modes, tail majorant

    For model residuals, if n_strip=n_full then q_rem≈0 and R_d→0.
    For arithmetic residuals, pass q_full from outside via
    ``remainder_diagnostic_from_q``.
    """
    if n_full < 1:
        raise ValueError("n_full >= 1")
    if n_strip < 0 or n_strip > n_full:
        raise ValueError("0 <= n_strip <= n_full")
    q_full, meta_full = truncated_mode_sum(u, T=T, n_zeros=n_full, form=form)
    r_full = float(energy_ratio(q_full, u, degree))
    q_rem, meta_peel = peel_via_remainder(
        q_full,
        u,
        T=T,
        n_strip=n_strip,
        form=form,
        fit_scale=fit_scale,
    )
    wts = _trapezoid_weights(u)
    l2_rem = float(np.sum(wts * q_rem * q_rem))
    if n_strip <= 0:
        r_rem = r_full
    elif l2_rem <= 1e-30:
        r_rem = 0.0  # full peel of model sum
    else:
        r_rem = float(energy_ratio(q_rem, u, degree))
    t = zeta_zero_ordinates(n_full)
    a = explicit_formula_amplitudes(t)
    # M5 bound on the *stripped* block (should decay in T)
    if n_strip > 0:
        b_kept = bound_R_d_finite_mode_sum(T, a[:n_strip], t[:n_strip], degree)
    else:
        b_kept = 0.0
    tail = bound_R_d_mode_tail(
        T, n_kept=n_strip, n_tail=n_tail_bound, d=degree
    )
    return {
        "T": float(T),
        "n_full": int(n_full),
        "n_strip": int(n_strip),
        "degree": int(degree),
        "R_d_full": r_full,
        "R_d_remainder": r_rem,
        "alpha": meta_peel.get("alpha", 0.0),
        "M5_bound_stripped_block": float(b_kept),
        "tail_majorant_R_d": tail["bound_R_d_tail"],
        "fit_scale": bool(fit_scale),
        "kind": "remainder_diagnostic_model",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": (
            "Model residual only. M5_bound is proved-style majorant for finite "
            "stripped block; tail_majorant is scaffolding for zeros beyond N."
        ),
    }


def remainder_diagnostic_from_q(
    q_full: np.ndarray,
    u: np.ndarray,
    *,
    T: float,
    n_strip: int,
    degree: int = 4,
    form: str = "cos",
    fit_scale: bool = True,
    n_tail_bound: int = 40,
) -> Dict[str, object]:
    """
    Same as ``remainder_diagnostic`` but for an externally supplied residual
    (e.g. arithmetic). Peels via shipped truncation path.
    """
    q_full = np.asarray(q_full, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    r_full = float(energy_ratio(q_full, u, degree))
    q_rem, meta_peel = peel_via_remainder(
        q_full,
        u,
        T=T,
        n_strip=n_strip,
        form=form,
        fit_scale=fit_scale,
    )
    l2_rem = float(np.sum(_trapezoid_weights(u) * q_rem * q_rem))
    r_rem = float(energy_ratio(q_rem, u, degree)) if l2_rem > 1e-30 else 0.0
    if n_strip > 0:
        t = zeta_zero_ordinates(n_strip)
        a = explicit_formula_amplitudes(t)
        b_kept = bound_R_d_finite_mode_sum(T, a, t, degree)
    else:
        b_kept = 0.0
    tail = bound_R_d_mode_tail(
        T, n_kept=n_strip, n_tail=n_tail_bound, d=degree
    )
    return {
        "T": float(T),
        "n_strip": int(n_strip),
        "degree": int(degree),
        "R_d_full": r_full,
        "R_d_remainder": r_rem,
        "alpha": meta_peel.get("alpha", 0.0),
        "M5_bound_stripped_block": float(b_kept),
        "tail_majorant_R_d": tail["bound_R_d_tail"],
        "fit_scale": bool(fit_scale),
        "kind": "remainder_diagnostic_external",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": (
            "External residual (e.g. arithmetic). Peeling model modes is a "
            "diagnostic; full Theorem A remains open."
        ),
    }


def multi_TN_remainder_scan(
    *,
    T_values: Sequence[float],
    n_full: int = 20,
    n_strips: Optional[Sequence[int]] = None,
    degree: int = 4,
    n_points: int = 4096,
    form: str = "cos",
) -> List[dict]:
    """
    Multi-(T, N_strip) remainder diagnostics on the model truncated sum.

    Returns a list of row dicts suitable for JSON export.
    """
    if n_strips is None:
        n_strips = [0, 1, 2, 5, 10, min(15, n_full), n_full]
    u = sample_grid(int(n_points))
    rows: List[dict] = []
    for T in T_values:
        for ns in n_strips:
            ns_i = int(ns)
            if ns_i > n_full:
                continue
            row = remainder_diagnostic(
                u,
                T=float(T),
                n_full=int(n_full),
                n_strip=ns_i,
                degree=degree,
                form=form,
            )
            rows.append(row)
    return rows
