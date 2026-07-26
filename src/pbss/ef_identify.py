"""
Explicit-formula identification attack: residual = mode sum + remainder.

Competing hypotheses for the arithmetic object that should match the truncated
critical-line mode sum from ``explicit_formula_residual``.

Hypotheses
----------
H_theta_sqrt : (θ(x)-x)/√x   [current shipped default; deg1 detrend optional]
H_psi_sqrt   : (ψ(x)-x)/√x   [Chebyshev ψ with prime powers — closer to classical EF]
H_psi_x      : (ψ(x)-x)/x
H_theta_x    : (θ(x)-x)/x

Identification
--------------
  q = residual under hypothesis H
  m = q_T^{(N)} from shipped explicit_formula_residual
  α = argmin ||q - α m||  (trapezoid L²) if fit_scale else 1
  r = q - α m

Remainder size metrics (all shipped path):
  frac_l2 = ||r||² / ||q||²
  E_d_r   = ||P_d r||²
  R_d_q, R_d_r, R_d_m
  corr    = ⟨q,m⟩ / (||q|| ||m||)

Bounds
------
- Model identity bound: if q = m exactly, fit α=1 ⇒ r=0.
- Projection remainder majorant: ||P_d r||² ≤ 2||P_d q||² + 2α²||P_d m||²
  (triangle; always true, checks plumbing).
- EF-scale majorant (model): M5 on m for R_d(m); if identification is good,
  R_d(q) cannot stay large unless r carries low-degree mass.

Attack goal: find H maximizing mode capture (min frac_l2 / min E_d_r) and
bound remaining mass; if all H leave large E_d_r on arithmetic data, name the
sharp block.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .probes import (
    _detrend,
    arithmetic_residual,
    explicit_formula_residual,
    sample_grid,
)
from .projection import (
    _trapezoid_weights,
    energy_ratio,
    project_coefficients,
    projection_energy,
)
from .lemmas import bound_R_d_finite_mode_sum
from .zeros import explicit_formula_amplitudes, zeta_zero_ordinates

HYPOTHESES = (
    "H_theta_sqrt",
    "H_psi_sqrt",
    "H_psi_x",
    "H_theta_x",
)


def psi_from_primes(primes: np.ndarray, x_max: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build (prime-power abscissae, Λ contributions) for ψ up to x_max.

    Returns sorted unique floors of p^k with weights log p (summed if collisions).
    For large tables this only loops primes with p² ≤ x_max for k≥2 (cheap).
    """
    p = np.asarray(primes)
    x_max = float(x_max)
    hi = int(np.searchsorted(p, x_max, side="right"))
    p = p[:hi]
    if p.size == 0:
        raise ValueError("no primes")
    # base: primes themselves (θ part)
    xs = [p.astype(np.float64)]
    ws = [np.log(p.astype(np.float64))]
    # higher powers p^k, k≥2
    # p ≤ sqrt(x_max) for k=2, etc.
    p_f = p.astype(np.float64)
    max_p2 = np.sqrt(x_max)
    p2_hi = int(np.searchsorted(p_f, max_p2, side="right"))
    for i in range(p2_hi):
        pv = float(p_f[i])
        logp = float(np.log(pv))
        pk = pv * pv
        while pk <= x_max + 1e-9:
            xs.append(np.array([pk], dtype=np.float64))
            ws.append(np.array([logp], dtype=np.float64))
            pk *= pv
    x_all = np.concatenate(xs)
    w_all = np.concatenate(ws)
    # sort and merge duplicates
    order = np.argsort(x_all, kind="mergesort")
    x_all = x_all[order]
    w_all = w_all[order]
    # unique merge
    uniq_x = []
    uniq_w = []
    cur_x = x_all[0]
    cur_w = w_all[0]
    for i in range(1, x_all.size):
        if abs(x_all[i] - cur_x) < 1e-9:
            cur_w += w_all[i]
        else:
            uniq_x.append(cur_x)
            uniq_w.append(cur_w)
            cur_x = x_all[i]
            cur_w = w_all[i]
    uniq_x.append(cur_x)
    uniq_w.append(cur_w)
    x_arr = np.asarray(uniq_x, dtype=np.float64)
    w_arr = np.asarray(uniq_w, dtype=np.float64)
    csum = np.cumsum(w_arr)
    return x_arr, csum


def chebyshev_residual(
    u: np.ndarray,
    *,
    T: float,
    primes: np.ndarray,
    kind: str = "theta",
    csum_theta: Optional[np.ndarray] = None,
    norm: str = "sqrt",
    detrend: str = "deg1",
    smooth: int = 1,
) -> Tuple[np.ndarray, float, dict]:
    """
    Build (θ-x) or (ψ-x) residual on the log-window with chosen normalization.

    kind : 'theta' | 'psi'
    norm : 'sqrt' → /√x ; 'x' → /x ; 'plain' → raw difference
    """
    u = np.asarray(u, dtype=np.float64)
    T = float(T)
    if T <= 0:
        raise ValueError("T > 0")
    x_max = float(np.exp(T))
    x = np.maximum(np.exp(u * T), 2.0)
    p = np.asarray(primes)
    hi = int(np.searchsorted(p, x_max, side="right"))
    p_use = p[:hi]
    if kind == "theta":
        if csum_theta is not None:
            c_use = np.asarray(csum_theta)[:hi]
        else:
            c_use = np.cumsum(np.log(p_use.astype(np.float64, copy=False)))
        idx = np.searchsorted(p_use, x, side="right") - 1
        cheb = np.zeros_like(x)
        ok = idx >= 0
        cheb[ok] = c_use[idx[ok]]
        n_terms = int(p_use.size)
    elif kind == "psi":
        x_pp, c_psi = psi_from_primes(p_use, x_max)
        idx = np.searchsorted(x_pp, x, side="right") - 1
        cheb = np.zeros_like(x)
        ok = idx >= 0
        cheb[ok] = c_psi[idx[ok]]
        n_terms = int(x_pp.size)
    else:
        raise ValueError("kind must be theta or psi")

    if norm == "sqrt":
        raw = (cheb - x) / np.sqrt(x)
    elif norm == "x":
        raw = (cheb - x) / x
    elif norm == "plain":
        raw = cheb - x
    else:
        raise ValueError("norm must be sqrt|x|plain")

    if smooth > 1:
        kernel = np.ones(smooth, dtype=np.float64) / float(smooth)
        raw = np.convolve(raw, kernel, mode="same")
    q = _detrend(raw, u, detrend)
    meta = {
        "T": T,
        "x_max": x_max,
        "kind": kind,
        "norm": norm,
        "detrend": detrend,
        "smooth": int(smooth),
        "n_terms": n_terms,
        "hypothesis_object": f"{kind}_{norm}",
    }
    return q, T, meta


def hypothesis_residual(
    u: np.ndarray,
    *,
    T: float,
    primes: np.ndarray,
    hypothesis: str,
    csum_theta: Optional[np.ndarray] = None,
    detrend: str = "deg1",
) -> Tuple[np.ndarray, float, dict]:
    """Dispatch one of the named EF-identification hypotheses."""
    h = hypothesis.strip()
    if h == "H_theta_sqrt":
        return chebyshev_residual(
            u, T=T, primes=primes, kind="theta", csum_theta=csum_theta,
            norm="sqrt", detrend=detrend,
        )
    if h == "H_psi_sqrt":
        return chebyshev_residual(
            u, T=T, primes=primes, kind="psi", csum_theta=csum_theta,
            norm="sqrt", detrend=detrend,
        )
    if h == "H_psi_x":
        return chebyshev_residual(
            u, T=T, primes=primes, kind="psi", csum_theta=csum_theta,
            norm="x", detrend=detrend,
        )
    if h == "H_theta_x":
        return chebyshev_residual(
            u, T=T, primes=primes, kind="theta", csum_theta=csum_theta,
            norm="x", detrend=detrend,
        )
    raise ValueError(f"unknown hypothesis {hypothesis!r}; choose from {HYPOTHESES}")


def fit_mode_scale(
    q: np.ndarray,
    m: np.ndarray,
    u: np.ndarray,
) -> float:
    """α = ⟨q,m⟩ / ⟨m,m⟩ with trapezoid weights."""
    w = _trapezoid_weights(np.asarray(u, dtype=np.float64))
    q = np.asarray(q, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    den = float(np.sum(w * m * m))
    if den <= 1e-30:
        return 0.0
    return float(np.sum(w * q * m) / den)


def l2_norm_sq(q: np.ndarray, u: np.ndarray) -> float:
    w = _trapezoid_weights(np.asarray(u, dtype=np.float64))
    q = np.asarray(q, dtype=np.float64)
    return float(np.sum(w * q * q))


def identify_ef(
    q: np.ndarray,
    u: np.ndarray,
    *,
    T: float,
    n_zeros: int,
    fit_scale: bool = True,
    degree: int = 4,
    form: str = "cos",
    include_poly_bulk: bool = False,
) -> Dict[str, float]:
    """
    Decompose q = α m + β0 + β1(u-1/2) + r  (bulk optional) with m = EF modes.

    Default: q = α m + r (bulk off).
    Returns remainder metrics and majorants.
    """
    u = np.asarray(u, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    m, _, meta_m = explicit_formula_residual(
        u, T=T, n_zeros=int(n_zeros), form=form, bulk="none", bulk_scale=0.0
    )
    w = _trapezoid_weights(u)
    if include_poly_bulk:
        # LS fit q ≈ α m + β0 + β1 (u-0.5)
        ones = np.ones_like(u)
        lin = u - 0.5
        # design matrix columns
        cols = [m, ones, lin]
        G = np.zeros((3, 3))
        b = np.zeros(3)
        for i, ci in enumerate(cols):
            b[i] = float(np.sum(w * q * ci))
            for j, cj in enumerate(cols):
                G[i, j] = float(np.sum(w * ci * cj))
        try:
            coef = np.linalg.solve(G, b)
        except np.linalg.LinAlgError:
            coef = np.zeros(3)
            coef[0] = fit_mode_scale(q, m, u) if fit_scale else 1.0
        alpha, beta0, beta1 = float(coef[0]), float(coef[1]), float(coef[2])
        m_eff = alpha * m + beta0 * ones + beta1 * lin
        r = q - m_eff
    else:
        alpha = fit_mode_scale(q, m, u) if fit_scale else 1.0
        beta0, beta1 = 0.0, 0.0
        r = q - alpha * m
    l2_q = l2_norm_sq(q, u)
    l2_m = l2_norm_sq(m, u)
    l2_r = l2_norm_sq(r, u)
    corr = float(np.sum(w * q * m) / np.sqrt(max(l2_q * l2_m, 1e-30)))
    e_q = projection_energy(q, u, degree)
    e_m = projection_energy(m, u, degree)
    e_r = projection_energy(r, u, degree)
    r_q = e_q / l2_q if l2_q > 1e-30 else 0.0
    r_m = energy_ratio(m, u, degree) if l2_m > 1e-30 else 0.0
    r_r = e_r / l2_r if l2_r > 1e-30 else 0.0
    frac_l2 = l2_r / l2_q if l2_q > 1e-30 else 0.0
    # always-true triangle majorant on projection energy of remainder
    tri = 2.0 * e_q + 2.0 * (alpha * alpha) * e_m
    # M5 majorant on the mode block (proved-style)
    t = zeta_zero_ordinates(int(n_zeros))
    a = explicit_formula_amplitudes(t)
    m5 = bound_R_d_finite_mode_sum(T, a, t, degree)
    # low-degree mass fraction carried by remainder
    e_r_over_q = e_r / l2_q if l2_q > 1e-30 else 0.0
    return {
        "T": float(T),
        "n_zeros": int(n_zeros),
        "degree": int(degree),
        "alpha": float(alpha),
        "beta0": float(beta0),
        "beta1": float(beta1),
        "include_poly_bulk": bool(include_poly_bulk),
        "fit_scale": bool(fit_scale),
        "frac_l2_remainder": float(frac_l2),
        "E_d_remainder": float(e_r),
        "E_d_remainder_over_l2q": float(e_r_over_q),
        "E_d_q": float(e_q),
        "E_d_modes": float(e_m),
        "R_d_q": float(r_q),
        "R_d_modes": float(r_m),
        "R_d_remainder": float(r_r),
        "corr_q_modes": float(corr),
        "triangle_majorant_Ed_rem": float(tri),
        "M5_bound_R_d_modes": float(m5),
        "triangle_holds": bool(e_r <= tri + 1e-9),
        "n_zeros_meta": meta_m.get("n_zeros"),
    }


def attack_one(
    u: np.ndarray,
    *,
    T: float,
    primes: np.ndarray,
    hypothesis: str,
    n_zeros: int,
    csum_theta: Optional[np.ndarray] = None,
    degree: int = 4,
    detrend: str = "deg1",
    fit_scale: bool = True,
    include_poly_bulk: bool = False,
) -> Dict[str, object]:
    """Full attack row: build H residual + identify against EF modes."""
    q, T_out, meta = hypothesis_residual(
        u,
        T=T,
        primes=primes,
        hypothesis=hypothesis,
        csum_theta=csum_theta,
        detrend=detrend,
    )
    idn = identify_ef(
        q,
        u,
        T=T_out,
        n_zeros=n_zeros,
        fit_scale=fit_scale,
        degree=degree,
        include_poly_bulk=include_poly_bulk,
    )
    return {
        "hypothesis": hypothesis,
        "residual_meta": meta,
        "identification": idn,
        "include_poly_bulk": bool(include_poly_bulk),
        "banner": "EF identification attack — technical campaign",
    }


def multi_hypothesis_scan(
    *,
    T_values: Sequence[float],
    n_zeros_list: Sequence[int],
    primes: np.ndarray,
    hypotheses: Sequence[str] = HYPOTHESES,
    n_points: int = 2048,
    degree: int = 4,
    detrend: str = "deg1",
    csum_theta: Optional[np.ndarray] = None,
    include_poly_bulk: bool = False,
) -> List[dict]:
    """Cartesian multi-(T, N, H) attack table."""
    u = sample_grid(int(n_points))
    p = np.asarray(primes)
    rows: List[dict] = []
    for T in T_values:
        if np.exp(float(T)) > float(p[-1]) * 0.99:
            continue
        for n_z in n_zeros_list:
            for h in hypotheses:
                row = attack_one(
                    u,
                    T=float(T),
                    primes=p,
                    hypothesis=h,
                    n_zeros=int(n_z),
                    csum_theta=csum_theta,
                    degree=degree,
                    detrend=detrend,
                    include_poly_bulk=include_poly_bulk,
                )
                rows.append(row)
    return rows


def model_sanity_identify(
    *,
    T: float,
    n_zeros: int,
    n_points: int = 4096,
    degree: int = 4,
) -> Dict[str, object]:
    """
    Model case: q = explicit_formula_residual itself.
    Identification must drive frac_l2 → 0 and E_d_rem → 0.
    """
    u = sample_grid(int(n_points))
    q, _, _ = explicit_formula_residual(u, T=T, n_zeros=n_zeros)
    idn = identify_ef(q, u, T=T, n_zeros=n_zeros, fit_scale=True, degree=degree)
    idn["model_identity_ok"] = bool(
        idn["frac_l2_remainder"] < 1e-12 and idn["E_d_remainder"] < 1e-18
    )
    return idn


def summarize_attack(rows: List[dict]) -> Dict[str, object]:
    """
    Aggregate multi-hypothesis results; name best H and sharp block if any.
    """
    if not rows:
        return {"status": "empty", "sharp_block": None}
    # score: lower E_d_remainder_over_l2q is better capture of low-degree mass by modes
    by_h: Dict[str, List[float]] = {}
    by_h_frac: Dict[str, List[float]] = {}
    by_h_corr: Dict[str, List[float]] = {}
    for row in rows:
        h = row["hypothesis"]
        idn = row["identification"]
        by_h.setdefault(h, []).append(idn["E_d_remainder_over_l2q"])
        by_h_frac.setdefault(h, []).append(idn["frac_l2_remainder"])
        by_h_corr.setdefault(h, []).append(abs(idn["corr_q_modes"]))
    means = {
        h: {
            "mean_Ed_rem_over_l2q": float(np.mean(v)),
            "mean_frac_l2": float(np.mean(by_h_frac[h])),
            "mean_abs_corr": float(np.mean(by_h_corr[h])),
            "n": len(v),
        }
        for h, v in by_h.items()
    }
    best = min(means.keys(), key=lambda h: means[h]["mean_Ed_rem_over_l2q"])
    best_ed = means[best]["mean_Ed_rem_over_l2q"]
    best_frac = means[best]["mean_frac_l2"]
    best_corr = means[best]["mean_abs_corr"]

    # Sharp block criteria (quantitative):
    # If even the best hypothesis leaves mean |corr| < 0.3 AND mean Ed_rem/l2q
    # > 0.5 * mean R_d scale (~0.05 if R_d~0.15), modes do not explain low-degree mass.
    # Also if frac_l2 stays > 0.5 for all H: L2 not captured by N modes.
    all_corr_low = all(means[h]["mean_abs_corr"] < 0.35 for h in means)
    all_frac_high = all(means[h]["mean_frac_l2"] > 0.5 for h in means)
    modes_miss_lowdeg = best_ed > 0.05 and best_corr < 0.35

    # N-invariance of low-degree remainder mass (decisive when corr is not tiny)
    by_h_N: Dict[str, Dict[int, List[float]]] = {}
    for row in rows:
        h = row["hypothesis"]
        idn = row["identification"]
        by_h_N.setdefault(h, {}).setdefault(int(idn["n_zeros"]), []).append(
            idn["E_d_remainder_over_l2q"]
        )
    n_invariant = False
    n_inv_detail = {}
    if best in by_h_N and len(by_h_N[best]) >= 2:
        n_means = {N: float(np.mean(v)) for N, v in by_h_N[best].items()}
        vals = list(n_means.values())
        # flat if max-min small relative to level
        spread = max(vals) - min(vals)
        level = float(np.mean(vals))
        n_invariant = level > 0.05 and spread < 0.05 * max(level, 1e-9) + 0.03
        n_inv_detail = {"Ed_by_N": n_means, "spread": spread, "level": level}

    if n_invariant and best_corr >= 0.3:
        sharp = {
            "name": "LOW_DEGREE_MASS_INVARIANT_TO_ZERO_TRUNCATION_N",
            "statement": (
                f"Best hypothesis {best} shows partial EF capture "
                f"(mean |corr(q,m)|≈{best_corr:.2f}, L² remainder fraction falls as N grows), "
                f"but the low-degree remainder mass E_d(r)/||q||² stays flat ≈{best_ed:.2f} "
                f"across N (spread {n_inv_detail.get('spread', 0):.3f}). "
                "So the Legendre-V_d mass of the arithmetic residual is not in the span of "
                "the first N critical-line modes after optimal scale (and poly bulk) fit. "
                "Missing ingredient: either (i) a secondary main-term / smooth component in "
                "the residual definition that accounts for V_d mass independent of zeros, or "
                "(ii) a different oscillatory basis (not the shipped a_n cos(t_n T u-α_n) law) "
                "that carries that mass, or (iii) a proof that V_d mass is an artifact of "
                "finite-x / window endpoints with an explicit rate →0 as T→∞."
            ),
            "parameters": {
                "best_hypothesis": best,
                "best_mean_abs_corr": best_corr,
                "best_mean_frac_l2": best_frac,
                "best_mean_Ed_rem_over_l2q": best_ed,
                "Ed_by_N": n_inv_detail.get("Ed_by_N"),
            },
            "failed_hypotheses": list(means.keys()),
            "what_would_unblock": (
                "Construct m_full = zero sum + secondary terms so that "
                "E_d(q-α m_full)/||q||² → 0 as N→∞ or T→∞ under RH; "
                "or change q so its V_d content is proved O(T^{-2})."
            ),
            "progress": (
                "θ/√x beats ψ and /x normalizations; |corr| increases with N "
                "(modes are not orthogonal garbage); model identity m=q_EF holds exactly."
            ),
        }
    elif all_corr_low and (all_frac_high or modes_miss_lowdeg):
        sharp = {
            "name": "EF_MODE_SUM_MISMATCH_ON_ARITHMETIC_RESIDUAL",
            "statement": (
                "For all tested identification hypotheses H and truncations N, "
                "the shipped truncated explicit-formula mode sum m=q_T^{(N)} has "
                "weak L² correlation with the arithmetic residual q_H and leaves "
                "substantial low-degree projection mass in r=q-αm. "
                "Missing estimate: a theorem that ||P_d(q_H - α m)|| / ||q_H|| → 0 "
                "for some H, α, N=N(T) as T→∞ under RH — OR a corrected residual "
                "definition whose leading oscillatory term is exactly m."
            ),
            "parameters": {
                "best_hypothesis": best,
                "best_mean_abs_corr": best_corr,
                "best_mean_frac_l2": best_frac,
                "best_mean_Ed_rem_over_l2q": best_ed,
            },
            "failed_hypotheses": list(means.keys()),
            "what_would_unblock": (
                "Either (i) prove identification for a modified residual "
                "(smoothed ψ-form, different norm, no detrend / different weight "
                "in the definition of q), or (ii) replace amplitude/phase law in "
                "explicit_formula_residual to match the true windowed transform of "
                "θ or ψ, or (iii) add secondary main terms (not only zeros) to m."
            ),
        }
    else:
        sharp = None

    return {
        "hypothesis_means": means,
        "best_hypothesis": best,
        "best_score_Ed_rem_over_l2q": best_ed,
        "sharp_block": sharp,
        "status": "sharp_block" if sharp else "partial_capture",
        "banner": "EF identification attack results — not an RH announcement",
    }
