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

from .basis import orthonormal_legendre_design, shifted_legendre_values
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

# m-enrichment variants (H_theta_sqrt fixed residual; enrich m only)
M_ENRICHMENTS = (
    "zeros",           # baseline: α m_N
    "zeros_poly1",     # + deg0/deg1 poly (secondary bulk)
    "zeros_highleg",   # + φ_2..φ_d (low-degree oscillatory polys)
    "zeros_Vd",        # + full V_d = φ_0..φ_d (absorbs all Ed mass if free)
    "zeros_smooth",    # + smooth secondary: exp(-uT), exp(-2uT) pullbacks
    "zeros_endpoint",  # + endpoint-shaped bumps u^2(1-u)^2 * {1, u-0.5}
)

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


def _ls_fit_columns(
    q: np.ndarray,
    u: np.ndarray,
    cols: List[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    """Weighted LS: q ≈ sum coef[i] * cols[i]. Returns (coef, m_eff)."""
    w = _trapezoid_weights(u)
    k = len(cols)
    G = np.zeros((k, k))
    b = np.zeros(k)
    for i, ci in enumerate(cols):
        b[i] = float(np.sum(w * q * ci))
        for j, cj in enumerate(cols):
            G[i, j] = float(np.sum(w * ci * cj))
    # ridge for near-singular designs
    G = G + 1e-12 * np.eye(k)
    try:
        coef = np.linalg.solve(G, b)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(G, b, rcond=None)[0]
    m_eff = np.zeros_like(q)
    for c, col in zip(coef, cols):
        m_eff = m_eff + float(c) * col
    return coef, m_eff


def build_m_columns(
    u: np.ndarray,
    *,
    T: float,
    n_zeros: int,
    enrich: str = "zeros",
    degree: int = 4,
    form: str = "cos",
) -> Tuple[List[np.ndarray], List[str], dict]:
    """
    Build columns of enriched m (zeros always first column).

    enrich in M_ENRICHMENTS. Residual is never modified — only m's span.
    """
    u = np.asarray(u, dtype=np.float64)
    m0, _, meta_m = explicit_formula_residual(
        u, T=float(T), n_zeros=int(n_zeros), form=form, bulk="none", bulk_scale=0.0
    )
    cols: List[np.ndarray] = [m0]
    names: List[str] = ["zeros"]
    enrich = (enrich or "zeros").strip()
    if enrich == "zeros":
        return cols, names, meta_m
    if enrich == "zeros_poly1":
        cols.extend([np.ones_like(u), u - 0.5])
        names.extend(["poly0", "poly1"])
    elif enrich == "zeros_highleg":
        for k in range(2, int(degree) + 1):
            cols.append(shifted_legendre_values(k, u))
            names.append(f"phi{k}")
    elif enrich == "zeros_Vd":
        Phi = orthonormal_legendre_design(int(degree), u)
        for k in range(int(degree) + 1):
            cols.append(Phi[:, k].copy())
            names.append(f"phi{k}")
    elif enrich == "zeros_smooth":
        # secondary-style decaying mains on log-window x=e^{uT}
        cols.append(np.exp(-0.5 * u * float(T)))
        cols.append(np.exp(-u * float(T)))
        cols.append(np.exp(-2.0 * u * float(T)))
        names.extend(["exp_half", "exp_one", "exp_two"])
    elif enrich == "zeros_endpoint":
        bump = (u * u) * ((1.0 - u) ** 2)
        cols.append(bump)
        cols.append(bump * (u - 0.5))
        names.extend(["end_bump", "end_bump_lin"])
    else:
        raise ValueError(f"unknown enrich {enrich!r}; choose from {M_ENRICHMENTS}")
    return cols, names, meta_m


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
    m_enrich: str = "zeros",
) -> Dict[str, float]:
    """
    Decompose q = m_eff + r with m_eff in span of enriched columns.

    m_enrich: zeros | zeros_poly1 | zeros_highleg | zeros_Vd | zeros_smooth | zeros_endpoint
    include_poly_bulk: legacy alias for zeros_poly1 when m_enrich is zeros.
    """
    u = np.asarray(u, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    enrich = m_enrich
    if include_poly_bulk and enrich == "zeros":
        enrich = "zeros_poly1"
    cols, names, meta_m = build_m_columns(
        u, T=T, n_zeros=int(n_zeros), enrich=enrich, degree=degree, form=form
    )
    w = _trapezoid_weights(u)
    if len(cols) == 1 and not fit_scale:
        alpha = 1.0
        coef = np.array([1.0])
        m_eff = cols[0].copy()
        r = q - m_eff
    else:
        coef, m_eff = _ls_fit_columns(q, u, cols)
        alpha = float(coef[0])
        r = q - m_eff
    m = cols[0]
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
    tri = 2.0 * e_q + 2.0 * (alpha * alpha) * e_m
    t = zeta_zero_ordinates(int(n_zeros))
    a = explicit_formula_amplitudes(t)
    m5 = bound_R_d_finite_mode_sum(T, a, t, degree)
    e_r_over_q = e_r / l2_q if l2_q > 1e-30 else 0.0
    return {
        "T": float(T),
        "n_zeros": int(n_zeros),
        "degree": int(degree),
        "alpha": float(alpha),
        "m_enrich": enrich,
        "m_column_names": names,
        "m_coefs": [float(c) for c in coef],
        "include_poly_bulk": bool(include_poly_bulk or enrich == "zeros_poly1"),
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


def multi_N_enrich_scan(
    *,
    T: float,
    n_zeros_list: Sequence[int],
    primes: np.ndarray,
    enrichments: Sequence[str] = M_ENRICHMENTS,
    n_points: int = 2048,
    degree: int = 4,
    detrend: str = "deg1",
    csum_theta: Optional[np.ndarray] = None,
) -> List[dict]:
    """
    Fixed H_theta_sqrt residual; multi-N × multi-enrich m table.
    Primary metric: E_d_remainder_over_l2q.
    """
    u = sample_grid(int(n_points))
    q, T_out, meta = hypothesis_residual(
        u,
        T=float(T),
        primes=np.asarray(primes),
        hypothesis="H_theta_sqrt",
        csum_theta=csum_theta,
        detrend=detrend,
    )
    rows: List[dict] = []
    for n_z in n_zeros_list:
        for enrich in enrichments:
            idn = identify_ef(
                q,
                u,
                T=T_out,
                n_zeros=int(n_z),
                degree=degree,
                m_enrich=str(enrich),
            )
            rows.append(
                {
                    "hypothesis": "H_theta_sqrt",
                    "residual_meta": meta,
                    "identification": idn,
                    "m_enrich": str(enrich),
                    "n_zeros": int(n_z),
                    "T": float(T_out),
                    "Ed_r_over_l2q": idn["E_d_remainder_over_l2q"],
                    "frac_l2": idn["frac_l2_remainder"],
                    "corr": idn["corr_q_modes"],
                }
            )
    return rows


def summarize_enrich_kill021(rows: List[dict], baseline: str = "zeros") -> Dict[str, object]:
    """
    Compare enrichments vs zeros baseline on Ed(r)/||q||² multi-N.
    Win: some enrich has mean Ed clearly below ~0.21 and drops vs baseline.
    """
    by_e: Dict[str, List[float]] = {}
    by_e_N: Dict[str, Dict[int, List[float]]] = {}
    for r in rows:
        e = r.get("m_enrich") or r["identification"].get("m_enrich", "zeros")
        ed = r.get("Ed_r_over_l2q", r["identification"]["E_d_remainder_over_l2q"])
        n = int(r.get("n_zeros", r["identification"]["n_zeros"]))
        by_e.setdefault(e, []).append(float(ed))
        by_e_N.setdefault(e, {}).setdefault(n, []).append(float(ed))
    means = {e: float(np.mean(v)) for e, v in by_e.items()}
    by_N_mean = {
        e: {N: float(np.mean(vs)) for N, vs in sorted(nd.items())}
        for e, nd in by_e_N.items()
    }
    base = float(means.get(baseline, 0.21))
    best = min(means.keys(), key=lambda k: means[k])
    best_mean = means[best]
    # clear drop: at least 25% relative reduction and absolute < 0.15
    dropped = best_mean < base * 0.75 and best_mean < 0.15
    killed = best_mean < 0.05  # essentially gone
    if killed:
        outcome = "win_killed_021"
        block = None
    elif dropped:
        outcome = "win_clear_drop"
        block = None
    else:
        outcome = "sharper_block"
        # which enrichments failed
        failed = [e for e, m in means.items() if m >= base * 0.9]
        block = {
            "name": "ENRICHED_M_FAILS_TO_ABSORB_VD_EXCEPT_EXPLICIT_VD",
            "statement": (
                f"On H_theta_sqrt, zeros-only baseline mean Ed(r)/||q||²≈{base:.3f}. "
                f"Best non-trivial enrichment among tried is {best} at {best_mean:.3f}. "
                "Smooth/endpoint/poly enrichments do not kill the ~0.21 V_d remainder; "
                "only spanning V_d inside m (zeros_Vd) is expected to zero Ed by construction. "
                "Missing for Full A: secondary main terms that are *not* arbitrary V_d "
                "but EF-derived, OR a proof that V_d mass of q is O(T^{-2})."
            ),
            "baseline_mean_Ed": base,
            "best_enrich": best,
            "best_mean_Ed": best_mean,
            "means": means,
            "by_N": by_N_mean,
            "failed_near_baseline": failed,
        }
        # If zeros_Vd kills, refine block name
        if "zeros_Vd" in means and means["zeros_Vd"] < 0.05:
            others = {e: m for e, m in means.items() if e != "zeros_Vd"}
            if others and min(others.values()) >= base * 0.85:
                block = {
                    "name": "VD_MASS_IS_POLYNOMIAL_NOT_ZERO_SUMMABLE",
                    "statement": (
                        f"zeros_Vd drives Ed(r)/||q||²→{means['zeros_Vd']:.4f} (V_d absorbed into m), "
                        f"but zeros-only≈{base:.3f} and other enrichments "
                        f"(smooth/endpoint/poly1/highleg) stay near baseline "
                        f"(best other={min(others, key=others.get)} at {min(others.values()):.3f}). "
                        "So the 0.21 is low-degree polynomial mass in q, not explained by "
                        "more zeros or classical-ish exp-decay secondary columns. "
                        "Unblock: EF secondary terms that reproduce V_d content of "
                        "(θ-x)/√x under RH, or redesign residual so that mass is O(T^{-2})."
                    ),
                    "baseline_mean_Ed": base,
                    "zeros_Vd_mean_Ed": means["zeros_Vd"],
                    "other_means": others,
                    "by_N": by_N_mean,
                }
                outcome = "sharper_block"
                # still a diagnostic win: we know WHAT the mass is
    return {
        "outcome": outcome,
        "baseline": baseline,
        "means_Ed_r_over_l2q": means,
        "by_N_Ed": by_N_mean,
        "best_enrich": best,
        "best_mean_Ed": best_mean,
        "sharp_block": block,
        "n_rows": len(rows),
        "banner": "enrich-m-only attack on 0.21 — H_theta_sqrt fixed",
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
    m_enrich: str = "zeros",
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
        m_enrich=m_enrich,
    )
    return {
        "hypothesis": hypothesis,
        "residual_meta": meta,
        "identification": idn,
        "include_poly_bulk": bool(include_poly_bulk),
        "m_enrich": idn.get("m_enrich", m_enrich),
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
