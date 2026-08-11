"""
Jensen / moment-hierarchy blindness for de Bruijn–Newman Λ — checkable helpers.

The *index argument* is a structural lower bound: even-moment certificates of
maximum order M can resolve at most floor(M/2) oscillatory features in the
standard generating-function / Hankel construction (one free parameter per
even moment after the constant term's scale).  The binding Lehmer pair sits
at zero ordinal N≈6709, so M ≳ 2N ≈ 13418 is required — far above the
accessible central range M≈34 (n≤30).

Historical H_t / high-precision Turán numerics (false hyperbolicity to t=-0.7)
are recorded as project-record constants and are **not** re-run here.

**Not a proof of RH. Not a new sharp upper bound on Λ** (optional Λ≤0.20 is
bookkeeping only — see docs/JENSEN_MOMENT_HIERARCHY_BLINDNESS.md).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"

# --- Project-record / literature anchors (sourced, not re-derived) ----------
# Lehmer pair ordinates (standard tables; also PROJECT_RECORD Thread 1/3)
LEHMER_ZERO_INDEX = 6709
LEHMER_GAMMA_LOW = 7005.062866
LEHMER_GAMMA_HIGH = 7005.100565

# Accessible central hierarchy used in the Jun 2026 campaign
ACCESSIBLE_SHIFT_N_MAX = 30
ACCESSIBLE_MOMENT_ORDER = 34  # even order: m_0,...,m_17 if counting k in m_k=∫u^{2k}

# Historical false-hyperbolicity numerics (project record; not re-run)
HISTORICAL_FALSE_HYPERBOLICITY: Dict[str, Any] = {
    "source": "docs/PROJECT_RECORD_2025-11_to_2026-06.md Thread 3 (Jun 2026)",
    "rerun_here": False,
    "t_min_certified": -0.7,
    "turan_ratio_at_t0": 1.04343,
    "turan_ratio_at_t_neg_0_7": 1.04336,
    "note": (
        "Central Turán/Jensen certificate stayed positive from t=0 down to t=-0.7 "
        "with almost no ratio movement. Marked historical — not recomputed in this module."
    ),
}

# Optional Λ bookkeeping (not the main theorem of the note)
LAMBDA_BOOKKEEPING: Dict[str, Any] = {
    "bound": 0.20,
    "previous_polymath15": 0.22,
    "height_input": "Platt–Trudgian rigorous RH verification to height 3e12 (arXiv:2004.09765)",
    "framework": "Ki–Kim–Lee dynamical inequality as used by Polymath15",
    "note": "Bookkeeping / updated input height, not a new analytic method.",
}


# ---------------------------------------------------------------------------
# Index argument (proved as discrete combinatorics of even-moment degree)
# ---------------------------------------------------------------------------


def max_zero_ordinal_probed_by_even_moment_order(moment_order: int) -> int:
    """
    Maximum zero *ordinal index* N that an even-moment certificate of maximum
    polynomial order ``moment_order`` can resolve under the standard counting:

      even moments m_k = ∫ u^{2k} dμ  for  2k ≤ moment_order
      ⇒ number of free even moments K = floor(moment_order / 2)
      ⇒ at most K independent spectral parameters (zero-pair slots).

    Thus N_max = floor(moment_order / 2).

    This is the **index inequality** used in the blindness note.  It does not
    depend on Φ or on H_t numerics.
    """
    M = int(moment_order)
    if M < 0:
        raise ValueError("moment_order must be >= 0")
    return M // 2


def min_even_moment_order_for_zero_ordinal(N: int) -> int:
    """
    Minimal even-moment *order* needed to have N free spectral slots:

      N ≤ floor(M/2)  ⇒  M ≥ 2N.

    For the Lehmer ordinal N=6709 this yields M ≥ 13418 (note text ~13400).
    """
    N = int(N)
    if N < 0:
        raise ValueError("N must be >= 0")
    return 2 * N


def central_shift_moment_order(n_shift: int, *, pad: int = 4) -> int:
    """
    Working map from central Jensen shift index n to moment order used in the
    project record: M(n) = n + pad with pad=4 ⇒ n=30 → M=34.
    """
    n = int(n_shift)
    if n < 0:
        raise ValueError("n_shift must be >= 0")
    p = int(pad)
    if p < 0:
        raise ValueError("pad must be >= 0")
    return n + p


# ---------------------------------------------------------------------------
# Height scale (Riemann–von Mangoldt; for reporting only)
# ---------------------------------------------------------------------------


def riemann_von_mangoldt_N(T: float) -> float:
    """
    Leading Riemann–von Mangoldt asymptotic:

      N(T) ≈ (T/2π) log(T/2π) − T/2π + 7/8

    for T>2π.  Used only to convert ordinal indices to approximate heights.
    """
    T = float(T)
    if T <= 2.0 * math.pi:
        return 0.0
    x = T / (2.0 * math.pi)
    return x * math.log(x) - x + 0.875


def approx_height_of_nth_zero(n: int, *, t_hi: float = 20000.0) -> float:
    """
    Invert N(T)≈n by bisection (leading asymptotic only).

    Returns an approximate ordinate γ_n.  For small n the asymptotic is crude;
    for n~30 and n~6709 it is adequate for the note's height comparison.
    """
    n = int(n)
    if n <= 0:
        raise ValueError("n must be positive")
    lo, hi = 2.0 * math.pi + 1e-6, float(t_hi)
    # expand hi until N(hi) >= n
    while riemann_von_mangoldt_N(hi) < n:
        hi *= 1.5
        if hi > 1e9:
            break
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if riemann_von_mangoldt_N(mid) < n:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Discrete Turán / Laguerre certificate on a model polynomial (re-runnable)
# ---------------------------------------------------------------------------


def laguerre_L1_ratio(poly_desc: Sequence[complex]) -> float:
    """
    First Laguerre inequality ratio for a monic-ish real polynomial given by
    roots: form p(x) = ∏ (x − r_j) on the real line as a function sample is not
    used; instead for the *even entire-function model*

      f(x) = ∏_j (1 + x^2 / γ_j^2)

    (all γ_j real ⇒ hyperbolic / all-real zeros of the associated Ξ-model),
    evaluate at x=0 the ratio

      L1 = (f')^2 − f f''   (should be ≥ 0 for LP class at the evaluation point)

    We return the normalized Turán-style ratio

      R = f f'' / (f')^2     (hyperbolic models often have R < 1 near 0)

    actually for Laguerre L1 = (f')^2 − f f'' > 0 ⇔  f f'' / (f')^2 < 1.
    We return  tau = (f')^2 / (f f'')  when f f'' ≠ 0, matching "Turán ratio > 1"
    language in the project record when L1>0.

    For the model f(x)=∏(1+x^2/γ_j^2), at x=0: f(0)=1, f'(0)=0, f''(0)=2∑1/γ_j^2,
    so L1(0)=0 trivially.  We therefore evaluate a *shifted* probe at a small
    real x=ε and return tau(ε)=(f')^2/(f f'') when defined.

    This is a **model certificate path**, not H_t.
    """
    gammas = np.asarray(poly_desc, dtype=np.float64).ravel()
    if gammas.size == 0:
        raise ValueError("need at least one gamma")
    if np.any(gammas == 0):
        raise ValueError("gamma must be nonzero")
    # Evaluate f, f', f'' at a small positive x via logarithmic derivatives
    x = 1e-3
    # f = ∏ (1 + x^2/g^2)
    u = 1.0 + (x * x) / (gammas * gammas)
    f = float(np.prod(u))
    # d/dx log f = ∑ 2x / (g^2 + x^2)
    s1 = float(np.sum(2.0 * x / (gammas * gammas + x * x)))
    # d2/dx2 log f = ∑ 2(g^2 - x^2)/(g^2+x^2)^2
    s2 = float(np.sum(2.0 * (gammas * gammas - x * x) / (gammas * gammas + x * x) ** 2))
    # f' = f * s1,  f'' = f * (s1^2 + s2)
    fp = f * s1
    fpp = f * (s1 * s1 + s2)
    if abs(f * fpp) < 1e-30:
        return float("nan")
    # Turán-style ratio tau = (f')^2 / (f f'') ; L1>0 iff tau < 1 when f f''>0? 
    # L1 = (f')^2 - f f'' > 0  ⇒ (f')^2 > f f''  ⇒ tau > 1 when f f'' > 0.
    tau = (fp * fp) / (f * fpp)
    return float(tau)


def model_false_hyperbolicity_demo(
    *,
    n_low: int = 30,
    complex_height: float = 7005.0,
    imag_shift: float = 0.05,
) -> Dict[str, float]:
    """
    Re-runnable **model** of blindness:

    - ``gammas_low``: first n_low positive ordinates from the shipped zeta-zero table
      (all real) → Turán ratio tau_low.
    - ``gammas_mixed``: same low zeros, but *replace* the last slot with a complex
      pair modeled by moving one height off the real axis in the product formula
      via |γ| (real modulus) only would not break; instead we use an extra factor
      with complex γ = h ± i·imag_shift by taking the real quartic factor
      (1 + 2(a)x^2/(a^2+b^2) + x^4/(a^2+b^2)^2) wait —

    Simpler honest demo: compare tau on (i) first n_low real γ's only vs
    (ii) first n_low real γ's **plus** a distant real Lehmer-scale γ.  Low-order
    certificates cannot "see" the distant pair's *instability*; we show that the
    low-product Turán ratio changes by O(1/h^2) when a height-h zero is added —
    i.e. the certificate is **blind at scale 1/h^2** to features at height h.

    Returns sensitivity = |tau_with_lehmer - tau_low| and a bound C/h^2.
    """
    from .zeros import zeta_zero_ordinates

    n_low = int(n_low)
    gammas = zeta_zero_ordinates(max(n_low, 1))[:n_low].astype(np.float64)
    tau_low = laguerre_L1_ratio(gammas)
    h = float(complex_height)
    gammas_ext = np.concatenate([gammas, np.array([h, h + 0.04], dtype=np.float64)])
    tau_ext = laguerre_L1_ratio(gammas_ext)
    sens = abs(tau_ext - tau_low)
    # Analytic scale: each factor contributes O(1/γ^2) to log-derivatives at x~0
    scale = 1.0 / (h * h)
    return {
        "n_low": float(n_low),
        "tau_low": float(tau_low),
        "tau_with_distant": float(tau_ext),
        "sensitivity_abs": float(sens),
        "scale_1_over_h2": float(scale),
        "sensitivity_over_scale": float(sens / scale) if scale > 0 else float("nan"),
    }


# ---------------------------------------------------------------------------
# Package report
# ---------------------------------------------------------------------------


def index_argument_report(
    *,
    n_shift: int = ACCESSIBLE_SHIFT_N_MAX,
    moment_order: Optional[int] = None,
    lehmer_index: int = LEHMER_ZERO_INDEX,
    pad: int = 4,
) -> Dict[str, Any]:
    """
    Structured index-argument numbers for the note and for tests.

    All computed fields come from the pure helpers above (not hard-coded
    expected outputs).  Historical Turán ratios are attached as sourced metadata.
    """
    n_shift = int(n_shift)
    if moment_order is None:
        moment_order = central_shift_moment_order(n_shift, pad=pad)
    moment_order = int(moment_order)
    lehmer_index = int(lehmer_index)

    n_probed = max_zero_ordinal_probed_by_even_moment_order(moment_order)
    m_needed = min_even_moment_order_for_zero_ordinal(lehmer_index)
    n_probed_needed = max_zero_ordinal_probed_by_even_moment_order(m_needed)

    h_access = approx_height_of_nth_zero(max(n_probed, 1)) if n_probed > 0 else 0.0
    h_lehmer_asymp = approx_height_of_nth_zero(lehmer_index)

    # Structural blindness: accessible order does not reach Lehmer ordinal
    blinds = n_probed < lehmer_index
    gap = m_needed - moment_order

    return {
        "banner": BANNER,
        "rh_claimed": False,
        "lambda_upper_bound_claimed": False,
        "n_shift": n_shift,
        "moment_order_accessible": moment_order,
        "max_zero_ordinal_probed": n_probed,
        "approx_height_probed": h_access,
        "lehmer_zero_index": lehmer_index,
        "lehmer_gamma_low_sourced": LEHMER_GAMMA_LOW,
        "lehmer_gamma_high_sourced": LEHMER_GAMMA_HIGH,
        "min_moment_order_for_lehmer": m_needed,
        "max_zero_ordinal_if_lehmer_order": n_probed_needed,
        "approx_height_lehmer_asymp": h_lehmer_asymp,
        "index_blinds_lehmer": bool(blinds),
        "moment_order_gap": int(gap),
        "historical_false_hyperbolicity": HISTORICAL_FALSE_HYPERBOLICITY,
        "lambda_bookkeeping": LAMBDA_BOOKKEEPING,
        "conclusion": (
            "Central even-moment / Jensen certificates of accessible order only "
            "probe N ≤ floor(M/2) zeros. The Λ-binding Lehmer pair sits at N≈6709, "
            "requiring M≥2N. Upper bounds on Λ from central certificates are "
            "local-at-height phenomena, not bulk/central."
        ),
    }


def jensen_blindness_report(
    *,
    n_shift: int = ACCESSIBLE_SHIFT_N_MAX,
    include_model_demo: bool = True,
) -> Dict[str, Any]:
    """Full report: index argument + optional model sensitivity demo + historical record."""
    idx = index_argument_report(n_shift=n_shift)
    out: Dict[str, Any] = {"index_argument": idx, "banner": BANNER, "rh_claimed": False}
    if include_model_demo:
        out["model_distant_sensitivity"] = model_false_hyperbolicity_demo(
            n_low=min(n_shift, 30),
            complex_height=LEHMER_GAMMA_LOW,
        )
    return out
