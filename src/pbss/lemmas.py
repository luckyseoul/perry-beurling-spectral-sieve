"""
Proved analytic facts about the PBSS diagnostic (continuous L²[0,1] model).

These are theorems about the *projection diagnostic*, not claims that RH is true.
Discrete trapezoid implementations approximate them; unit tests check the shipped
``energy_ratio`` path against these identities within quadrature tolerance.

Lemma M1 (pure modes)
---------------------
If {φ_k} is orthonormal on [0,1] and q = φ_m, then
  R_d(q) = 1  if m ≤ d,   and   R_d(q) = 0  if m > d.

Lemma M2 (orthogonal defect formula)
------------------------------------
If j ≤ d, ε ∈ [0,1], f ⊥ V_d := span{φ_0,…,φ_d}, ‖f‖₂=1, and
  q = √(1-ε²) f + ε φ_j ,
then R_d(q) = ε².

Lemma M3 (critical-line pure mode decay)
----------------------------------------
Let q_ω(u) = sin(ω u) for ω > 0. Each coefficient c_k = ⟨q_ω, φ_k⟩ satisfies
  |c_k| ≤ C_k / ω
because φ_k is C¹ (in fact a polynomial) and
  ∫_0^1 sin(ωu) φ_k(u) du = O(ω^{-1})
by integration by parts. Hence with ‖q_ω‖₂² → 1/2,
  R_d(q_ω) = O_d(ω^{-2}) as ω → ∞.
In particular for the critical-line model q_T(u)=sin(t T u) (ω = t T),
  R_d(q_T) = O(T^{-2}).

(The archive's sharper heuristic R_d = O(T^{-2(d+1)}) is *not* claimed as proved here.)

Lemma M4 (persistent defect blocks vanishing)
---------------------------------------------
Under the hypotheses of M2 with ε fixed > 0, R_d(q) = ε² does not tend to 0
as any auxiliary parameter (e.g. frequency of f) varies. Thus R_d → 0 along a
family forces the low-degree mass ε → 0 for that family.

Lemma M5 (finite critical-line superposition — finite-mode A₀)
--------------------------------------------------------------
Let N < ∞, amplitudes a_n ∈ ℝ, distinct ordinates t_n > 0, phases φ_n ∈ ℝ, and
  q_T(u) = ∑_{n=1}^N a_n sin(t_n T u + φ_n)   (or the cosine form).
Integration by parts on each mode gives |⟨q_T, φ_k⟩| ≤ C_k ∑_n |a_n|/(t_n T).
Hence ‖P_d q_T‖₂² = O_{d,N,{a,t}}(T^{-2}). For large T, ‖q_T‖₂² is bounded
below by a positive constant depending only on {a_n} (diagonal sine terms →
½∑a_n²; cross terms O(T^{-1})), so when a ≢ 0,
  R_d(q_T) = O_d(T^{-2}) as T → ∞
at the **same order** as pure-mode M3. This is the finite-mode extension of A₀
for truncated explicit-formula residuals. It does **not** prove full Theorem A
for the arithmetic prime residual or RH.

Lemma M6 (admissible weight preserves model-mode decay)
-------------------------------------------------------
Let w be fixed and absolutely continuous with ‖w‖_∞ ≤ W,
‖w'‖₁ ≤ V < ∞, and ‖w q‖₂² ≥ c_w > 0 for the residual
family under consideration (true for Hanning/Tukey on pure CL modes at large T).
If q_ω(u)=sin(ω u), then
  |⟨w q_ω, φ_k⟩| ≤ (|(wφ_k)(0)| + |(wφ_k)(1)| + ‖(wφ_k)'‖₁) / ω
by integration by parts on wφ_k. The constant depends on W and V; boundedness
alone does not imply this rate. Hence R_d(w q_ω) = O_{d,w}(ω^{-2}).
The same argument with a finite CL sum gives weighted finite-mode A₀:
  R_d(w q_T^{(N)}) = O_d(T^{-2}).
This is a model lemma about the diagnostic under fixed admissible weights.
It does **not** by itself prove R_d(w q_T^{arith})→0.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

from .basis import shifted_legendre_values
from .basis import _validate_degree
from .probes import probe_high_frequency, sample_grid
from .projection import energy_ratio


def continuous_R_d_pure_mode(m: int, d: int) -> float:
    """Exact continuous value from Lemma M1."""
    if m < 0 or d < 0:
        raise ValueError("m,d >= 0")
    return 1.0 if m <= d else 0.0


def continuous_R_d_orthogonal_defect(eps: float) -> float:
    """Exact continuous value from Lemma M2 (any j ≤ d)."""
    if not 0.0 <= eps <= 1.0:
        raise ValueError("eps in [0,1]")
    return float(eps * eps)


def synthetic_orthogonal_defect(
    u: np.ndarray,
    eps: float,
    j: int = 0,
    waves: int = 60,
) -> np.ndarray:
    """
    Build a discrete surrogate for Lemma M2:
      q = √(1-ε²) f + ε φ_j
    with f a high-frequency sinusoid projected orthogonal to V_d is expensive;
    we use a pure high-frequency mode and then *remove* its empirical projection
    onto φ_0..φ_j by Gram-Schmidt in the discrete inner product — simpler path:

    For tests of M2 we construct f already nearly orthogonal by using high waves
    and then form the mixture; the identity R≈ε² holds asymptotically as waves→∞
    and n→∞. For a sharp test, we build f by taking high-frequency and subtracting
    its components in V_d via the shipped design matrix.
    """
    from .basis import orthonormal_legendre_design
    from .projection import _trapezoid_weights, project_coefficients

    u = np.asarray(u, dtype=np.float64)
    if not 0.0 <= eps <= 1.0:
        raise ValueError("eps in [0,1]")
    # start with HF, remove discrete components in degrees 0..max(j,4) so residual
    # is nearly orthogonal to V_d for d>=j used in tests
    d_kill = max(j, 4)
    f0 = probe_high_frequency(u, waves=waves)
    w = _trapezoid_weights(u)
    c = project_coefficients(f0, u, d_kill, weights=w)
    Phi = orthonormal_legendre_design(d_kill, u)
    f = f0 - Phi @ c
    # renorm f
    nrm = np.sqrt(float(np.sum(w * f * f)))
    if nrm <= 1e-15:
        raise RuntimeError("orthogonalization wiped f; increase waves")
    f = f / nrm
    phi_j = shifted_legendre_values(j, u)
    # φ_j already unit in continuous; discrete renorm lightly
    nrm_j = np.sqrt(float(np.sum(w * phi_j * phi_j)))
    phi_j = phi_j / nrm_j
    return np.sqrt(1.0 - eps * eps) * f + eps * phi_j


def critical_line_omega(T: float, t: float = 14.134725) -> float:
    """Angular frequency ω = t T for the critical-line pure mode on [0,1]."""
    return float(t * T)


def bound_R_d_sine_order(omega: float, d: int) -> float:
    """
    Explicit continuous M3 bound for omega >= 1.

    |c_k| <= C_k/omega, C_k = 2(k+1)sqrt(2k+1), and
    ||sin(omega .)||² >= 1/4. For 0 < omega < 1 return the contraction
    bound 1. Constants account for normalized Legendre endpoints/variation.
    """
    omega = _positive_finite(omega, "omega")
    constants = _coefficient_bounds(d)
    if omega < 1.0:
        return 1.0
    return float(4.0 * np.sum((constants / omega) ** 2))


def predicted_R_d_critical_scaling(T: float, t: float = 14.134725) -> float:
    """
    Exact continuous R_0(sin(t T u)), with its continuous limit at zero:
      c_0 = (1 - cos(tT))/(tT),  ‖q‖² → 1/2,
      R_0 ∼ 2 (1-cos(ω))² / ω²   with ω=tT.
    This is R_0 only. Other coefficients can have the same order or dominate;
    in particular c_0 vanishes at nonzero even multiples of pi.
    """
    if not np.isfinite(T) or not np.isfinite(t):
        raise ValueError("T and t must be finite")
    omega = float(t * T)
    if not np.isfinite(omega):
        raise ValueError("t*T must be finite")
    if abs(omega) < 1e-3:
        # sin(omega*u)/omega -> u, so R_0 -> (1/2)^2/(1/3)=3/4.
        z = omega * omega
        return float(0.75 + z / 40.0 + z * z / 11200.0)
    c0 = 2.0 * np.sin(0.5 * omega) ** 2 / omega
    # ‖sin‖² = ∫ sin² = 1/2 - sin(2ω)/(4ω) → 1/2
    l2 = 0.5 - np.sin(2.0 * omega) / (4.0 * omega)
    return float((c0 * c0) / l2)


def _positive_finite(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _coefficient_bounds(d: int, w_linf: float = 1.0, w_variation: float = 0.0) -> np.ndarray:
    """IBP constants for absolutely continuous weights of bounded variation.

    |phi_k| <= sqrt(2k+1). A degree-k Legendre polynomial has at most k
    monotone pieces, hence TV(phi_k) <= 2k sqrt(2k+1). Apply the product
    rule and add the two endpoint values of w*phi_k.
    """
    d = _validate_degree(d)
    W, V = float(w_linf), float(w_variation)
    if not np.isfinite(W) or W < 0.0 or not np.isfinite(V) or V < 0.0:
        raise ValueError("w_linf and w_variation must be finite and nonnegative")
    k = np.arange(d + 1, dtype=np.float64)
    return np.sqrt(2.0 * k + 1.0) * (2.0 * (k + 1.0) * W + V)


def _finite_mode_inputs(amplitudes, ordinates) -> tuple[np.ndarray, np.ndarray]:
    a = np.asarray(amplitudes, dtype=np.float64)
    t = np.asarray(ordinates, dtype=np.float64)
    if a.ndim != 1 or t.ndim != 1 or a.size == 0 or a.size != t.size:
        raise ValueError("amplitudes and ordinates must be nonempty 1D arrays of equal length")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(t)) or np.any(t <= 0.0):
        raise ValueError("amplitudes must be finite and ordinates positive and finite")
    return a, t


def _default_finite_mode_floor(T: float, a: np.ndarray, t: np.ndarray) -> float:
    """Certify the conventional quarter-amplitude floor uniformly in phases.

    Integrating each product of sines bounds the deviation from half the
    amplitude energy by error/T. A nonpositive return means the asymptotic
    floor is not certified at this T and the caller uses R_d <= 1 instead.
    """
    active = a != 0.0
    a, t = a[active], t[active]
    if not a.size:
        raise ValueError("nontrivial amplitudes required")
    if np.unique(t).size != t.size:
        raise ValueError("distinct active ordinates required without an explicit l2_floor")
    mass = float(np.dot(a, a))
    error = float(np.sum(a * a / (2.0 * t)))
    for i in range(a.size):
        error += float(np.sum(
            2.0 * np.abs(a[i] * a[i + 1:])
            * (1.0 / np.abs(t[i] - t[i + 1:]) + 1.0 / (t[i] + t[i + 1:]))
        ))
    return 0.25 * mass if error / T <= 0.25 * mass else 0.0


def bound_R_d_finite_mode_sum(
    T: float,
    amplitudes: Union[Sequence[float], np.ndarray],
    ordinates: Union[Sequence[float], np.ndarray],
    d: int,
    *,
    l2_floor: Optional[float] = None,
) -> float:
    """
    Explicit majorant consistent with Lemma M5 (finite-mode A₀).

    With C_k=2(k+1)sqrt(2k+1), the numerator bound is
      sum_k C_k² (sum_n |a_n|/(t_n T))².
    An explicit l2_floor is a caller-certified lower bound on ||q||².
    Otherwise a phase-uniform cross-term estimate certifies ¼ sum a_n²;
    when T is too small to certify that floor, return R_d <= 1.
    Active ordinates must be distinct when using the default floor.
    """
    T = _positive_finite(T, "T")
    constants = _coefficient_bounds(d)
    a, t = _finite_mode_inputs(amplitudes, ordinates)
    sum_term = float(np.sum(np.abs(a) / (t * T)))
    proj_energy_bound = float(np.sum(constants * constants)) * sum_term**2
    if l2_floor is None:
        l2_floor = _default_finite_mode_floor(T, a, t)
        if l2_floor == 0.0:
            return 1.0
    l2_floor = _positive_finite(l2_floor, "l2_floor")
    return float(proj_energy_bound / l2_floor)


def finite_mode_R_d_order_T(
    T: float,
    amplitudes: Union[Sequence[float], np.ndarray],
    ordinates: Union[Sequence[float], np.ndarray],
    d: int,
) -> float:
    """
    Leading O(T^{-2}) factor times T² from the M5 majorant:
      T² · bound_R_d_finite_mode_sum(T, ...) is T-independent once the
      default denominator floor is certified at that T.
    Useful for scaling checks (value stable as T grows).
    """
    return float(T * T * bound_R_d_finite_mode_sum(T, amplitudes, ordinates, d))


def bound_R_d_weighted_sine_order(
    omega: float,
    d: int,
    *,
    w_linf: float = 1.0,
    w_variation: float = 2.0,
    wq_l2_floor: float = 0.05,
) -> float:
    """
    Explicit majorant consistent with Lemma M6 for q=sin(ω·) under a weight
    with an absolutely continuous fixed weight satisfying ||w||_inf <=
    w_linf, ||w'||_1 <= w_variation, and ||w q||² >= wq_l2_floor.

    The coefficient constant is sqrt(2k+1)[2(k+1)w_linf+w_variation].
    Default variation 2 covers unit Hann/Tukey tapers; it is not inferred
    from w_linf. The denominator floor is an explicit assumption, including
    its default 0.05, and is not certified from these weight bounds alone.
    This is a continuous model bound, not an arithmetic Theorem A claim.
    """
    omega = _positive_finite(omega, "omega")
    floor = _positive_finite(wq_l2_floor, "wq_l2_floor")
    constants = _coefficient_bounds(d, w_linf, w_variation)
    return float(np.sum((constants / omega) ** 2) / floor)


def bound_R_d_weighted_finite_mode_sum(
    T: float,
    amplitudes: Union[Sequence[float], np.ndarray],
    ordinates: Union[Sequence[float], np.ndarray],
    d: int,
    *,
    w_linf: float = 1.0,
    w_variation: float = 2.0,
    wq_l2_floor: Optional[float] = None,
) -> float:
    """
    Weighted finite-mode A₀ majorant (Lemma M6 + M5).

    Use the same absolute-continuity and variation assumptions as the
    weighted sine bound, and sum |a_n|/(t_n*T) modewise.
    If wq_l2_floor is omitted, 0.05*sum(a_n²) is an assumed floor, suitable
    for distinct fixed modes under Hann/Tukey only for sufficiently large T.
    This default is not a certified finite-T lower bound for an arbitrary w.
    """
    T = _positive_finite(T, "T")
    constants = _coefficient_bounds(d, w_linf, w_variation)
    a, t = _finite_mode_inputs(amplitudes, ordinates)
    sum_term = float(np.sum(np.abs(a) / (t * T)))
    proj_energy_bound = float(np.sum(constants * constants)) * sum_term**2
    if wq_l2_floor is None:
        if np.unique(t[a != 0.0]).size != np.count_nonzero(a):
            raise ValueError("distinct active ordinates required without an explicit wq_l2_floor")
        wq_l2_floor = 0.05 * float(np.sum(a * a))
    wq_l2_floor = _positive_finite(wq_l2_floor, "wq_l2_floor")
    return float(proj_energy_bound / wq_l2_floor)


def weighted_cl_R_d_order_T(
    T: float,
    t: float = 14.134725,
    d: int = 4,
    *,
    w_linf: float = 1.0,
    w_variation: float = 2.0,
    wq_l2_floor: float = 0.05,
) -> float:
    """T² · bound_R_d_weighted_sine_order(t T, d) — T-independent scaling check."""
    omega = float(t * T)
    return float(
        T
        * T
        * bound_R_d_weighted_sine_order(
            omega, d, w_linf=w_linf, w_variation=w_variation,
            wq_l2_floor=wq_l2_floor
        )
    )
