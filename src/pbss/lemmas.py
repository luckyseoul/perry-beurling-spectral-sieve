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
Let N < ∞, amplitudes a_n ∈ ℝ, ordinates t_n > 0, phases φ_n ∈ ℝ, and
  q_T(u) = ∑_{n=1}^N a_n sin(t_n T u + φ_n)   (or the cosine form).
Integration by parts on each mode gives |⟨q_T, φ_k⟩| ≤ C_k ∑_n |a_n|/(t_n T).
Hence ‖P_d q_T‖₂² = O_{d,N,{a,t}}(T^{-2}). For large T, ‖q_T‖₂² is bounded
below by a positive constant depending only on {a_n} (diagonal sine terms →
½∑a_n²; cross terms O(T^{-1})), so when a ≢ 0,
  R_d(q_T) = O_d(T^{-2}) as T → ∞
at the **same order** as pure-mode M3. This is the finite-mode extension of A₀
for truncated explicit-formula residuals. It does **not** prove full Theorem A
for the arithmetic prime residual or RH.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

from .basis import shifted_legendre_values
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
    Crude explicit majorant consistent with Lemma M3:
    |c_k| ≤ 2/ω for each k (loose), ‖q‖² ≥ 1/4 for ω≥2π, so
    R_d ≤ 4 d' (2/ω)² with d'=d+1 — not sharp, used only as a scaling check.
    """
    if omega <= 0:
        raise ValueError("omega > 0")
    return float(4.0 * (d + 1) * (2.0 / omega) ** 2)


def predicted_R_d_critical_scaling(T: float, t: float = 14.134725) -> float:
    """
    Leading continuous asymptotics for R_0(sin(t T u)):
      c_0 = (1 - cos(tT))/(tT),  ‖q‖² → 1/2,
      R_0 ∼ 2 (1-cos(ω))² / ω²   with ω=tT.
    Higher d only adds smaller terms; this is the dominant piece of R_d.
    """
    omega = t * T
    if abs(omega) < 1e-12:
        return 1.0
    c0 = (1.0 - np.cos(omega)) / omega
    # ‖sin‖² = ∫ sin² = 1/2 - sin(2ω)/(4ω) → 1/2
    l2 = 0.5 - np.sin(2.0 * omega) / (4.0 * omega)
    return float((c0 * c0) / l2)


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

    Unnormalized crude bound (same style as ``bound_R_d_sine_order``):
      |c_k| ≤ ∑_n |a_n| · 2 / (t_n T)   (loose integration-by-parts),
      ‖P_d q‖² ≤ (d+1) (∑ |a_n|·2/(t_n T))²,
      ‖q‖² ≥ l2_floor  (default: ¼ ∑ a_n², valid for large T when a≢0).

    Hence R_d ≤ 4(d+1) (∑ |a_n|/(t_n T))² / l2_floor = O_d(T^{-2}).

    This is a scaling majorant for tests, not a sharp constant.
    """
    if T <= 0:
        raise ValueError("T > 0")
    if d < 0:
        raise ValueError("d >= 0")
    a = np.asarray(amplitudes, dtype=np.float64).ravel()
    t = np.asarray(ordinates, dtype=np.float64).ravel()
    if a.size == 0 or a.size != t.size:
        raise ValueError("amplitudes and ordinates must be nonempty and same length")
    if np.any(t <= 0):
        raise ValueError("ordinates must be positive")
    sum_term = float(np.sum(np.abs(a) * 2.0 / (t * T)))
    proj_energy_bound = (d + 1) * (sum_term**2)
    if l2_floor is None:
        l2_floor = 0.25 * float(np.sum(a * a))
    if l2_floor <= 0:
        raise ValueError("l2_floor must be positive (nontrivial amplitudes)")
    return float(proj_energy_bound / l2_floor)


def finite_mode_R_d_order_T(
    T: float,
    amplitudes: Union[Sequence[float], np.ndarray],
    ordinates: Union[Sequence[float], np.ndarray],
    d: int,
) -> float:
    """
    Leading O(T^{-2}) factor times T² from the M5 majorant:
      T² · bound_R_d_finite_mode_sum(T, ...)  is T-independent.
    Useful for scaling checks (value stable as T grows).
    """
    return float(T * T * bound_R_d_finite_mode_sum(T, amplitudes, ordinates, d))
