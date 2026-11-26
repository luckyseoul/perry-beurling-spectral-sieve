"""
Ramp statistic of the 2026-11 RH proof, exact prime-sum moments, and the
RH asymptotic expansion of PBSS's continuous Legendre moments.

Conventions (as in ``docs/PSI_THETA_CORRECTION.md``)::

    T = log X,   y in [0, T],   phi_k(u) = sqrt(2k+1) P_k(2u - 1),
    c_k[R] = (1/T) int_0^T R(y) phi_k(y/T) dy,
    F(y) = e^{-y/2} (psi(e^y) - e^y),     Q(y) = e^{-y/2} (theta(e^y) - e^y).

Three tools, all documented in ``docs/RH_PROOF_SIMPLIFICATIONS.md``:

1. **Ramp statistic** (step 1 of the proof, Rubinstein's exposition section I.1)::

       S(x) = sum_{n <= x} Lambda(n) (x/n - 1),   G(x) = S(x) - x log x + (1 + gamma) x.

   Explicit formula: G(x) = log 2pi - sum_rho x^rho / (rho (rho - 1))
   - sum_{m >= 1} x^{-2m} / (2m (2m + 1)).  Under RH,
   |G(x) - log 2pi| <= A2 sqrt(x) + (1 - log 2) x^{-2} with
   A2 = sum_rho 1/|rho|^2 = 2 + gamma - log 4pi.  In PBSS coordinates
   g(y) = e^{-y/2} G(e^y) satisfies F = g' - g/2 exactly.

2. **Prime-sum moments.**  For residuals built from jumps w_n at log n::

       T c_k = sum_{n <= X} w_n W_k(log n) - B_k(T),
       W_k(y) = int_y^T e^{-u/2} phi_k(u/T) du,   B_k(T) = int_0^T e^{u/2} phi_k(u/T) du,

   with the B term absent for psi - theta.  This is an independent route to
   the interval-integrated moments of ``pbss.arithmetic_moments``.

3. **RH expansion.**  With Phi(s) = -zeta'(s)/(s zeta(s)) - 1/(s - 1), the
   Mellin transform of F at s = 1/2 + eps, put mu_j = (-1)^j Phi^{(j)}(1/2).
   Under RH, for every T > 0 and k >= 0::

       T c_k[F] = S_k(T) + O_k(T),
       S_k(T)   = sum_{j<=k} (phi_k^{(j)}(0)/j!) T^{-j} [mu_j - int_T^inf y^j e_F(y) dy],
       O_k(T)   = -sum_rho (e^{i gamma T}/rho) sum_{i<=k} (-1)^i phi_k^{(i)}(1) / ((i gamma)^{i+1} T^i),

   where e_F(y) = -log(2pi) e^{-y/2} - (1/2) e^{-y/2} log(1 - e^{-2y}).  The
   zero sum converges absolutely and |O_k(T)| <= C1 sum_i |phi_k^{(i)}(1)| (gamma_1 T)^{-i}
   with C1 = A2 (1 + 1/(4 gamma_1^2)).  For theta, c_k[Q] = c_k[F] - Delta_k exactly,
   with Delta_k from ``kind="psi_minus_theta"`` (no RH input).

None of this proves RH.  The expansion and envelope are consequences of RH;
RH has a machine-checked Lean proof (Alpoge, 2026-11-20) whose human review is
pending.  See the doc for status and non-claims.
"""

from __future__ import annotations

import math
from functools import lru_cache
from math import comb, factorial
from typing import Any, Iterable, Sequence

import mpmath as mp
import numpy as np

from .probes import primes_upto
from .zeros import zeta_zero_ordinates

EULER_GAMMA = 0.57721566490153286061
LOG_2PI = math.log(2.0 * math.pi)
#: sum over nontrivial zeros (with multiplicity) of 1/(rho (1 - rho)); equals sum 1/|rho|^2 under RH.
RH_A2 = 2.0 + EULER_GAMMA - math.log(4.0 * math.pi)
#: sum_{m>=1} 1/(2m(2m+1)): the trivial-zero term of the ramp explicit formula at x = 1.
TRIVIAL_ZERO_SUM = 1.0 - math.log(2.0)
#: ordinate of the first zero; no zero has 0 < |gamma| < GAMMA_1.
GAMMA_1 = 14.134725141734693
#: upper bound for sum_rho 1/gamma^2 under RH.
RH_C1 = RH_A2 * (1.0 + 1.0 / (4.0 * GAMMA_1 ** 2))

KINDS = ("theta", "psi", "psi_minus_theta")


# --------------------------------------------------------------------------------------------
# Legendre basis
# --------------------------------------------------------------------------------------------
def phi_coefficients(k: int) -> np.ndarray:
    """Coefficients a_{k,j} (j = 0..k) of phi_k(u) = sqrt(2k+1) P_k(2u-1) in powers of u."""
    if k < 0:
        raise ValueError("k >= 0")
    s = math.sqrt(2 * k + 1)
    return np.array([s * (-1) ** (k + j) * comb(k, j) * comb(k + j, j) for j in range(k + 1)], dtype=float)


def phi_derivatives(k: int, u: float) -> np.ndarray:
    """phi_k^{(i)}(u) for i = 0..k."""
    a = phi_coefficients(k)
    out = np.zeros(k + 1)
    for i in range(k + 1):
        out[i] = sum(a[j] * factorial(j) / factorial(j - i) * u ** (j - i) for j in range(i, k + 1))
    return out


# --------------------------------------------------------------------------------------------
# Prime-power events
# --------------------------------------------------------------------------------------------
def _higher_powers(primes: np.ndarray, n_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Prime powers p^j <= n_max with j >= 2: values (int64) and log p."""
    vals, wts = [], []
    j = 2
    while 2 ** j <= n_max:
        root = int(round(n_max ** (1.0 / j))) + 1
        cand = primes[primes <= root]
        pw = cand ** j
        keep = pw <= n_max
        if keep.any():
            vals.append(pw[keep])
            wts.append(np.log(cand[keep].astype(np.float64)))
        j += 1
    if not vals:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    return np.concatenate(vals), np.concatenate(wts)


def prime_power_events(X: float, kind: str = "psi") -> tuple[np.ndarray, np.ndarray]:
    """
    Jump values n <= X (int64) and jump heights for the residual ``kind``:
    theta (primes, height log p), psi (all prime powers, height Lambda(n)),
    psi_minus_theta (p^j with j >= 2, height log p).  Sorted by n.
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    n_max = int(math.floor(X))
    if n_max < 2:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    primes = primes_upto(n_max if kind != "psi_minus_theta" else math.isqrt(n_max)).astype(np.int64)
    parts_n, parts_w = [], []
    if kind in ("theta", "psi"):
        parts_n.append(primes)
        parts_w.append(np.log(primes.astype(np.float64)))
    if kind in ("psi", "psi_minus_theta"):
        hn, hw = _higher_powers(primes[primes <= math.isqrt(n_max)], n_max)
        parts_n.append(hn)
        parts_w.append(hw)
    n = np.concatenate(parts_n)
    w = np.concatenate(parts_w)
    order = np.argsort(n, kind="stable")
    return n[order], w[order]


# --------------------------------------------------------------------------------------------
# Exact moments as prime sums
# --------------------------------------------------------------------------------------------
def _upper_gamma_int(j: int, z):
    """Gamma(j+1, z) = j! e^{-z} sum_{i<=j} z^i / i!  (mpmath)."""
    return mp.factorial(j) * mp.exp(-z) * mp.fsum(z ** i / mp.factorial(i) for i in range(j + 1))


def prime_sum_moments(X: float, kmax: int = 4, kind: str = "psi", chunk: int = 2_000_000,
                      events: tuple[np.ndarray, np.ndarray] | None = None) -> dict[str, Any]:
    """
    Exact continuous moments c_0..c_kmax of the residual ``kind`` on [0, log X].

    For theta/psi this is c_k[Q] / c_k[F]; for psi_minus_theta it is
    Delta_k = c_k[F] - c_k[Q] (which includes the mean-one part in Delta_0).
    Event sums are accumulated in long double; the closed-form terms in mpmath.
    """
    n, w = events if events is not None else prime_power_events(X, kind)
    with mp.workdps(40):
        T = mp.log(mp.mpf(int(math.floor(X))) if float(X).is_integer() else mp.mpf(X))
        Tl = np.longdouble(str(T))
        M = [mp.mpf(0)] * (kmax + 1)
        sw = mp.mpf(0)
        for s in range(0, n.size, chunk):
            y = np.log(n[s:s + chunk].astype(np.longdouble))
            ww = w[s:s + chunk].astype(np.longdouble)
            e = np.exp(-y / 2) * ww
            term = np.ones_like(y)
            poly = np.zeros_like(y)
            sw += mp.mpf(str(np.sum(ww)))
            for j in range(kmax + 1):
                if j > 0:
                    term = term * (y / 2) / j
                poly = poly + term
                # E_j(y) = 2^{j+1} j! e^{-y/2} sum_{i<=j} (y/2)^i / i!
                M[j] += mp.mpf(str(np.sum(e * poly))) * 2 ** (j + 1) * mp.factorial(j)
        c, tc = [], []
        for k in range(kmax + 1):
            a = [mp.sqrt(2 * k + 1) * (-1) ** (k + j) * comb(k, j) * comb(k + j, j) for j in range(k + 1)]
            val = mp.mpf(0)
            for j in range(k + 1):
                E_T = 2 ** (j + 1) * _upper_gamma_int(j, T / 2)
                val += a[j] / T ** j * (M[j] - E_T * sw)
                if kind != "psi_minus_theta":
                    # B_k: int_0^T e^{u/2} u^j du = 2^{j+1} j! [e^{T/2} sum_i (-1)^{j-i} (T/2)^i/i! - (-1)^j]
                    J = 2 ** (j + 1) * mp.factorial(j) * (
                        mp.exp(T / 2) * mp.fsum((-1) ** (j - i) * (T / 2) ** i / mp.factorial(i) for i in range(j + 1))
                        - (-1) ** j)
                    val -= a[j] / T ** j * J
            tc.append(float(val))
            c.append(float(val / T))
        del Tl
    return {"X": float(X), "T": float(T), "kind": kind, "events": int(n.size), "c": c, "Tc": tc}


# --------------------------------------------------------------------------------------------
# Ramp statistic
# --------------------------------------------------------------------------------------------
def ramp_statistic(x: float | Sequence[float] | np.ndarray,
                   events: tuple[np.ndarray, np.ndarray] | None = None) -> np.ndarray:
    """G(x) = S(x) - x log x + (1 + gamma) x, exact from prime powers (long double)."""
    xs = np.atleast_1d(np.asarray(x, dtype=np.float64))
    if np.any(xs < 1):
        raise ValueError("x >= 1")
    n, w = events if events is not None else prime_power_events(float(xs.max()), "psi")
    nl = n.astype(np.longdouble)
    wl = w.astype(np.longdouble)
    cum_psi = np.concatenate([[np.longdouble(0)], np.cumsum(wl)])
    cum_a = np.concatenate([[np.longdouble(0)], np.cumsum(wl / nl)])
    idx = np.searchsorted(n, xs, side="right")
    xl = xs.astype(np.longdouble)
    G = xl * (cum_a[idx] + np.longdouble(1) + np.longdouble(EULER_GAMMA) - np.log(xl)) - cum_psi[idx]
    return G.astype(np.float64)


def ramp_residual(x, events=None) -> np.ndarray:
    """r(x) = (G(x) - log 2pi) / sqrt(x); under RH |r| <= A2 + (1 - log 2) x^{-5/2}."""
    xs = np.atleast_1d(np.asarray(x, dtype=np.float64))
    return (ramp_statistic(xs, events) - LOG_2PI) / np.sqrt(xs)


def ramp_window(Y: float, power: float = 1.5, events=None) -> dict[str, Any]:
    """
    One-sided window statistic of the proof's step 1 on [Y, Y^power].

    G is continuous and, between consecutive prime powers, G'(x) = sum_{n<=x} Lambda(n)/n
    + gamma - log x, so the extremes of r over the window are attained at prime
    powers up to O(log^3 x / x); the window endpoints are included.
    """
    hi = Y ** power
    n, w = events if events is not None else prime_power_events(hi, "psi")
    inside = n[(n >= Y) & (n <= hi)].astype(np.float64)
    pts = np.unique(np.concatenate([[float(Y), float(hi)], inside]))
    r = ramp_residual(pts, (n, w))
    i_max, i_min = int(np.argmax(r)), int(np.argmin(r))
    lower = -RH_A2 - TRIVIAL_ZERO_SUM * float(Y) ** -2.5
    return {"Y": float(Y), "power": power, "points": int(pts.size),
            "max_r": float(r[i_max]), "argmax_x": float(pts[i_max]),
            "min_r": float(r[i_min]), "argmin_x": float(pts[i_min]),
            "rh_upper": RH_A2, "rh_lower": lower,
            "within_rh_envelope": bool(r[i_max] <= RH_A2 and r[i_min] >= lower)}


def ramp_explicit_formula(x: float, ordinates: Iterable[float] | None = None) -> float:
    """log 2pi - sum_{|gamma|<=gamma_N} x^rho/(rho(rho-1)) - trivial-zero sum (RH form)."""
    gam = np.asarray(zeta_zero_ordinates() if ordinates is None else list(ordinates), dtype=float)
    rho = 0.5 + 1j * gam
    zs = 2.0 * np.real(np.exp(rho * math.log(x)) / (rho * (rho - 1)))
    z = 1.0 / x
    trivial = sum(z ** (2 * m) / (2 * m * (2 * m + 1)) for m in range(1, 200) if z ** (2 * m) > 1e-30)
    return LOG_2PI - float(np.sum(zs)) - trivial


def ramp_zero_tail_envelope(x: float, ordinates: Iterable[float] | None = None) -> float:
    """Under RH, |G(x) - ramp_explicit_formula(x)| <= sqrt(x) (A2 - sum_{n<=N} 2/|rho_n|^2)."""
    gam = np.asarray(zeta_zero_ordinates() if ordinates is None else list(ordinates), dtype=float)
    return math.sqrt(x) * (RH_A2 - float(np.sum(2.0 / (0.25 + gam ** 2))))


# --------------------------------------------------------------------------------------------
# RH asymptotic expansion
# --------------------------------------------------------------------------------------------
def _Phi(s):
    return -mp.zeta(s, 1, 1) / (s * mp.zeta(s)) - 1 / (s - 1)


@lru_cache(maxsize=8)
def psi_mellin_constants(jmax: int = 8, dps: int = 30, radius: str = "0.25", points: int = 96) -> tuple:
    """
    mu_j = (-1)^j Phi^{(j)}(1/2) = Abel-regularized int_0^inf y^j F(y) dy, j = 0..jmax.

    Cauchy formula on |s - 1/2| = radius (Phi is analytic for |s - 1/2| < 1/2).
    mu_0 = 2 - 2 zeta'(1/2)/zeta(1/2).
    """
    with mp.workdps(dps):
        r, half = mp.mpf(radius), mp.mpf(1) / 2
        acc = [mp.mpc(0)] * (jmax + 1)
        for m in range(points):
            e = r * mp.expjpi(mp.mpf(2 * m) / points)
            v = _Phi(half + e)
            for j in range(jmax + 1):
                acc[j] += v * e ** (-j)
        return tuple(float((-1) ** j * mp.factorial(j) * mp.re(a) / points) for j, a in enumerate(acc))


def _eF_upper_moment(j: int, T: float) -> float:
    """int_T^inf y^j e_F(y) dy with e_F = -log(2pi) e^{-y/2} + (1/2) sum_m e^{-(2m+1/2) y}/m."""
    with mp.workdps(30):
        T = mp.mpf(T)
        val = -mp.log(2 * mp.pi) * 2 ** (j + 1) * _upper_gamma_int(j, T / 2)
        m = 1
        while True:
            a = 2 * m + mp.mpf(1) / 2
            t = _upper_gamma_int(j, a * T) / a ** (j + 1) / (2 * m)
            val += t
            if abs(t) < mp.mpf(10) ** -25 or m > 400:
                break
            m += 1
        return float(val)


def psi_smooth_part(T: float, k: int) -> float:
    """S_k(T): the zero-free part of T c_k[F] under RH (see module docstring)."""
    mu = psi_mellin_constants(max(8, k))
    a = phi_coefficients(k)
    return float(sum(a[j] / T ** j * (mu[j] - _eF_upper_moment(j, T)) for j in range(k + 1)))


def psi_limit(k: int) -> float:
    """lim_{T->inf} S_k(T) = phi_k(0) mu_0 = (-1)^k sqrt(2k+1) (2 - 2 zeta'/zeta(1/2))."""
    return float(phi_coefficients(k)[0] * psi_mellin_constants()[0])


def t2_rd_limit_psi(d: int = 4) -> float:
    """sum_{k=2}^d (2k+1) mu_0^2 / A2: the limit of T^2 R_d[F] if all zeros are simple
    (otherwise A2 in the denominator becomes sum m_rho^2/|rho|^2 >= A2, so this is an upper bound
    for the non-oscillating part)."""
    return sum(2 * k + 1 for k in range(2, d + 1)) * psi_mellin_constants()[0] ** 2 / RH_A2


def rh_oscillation_envelope(T: float, k: int) -> float:
    """Under RH: |T c_k[F] - S_k(T)| <= C1 sum_i |phi_k^{(i)}(1)| / (gamma_1 T)^i."""
    d = np.abs(phi_derivatives(k, 1.0))
    return float(RH_C1 * sum(d[i] / (GAMMA_1 * T) ** i for i in range(k + 1)))


def zero_sum_oscillation(T: float, k: int, ordinates: Iterable[float] | None = None) -> float:
    """O_k(T) summed over the given positive ordinates (and their conjugates)."""
    gam = np.asarray(zeta_zero_ordinates() if ordinates is None else list(ordinates), dtype=float)
    d = phi_derivatives(k, 1.0)
    rho = 0.5 + 1j * gam
    inner = sum((-1) ** i * d[i] / ((1j * gam) ** (i + 1) * T ** i) for i in range(k + 1))
    return float(-2.0 * np.sum(np.real(np.exp(1j * gam * T) / rho * inner)))


def zero_tail_envelope(T: float, k: int, ordinates: Iterable[float] | None = None) -> float:
    """Under RH: |O_k(T) - zero_sum_oscillation(T, k, ordinates)| <= this bound."""
    gam = np.asarray(zeta_zero_ordinates() if ordinates is None else list(ordinates), dtype=float)
    tail = max(RH_C1 - float(np.sum(2.0 / gam ** 2)), 0.0)
    d = np.abs(phi_derivatives(k, 1.0))
    g_n = float(gam.max())
    return float(tail * sum(d[i] / (g_n * T) ** i for i in range(k + 1)))


def predict_moments(X: float, kmax: int = 4, ordinates: Iterable[float] | None = None) -> dict[str, Any]:
    """
    RH predictions for T c_k[F] and T c_k[Q] at X, using the zeros in ``ordinates``.
    The theta prediction subtracts the exact T Delta_k (primes <= sqrt X only).
    """
    gam = list(zeta_zero_ordinates() if ordinates is None else ordinates)
    delta = prime_sum_moments(X, kmax, "psi_minus_theta")
    T = delta["T"]
    rows = []
    for k in range(kmax + 1):
        s = psi_smooth_part(T, k)
        o = zero_sum_oscillation(T, k, gam)
        rows.append({"k": k, "S_k": s, "O_k_zeros": o, "zero_tail_bound": zero_tail_envelope(T, k, gam),
                     "rh_envelope": rh_oscillation_envelope(T, k),
                     "Tc_psi_pred": s + o, "T_delta": delta["Tc"][k], "Tc_theta_pred": s + o - delta["Tc"][k]})
    return {"X": float(X), "T": T, "zeros": len(gam), "rows": rows}
