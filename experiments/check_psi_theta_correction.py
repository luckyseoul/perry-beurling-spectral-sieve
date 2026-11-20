#!/usr/bin/env python3
"""psi - theta correction to the continuous Legendre moments, and the R_4 plateau test (2026-11-20).

Notation as in docs/CONTINUOUS_THEOREM_A.md:
    F(y) = e^{-y/2} (psi(e^y) - e^y),   Q(y) = e^{-y/2} (theta(e^y) - e^y),   0 <= y <= T = log X,
    c_k[R] = (1/T) int_0^T R(y) phi_k(y/T) dy,   phi_k(u) = sqrt(2k+1) L_k(2u-1).
    g(y) = F(y) - Q(y) = e^{-y/2} (psi - theta)(e^y) >= 0,   h = 1 - g  (the h of Theorem A, section 2).

Parts (see docs/PSI_THETA_CORRECTION.md):
  A. Delta_k(X) = c_k[F] - c_k[Q] = <g_T, phi_k>, exact: only prime powers p^j <= X with j >= 2 enter,
     so p <= sqrt(X). Interval-by-interval closed-form integrals in mpmath. For k = 2 also the exact
     three-term split Delta_2 = (sqrt5/T) [K0(T) - 6 K1(T)/T + 6 K2(T)/T^2], K_m(T) = int_0^T (g-1) y^m dy.
  B. Limits: K0(inf) = 2 (C3 - B - gamma - 1) in closed form, with
     B = sum_p log p / (p(p-1)) and C3 = sum_p log p / (p^{3/2} - p); K0, K1, K2 at infinity also by
     Cauchy's formula on G(eps) = int_0^inf (g-1) e^{-eps y} dy (independent check of the closed form).
  C. Continuous affine-detrended R_4 of the theta and psi residuals from exact prime-interval integrals
     (local stable form of docs/EXACT_ARITHMETIC_MOMENTS.md, NumPy longdouble, chunked). The theta path
     is checked against pbss.arithmetic_moments at X = 1e5.

Floating-point evaluations of exact identities; not interval certificates; no RH claim.
Compute: serial by design. Part A touches only primes <= sqrt(X); part C is one vectorized,
memory-bound pass per window (chunked), so worker fan-out buys little at these sizes.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from math import comb, factorial

import numpy as np

try:
    import mpmath as mpm
except ImportError:  # pragma: no cover
    sys.exit("mpmath is required for this experiment: pip install mpmath")

LD = np.longdouble
DEGREE = 4
NSER = 30


def sieve(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=np.uint8)
    s[:2] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = 0
    return np.nonzero(s)[0].astype(np.int64)


def mobius(n: int) -> int:
    r, p, m = 1, 2, n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r


def legendre_coeffs_in_y(k: int, T):
    """phi_k(y/T) = sum_j a_j y^j (mpmath)."""
    return [mpm.sqrt(2 * k + 1) * (-1) ** (k + j) * comb(k, j) * comb(k + j, j) / T ** j for j in range(k + 1)]


# ----------------------------------------------------------------------------------------- part A
def prime_power_events(X):
    """(log p^j, log p) for p^j <= X, j >= 2, sorted by position."""
    ev = []
    for p in sieve(int(mpm.sqrt(X)) + 2):
        p = int(p)
        lp = mpm.log(p)
        pk, j = p * p, 2
        while pk <= X:
            ev.append((j * lp, lp))
            pk *= p
            j += 1
    ev.sort(key=lambda t: t[0])
    return ev


def _antiderivative(n: int):
    """d/dy of the result is y^n e^{-y/2}."""
    return lambda y: -2 * mpm.exp(-y / 2) * sum(
        mpm.mpf(factorial(n)) / factorial(n - r) * 2 ** r * y ** (n - r) for r in range(n + 1))


def integral_g_times_poly(ev, T, coeffs):
    """int_0^T e^{-y/2} (psi - theta)(e^y) * sum_n coeffs[n] y^n dy, exact on each constant piece."""
    anti = [_antiderivative(n) for n in range(len(coeffs))]
    F = lambda y: sum(c * anti[n](y) for n, c in enumerate(coeffs) if c != 0)
    total, C = mpm.mpf(0), mpm.mpf(0)
    for i, (y0, lp) in enumerate(ev):
        C += lp
        y1 = ev[i + 1][0] if i + 1 < len(ev) else T
        total += C * (F(y1) - F(y0))
    return total


def part_a(xs, dps):
    mpm.mp.dps = dps
    rows = []
    for X in xs:
        X = mpm.mpf(X)
        T = mpm.log(X)
        ev = prime_power_events(X)
        delta = {k: integral_g_times_poly(ev, T, legendre_coeffs_in_y(k, T)) / T for k in range(1, DEGREE + 1)}
        K = [integral_g_times_poly(ev, T, [0] * m + [1]) - T ** (m + 1) / (m + 1) for m in range(3)]
        split = [mpm.sqrt(5) * K[0] / T, -6 * mpm.sqrt(5) * K[1] / T ** 2, 6 * mpm.sqrt(5) * K[2] / T ** 3]
        rows.append({
            "X": mpm.nstr(X, 6), "T": mpm.nstr(T, 15), "prime_power_jumps": len(ev),
            "delta_k": {str(k): mpm.nstr(v, 15) for k, v in delta.items()},
            "K_m_T": [mpm.nstr(v, 15) for v in K],
            "delta_2_three_terms": [mpm.nstr(v, 12) for v in split],
            "delta_2_from_terms": mpm.nstr(sum(split), 15),
        })
        print(f"[A] X={mpm.nstr(X, 4):>8} T={mpm.nstr(T, 8):>10}  Delta_2={mpm.nstr(delta[2], 13):>17}"
              f"  Delta_3={mpm.nstr(delta[3], 9):>12}  Delta_4={mpm.nstr(delta[4], 9):>12}"
              f"  terms={[mpm.nstr(v, 6) for v in split]}", flush=True)
    return rows


# ----------------------------------------------------------------------------------------- part B
def prime_log_zeta(w):
    """P(w) = sum_p log p p^{-w} = -sum_j mu(j) zeta'(jw)/zeta(jw), Re w > 1."""
    tot = 0
    for j in range(1, 2000):
        if mpm.mpf(2) ** (-j * mpm.re(w)) < mpm.mpf(10) ** (-mpm.mp.dps - 5):
            break
        m = mobius(j)
        if m:
            tot += -m * mpm.zeta(j * w, derivative=1) / mpm.zeta(j * w)
    return tot


def _converging(w):
    return mpm.mpf(2) ** (-mpm.re(w)) > mpm.mpf(10) ** (-mpm.mp.dps - 3)


def G_eps(eps):
    """G(eps) = int_0^inf (g - 1) e^{-eps y} dy, analytic near 0 (radius 1/6)."""
    s = mpm.mpf(1) / 2 + eps
    w = 2 * s
    core = -mpm.zeta(w, derivative=1) / mpm.zeta(w) - 1 / (w - 1)
    k3 = sum(prime_log_zeta(k * s) for k in range(3, 400) if _converging(k * s))
    j2 = sum(prime_log_zeta(2 * j * s) for j in range(2, 200) if _converging(2 * j * s))
    return (core - 1 + k3 - j2) / s


def part_b(dps_closed, dps_cauchy):
    mpm.mp.dps = dps_closed
    B = sum(prime_log_zeta(mpm.mpf(m)) for m in range(2, 400) if _converging(mpm.mpf(m)))
    C3 = sum(prime_log_zeta(mpm.mpf(m) / 2) for m in range(3, 800) if _converging(mpm.mpf(m) / 2))
    K0_closed = 2 * (C3 - B - mpm.euler - 1)
    out = {"B": mpm.nstr(B, 25), "C3": mpm.nstr(C3, 25), "K0_inf_closed_form": mpm.nstr(K0_closed, 25),
           "large_window_constant_sqrt5_K0": mpm.nstr(mpm.sqrt(5) * K0_closed, 25)}
    print(f"[B] B={mpm.nstr(B, 18)}  C3={mpm.nstr(C3, 18)}  K0(inf)={mpm.nstr(K0_closed, 18)}"
          f"  sqrt5*K0={mpm.nstr(mpm.sqrt(5) * K0_closed, 15)}", flush=True)
    mpm.mp.dps = dps_cauchy
    r, M = mpm.mpf("0.05"), 48
    pts = [r * mpm.exp(2j * mpm.pi * k / M) for k in range(M)]
    vals = [G_eps(z) for z in pts]
    K = [mpm.re((-1) ** m * factorial(m) * sum(v / z ** m for v, z in zip(vals, pts)) / M) for m in range(3)]
    roots = mpm.polyroots([6 * K[2], -6 * K[1], K[0]])
    Ts = sorted(1 / mpm.re(z) for z in roots if abs(mpm.im(z)) < mpm.mpf(10) ** (-10))
    out.update({
        "K_m_inf_cauchy": [mpm.nstr(v, 15) for v in K],
        "closed_minus_cauchy_K0": mpm.nstr(K0_closed - K[0], 3),
        "delta_2_sign_changes_T_limiting_moments": [mpm.nstr(t, 8) for t in Ts],
        "delta_2_sign_changes_log10X": [mpm.nstr(t / mpm.log(10), 6) for t in Ts],
        "leading_term_exceeds_next_for_T_above": mpm.nstr(6 * K[1] / K[0], 8),
    })
    print(f"[B] Cauchy K0,K1,K2 = {[mpm.nstr(v, 12) for v in K]}  sign changes at T = {[mpm.nstr(t, 6) for t in Ts]}"
          f"  (6 K1/K0 = {mpm.nstr(6 * K[1] / K[0], 6)})", flush=True)
    return out


# ----------------------------------------------------------------------------------------- part C
def step_events(X: int, kind: str):
    P = sieve(int(X))
    vals, wts = [P], [np.log(P.astype(LD))]
    if kind == "psi":
        base = P[P.astype(np.float64) ** 2 <= X]
        cur, lb = base * base, np.log(base.astype(LD))
        while cur.size:
            keep = cur <= X
            vals.append(cur[keep])
            wts.append(lb[keep])
            base, lb = base[keep], lb[keep]
            cur = cur[keep] * base
    v, w = np.concatenate(vals), np.concatenate(wts)
    order = np.argsort(v, kind="stable")
    return v[order], w[order]


def continuous_moments(X: int, kind: str, chunk: int = 400_000) -> dict:
    """Exact continuous c_0..c_4, N, N_aff and affine-detrended R_4 for the theta or psi residual."""
    T = np.log(LD(X))
    v, w = step_events(X, kind)
    C = np.concatenate([[LD(0)], np.cumsum(w)])
    start = np.concatenate([[LD(1)], v.astype(LD)])
    a_all = np.log(start)
    b_all = np.concatenate([a_all[1:], [T]])
    n = np.arange(NSER, dtype=np.int64)
    fact = np.array([LD(factorial(int(i))) for i in n])
    odd = n % 2 == 1
    sgn, sgn1 = ((-1) ** n).astype(LD), ((-1) ** (n + 1)).astype(LD)
    acoef = [[LD(np.sqrt(LD(2 * k + 1))) * (-1) ** (k + m) * comb(k, m) * comb(k + m, m) for m in range(k + 1)]
             for k in range(DEGREE + 1)]
    n_sum, c_sum = LD(0), [LD(0)] * (DEGREE + 1)
    for s0 in range(0, a_all.size, chunk):
        a = a_all[s0:s0 + chunk]
        h = b_all[s0:s0 + chunk] - a
        x = start[s0:s0 + chunk]
        delta = C[s0:s0 + chunk] / x - 1
        pw = (h[:, None] / 2) ** n[None, :]
        f = [np.sum(sgn[None, :] * pw / (fact[None, :] * (j + n + 1)), axis=1) for j in range(DEGREE + 1)]
        g = [2 * np.sum(np.where(odd[None, :], pw / (fact[None, :] * (j + n + 1)), 0), axis=1)
             for j in range(DEGREE + 1)]
        H = h[:, None] ** n[None, :]
        a0 = np.sum(np.where(n[None, :] >= 1, sgn1[None, :] * H / fact[None, :], 0), axis=1)
        b0 = np.sum(np.where(n[None, :] >= 2, sgn[None, :] * H / fact[None, :], 0), axis=1)
        c0 = 2 * np.sum(np.where((n[None, :] >= 3) & odd[None, :], H / fact[None, :], 0), axis=1)
        n_sum += np.sum(x * (a0 * delta ** 2 - 2 * b0 * delta + c0))
        u0, wv = a / T, h / T
        for k in range(DEGREE + 1):
            tot = np.zeros_like(h)
            for j in range(k + 1):
                bkj = wv ** j * sum(acoef[k][m] * comb(m, j) * u0 ** (m - j) for m in range(j, k + 1))
                tot = tot + bkj * (delta * f[j] - g[j])
            c_sum[k] += np.sum(np.sqrt(x) * h * tot)
    N = n_sum / T
    c = [cs / T for cs in c_sum]
    n_aff = N - c[0] ** 2 - c[1] ** 2
    num = sum(ck ** 2 for ck in c[2:])
    return {"X": X, "T": float(T), "intervals": int(a_all.size), "c": [float(ck) for ck in c],
            "N": float(N), "N_aff": float(n_aff), "R4": float(num / n_aff), "T2_R4": float(T * T * num / n_aff)}


def self_check_theta(X: int = 100_000) -> dict:
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
        from pbss.arithmetic_moments import arithmetic_moments
    except Exception as exc:  # pragma: no cover
        return {"status": "skipped", "reason": repr(exc)}
    ref = arithmetic_moments(T=float(np.log(LD(X))), degree=DEGREE, max_prime_limit=X).to_dict()
    mine = continuous_moments(X, "theta")
    ref_c = [float(s) for s in ref["coefficients_decimal"]]
    dc = max(abs(a - b) for a, b in zip(mine["c"], ref_c))
    dr = abs(mine["R4"] - ref["affine_detrended_energy_ratio"])
    ok = dc < 1e-12 and dr < 1e-12
    print(f"[C] self-check vs pbss.arithmetic_moments at X={X:.0e}: max|dc|={dc:.2e}  |dR4|={dr:.2e}  ->"
          f" {'OK' if ok else 'MISMATCH'}", flush=True)
    return {"status": "ok" if ok else "mismatch", "X": X, "max_abs_coefficient_diff": dc, "abs_R4_diff": dr}


def part_c(xs):
    check = self_check_theta()
    rows = []
    for X in xs:
        t0 = time.time()
        th, ps = continuous_moments(X, "theta"), continuous_moments(X, "psi")
        rows.append({"theta": th, "psi": ps})
        print(f"[C] X={X:.0e} T={th['T']:.4f}  R4 theta={th['R4']:.6f} (T^2 R4={th['T2_R4']:.2f})"
              f"  R4 psi={ps['R4']:.6f} (T^2 R4={ps['T2_R4']:.2f})  [{time.time() - t0:.1f}s]", flush=True)
    return check, rows


# ----------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--corr-x", default="1e6,1e8,1e10,5e10,1e12", help="endpoints X for part A")
    ap.add_argument("--r4-x", default="1e3,1e4,1e5,1e6,1e7", help="endpoints X for part C (1e8 takes minutes)")
    ap.add_argument("--dps", type=int, default=40, help="mpmath digits for parts A and B (closed form)")
    ap.add_argument("--skip-b", action="store_true")
    ap.add_argument("--skip-c", action="store_true")
    ap.add_argument("--out", default="results/psi_theta_correction_20261120/psi_theta_correction.json")
    args = ap.parse_args(argv)
    t_start = time.time()
    to_x = lambda s: [int(float(v)) for v in s.split(",") if v.strip()]
    report = {
        "experiment": "psi_theta_correction_20261120",
        "doc": "docs/PSI_THETA_CORRECTION.md",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(), "python": platform.python_version(),
        "numpy": np.__version__, "mpmath": mpm.__version__,
        "longdouble_mantissa_bits": int(np.finfo(LD).nmant),
        "certified_interval_result": False, "rh_claim": False,
        "conventions": "y in [0, T], T = log X; phi_k orthonormal shifted Legendre; c_k = (1/T) int R phi_k(y/T) dy",
    }
    try:
        report["source_git_head"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__)), text=True).strip()
    except Exception:
        report["source_git_head"] = None
    report["part_a_exact_corrections"] = part_a(to_x(args.corr_x), args.dps)
    if not args.skip_b:
        report["part_b_limits"] = part_b(args.dps, 30)
    if not args.skip_c:
        check, rows = part_c(to_x(args.r4_x))
        report["part_c_self_check"] = check
        report["part_c_continuous_R4_theta_vs_psi"] = rows
    report["runtime_seconds"] = round(time.time() - t_start, 1)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"wrote {args.out} ({report['runtime_seconds']} s)")


if __name__ == "__main__":
    main()
