#!/usr/bin/env python3
"""Independent check of the pure-mode phase-trough constant.

q_T(u) = e^{A T u} sin(t T u), R_d = sum_{k=0..d} <q,phi_k>^2 / ||q||^2.
Claim (docs/OFF_CRITICAL_PHASE_OBSTRUCTION.md eq. (3)):
    liminf_{T->inf} T^3 R_d(T) = A d^2 (d+1)^2 (d+2)^2 / (3 (A^2+t^2)^2)   (d>=1)
Test point here: A=0.5, t=3  ->  d=2: 1536/1369, d=4: 38400/1369.
"""
from __future__ import annotations

import json
import sys

import mpmath as mp

sys.path.insert(0, "/home/nick/scratch/pbss_attempt")
from model_family import shifted_legendre_powers  # exact rational power coeffs


def ck_of_T(T, polys, A, t, dmax):
    lam = mp.mpc(A, t) * T
    # int_0^1 e^{lam u} u^j du via recursion
    emu = mp.e ** lam
    I = [(emu - 1) / lam]
    for j in range(1, dmax + 1):
        I.append((emu - j * I[-1]) / lam)
    out = []
    for k in range(dmax + 1):
        pw, sc = polys[k]
        s = mp.mpf(0)
        for c, v in zip(pw, I):
            s += c * v
        out.append(mp.im(s) * sc)
    return out


def norm_of_T(T, A, t):
    a = 2 * A * T
    b = 2 * t * T
    main = (mp.e ** a - 1) / a
    osc = mp.re((mp.e ** mp.mpc(a, b) - 1) / mp.mpc(a, b))
    return mp.mpf("0.5") * (main - osc)


def Rd(T, polys, A, t, d):
    c = ck_of_T(T, polys, A, t, d)
    n = norm_of_T(T, A, t)
    return sum(v ** 2 for v in c) / n


def main():
    mp.mp.dps = 40
    A = mp.mpf("0.5")
    t = mp.mpf(3)
    dmax = 4
    polys = {k: shifted_legendre_powers(k) for k in range(dmax + 1)}
    r2 = A * A + t * t
    pred = {d: A * d ** 2 * (d + 1) ** 2 * (d + 2) ** 2 / (3 * r2 ** 2) for d in (1, 2, 3, 4)}

    beta = mp.atan2(t, A)
    rows = []
    best = {d: (None, None) for d in (1, 2, 3, 4)}
    n0 = int((1000 * t - beta) / mp.pi) + 1
    n1 = int((1300 * t - beta) / mp.pi)
    for n in range(n0, n1 + 1):
        theta_n = beta + n * mp.pi
        for d in (1, 2, 3, 4):
            K = (d + 1) ** 2
            L = mp.mpf(d) * (d + 1) ** 2 * (d + 2) / 2
            D = t * L / (r2 * K)
            Tstar = (theta_n + mp.sqrt(theta_n ** 2 - 4 * t * D)) / (2 * t)
            f = lambda T: T ** 3 * Rd(T, polys, A, t, d)
            lo, hi = Tstar - mp.mpf("0.02"), Tstar + mp.mpf("0.02")
            gr = (mp.sqrt(5) - 1) / 2
            c_, d_ = hi - gr * (hi - lo), lo + gr * (hi - lo)
            fc, fd = f(c_), f(d_)
            for _ in range(80):
                if fc < fd:
                    hi, d_, fd = d_, c_, fc
                    c_ = hi - gr * (hi - lo)
                    fc = f(c_)
                else:
                    lo, c_, fc = c_, d_, fd
                    d_ = lo + gr * (hi - lo)
                    fd = f(d_)
            Tb = (lo + hi) / 2
            val = Tb ** 3 * Rd(Tb, polys, A, t, d)
            if best[d][0] is None or val < best[d][0]:
                best[d] = (val, Tb)
        rows.append({"n": n, "T": float(theta_n / t)})
    out = {"A": float(A), "t": float(t), "scan": f"T in [{float(rows[0]['T']):.1f}, {float(rows[-1]['T']):.1f}] ({(len(rows))} troughs)"}
    print(f"scan: T from {float(rows[0]['T']):.2f} to {float(rows[-1]['T']):.2f}, {len(rows)} trough windows")
    for d in (1, 2, 3, 4):
        val, Tb = best[d]
        p = pred[d]
        print(f"d={d}: observed min T^3 R_d = {mp.nstr(val, 15)}  at T={mp.nstr(Tb, 10)}   "
              f"predicted = {mp.nstr(p, 15)}   ratio = {mp.nstr(val/p, 12)}")
        out[f"d{d}"] = {"observed_min_T3Rd": mp.nstr(val, 18), "at_T": mp.nstr(Tb, 12),
                        "predicted": mp.nstr(p, 18), "ratio": mp.nstr(val / p, 12)}
    json.dump(out, open("/home/nick/scratch/pbss_attempt/phase_trough_check.json", "w"), indent=1)
    print("WROTE /home/nick/scratch/pbss_attempt/phase_trough_check.json")


if __name__ == "__main__":
    main()
