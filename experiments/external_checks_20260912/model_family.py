#!/usr/bin/env python3
"""Exact R_d(T) for explicit exponential-spectrum families (obstruction robustness).

Model
-----
    Q(y) = sum_{n=n0}^{N} A_n * e^{a_n y} * cos(gamma_n y)

    gamma_n = exp(n**beta)          # ordinate scale
    A_n     = gamma_n**(-q)         # q=1  ==  arithmetic residue size 1/|rho|
    a_n     = A - delta_n           # A = 1/4, unattained boundary
        delta_n = 1/log(gamma_n)            (deficit mode 'log'; (A-a_n)log gamma_n = 1)
        delta_n = D * n**(-p)               (deficit mode 'pow')

Diagnostics (continuous L^2([0,1],du), shifted-Legendre phi_k):
    c_k(T) = int_0^1 Q(Tu) phi_k(u) du
    V(T)   = int_0^T Q(y)^2 dy           (self terms exactly; top-K cross terms exactly)
    R_d    = sum_{k=2..d} c_k(T)^2 / ( V(T)/T - c_0^2 - c_1^2 )

All closed-form; mpmath high precision.  Emits one JSON per case.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from fractions import Fraction

import mpmath as mp


def shifted_legendre_powers(k: int):
    """phi_k(u) = sqrt(2k+1) L_k(2u-1) as [float coeffs] in powers of u (exact rationals then float)."""
    # L_k(x) in powers of x via recurrence on rationals
    # L_0 = 1, L_1 = x, (n+1)L_{n+1} = (2n+1)x L_n - n L_{n-1}
    polys = {0: [Fraction(1)], 1: [Fraction(0), Fraction(1)]}
    for n in range(1, k):
        prev, cur = polys[n - 1], polys[n]
        m = len(cur) + 1
        nxt = [Fraction(0)] * m
        for i, c in enumerate(cur):
            nxt[i + 1] += Fraction(2 * n + 1, n + 1) * c
        for i, c in enumerate(prev):
            nxt[i] -= Fraction(n, n + 1) * c
        polys[n + 1] = nxt
    base = polys[k]
    # compose x = 2u-1
    out = [Fraction(0)] * (len(base))
    for j, c in enumerate(base):
        if c == 0:
            continue
        # (2u-1)^j
        term = [Fraction(1)]
        for _ in range(j):
            new = [Fraction(0)] * (len(term) + 1)
            for i, d in enumerate(term):
                new[i] += -d
                new[i + 1] += 2 * d
            term = new
        for i, d in enumerate(term):
            out[i] = out[i] + c * d if i < len(out) else c * d
    scale = mp.sqrt(2 * k + 1)
    return [mp.mpf(f.numerator) / mp.mpf(f.denominator) for f in out], scale


def moment_term(lam, T, powers, scale):
    """int_0^1 exp(lam*T*u) * phi_k(u) du  for lam = a + i*gamma."""
    mu = lam * T
    # I_j = int_0^1 u^j e^{mu u} du ;  I_j = (e^mu - j I_{j-1})/mu
    emu = mp.e ** mu
    I = []
    prev = (emu - 1) / mu
    I.append(prev)
    for j in range(1, len(powers)):
        prev = (emu - j * prev) / mu
        I.append(prev)
    s = mp.mpf(0)
    for c, v in zip(powers, I):
        s += c * v
    return s * scale


def model_params(beta, q, A, n0, nmax, deficit, D, p):
    rows = []
    for n in range(n0, nmax + 1):
        g = mp.e ** (mp.mpf(n) ** beta)
        amp = g ** (-q)
        if deficit == "log":
            d = 1 / mp.log(g)
        elif deficit == "zero":
            d = mp.mpf(0)
        elif deficit == "expn":
            d = mp.e ** (-D * n)
        elif deficit == "shift":          # delta_n = D + p/n  (uniformly positive, sup still unattained)
            d = D + p / n
        elif deficit == "attained":
            d = mp.mpf(0) if n <= int(D) else mp.mpf(p) * (n - int(D))
        elif deficit == "expn2":
            d = mp.e ** (-D * n * n)
        else:
            d = D * mp.mpf(n) ** (-p)
        rows.append((n, g, amp, A - d, d))
    return rows


def term_logs(rows, T, A, kind):
    """log10 of |contribution| proxies to pick truncation / dominant indices."""
    out = []
    for (n, g, amp, a, d) in rows:
        if kind == "energy":
            # self term ~ amp^2 * e^{2 a T} / (2 d)  (if dT large); use e^{2 a T} min(T, 1/d)
            damp = min(T, 1 / d) if d > 0 else T
            val = 2 * a * T + 2 * mp.log(amp) + mp.log(damp / 2)
        else:
            # moment ~ amp * e^{aT} / (g T)
            val = a * T + mp.log(amp) - mp.log(g * T)
        out.append(float(val / mp.log(10)))
    return out


def run_case(case):
    A = case.get("A", 0.25)
    beta = case["beta"]
    q = case.get("q", 1.0)
    n0 = case.get("n0", 8)
    deficit = case.get("deficit", "log")
    D = case.get("D", 1.0)
    p = case.get("p", 1.0)
    Ts = case["Ts"]
    dmax = case.get("dmax", 4)
    topk = case.get("topk", 48)
    extra_dps = case.get("dps", 60)

    out = {"case": case, "per_T": {}, "meta": {}}

    # precision: phases gamma_n * T need enough digits
    res = []
    for T in Ts:
        T = mp.mpf(T)
        # choose truncation by a first pass over a wide range using log proxies
        W = int(case.get("scan_nmax", 20000))
        rows_scan = model_params(beta, q, A, n0, min(W, 6000), deficit, D, p)
        el = term_logs(rows_scan, T, A, "energy")
        ml = term_logs(rows_scan, T, A, "moment")
        # find dominant index and cutoff where both profiles fall 45 orders below peak
        def cutoff(logs):
            peak = max(logs)
            idx = max(range(len(logs)), key=lambda i: logs[i])
            for i in range(idx + 1, len(logs)):
                if logs[i] < peak - 45:
                    return i + 1
            return len(logs)
        N = max(cutoff(el), cutoff(ml))
        N = max(N, n0 + 4)
        N = min(N, case.get("nmax_cap", 8000))
        # digits needed for gamma_N * T
        digits = int(mp.log10(rows_scan[min(N - 1, len(rows_scan) - 1)][1] * T)) + 1
        dps = max(extra_dps, digits + 40)
        mp.mp.dps = dps
        rows = model_params(beta, q, A, n0, N, deficit, D, p)
        el_rows = term_logs(rows, T, A, "energy")
        ml_rows = term_logs(rows, T, A, "moment")
        nstar_energy = max(range(len(rows)), key=lambda i: el_rows[i])
        nstar_moment = max(range(len(rows)), key=lambda i: ml_rows[i])

        powers = {}
        scales = {}
        for k in range(dmax + 1):
            pw, sc = shifted_legendre_powers(k)
            powers[k] = pw
            scales[k] = sc

        ck = [mp.mpf(0)] * (dmax + 1)
        for (n, g, amp, a, d) in rows:
            lam = mp.mpc(a, g)
            for k in range(dmax + 1):
                ck[k] += amp * mp.re(moment_term(lam, T, powers[k], scales[k]))

        # energy: diagonal self terms exactly
        V_self = mp.mpf(0)
        self_parts = []
        for (n, g, amp, a, d) in rows:
            tail = (mp.e ** ((2 * a + 2j * g) * T) - 1) / (2 * a + 2j * g)
            if abs(a) < mp.mpf('1e-12'):
                # limiting case of a pole exactly on the critical line: (e^{2aT}-1)/(2a) -> T
                main = T
            else:
                main = (mp.e ** (2 * a * T) - 1) / (2 * a)
            val = amp ** 2 * mp.mpf(0.5) * (main + mp.re(tail))
            V_self += val
            self_parts.append((float(val), n))

        # cross terms among top-K amplitudes-of-contribution indices for the energy
        order = sorted(range(len(rows)), key=lambda i: -self_parts[i][0])
        keep = order[:topk]
        V_cross = mp.mpf(0)
        for ii in range(len(keep)):
            for jj in range(ii + 1, len(keep)):
                n1, g1, A1, a1, d1 = rows[keep[ii]]
                n2, g2, A2, a2, d2 = rows[keep[jj]]
                sig = a1 + a2
                tot = mp.mpf(0)
                for om in (abs(g1 - g2), g1 + g2):
                    intg = (mp.e ** ((sig + 1j * om) * T) - 1) / (sig + 1j * om)
                    tot += mp.re(intg)
                V_cross += 2 * A1 * A2 * mp.mpf(0.5) * tot  # both orders of the double sum

        V = V_self + V_cross
        norm2 = V / T
        denom = norm2 - ck[0] ** 2 - ck[1] ** 2
        Rd = {f"d{k}": float(sum(ck[j] ** 2 for j in range(2, k + 1)) / denom) for k in range(2, dmax + 1)}
        out["per_T"][str(T)] = {
            "T": float(T),
            "N": N,
            "dps": dps,
            "nstar_energy": rows[nstar_energy][0],
            "nstar_energy_gamma": float(rows[nstar_energy][1]),
            "nstar_moment": rows[nstar_moment][0],
            "ck": [float(v) for v in ck],
            "V_self": float(V_self),
            "V_cross_topk": float(V_cross),
            "V_total": float(V),
            "norm2": float(norm2),
            "detrend_frac": float((ck[0] ** 2 + ck[1] ** 2) / norm2),
            "Rd": Rd,
            "log10_Rd_d2": float(mp.log10(abs(Rd["d2"]))) if Rd["d2"] > 0 else None,
        }
    return out


def selftest():
    """Brute-force check of moment/V formulas on a small low-frequency model."""
    mp.mp.dps = 40
    beta, q, A, n0, N = 1.0, 1.0, 0.25, 1, 3
    case = {"beta": beta, "q": q, "A": A, "n0": n0, "deficit": "log", "Ts": [3.0], "dmax": 4}
    # override ordinates to small values for brute force
    def Qfun(y):
        s = mp.mpf(0)
        for n in range(n0, N + 1):
            g = mp.mpf(3 + n)  # small ordinates
            amp = mp.mpf(0.1) * n
            a = A - mp.mpf(1) / n
            s += amp * mp.e ** (a * y) * mp.cos(g * y)
        return s

    T = mp.mpf(3)
    # brute force c_k by quadrature
    for k in range(5):
        pw, sc = shifted_legendre_powers(k)
        f = lambda u: Qfun(T * u) * sum(c * u ** j for j, c in enumerate(pw)) * sc
        ck_q = mp.quad(f, [0, 0.5, 1])
        # closed form
        ck_cf = mp.mpf(0)
        for n in range(n0, N + 1):
            g = mp.mpf(3 + n)
            amp = mp.mpf(0.1) * n
            a = A - mp.mpf(1) / n
            ck_cf += amp * mp.re(moment_term(mp.mpc(a, g), T, pw, sc))
        print(f"k={k} quad={mp.nstr(ck_q, 20)} closed={mp.nstr(ck_cf, 20)} diff={mp.nstr(abs(ck_q - ck_cf), 5)}")
    Vq = mp.quad(lambda y: Qfun(y) ** 2, [0, T])
    Vs = mp.mpf(0)
    for n in range(n0, N + 1):
        g = mp.mpf(3 + n)
        amp = mp.mpf(0.1) * n
        a = A - mp.mpf(1) / n
        emu2 = mp.e ** ((2 * a + 2j * g) * T)
        Vs += amp ** 2 * mp.mpf(0.5) * ((mp.e ** (2 * a * T) - 1) / (2 * a) + mp.re((emu2 - 1) / (2 * a + 2j * g)))
    Vc = mp.mpf(0)
    for n1 in range(n0, N + 1):
        for n2 in range(n1 + 1, N + 1):
            g1, g2 = mp.mpf(3 + n1), mp.mpf(3 + n2)
            A1, A2 = mp.mpf(0.1) * n1, mp.mpf(0.1) * n2
            a1, a2 = A - mp.mpf(1) / n1, A - mp.mpf(1) / n2
            sig = a1 + a2
            tot = sum(mp.re((mp.e ** ((sig + 1j * om) * T) - 1) / (sig + 1j * om)) for om in (abs(g1 - g2), g1 + g2))
            Vc += 2 * A1 * A2 * mp.mpf(0.5) * tot
    print("V quad =", mp.nstr(Vq, 20), " self+cross =", mp.nstr(Vs + Vc, 20), " diff=", mp.nstr(abs(Vq - (Vs + Vc)), 5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--cases", required=False, help="JSON file with list of case dicts")
    ap.add_argument("--out", required=False)
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    with open(args.cases) as fh:
        cases = json.load(fh)
    results = []
    t0 = time.time()
    for case in cases:
        r = run_case(case)
        results.append(r)
        print(json.dumps({"done": case.get("name", "?"), "elapsed_s": time.time() - t0}), flush=True)
    with open(args.out, "w") as fh:
        json.dump({"results": results, "elapsed_s": time.time() - t0}, fh, indent=1)
    print("WROTE", args.out)


if __name__ == "__main__":
    main()
