#!/usr/bin/env python3
"""Independent (parent-authored) high-precision check of the shipped arithmetic
moment evaluator, using the repo's prime checkpoints but our own formulas.

Q(y) = e^{-y/2}(theta(e^y) - e^y);  c_k(T) = (1/T) int_0^T Q(y) phi_k(y/T) dy
V(T)/T = (1/T) int_0^T Q(y)^2 dy ;  R_d = sum_{k=2..d} c_k^2 / (V/T - c_0^2 - c_1^2)

Interval identities (theta = Theta on [a,b]):
  int_a^b e^{-y/2}(Theta - e^y) phi_k(y/T) dy = Theta*I^{-}_k(a,b) - I^{+}_k(a,b)
      I^{s}_k(a,b) = int_a^b e^{s y} phi_k(y/T) dy,  s=±1/2
  int_a^b e^{-y}(Theta - e^y)^2 dy = Theta^2(e^{-a}-e^{-b}) - 2 Theta (b-a) + (e^b-e^a)
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

CKPT = "/home/nick/perry-beurling-spectral-sieve/results/prime_checkpoints/primes_le_10000000000.npy"


def legendre_power_coeffs(k):
    from numpy.polynomial import legendre as Lg, polynomial as Pl
    c = Lg.Legendre.basis(k).convert(kind=Pl.Polynomial).coef
    out = np.zeros(len(c))
    for j, cj in enumerate(c):
        if cj == 0:
            continue
        term = np.array([1.0])
        for _ in range(j):
            new = np.zeros(len(term) + 1)
            new[:-1] += -term
            new[1:] += 2 * term
            term = new
        out[: len(term)] += cj * term
    return out * np.sqrt(2 * k + 1)


def _partial_log_sums(args):
    import mpmath as mp
    i, lo, hi, primes_np, digits = args
    with mp.workdps(digits):
        s = mp.fsum([mp.log(int(p)) for p in primes_np[lo:hi]])
    return i, mp.nstr(s, digits)


def _chunk_contrib(args):
    import mpmath as mp
    i, lo, hi, primes_np, T_str, theta0_str, degree, digits, next_prime = args
    with mp.workdps(digits):
        T = mp.mpf(T_str)
        primes = [int(p) for p in primes_np[lo:hi]]
        logs = [mp.log(p) for p in primes]
        q_polys = [[mp.mpf(float(c)) for c in legendre_power_coeffs(k)] for k in range(degree + 1)]
        ck = [mp.mpf(0)] * (degree + 1)
        norm = mp.mpf(0)
        cum = mp.mpf(theta0_str)
        half = mp.mpf("0.5")
        b_end = mp.log(next_prime) if next_prime else T
        intervals = []
        if i == 0:
            intervals.append((mp.mpf(0), logs[0], mp.mpf(0)))  # leading segment [0, log 2), theta = 0
        running = mp.mpf(theta0_str)
        for j in range(len(logs)):
            running += logs[j]
            intervals.append((logs[j], logs[j + 1] if j + 1 < len(logs) else b_end, running))
        for a, b, cum in intervals:
            if b <= a:
                continue
            ea2, eb2 = mp.e ** (-a / 2), mp.e ** (-b / 2)
            eA2, eB2 = mp.e ** (a / 2), mp.e ** (b / 2)
            h = b - a
            # powers
            ap = [mp.mpf(1)] * 6
            bp = [mp.mpf(1)] * 6
            for j in range(1, 6):
                ap[j] = ap[j - 1] * a
                bp[j] = bp[j - 1] * b
            # I^-_j : s=-1/2 ;  I^+_j : s=+1/2   (j=0..5)
            Im = []
            Ip = []
            for j in range(6):
                if j == 0:
                    Im0 = (eb2 - ea2) / (-half)
                    Ip0 = (eB2 - eA2) / half
                    Im.append(Im0)
                    Ip.append(Ip0)
                else:
                    Im.append((bp[j] * eb2 - ap[j] * ea2) / (-half) + (2 * j) * Im[j - 1])
                    Ip.append((bp[j] * eB2 - ap[j] * eA2) / half - (2 * j) * Ip[j - 1])
            for k in range(degree + 1):
                q = q_polys[k]
                acc_m = mp.mpf(0)
                acc_p = mp.mpf(0)
                scale = mp.mpf(1)
                for j in range(len(q)):
                    cj = q[j] / scale
                    acc_m += cj * Im[j]
                    acc_p += cj * Ip[j]
                    scale *= T
                ck[k] += cum * acc_m - acc_p
            norm += cum ** 2 * (mp.e ** (-a) - mp.e ** (-b)) - 2 * cum * h + (mp.e ** b - mp.e ** a)
        return i, [v / T for v in ck], norm / T


def run(T_target, workers, digits):
    primes = np.load(CKPT, mmap_mode="r")
    n_use = int(np.searchsorted(primes, int(np.exp(T_target)) + 1))
    primes_np = np.asarray(primes[:n_use])
    print(f"T={T_target}: primes used = {n_use}", flush=True)
    bounds = np.linspace(0, n_use, workers + 1).astype(int)
    degree = 4
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        sums = dict(ex.map(_partial_log_sums,
                           [(i, bounds[i], bounds[i + 1], primes_np, digits) for i in range(workers)]))
    offsets = {}
    import mpmath as mp
    with mp.workdps(digits):
        acc = mp.mpf(0)
        for i in range(workers):
            offsets[i] = mp.nstr(acc, digits)
            acc += mp.mpf(sums[i])
        theta_total = mp.nstr(acc, 25)
    print(f"  pass1 {time.time()-t0:.1f}s theta_total={theta_total}", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_chunk_contrib,
                          [(i, bounds[i], bounds[i + 1], primes_np, repr(float(T_target)), offsets[i], degree, digits,
                            (int(primes_np[bounds[i + 1]]) if i + 1 < workers and bounds[i + 1] < n_use else 0))
                           for i in range(workers) if bounds[i + 1] > bounds[i]]))
    with mp.workdps(digits):
        ck = [mp.mpf(0)] * (degree + 1)
        norm = mp.mpf(0)
        for _, cks, ns in res:
            for k in range(degree + 1):
                ck[k] += cks[k]
            norm += ns
        return {
            "T": T_target, "n_primes": n_use, "digits": digits, "workers": workers,
            "c": [mp.nstr(v, 25) for v in ck],
            "l2_norm_sq": mp.nstr(norm, 25),
            "affine_detrended_l2_norm_sq": mp.nstr(norm - ck[0] ** 2 - ck[1] ** 2, 25),
            "affine_Rd": mp.nstr(sum(ck[k] ** 2 for k in range(2, degree + 1)) / (norm - ck[0] ** 2 - ck[1] ** 2), 25),
            "elapsed_s": time.time() - t0,
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, nargs="+", default=[12.0, 14.0])
    ap.add_argument("--workers", type=int, default=40)
    ap.add_argument("--digits", type=int, default=45)
    ap.add_argument("--out", default="/home/nick/scratch/pbss_attempt/my_check.json")
    a = ap.parse_args()
    out = []
    for T in a.T:
        out.append(run(T, a.workers, a.digits))
        print(json.dumps(out[-1], indent=1), flush=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("WROTE", a.out)


if __name__ == "__main__":
    main()
