#!/usr/bin/env python3
"""
Ramp statistic and RH moment asymptotics for PBSS (2026-11-26).

Doc: docs/RH_PROOF_SIMPLIFICATIONS.md.  Library: src/pbss/rh_moment_asymptotics.py.

Part A  Ramp statistic G(x) of the proof's step 1: reproduce the published sanity values,
        compare with the explicit formula over the first N zeros, scan r(x) = (G - log 2pi)/sqrt x
        against the RH envelope A2 = 2 + gamma - log 4pi, and evaluate the one-sided windows
        [Y, Y^{3/2}] that the proof's step 1 consumes.
Part B  Mellin constants mu_j = (-1)^j Phi^{(j)}(1/2) and their small-x (s = 0 pole) structure.
Part C  Exact continuous moments T c_k (psi and theta) versus the RH expansion S_k(T) + O_k(T),
        X = 10^3 .. 10^9, with the rigorous zero-tail bound; predictions at 5e10 and 1e12.
Part D  The R4 plateau: predicted R4 from the expansion and the measured normalizer (PR #2
        artifact), the fraction of each moment's limit reached, and the asymptotic T^2 R4.

Compute: worker count follows full_workers = max(1, nproc - 2).  On a 2-core host that is 1,
so the X loop runs serially; each X is one vectorized long-double pass.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from pbss import rh_moment_asymptotics as R  # noqa: E402

PR2_JSON = os.path.join(HERE, "..", "results", "psi_theta_correction_20261120", "psi_theta_correction.json")
PUBLISHED_G = {10_000: 1.749, 100_000: 3.671, 1_000_000: 11.040, 2_000_000: -12.268}


def zeros(n: int) -> list[float]:
    with mp.workdps(20):
        return [float(mp.im(mp.zetazero(i))) for i in range(1, n + 1)]


def part_a(gam, scan_max: float, window_ys) -> dict:
    t0 = time.time()
    xs = sorted(PUBLISHED_G)
    G = R.ramp_statistic(xs)
    published = [{"x": x, "G": float(g), "published": PUBLISHED_G[x], "abs_diff": abs(float(g) - PUBLISHED_G[x])}
                 for x, g in zip(xs, G)]
    explicit = []
    for x in [1e3, 1e4, 1e5, 1e6, 2e6, 1e7]:
        g = float(R.ramp_statistic([x])[0])
        gz = R.ramp_explicit_formula(x, gam)
        explicit.append({"x": x, "G": g, "G_zeros": gz, "diff": g - gz, "rh_tail_bound": R.ramp_zero_tail_envelope(x, gam)})
    n, w = R.prime_power_events(scan_max, "psi")
    r = R.ramp_residual(n.astype(np.float64), (n, w))
    big = n >= 100
    scan = {"x_max": scan_max, "prime_powers": int(n.size), "max_r": float(r[big].max()),
            "argmax": float(n[big][np.argmax(r[big])]), "min_r": float(r[big].min()),
            "argmin": float(n[big][np.argmin(r[big])]), "max_abs_r_over_A2": float(np.abs(r[big]).max() / R.RH_A2)}
    windows = [R.ramp_window(Y, 1.5) for Y in window_ys]
    print(f"[A] published G reproduced (max diff {max(p['abs_diff'] for p in published):.1e}); "
          f"scan to {scan_max:.0e}: max|r|/A2 = {scan['max_abs_r_over_A2']:.3f}; "
          f"windows all within envelope: {all(wd['within_rh_envelope'] for wd in windows)} ({time.time() - t0:.1f}s)",
          flush=True)
    return {"published": published, "explicit_formula": explicit, "scan": scan, "windows": windows,
            "A2": R.RH_A2, "trivial_zero_sum": R.TRIVIAL_ZERO_SUM}


def part_b() -> dict:
    mu = R.psi_mellin_constants(8)
    with mp.workdps(30):
        mu0 = 2 - 2 * mp.zeta(mp.mpf(1) / 2, 1, 1) / mp.zeta(mp.mpf(1) / 2)
    pole = [-2 * R.LOG_2PI * 2 ** j * math.factorial(j) for j in range(9)]
    print(f"[B] mu_0 = {mu[0]:.15f} (2 - 2 zeta'/zeta(1/2) = {float(mu0):.15f}); "
          f"mu_8 / pole term = {mu[8] / pole[8]:.6f}", flush=True)
    return {"mu": list(mu), "mu0_closed_form": float(mu0), "s0_pole_term": pole,
            "ratio_to_pole_term": [m / p for m, p in zip(mu, pole)]}


def moments_at(X: float, kmax: int, gam) -> dict:
    t0 = time.time()
    psi = R.prime_sum_moments(X, kmax, "psi")
    dlt = R.prime_sum_moments(X, kmax, "psi_minus_theta")
    T = psi["T"]
    rows = []
    for k in range(kmax + 1):
        s, o = R.psi_smooth_part(T, k), R.zero_sum_oscillation(T, k, gam)
        tc_t = psi["Tc"][k] - dlt["Tc"][k]
        rows.append({"k": k, "Tc_psi": psi["Tc"][k], "Tc_theta": tc_t, "T_delta": dlt["Tc"][k],
                     "S_k": s, "O_k_zeros": o, "residual": psi["Tc"][k] - s,
                     "after_zeros": psi["Tc"][k] - s - o, "zero_tail_bound": R.zero_tail_envelope(T, k, gam),
                     "rh_envelope": R.rh_oscillation_envelope(T, k), "fraction_of_limit": s / R.psi_limit(k)})
    return {"X": X, "T": T, "prime_powers": psi["events"], "rows": rows, "seconds": time.time() - t0}


def part_c(xs, kmax, gam, workers) -> list:
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            out = list(ex.map(moments_at, xs, [kmax] * len(xs), [gam] * len(xs)))
    else:
        out = [moments_at(X, kmax, gam) for X in xs]
    for o in out:
        worst = max(abs(r["after_zeros"]) / r["zero_tail_bound"] for r in o["rows"])
        inside = all(abs(r["residual"]) <= r["rh_envelope"] for r in o["rows"])
        print(f"[C] X={o['X']:.0e}  T={o['T']:.3f}  max|after zeros|/tail bound = {worst:.3f}  "
              f"within RH envelope: {inside}  ({o['seconds']:.1f}s)", flush=True)
    return out


def crossover(k: int, frac: float) -> float:
    """Smallest T (bisection on a grid) with S_k(T)/limit >= frac thereafter (monotone here)."""
    lo, hi = 2.0, 2000.0
    f = lambda T: R.psi_smooth_part(T, k) / R.psi_limit(k) - frac  # noqa: E731
    if f(hi) < 0:
        return float("nan")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return hi


def part_d(c_rows, d: int = 4) -> dict:
    with open(PR2_JSON) as fh:
        stored = json.load(fh)["part_c_continuous_R4_theta_vs_psi"]
    by_x = {round(o["X"]): o for o in c_rows}
    table = []
    for row in stored:
        X = row["psi"]["X"]
        if X not in by_x:
            continue
        o, T = by_x[X], by_x[X]["T"]
        for kind, key in (("psi", "Tc_psi"), ("theta", "Tc_theta")):
            meas = row[kind]
            pred_tc = [r["S_k"] + r["O_k_zeros"] - (r["T_delta"] if kind == "theta" else 0.0) for r in o["rows"]]
            num_pred = sum(pred_tc[k] ** 2 for k in range(2, d + 1)) / T ** 2
            num_meas = sum(r[key] ** 2 for r in o["rows"][2:d + 1]) / T ** 2
            table.append({"X": X, "kind": kind, "T": T, "N": meas["N"], "N_aff": meas["N_aff"],
                          "R4_measured": meas["R4"], "R4_predicted": num_pred / meas["N_aff"],
                          "numerator_measured": num_meas, "numerator_predicted": num_pred,
                          "T2_R4_measured": meas["T2_R4"], "NT_minus_A2T": (meas["N"] - R.RH_A2) * T})
    cross = {k: {"T_50pct": crossover(k, 0.5), "T_90pct": crossover(k, 0.9), "2k(k+1)": 2 * k * (k + 1)}
             for k in range(2, d + 1)}
    n1 = [t["NT_minus_A2T"] for t in table if t["kind"] == "psi"]
    out = {"table": table, "crossover_T": cross, "t2_r4_limit_psi_simple_zeros": R.t2_rd_limit_psi(d),
           "psi_N_T_minus_A2_T": n1, "n1_over_A2": (sum(n1) / len(n1)) / R.RH_A2}
    worst = max(abs(t["R4_predicted"] - t["R4_measured"]) / t["R4_measured"] for t in table)
    print(f"[D] R4 predicted from the expansion vs measured: max relative diff {worst:.2e}; "
          f"T^2 R4 limit (psi, simple zeros) = {out['t2_r4_limit_psi_simple_zeros']:.1f}; "
          f"T_90% for k=2,3,4 = {[round(cross[k]['T_90pct'], 1) for k in cross]}", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zeros", type=int, default=200, help="zeta zeros for the oscillatory term (mpmath)")
    ap.add_argument("--x-max-exp", type=int, default=9, help="largest exponent in the X = 10^e grid")
    ap.add_argument("--scan-max", type=float, default=1e8, help="ramp scan upper limit")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "rh_moment_asymptotics_20261126",
                                                  "rh_moment_asymptotics.json"))
    args = ap.parse_args(argv)
    workers = max(1, (os.cpu_count() or 1) - 2)
    t0 = time.time()
    gam = zeros(args.zeros)
    print(f"[0] {len(gam)} zeros (gamma_N = {gam[-1]:.3f}); workers = {workers}", flush=True)
    a = part_a(gam, args.scan_max, [10 ** e for e in range(2, 6)])
    b = part_b()
    xs = [float(10 ** e) for e in range(3, args.x_max_exp + 1)]
    c = part_c(xs, 4, gam, workers)
    dd = part_d(c)
    pred = [R.predict_moments(X, 4, gam) for X in (5e10, 1e12)]
    for p in pred:
        print(f"[C] prediction X={p['X']:.0e}: T c_k[theta] = "
              + ", ".join(f"{r['Tc_theta_pred']:+.4f}±{r['zero_tail_bound']:.4f}" for r in p["rows"]), flush=True)
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    except Exception:  # pragma: no cover
        head = ""
    out = {"experiment": "rh_moment_asymptotics", "doc": "docs/RH_PROOF_SIMPLIFICATIONS.md",
           "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": platform.node(),
           "python": platform.python_version(), "numpy": np.__version__, "mpmath": mp.__version__,
           "source_git_head": head, "workers": workers,
           "rh_status": "RH is used as a hypothesis for S_k/O_k and the envelopes; it has a machine-checked "
                        "Lean proof (2026-11-20/26) pending human review. Nothing here proves RH.",
           "zeros": {"count": len(gam), "gamma_N": gam[-1], "ordinates": gam},
           "part_a_ramp": a, "part_b_mellin": b, "part_c_moments": c, "part_d_r4_plateau": dd,
           "predictions": pred, "runtime_seconds": time.time() - t0}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[done] {args.out} ({out['runtime_seconds']:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
