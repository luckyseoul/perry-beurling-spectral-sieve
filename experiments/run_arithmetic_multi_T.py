#!/usr/bin/env python3
"""
Arithmetic residual multi-T campaign for PBSS.

For each logarithmic window length T:
  - build q_T from real primes via shipped arithmetic_residual
  - project with shipped project() → R_d(q_T)
  - compare to critical-line model mode and persistent-defect control

Primes are sieved once up to max x_max = exp(T_max), then shared with workers.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _worker(payload: dict) -> dict:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"

    from pbss.probes import (
        arithmetic_residual,
        probe_critical_line_mode,
        probe_persistent_defect,
        sample_grid,
    )
    from pbss.projection import project

    T = float(payload["T"])
    degree = int(payload["degree"])
    n_points = int(payload["n_points"])
    eps = float(payload["eps"])
    detrend = str(payload["detrend"])
    smooth = int(payload["smooth"])
    primes = np.asarray(payload["primes"], dtype=np.int64)

    u = sample_grid(n_points)
    q_ar, T_ar, meta = arithmetic_residual(
        u, T=T, primes=primes, detrend=detrend, smooth=smooth
    )
    r_ar = project(q_ar, u, degree=degree, T=T_ar)

    r_cl = project(probe_critical_line_mode(u, T=T), u, degree=degree, T=T)
    r_def = project(
        probe_persistent_defect(u, eps=eps, j=0, waves=80),
        u,
        degree=degree,
        T=T,
    )

    return {
        "T": T_ar,
        "x_max": meta["x_max"],
        "n_primes_used": meta["n_primes"],
        "detrend": meta["detrend"],
        "arithmetic": {
            "R_d": r_ar.energy_ratio,
            "S_d": r_ar.scaled_strength,
            "energy": r_ar.energy,
            "l2_norm_sq": r_ar.l2_norm_sq,
        },
        "critical_line": {
            "R_d": r_cl.energy_ratio,
            "S_d": r_cl.scaled_strength,
        },
        "persistent_defect": {
            "R_d": r_def.energy_ratio,
            "S_d": r_def.scaled_strength,
            "eps2": eps * eps,
        },
        "gap_def_minus_arith": r_def.energy_ratio - r_ar.energy_ratio,
        "gap_arith_minus_cl": r_ar.energy_ratio - r_cl.energy_ratio,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument(
        "--T-list",
        type=str,
        default="8,9,10,11,12,13,14,15,16",
        help="comma-separated T = log(x_max) values",
    )
    ap.add_argument("--detrend", type=str, default="deg1")
    ap.add_argument("--smooth", type=int, default=1)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--scratch", type=str, default="")
    args = ap.parse_args()

    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    T_values = sorted(T_values)
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)

    from pbss.probes import primes_upto, xmax_from_T

    x_max_global = xmax_from_T(max(T_values))
    print(
        f"Sieving primes up to x_max={x_max_global:.3e} (T_max={max(T_values)}) …",
        flush=True,
    )
    t0 = time.time()
    primes = primes_upto(int(x_max_global))
    print(
        f"  n_primes={primes.size} in {time.time()-t0:.2f}s; "
        f"workers={workers} d={args.degree} detrend={args.detrend}",
        flush=True,
    )

    # primes as list for pickling (numpy ok too)
    payloads = [
        {
            "T": T,
            "degree": args.degree,
            "n_points": args.n_points,
            "eps": args.eps,
            "detrend": args.detrend,
            "smooth": args.smooth,
            "primes": primes,
        }
        for T in T_values
    ]

    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for row in ex.map(_worker, payloads):
            rows.append(row)
            print(
                f"T={row['T']:6.2f}  x_max={row['x_max']:.3e}  "
                f"R_ar={row['arithmetic']['R_d']:.4e}  "
                f"R_cl={row['critical_line']['R_d']:.4e}  "
                f"R_def={row['persistent_defect']['R_d']:.4e}  "
                f"gap_def-ar={row['gap_def_minus_arith']:.4e}",
                flush=True,
            )

    R_ar = np.array([r["arithmetic"]["R_d"] for r in rows])
    R_cl = np.array([r["critical_line"]["R_d"] for r in rows])
    R_def = np.array([r["persistent_defect"]["R_d"] for r in rows])
    T_arr = np.array([r["T"] for r in rows])

    # qualitative reading (honest — may not decay)
    if len(R_ar) >= 3:
        # linear slope of log R vs log T if positive R
        logT = np.log(T_arr)
        logR = np.log(np.maximum(R_ar, 1e-30))
        slope = float(np.polyfit(logT, logR, 1)[0])
    else:
        slope = float("nan")

    mid = len(rows) // 2
    reading = {
        "R_ar_first": float(R_ar[0]),
        "R_ar_last": float(R_ar[-1]),
        "R_ar_min": float(np.min(R_ar)),
        "R_ar_max": float(np.max(R_ar)),
        "loglog_slope_R_ar_vs_T": slope,
        "arith_decreases_end_vs_start": bool(R_ar[-1] < R_ar[0]),
        "arith_below_defect_all": bool(np.all(R_ar < R_def - 1e-6)),
        "arith_above_critical_line_all": bool(np.all(R_ar > R_cl)),
        "mean_gap_def_minus_ar": float(np.mean(R_def - R_ar)),
        "mean_gap_ar_minus_cl": float(np.mean(R_ar - R_cl)),
        "max_T": float(T_arr[-1]),
        "max_x_max": float(rows[-1]["x_max"]),
    }

    if reading["arith_decreases_end_vs_start"] and slope < -0.2:
        narrative = (
            "Arithmetic R_d trends downward over the scanned T range "
            f"(log-log slope ≈ {slope:.2f}), remaining below the persistent-defect "
            "floor and above the pure critical-line mode — consistent with partial "
            "A0-like behavior but not a proof of full Theorem A."
        )
    elif reading["arith_below_defect_all"]:
        narrative = (
            "Arithmetic R_d stays below the persistent-defect control but does not "
            "show clean T^{-2}-style decay over this window (plateau / slow drift). "
            "Staircase and intermediate spectral mass still dominate at max T="
            f"{reading['max_T']:.1f}."
        )
    else:
        narrative = (
            "Arithmetic R_d does not cleanly separate from controls at the T reached; "
            "see table. Do not claim residual decay."
        )

    summary = {
        "degree": args.degree,
        "n_points": args.n_points,
        "eps": args.eps,
        "detrend": args.detrend,
        "smooth": args.smooth,
        "T_values": T_values,
        "workers": workers,
        "n_primes_global": int(primes.size),
        "x_max_global": float(x_max_global),
        "rows": rows,
        "reading": reading,
        "narrative": narrative,
        "not_an_RH_proof": True,
        "note": (
            "Arithmetic residual uses θ(x)-x over √x on u=log(x)/T with detrend. "
            "A0 pure-mode decay is separate (proved). Full Theorem A under RH remains open."
        ),
    }

    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / "arithmetic_multi_T.json"
    txt_path = results_dir / "arithmetic_multi_T.txt"
    json_path.write_text(json.dumps(summary, indent=2))

    lines = [
        "PBSS arithmetic residual multi-T scan",
        f"d={args.degree} n_points={args.n_points} detrend={args.detrend} "
        f"eps={args.eps} workers={workers}",
        f"primes ≤ {x_max_global:.3e}: n={primes.size}",
        "",
        f"{'T':>6} {'x_max':>12} {'R_ar':>12} {'R_cl':>12} {'R_def':>12} "
        f"{'def-ar':>12} {'ar-cl':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['T']:6.2f} {r['x_max']:12.3e} "
            f"{r['arithmetic']['R_d']:12.4e} "
            f"{r['critical_line']['R_d']:12.4e} "
            f"{r['persistent_defect']['R_d']:12.4e} "
            f"{r['gap_def_minus_arith']:12.4e} "
            f"{r['gap_arith_minus_cl']:12.4e}"
        )
    lines += [
        "",
        f"R_ar: {reading['R_ar_first']:.4e} → {reading['R_ar_last']:.4e}  "
        f"loglog_slope={reading['loglog_slope_R_ar_vs_T']:.3f}",
        f"arith_below_defect_all={reading['arith_below_defect_all']}  "
        f"arith_decreases={reading['arith_decreases_end_vs_start']}",
        "",
        narrative,
        "",
        "Not an unconditional proof of RH or of full Theorem A.",
    ]
    text = "\n".join(lines) + "\n"
    txt_path.write_text(text)
    print(text)

    if args.scratch:
        scratch = Path(args.scratch)
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "arithmetic_multi_T.json").write_text(json.dumps(summary, indent=2))
        (scratch / "arithmetic_multi_T.txt").write_text(text)
        print(f"Wrote scratch copies under {scratch}", flush=True)


if __name__ == "__main__":
    main()
