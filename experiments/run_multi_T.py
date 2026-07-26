#!/usr/bin/env python3
"""
Multi-T campaign: RH-like critical-line mode vs persistent defect vs off-critical.

Uses shipped projection only. Parallel over T with ProcessPool (OMP=1 per worker).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
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

    from pbss.lemmas import predicted_R_d_critical_scaling
    from pbss.probes import (
        probe_critical_line_mode,
        probe_off_critical_mode,
        probe_persistent_defect,
        sample_grid,
    )
    from pbss.projection import project

    T = float(payload["T"])
    degree = int(payload["degree"])
    n_points = int(payload["n_points"])
    eps = float(payload["eps"])
    sigma = float(payload["sigma"])

    u = sample_grid(n_points)
    r_cl = project(probe_critical_line_mode(u, T=T), u, degree=degree, T=T)
    r_def = project(
        probe_persistent_defect(u, eps=eps, j=0, waves=100),
        u,
        degree=degree,
        T=T,
    )
    r_off = project(
        probe_off_critical_mode(u, T=T, sigma=sigma),
        u,
        degree=degree,
        T=T,
    )
    return {
        "T": T,
        "critical_line": {
            "R_d": r_cl.energy_ratio,
            "S_d": r_cl.scaled_strength,
            "R0_theory": predicted_R_d_critical_scaling(T),
        },
        "persistent_defect": {
            "R_d": r_def.energy_ratio,
            "S_d": r_def.scaled_strength,
            "eps2": eps * eps,
        },
        "off_critical": {
            "R_d": r_off.energy_ratio,
            "S_d": r_off.scaled_strength,
            "sigma": sigma,
        },
        "gap_defect_minus_cl": r_def.energy_ratio - r_cl.energy_ratio,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=8192)
    ap.add_argument("--eps", type=float, default=0.5, help="persistent defect mass")
    ap.add_argument("--sigma", type=float, default=0.75, help="off-line real part")
    ap.add_argument(
        "--T-list",
        type=str,
        default="3,5,8,12,16,20,28,40,56,80",
        help="comma-separated T values",
    )
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()

    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)
    payloads = [
        {
            "T": T,
            "degree": args.degree,
            "n_points": args.n_points,
            "eps": args.eps,
            "sigma": args.sigma,
        }
        for T in T_values
    ]

    print(
        f"multi-T scan d={args.degree} n={args.n_points} workers={workers} "
        f"T={T_values}",
        flush=True,
    )
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for row in ex.map(_worker, payloads):
            rows.append(row)
            print(
                f"T={row['T']:7.2f}  R_cl={row['critical_line']['R_d']:.4e}  "
                f"R_def={row['persistent_defect']['R_d']:.4e}  "
                f"R_off={row['off_critical']['R_d']:.4e}  "
                f"gap={row['gap_defect_minus_cl']:.4e}",
                flush=True,
            )

    # qualitative checks for the writeup
    r_cl = np.array([r["critical_line"]["R_d"] for r in rows])
    r_def = np.array([r["persistent_defect"]["R_d"] for r in rows])
    T_arr = np.array([r["T"] for r in rows])
    # last half: cl should be below def and cl should trend down vs first points
    mid = len(rows) // 2
    gap_ok = bool(np.all(r_def[mid:] - r_cl[mid:] > 0.15))
    decay_ok = bool(r_cl[-1] < r_cl[0] * 0.5)
    flat_def_ok = bool(np.std(r_def) < 0.08)

    summary = {
        "degree": args.degree,
        "n_points": args.n_points,
        "eps": args.eps,
        "sigma": args.sigma,
        "T_values": T_values,
        "workers": workers,
        "rows": rows,
        "qualitative": {
            "gap_defect_minus_cl_large_T_ok": gap_ok,
            "critical_line_decays_ok": decay_ok,
            "persistent_defect_flat_ok": flat_def_ok,
            "R_cl_first": float(r_cl[0]),
            "R_cl_last": float(r_cl[-1]),
            "R_def_mean": float(np.mean(r_def)),
            "min_gap_large_T": float(np.min(r_def[mid:] - r_cl[mid:])),
        },
        "interpretation": (
            "Supports model Theorem A0 (R_d of critical-line mode → 0) and "
            "model B0 (persistent defect keeps R_d≈ε²). Not an RH proof."
        ),
    }

    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / "multi_T_scan.json"
    txt_path = results_dir / "multi_T_scan.txt"
    json_path.write_text(json.dumps(summary, indent=2))

    lines = [
        "PBSS multi-T scan (shipped projection)",
        f"d={args.degree} n_points={args.n_points} eps={args.eps} sigma={args.sigma}",
        "",
        f"{'T':>8}  {'R_cl':>12}  {'R0_thy':>12}  {'R_def':>12}  {'R_off':>12}  {'gap':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['T']:8.2f}  {r['critical_line']['R_d']:12.4e}  "
            f"{r['critical_line']['R0_theory']:12.4e}  "
            f"{r['persistent_defect']['R_d']:12.4e}  "
            f"{r['off_critical']['R_d']:12.4e}  "
            f"{r['gap_defect_minus_cl']:12.4e}"
        )
    q = summary["qualitative"]
    lines += [
        "",
        f"critical_line_decays_ok={q['critical_line_decays_ok']}  "
        f"R_cl {q['R_cl_first']:.4e} → {q['R_cl_last']:.4e}",
        f"persistent_defect_flat_ok={q['persistent_defect_flat_ok']}  "
        f"mean R_def={q['R_def_mean']:.4e}",
        f"gap_ok={q['gap_defect_minus_cl_large_T_ok']}  "
        f"min_gap_large_T={q['min_gap_large_T']:.4e}",
        "",
        "Not an unconditional proof of RH. Model A0/B0 support only.",
    ]
    text = "\n".join(lines) + "\n"
    txt_path.write_text(text)
    print(text)

    if args.scratch:
        scratch = Path(args.scratch)
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "multi_T_scan.json").write_text(json.dumps(summary, indent=2))
        (scratch / "multi_T_scan.txt").write_text(text)
        print(f"Wrote scratch copies under {scratch}", flush=True)

    if not (gap_ok and decay_ok):
        raise SystemExit("multi-T qualitative separation failed")


if __name__ == "__main__":
    main()
