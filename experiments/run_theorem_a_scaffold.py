#!/usr/bin/env python3
"""
Theorem-A scaffolding campaigns: weight-class endpoints + truncation remainder.

Phases (resume via PHASE_* stamps under --out-dir):
  weight     — multi-T bulk vs weighted R_d on model CL / truncated EF residuals
  remainder  — multi-(T,N) peel-via-remainder + M5 / tail majorants
  arithmetic — optional light arithmetic residual weight/remainder if primes given

Not an RH proof. Not full Theorem A. Smoke: separate --out-dir.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "src"))


def _env1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _stamp(out: Path, name: str, payload: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / f"PHASE_{name}_COMPLETE").write_text(json.dumps(payload, indent=2) + "\n")
    sp = out / "scaffold_state.json"
    state = json.loads(sp.read_text()) if sp.exists() else {"phases": {}}
    state.setdefault("phases", {})[name] = payload
    state["updated_unix"] = time.time()
    state["banner"] = "NOT AN UNCONDITIONAL PROOF OF RH"
    sp.write_text(json.dumps(state, indent=2))


def _done(out: Path, name: str) -> bool:
    return (out / f"PHASE_{name}_COMPLETE").exists()


def _weight_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import (
        explicit_formula_residual,
        probe_critical_line_mode,
        sample_grid,
    )
    from pbss.weights import bulk_vs_weighted_report

    u = sample_grid(int(payload["n_points"]))
    T = float(payload["T"])
    kind = payload["kind"]
    if kind == "cl":
        q = probe_critical_line_mode(u, T=T)
    elif kind == "ef":
        q, _, _ = explicit_formula_residual(
            u, T=T, n_zeros=int(payload["n_zeros"]), form="cos"
        )
    else:
        raise ValueError(kind)
    rep = bulk_vs_weighted_report(
        q,
        u,
        degree=int(payload["degree"]),
        alpha=float(payload["alpha"]),
        weight_name=str(payload["weight_name"]),
    )
    rep["T"] = T
    rep["kind"] = kind
    rep["n_zeros"] = payload.get("n_zeros")
    return rep


def phase_weight(out: Path, args) -> dict:
    name = "WEIGHT"
    phase_dir = out / "weight"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    payloads = []
    for T in T_values:
        for kind in ("cl", "ef"):
            for wname in ("flat", "tukey", "hanning"):
                payloads.append(
                    {
                        "T": T,
                        "kind": kind,
                        "weight_name": wname,
                        "alpha": args.alpha,
                        "degree": args.degree,
                        "n_points": args.n_points,
                        "n_zeros": args.n_zeros,
                    }
                )
    workers = max(1, min(int(args.workers) or 1, len(payloads), 32))
    print(f"[weight] payloads={len(payloads)} workers={workers}", flush=True)
    rows = []
    if workers <= 1:
        for pl in payloads:
            rows.append(_weight_job(pl))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(_weight_job, payloads, chunksize=1))
    elapsed = time.time() - t0
    summary = {
        "class": "weight",
        "hypothesis": (
            "Endpoint zones inflate low-degree mass; admissible Tukey/Hanning weights "
            "should lower R_d vs flat on bulk-sensitive residuals while preserving "
            "M3-style decay for pure CL modes."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "T_values": T_values,
        "alpha": args.alpha,
        "degree": args.degree,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    (phase_dir / "weight.json").write_text(json.dumps(summary, indent=2))
    (phase_dir / "weight.txt").write_text(
        f"class=weight\nelapsed_s={elapsed}\nn_rows={len(rows)}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    _stamp(out, name, {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows)})
    print(f"[weight] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def phase_remainder(out: Path, args) -> dict:
    name = "REMAINDER"
    phase_dir = out / "remainder"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    from pbss.remainder import multi_TN_remainder_scan

    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    n_strips = [int(x) for x in args.n_strips.split(",") if x.strip()]
    print(
        f"[remainder] T={len(T_values)} strips={n_strips} n_full={args.n_zeros}",
        flush=True,
    )
    rows = multi_TN_remainder_scan(
        T_values=T_values,
        n_full=int(args.n_zeros),
        n_strips=n_strips,
        degree=int(args.degree),
        n_points=int(args.n_points),
    )
    elapsed = time.time() - t0
    summary = {
        "class": "remainder",
        "hypothesis": (
            "Finite truncation remainder after stripping N modes obeys M5 O(T^{-2}) "
            "on the stripped block; full strip of q^{(N)} leaves ~0 residual; tail "
            "majorant is scaffolding only for zeros beyond N."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "T_values": T_values,
        "n_strips": n_strips,
        "n_full": int(args.n_zeros),
        "degree": int(args.degree),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    (phase_dir / "remainder.json").write_text(json.dumps(summary, indent=2))
    (phase_dir / "remainder.txt").write_text(
        f"class=remainder\nelapsed_s={elapsed}\nn_rows={len(rows)}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    _stamp(out, name, {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows)})
    print(f"[remainder] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def phase_arithmetic_light(out: Path, args) -> dict:
    """Optional light arithmetic weight/remainder if primes_path set."""
    name = "ARITHMETIC"
    phase_dir = out / "arithmetic"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if not args.primes_path or not Path(args.primes_path).exists():
        summary = {
            "class": "arithmetic",
            "status": "skipped",
            "reason": "no primes_path",
            "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
            "elapsed_s": 0.0,
            "n_rows": 0,
            "rows": [],
        }
        (phase_dir / "arithmetic.json").write_text(json.dumps(summary, indent=2))
        _stamp(out, name, {"status": "skipped", "elapsed_s": 0.0, "n_rows": 0})
        print("[arithmetic] skipped (no primes)", flush=True)
        return summary

    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.remainder import remainder_diagnostic_from_q
    from pbss.weights import bulk_vs_weighted_report

    primes = np.load(args.primes_path, mmap_mode="r")
    csum = (
        np.load(args.csum_path, mmap_mode="r")
        if args.csum_path and Path(args.csum_path).exists()
        else None
    )
    u = sample_grid(int(args.n_points))
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    # keep T such that exp(T) is within primes
    p_max = float(primes[-1])
    rows = []
    for T in T_values:
        if np.exp(T) > p_max * 0.99:
            continue
        q, T_out, meta = arithmetic_residual(
            u, T=T, primes=primes, csum=csum, detrend="deg1"
        )
        wrep = bulk_vs_weighted_report(
            q, u, degree=args.degree, alpha=args.alpha, weight_name="tukey"
        )
        wrep["T"] = float(T_out)
        wrep["source"] = "arithmetic_deg1"
        rem = remainder_diagnostic_from_q(
            q, u, T=float(T_out), n_strip=min(10, args.n_zeros), degree=args.degree
        )
        rem["source"] = "arithmetic_deg1"
        rows.append({"weight": wrep, "remainder": rem})
    elapsed = time.time() - t0
    summary = {
        "class": "arithmetic",
        "status": "completed",
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": "Light arithmetic diagnostic only; not full Theorem A.",
    }
    (phase_dir / "arithmetic.json").write_text(json.dumps(summary, indent=2))
    (phase_dir / "arithmetic.txt").write_text(
        f"class=arithmetic\nelapsed_s={elapsed}\nn_rows={len(rows)}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    _stamp(
        out, name, {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows)}
    )
    print(f"[arithmetic] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default="results/theorem_a_scaffold")
    ap.add_argument(
        "--phases",
        type=str,
        default="weight,remainder,arithmetic",
        help="comma list: weight,remainder,arithmetic",
    )
    ap.add_argument("--T-list", type=str, default="12,16,20,24,28,32")
    ap.add_argument("--n-strips", type=str, default="0,1,2,5,10,15,20")
    ap.add_argument("--n-zeros", type=int, default=20)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--primes-path", type=str, default="")
    ap.add_argument("--csum-path", type=str, default="")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if args.workers <= 0:
        args.workers = max(1, (os.cpu_count() or 4) - 2)

    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)
    print(f"theorem_a_scaffold out={out} workers~{args.workers}", flush=True)

    wanted = [p.strip() for p in args.phases.split(",") if p.strip()]
    results = {}
    runners = {
        "weight": phase_weight,
        "remainder": phase_remainder,
        "arithmetic": phase_arithmetic_light,
    }
    stamp_map = {
        "weight": "WEIGHT",
        "remainder": "REMAINDER",
        "arithmetic": "ARITHMETIC",
    }
    t0 = time.time()
    for p in wanted:
        st = stamp_map[p]
        if _done(out, st):
            print(f"[{p}] skip stamp {st}", flush=True)
            results[p] = json.loads((out / f"PHASE_{st}_COMPLETE").read_text())
            continue
        results[p] = runners[p](out, args)

    elapsed = time.time() - t0
    final = {
        "status": "completed",
        "elapsed_s": elapsed,
        "phases": wanted,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": "Scaffolding for Theorem A — weight class + truncation remainder. Full A open.",
    }
    (out / "SCAFFOLD_COMPLETE").write_text(json.dumps(final, indent=2) + "\n")
    (out / "scaffold_summary.json").write_text(json.dumps(final, indent=2) + "\n")
    print(json.dumps(final, indent=2), flush=True)


if __name__ == "__main__":
    main()
