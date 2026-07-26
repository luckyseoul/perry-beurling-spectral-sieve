#!/usr/bin/env python3
"""
Multi-T arithmetic residual under admissible weights (goal 4).

Compares raw deg1 R_d, E_end, and weighted R_d (tukey/hanning) on the largest
available prime table. Resume via PHASE stamps. Multi-core over T jobs.

Not an RH proof. Not A0 collapse claim.
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
    sp = out / "arith_weight_state.json"
    state = json.loads(sp.read_text()) if sp.exists() else {"phases": {}}
    state.setdefault("phases", {})[name] = payload
    state["updated_unix"] = time.time()
    state["banner"] = "NOT AN UNCONDITIONAL PROOF OF RH"
    sp.write_text(json.dumps(state, indent=2))


def _done(out: Path, name: str) -> bool:
    return (out / f"PHASE_{name}_COMPLETE").exists()


def _arith_weight_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection import energy_ratio
    from pbss.remainder import remainder_diagnostic_from_q
    from pbss.weights import (
        bulk_vs_weighted_report,
        endpoint_contribution,
    )

    p = np.load(payload["primes_path"], mmap_mode="r")
    c = (
        np.load(payload["csum_path"], mmap_mode="r")
        if payload.get("csum_path")
        else None
    )
    T = float(payload["T"])
    u = sample_grid(int(payload["n_points"]))
    hi = int(np.searchsorted(p, min(float(np.exp(T)), float(p[-1])), side="right"))
    if hi > 350_000_000:
        p_use, c_use = p[:hi], (c[:hi] if c is not None else None)
    else:
        p_use = np.array(p[:hi], dtype=np.int64, copy=True)
        c_use = (
            np.array(c[:hi], dtype=np.float64, copy=True) if c is not None else None
        )
    q, T_out, meta = arithmetic_residual(
        u,
        T=T,
        primes=p_use,
        csum=c_use,
        detrend=str(payload.get("detrend", "deg1")),
        smooth=int(payload.get("smooth", 1)),
    )
    degree = int(payload["degree"])
    alpha = float(payload["alpha"])
    end = endpoint_contribution(q, u, degree=degree, alpha=alpha)
    rows_w = {}
    for wname in payload["weights"]:
        rep = bulk_vs_weighted_report(
            q, u, degree=degree, alpha=alpha, weight_name=wname
        )
        rows_w[wname] = {
            "R_d_weighted": rep["R_d_weighted"],
            "R_d_raw": rep["R_d_raw"],
            "E_end": rep["E_end"],
            "E_bulk": rep["E_bulk"],
        }
    rem = remainder_diagnostic_from_q(
        q,
        u,
        T=float(T_out),
        n_strip=int(payload.get("n_strip", 10)),
        degree=degree,
        fit_scale=True,
    )
    return {
        "T": float(T_out),
        "degree": degree,
        "detrend": payload.get("detrend", "deg1"),
        "R_d_raw": end["R_d"],
        "E_end": end["E_end"],
        "E_bulk": end["E_bulk"],
        "R_d_bulk_only": end["R_d_bulk"],
        "weights": rows_w,
        "remainder_peel": {
            "n_strip": rem["n_strip"],
            "R_d_remainder": rem["R_d_remainder"],
            "alpha": rem["alpha"],
            "tail_majorant_R_d": rem["tail_majorant_R_d"],
        },
        "n_primes": int(meta.get("n_primes", hi)),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default="results/arithmetic_weights")
    ap.add_argument(
        "--primes-path",
        type=str,
        default="results/prime_checkpoints/primes_le_50000000000.npy",
    )
    ap.add_argument(
        "--csum-path",
        type=str,
        default="results/extend_x_scan/csum_le_50000000000.npy",
    )
    ap.add_argument(
        "--T-list",
        type=str,
        default="14,16,18,20,22,24",
        help="T values with exp(T) within prime table",
    )
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--weights", type=str, default="flat,tukey,hanning")
    ap.add_argument("--n-strip", type=int, default=10)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--detrend", type=str, default="deg1")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    primes_path = Path(args.primes_path)
    if not primes_path.exists():
        # fall back to 1e10
        alt = ROOT / "results/prime_checkpoints/primes_le_10000000000.npy"
        if alt.exists():
            primes_path = alt
            print(f"[arith_weights] fallback primes {primes_path}", flush=True)
        else:
            raise SystemExit(f"no primes at {args.primes_path}")

    csum_path = args.csum_path if Path(args.csum_path).exists() else ""
    p = np.load(primes_path, mmap_mode="r")
    p_max = float(p[-1])
    T_values = []
    for x in args.T_list.split(","):
        if not x.strip():
            continue
        T = float(x)
        if np.exp(T) <= p_max * 0.99:
            T_values.append(T)
        else:
            print(f"[arith_weights] skip T={T} (exp(T)>{p_max:.3e})", flush=True)
    if not T_values:
        raise SystemExit("no feasible T values for prime table")

    if _done(out, "ARITH_WEIGHTS"):
        print("[arith_weights] skip stamp ARITH_WEIGHTS", flush=True)
        print(json.dumps(json.loads((out / "PHASE_ARITH_WEIGHTS_COMPLETE").read_text()), indent=2))
        return

    weights = [w.strip() for w in args.weights.split(",") if w.strip()]
    workers = int(args.workers) if args.workers > 0 else max(1, (os.cpu_count() or 4) - 2)
    # mmap safety: cap workers for huge tables
    if p_max >= 1e10:
        workers = min(workers, 6)
    payloads = [
        {
            "T": T,
            "primes_path": str(primes_path),
            "csum_path": csum_path,
            "n_points": args.n_points,
            "degree": args.degree,
            "alpha": args.alpha,
            "weights": weights,
            "n_strip": args.n_strip,
            "detrend": args.detrend,
            "smooth": 1,
        }
        for T in T_values
    ]
    print(
        f"[arith_weights] T={T_values} workers={workers} primes<={p_max:.3e}",
        flush=True,
    )
    t0 = time.time()
    if workers <= 1:
        rows = [_arith_weight_job(pl) for pl in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(_arith_weight_job, payloads, chunksize=1))
    elapsed = time.time() - t0

    summary = {
        "campaign": "arithmetic_weights",
        "hypothesis": (
            "Admissible weights reduce E_end and may lower R_d vs flat on arithmetic "
            "deg1 residual; multi-T scan tests whether weighted R_d decays (A0) or plateaus."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "T_values": T_values,
        "x_max_primes": p_max,
        "primes_path": str(primes_path),
        "weights": weights,
        "degree": args.degree,
        "alpha": args.alpha,
        "workers": workers,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": "Numeric diagnostic only; not full Theorem A or RH.",
    }
    (out / "arithmetic_weights.json").write_text(json.dumps(summary, indent=2))
    # compact txt
    lines = [
        f"campaign=arithmetic_weights n_rows={len(rows)} elapsed_s={elapsed:.1f}",
        f"T={T_values}",
        f"x_max_primes={p_max:.6e}",
    ]
    for r in rows:
        tw = r["weights"].get("tukey", {})
        lines.append(
            f"T={r['T']:.3f} R_raw={r['R_d_raw']:.4f} E_end={r['E_end']:.4f} "
            f"R_tukey={tw.get('R_d_weighted', float('nan')):.4f} "
            f"R_peel={r['remainder_peel']['R_d_remainder']:.4f}"
        )
    lines.append("NOT AN UNCONDITIONAL PROOF OF RH")
    (out / "arithmetic_weights.txt").write_text("\n".join(lines) + "\n")
    _stamp(
        out,
        "ARITH_WEIGHTS",
        {
            "status": "completed",
            "elapsed_s": elapsed,
            "n_rows": len(rows),
            "T_values": T_values,
        },
    )
    (out / "ARITH_WEIGHTS_COMPLETE").write_text(
        json.dumps(
            {
                "status": "completed",
                "elapsed_s": elapsed,
                "n_rows": len(rows),
                "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"[arith_weights] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    print(json.dumps({"status": "completed", "n_rows": len(rows), "elapsed_s": elapsed}, indent=2))


if __name__ == "__main__":
    main()
