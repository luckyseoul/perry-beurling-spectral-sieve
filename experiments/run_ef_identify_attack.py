#!/usr/bin/env python3
"""
Multi-hypothesis EF identification attack on arithmetic residual.

Compares H_theta_sqrt, H_psi_sqrt, H_psi_x, H_theta_x across multi-T × multi-N.
Outputs remainder metrics and a sharp block if all hypotheses fail mode capture.

Not a public RH announcement — technical attack log.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.ef_identify import (  # noqa: E402
    HYPOTHESES,
    model_sanity_identify,
    multi_hypothesis_scan,
    summarize_attack,
)
from pbss.probes import prime_log_cumsum  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default="results/ef_identify_attack")
    ap.add_argument(
        "--primes-path",
        type=str,
        default="results/prime_checkpoints/primes_le_10000000000.npy",
    )
    ap.add_argument(
        "--T-list",
        type=str,
        default="12,14,16,18,20",
    )
    ap.add_argument("--n-zeros-list", type=str, default="5,10,20,30")
    ap.add_argument("--n-points", type=int, default=2048)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--detrend", type=str, default="deg1")
    ap.add_argument(
        "--hypotheses",
        type=str,
        default="H_theta_sqrt,H_psi_sqrt,H_psi_x,H_theta_x",
    )
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("EF identification attack — technical campaign", flush=True)

    # Model sanity first
    model_rows = [
        model_sanity_identify(T=T, n_zeros=10, n_points=args.n_points, degree=args.degree)
        for T in (15.0, 25.0, 40.0)
    ]
    assert all(r["model_identity_ok"] for r in model_rows)

    primes_path = Path(args.primes_path)
    if not primes_path.exists():
        alt = ROOT / "results/prime_checkpoints/primes_le_50000000000.npy"
        if alt.exists():
            primes_path = alt
        else:
            raise SystemExit(f"missing primes {args.primes_path}")

    primes = np.load(primes_path, mmap_mode="r")
    p_max = float(primes[-1])
    print(f"primes path={primes_path} n={primes.size} p_max={p_max:.3e}", flush=True)

    # precompute theta csum for speed (may be large but once)
    t0 = time.time()
    # only need prefix for max T
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    T_values = [T for T in T_values if np.exp(T) <= p_max * 0.99]
    max_x = int(np.floor(np.exp(max(T_values))))
    hi = int(np.searchsorted(primes, max_x, side="right"))
    print(f"building theta csum prefix hi={hi} for max_x={max_x}", flush=True)
    p_pref = np.array(primes[:hi], dtype=np.int64, copy=True)
    csum = prime_log_cumsum(p_pref)
    print(f"csum ready in {time.time()-t0:.1f}s", flush=True)

    n_zeros_list = [int(x) for x in args.n_zeros_list.split(",") if x.strip()]
    hyps = [h.strip() for h in args.hypotheses.split(",") if h.strip()]

    t1 = time.time()
    rows = multi_hypothesis_scan(
        T_values=T_values,
        n_zeros_list=n_zeros_list,
        primes=p_pref,
        hypotheses=hyps,
        n_points=args.n_points,
        degree=args.degree,
        detrend=args.detrend,
        csum_theta=csum,
    )
    elapsed = time.time() - t1
    summary = summarize_attack(rows)
    summary["elapsed_s"] = elapsed
    summary["n_rows"] = len(rows)
    summary["T_values"] = T_values
    summary["n_zeros_list"] = n_zeros_list
    summary["hypotheses"] = hyps
    summary["primes_path"] = str(primes_path)
    summary["x_max_used"] = float(max_x)
    summary["model_sanity"] = model_rows
    summary["detrend"] = args.detrend

    # compact rows for JSON (drop huge nothing — already compact)
    payload = {
        **summary,
        "rows": rows,
    }
    (out / "ef_identify_attack.json").write_text(json.dumps(payload, indent=2))

    lines = [
        f"n_rows={len(rows)} elapsed_s={elapsed:.1f}",
        f"T={T_values} N={n_zeros_list}",
        f"best={summary.get('best_hypothesis')} score_Ed={summary.get('best_score_Ed_rem_over_l2q')}",
        f"status={summary.get('status')}",
    ]
    for h, m in summary.get("hypothesis_means", {}).items():
        lines.append(
            f"  {h}: corr={m['mean_abs_corr']:.4f} frac_l2={m['mean_frac_l2']:.4f} "
            f"Ed_rem/l2q={m['mean_Ed_rem_over_l2q']:.4f}"
        )
    if summary.get("sharp_block"):
        sb = summary["sharp_block"]
        lines.append(f"SHARP_BLOCK={sb['name']}")
        lines.append(sb["statement"][:500])
        lines.append("UNBLOCK: " + sb["what_would_unblock"][:400])
    lines.append("technical attack log — not an RH announcement")
    (out / "ef_identify_attack.txt").write_text("\n".join(lines) + "\n")
    (out / "EF_IDENTIFY_ATTACK_COMPLETE").write_text(
        json.dumps(
            {
                "status": summary.get("status"),
                "n_rows": len(rows),
                "elapsed_s": elapsed,
                "sharp_block": summary.get("sharp_block", {}).get("name")
                if summary.get("sharp_block")
                else None,
                "best_hypothesis": summary.get("best_hypothesis"),
            },
            indent=2,
        )
        + "\n"
    )
    print("\n".join(lines), flush=True)
    print(json.dumps({"status": summary.get("status"), "n_rows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
