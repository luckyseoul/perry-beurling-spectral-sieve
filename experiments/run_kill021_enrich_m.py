#!/usr/bin/env python3
"""
Kill 0.21: H_theta_sqrt fixed; enrich m only; multi-N Ed(r)/||q||² table.
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
    M_ENRICHMENTS,
    multi_N_enrich_scan,
    summarize_enrich_kill021,
)
from pbss.probes import prime_log_cumsum  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default="results/kill021_enrich_m")
    ap.add_argument(
        "--primes-path",
        type=str,
        default="results/prime_checkpoints/primes_le_10000000000.npy",
    )
    ap.add_argument("--T-list", type=str, default="14,16,18,20")
    ap.add_argument("--n-zeros-list", type=str, default="5,10,20,40")
    ap.add_argument(
        "--enrichments",
        type=str,
        default=",".join(M_ENRICHMENTS),
    )
    ap.add_argument("--n-points", type=int, default=2048)
    ap.add_argument("--degree", type=int, default=4)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("kill 0.21 — enrich m only, H_theta_sqrt", flush=True)

    primes_path = Path(args.primes_path)
    if not primes_path.exists():
        raise SystemExit(f"missing {primes_path}")
    primes = np.load(primes_path, mmap_mode="r")
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    T_values = [T for T in T_values if np.exp(T) <= float(primes[-1]) * 0.99]
    n_list = [int(x) for x in args.n_zeros_list.split(",") if x.strip()]
    enrichments = [e.strip() for e in args.enrichments.split(",") if e.strip()]

    max_x = int(np.floor(np.exp(max(T_values))))
    hi = int(np.searchsorted(primes, max_x, side="right"))
    p = np.array(primes[:hi], dtype=np.int64, copy=True)
    csum = prime_log_cumsum(p)
    print(f"primes hi={hi} max_x={max_x} T={T_values} N={n_list}", flush=True)

    t0 = time.time()
    all_rows = []
    for T in T_values:
        rows = multi_N_enrich_scan(
            T=T,
            n_zeros_list=n_list,
            primes=p,
            enrichments=enrichments,
            n_points=args.n_points,
            degree=args.degree,
            detrend="deg1",
            csum_theta=csum,
        )
        all_rows.extend(rows)
        print(f"  T={T} rows={len(rows)}", flush=True)
    elapsed = time.time() - t0
    summary = summarize_enrich_kill021(all_rows)
    summary["elapsed_s"] = elapsed
    summary["T_values"] = T_values
    summary["n_zeros_list"] = n_list
    summary["enrichments"] = enrichments
    summary["primes_path"] = str(primes_path)
    summary["degree"] = args.degree
    summary["rows"] = all_rows

    (out / "kill021_enrich_m.json").write_text(json.dumps(summary, indent=2))
    lines = [
        f"outcome={summary['outcome']} elapsed_s={elapsed:.1f} n_rows={len(all_rows)}",
        f"baseline zeros mean Ed={summary['means_Ed_r_over_l2q'].get('zeros', float('nan')):.4f}",
        f"best={summary['best_enrich']} mean Ed={summary['best_mean_Ed']:.4f}",
        "by enrichment mean Ed(r)/||q||²:",
    ]
    for e, m in sorted(summary["means_Ed_r_over_l2q"].items(), key=lambda x: x[1]):
        lines.append(f"  {e}: {m:.4f}")
    lines.append("by N (zeros vs best):")
    zN = summary["by_N_Ed"].get("zeros", {})
    bN = summary["by_N_Ed"].get(summary["best_enrich"], {})
    for N in sorted(set(zN) | set(bN)):
        lines.append(
            f"  N={N}: zeros={zN.get(N, float('nan')):.4f} "
            f"{summary['best_enrich']}={bN.get(N, float('nan')):.4f}"
        )
    if summary.get("sharp_block"):
        sb = summary["sharp_block"]
        lines.append(f"SHARP_BLOCK={sb['name']}")
        lines.append(sb["statement"][:600])
    lines.append("enrich-m-only — H_theta_sqrt fixed residual")
    (out / "kill021_enrich_m.txt").write_text("\n".join(lines) + "\n")
    (out / "KILL021_COMPLETE").write_text(
        json.dumps(
            {
                "outcome": summary["outcome"],
                "best_enrich": summary["best_enrich"],
                "best_mean_Ed": summary["best_mean_Ed"],
                "baseline_zeros": summary["means_Ed_r_over_l2q"].get("zeros"),
                "sharp_block": (summary.get("sharp_block") or {}).get("name"),
                "n_rows": len(all_rows),
                "elapsed_s": elapsed,
            },
            indent=2,
        )
        + "\n"
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
