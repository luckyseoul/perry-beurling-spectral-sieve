#!/usr/bin/env python3
"""
Arithmetic residual × zero-peel multi-(T, N) campaign.

For each T and strip depth N:
  build q_arith = (θ-x)/√x (detrend), strip first N EF modes (optional fit α),
  record R_d.

Uses existing prime checkpoint (default 1e10). Not an RH proof.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

AZP_PRIMES = None
AZP_CSUM = None
AZP_N_POINTS = 16384
AZP_FIT = True


def _peel_job(payload: dict) -> dict:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"
    from pbss.probes import arithmetic_zero_peel, sample_grid
    from pbss.projection import energy_ratio

    u = sample_grid(AZP_N_POINTS)
    N = int(payload["N"])
    q, T_out, meta = arithmetic_zero_peel(
        u,
        T=float(payload["T"]),
        primes=AZP_PRIMES,
        csum=AZP_CSUM,
        n_strip=N,
        detrend=str(payload["detrend"]),
        fit_scale=AZP_FIT if N > 0 else False,
    )
    r = float(energy_ratio(q, u, degree=int(payload["degree"])))
    return {
        "T": float(T_out),
        "N": N,
        "degree": int(payload["degree"]),
        "detrend": payload["detrend"],
        "R_d": r,
        "mode_scale": meta.get("mode_scale", 0.0),
        "x_max": meta.get("x_max"),
        "n_primes": meta.get("n_primes"),
        "fit_scale": meta.get("fit_scale", False),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prime-dir", type=str, default=str(ROOT / "results" / "prime_checkpoints"))
    ap.add_argument("--x-max", type=float, default=1e10)
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "arith_zero_peel"))
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--n-points", type=int, default=16384)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--degrees", type=str, default="2,4,6")
    ap.add_argument("--detrends", type=str, default="deg1,none")
    ap.add_argument("--T-list", type=str, default="10,12,14,16,18,20,22,23")
    ap.add_argument("--N-list", type=str, default="0,1,2,5,10,20,30,50")
    ap.add_argument("--fit-scale", action="store_true", default=True)
    ap.add_argument("--no-fit-scale", action="store_true")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    fit_scale = bool(args.fit_scale) and not bool(args.no_fit_scale)
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    N_values = [int(x) for x in args.N_list.split(",") if x.strip()]
    degrees = [int(x) for x in args.degrees.split(",") if x.strip()] or [args.degree]
    detrends = [x.strip() for x in args.detrends.split(",") if x.strip()]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    from pbss.primes_io import ensure_primes
    from pbss.probes import prime_log_cumsum

    global AZP_PRIMES, AZP_CSUM, AZP_N_POINTS, AZP_FIT
    t0 = time.time()
    workers = max(1, (os.cpu_count() or 4) - 2)
    print(
        f"arith zero-peel x_max={args.x_max:.3e} T={T_values} N={N_values} "
        f"d={degrees} fit_scale={fit_scale} workers={workers}",
        flush=True,
    )
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    primes, pmeta = ensure_primes(int(args.x_max), args.prime_dir, workers=workers)
    AZP_PRIMES = primes
    AZP_CSUM = prime_log_cumsum(primes)
    AZP_N_POINTS = int(args.n_points)
    AZP_FIT = fit_scale

    payloads = []
    for T in T_values:
        if T > np.log(args.x_max) + 1e-9:
            continue
        for N in N_values:
            for det in detrends:
                for d in degrees:
                    payloads.append(
                        {
                            "T": float(T),
                            "N": int(N),
                            "detrend": det,
                            "degree": int(d),
                        }
                    )

    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(
        max_workers=min(workers, max(1, len(payloads))), mp_context=ctx
    ) as ex:
        rows = list(ex.map(_peel_job, payloads, chunksize=1))
    for row in sorted(rows, key=lambda z: (z["T"], z["N"], z["detrend"], z["degree"])):
        print(
            f"T={row['T']:6.2f} N={row['N']:3d} d={row['degree']} {row['detrend']:5s} "
            f"R_d={row['R_d']:.4e} α={row['mode_scale']:.4e}",
            flush=True,
        )

    elapsed = time.time() - t0
    # focus table d=4 deg1
    focus = [
        r
        for r in rows
        if r["degree"] == (4 if 4 in degrees else degrees[0]) and r["detrend"] == "deg1"
    ]
    summary = {
        "status": "completed",
        "elapsed_s": elapsed,
        "x_max": float(args.x_max),
        "n_primes": int(pmeta.get("n_primes", 0)),
        "T_values": T_values,
        "N_values": N_values,
        "degrees": degrees,
        "detrends": detrends,
        "fit_scale": fit_scale,
        "n_points": args.n_points,
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "interpretation": (
            "Arithmetic residual with truncated CL modes peeled. "
            "Whether R_d drops with N is diagnostic only; full Theorem A / RH open."
        ),
        "prime_meta": {k: pmeta[k] for k in pmeta if k in ("x_max", "n_primes", "path", "method")},
    }
    (out_dir / "arith_zero_peel.json").write_text(json.dumps(summary, indent=2))
    lines = [
        "PBSS arithmetic zero-peel campaign",
        f"status=completed elapsed_s={elapsed:.1f} x_max={args.x_max:.3e}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
        "",
        f"{'T':>8} {'N':>4} {'d':>3} {'det':>6} {'R_d':>12} {'alpha':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['T']:8.2f} {r['N']:4d} {r['degree']:3d} {r['detrend']:>6} "
            f"{r['R_d']:12.4e} {r['mode_scale']:12.4e}"
        )
    lines.append("")
    lines.append(summary["interpretation"])
    text = "\n".join(lines) + "\n"
    (out_dir / "arith_zero_peel.txt").write_text(text)
    print(text, flush=True)

    if not args.no_plot and focus:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
            for T in sorted(set(r["T"] for r in focus)):
                rs = sorted([r for r in focus if r["T"] == T], key=lambda z: z["N"])
                ax.plot(
                    [r["N"] for r in rs],
                    [r["R_d"] for r in rs],
                    "o-",
                    label=f"T={T:g}",
                    ms=4,
                )
            ax.set_xlabel("N (modes stripped)")
            ax.set_ylabel(r"$R_d$ arithmetic peel")
            ax.set_title("Arithmetic zero-peel (not an RH proof)")
            ax.grid(True, alpha=0.35)
            ax.legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(out_dir / "arith_zero_peel_Rd_vs_N.png", bbox_inches="tight")
            plt.close()
        except Exception as exc:
            print(f"plot skipped: {exc}", flush=True)

    if args.scratch:
        import shutil

        sc = Path(args.scratch)
        sc.mkdir(parents=True, exist_ok=True)
        for name in (
            "arith_zero_peel.json",
            "arith_zero_peel.txt",
            "arith_zero_peel_Rd_vs_N.png",
        ):
            p = out_dir / name
            if p.exists():
                shutil.copy(p, sc / name)

    print(f"GRAND_ARITH_PEEL_COMPLETE elapsed_s={elapsed:.1f}", flush=True)


if __name__ == "__main__":
    main()
