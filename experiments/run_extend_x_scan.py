#!/usr/bin/env python3
"""
Extend prime checkpoint beyond 1e10 and multi-T arithmetic residual scan.

Target largest feasible in {1e11, 1e12} given RAM/disk; document stop reason.
Not an RH proof.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

# Fork-inherited by multi-T residual workers (set in main before pool).
WORKER_PRIMES = None
WORKER_CSUM = None
WORKER_N_POINTS = 16384
WORKER_DEGREES: list = []


def _env1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _arith_group_job(payload: dict) -> list:
    """
    One residual build per (T, detrend); then energy_ratio for all degrees.
    Parallel unit is (T, detrend) — not (T, detrend, d) — to cut work 3×.
    """
    _env1()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection import energy_ratio

    u = sample_grid(WORKER_N_POINTS)
    q, T_out, meta = arithmetic_residual(
        u,
        T=float(payload["T"]),
        primes=WORKER_PRIMES,
        csum=WORKER_CSUM,
        detrend=str(payload["detrend"]),
    )
    rows = []
    for d in WORKER_DEGREES:
        r = float(energy_ratio(q, u, degree=int(d)))
        rows.append(
            {
                "T": float(T_out),
                "x_max_window": float(meta["x_max"]),
                "degree": int(d),
                "detrend": payload["detrend"],
                "R_d": r,
                "n_primes": meta["n_primes"],
            }
        )
    return rows


def avail_ram() -> int:
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except Exception:
        pass
    return 40 * 1024**3


def avail_disk(path: Path) -> int:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def choose_target(candidates: list[float], ram: int, disk: int) -> tuple[float, str]:
    """Pick largest x with est prime storage < 55% RAM and < 45% free disk."""
    notes = []
    for x in sorted(candidates, reverse=True):
        n = float(x)
        n_primes = n / max(math.log(n), 2.0)
        need = n_primes * 8 * 1.25 + 3e9
        ok_ram = need < 0.55 * ram
        ok_disk = need < 0.45 * disk
        notes.append(
            f"x={n:.3e} est_primes~{n_primes:.3e} need~{need/1e9:.1f}GiB "
            f"ram_ok={ok_ram} disk_ok={ok_disk}"
        )
        if ok_ram and ok_disk:
            return n, "; ".join(notes)
    # fallback: smallest candidate still > 1e10 if any
    for x in sorted(candidates):
        if x > 1e10:
            return float(x), "forced_min_above_1e10; " + "; ".join(notes)
    return float(candidates[0]), "fallback; " + "; ".join(notes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prime-dir", type=str, default=str(ROOT / "results" / "prime_checkpoints"))
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "extend_x_scan"))
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument(
        "--xmax-candidates",
        type=str,
        default="1e12,1e11,5e10,2e10",
    )
    ap.add_argument("--force-x-max", type=float, default=0.0, help="override chooser if >0")
    ap.add_argument("--n-points", type=int, default=16384)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--degrees", type=str, default="2,4,6")
    ap.add_argument("--detrends", type=str, default="deg1,none")
    ap.add_argument(
        "--segment-size",
        type=int,
        default=20_000_000,
        help="segment length; smaller → more parallel segments",
    )
    ap.add_argument("--T-list", type=str, default="", help="default auto from log(x_max)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--skip-sieve", action="store_true", help="only scan if checkpoint exists")
    ap.add_argument(
        "--workers",
        type=int,
        default=0,
        help="parallel sieve segments + multi-T residual jobs (default: nproc-2)",
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prime_dir = Path(args.prime_dir)
    candidates = [float(x) for x in args.xmax_candidates.split(",") if x.strip()]
    ram = avail_ram()
    disk = avail_disk(prime_dir if prime_dir.exists() else ROOT)
    if args.force_x_max > 0:
        x_max = float(args.force_x_max)
        choose_note = f"forced x_max={x_max:.3e}"
    else:
        x_max, choose_note = choose_target(candidates, ram, disk)

    print(f"extend-x: chosen x_max={x_max:.3e}", flush=True)
    print(f"chooser: {choose_note}", flush=True)
    print(f"ram={ram/1e9:.1f}GiB disk_free={disk/1e9:.1f}GiB", flush=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    from pbss.primes_io import ensure_primes
    from pbss.probes import prime_log_cumsum, sample_grid

    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)
    print(f"workers={workers} (parallel sieve + multi-T)", flush=True)

    t0 = time.time()
    if args.skip_sieve:
        from pbss.primes_io import load_primes_checkpoint

        primes, pmeta = load_primes_checkpoint(prime_dir, int(x_max), mmap=True)
    else:
        primes, pmeta = ensure_primes(
            int(x_max),
            prime_dir,
            segment_size=args.segment_size,
            extend_from_existing=True,
            workers=workers,
        )
    sieve_s = float(pmeta.get("sieve_seconds", 0.0) or 0.0)
    print(
        f"primes ready n={pmeta.get('n_primes')} method={pmeta.get('method')} "
        f"sieve_s={sieve_s:.1f}",
        flush=True,
    )

    # Residual phase strategy for multi-GB prime tables:
    # 1) Ensure csum on disk
    # 2) Load primes+csum into RAM once (sequential read → page cache)
    # 3) Parallelize over (T, detrend) groups with modest workers (COW-friendly)
    global WORKER_PRIMES, WORKER_CSUM, WORKER_N_POINTS, WORKER_DEGREES
    print("building/loading θ-prefix csum…", flush=True)
    csum_path = out_dir / f"csum_le_{int(x_max)}.npy"
    if not csum_path.exists():
        # stream log+cumsum without full float copy of int64 if possible
        p_arr = np.asarray(primes, dtype=np.int64)
        csum_dense = prime_log_cumsum(p_arr)
        np.save(csum_path, csum_dense)
        del csum_dense, p_arr
    n_primes = int(pmeta.get("n_primes") or getattr(primes, "size", 0))
    # Load into RAM for residual (avoids multi-worker memmap thrash / OOM)
    need_gib = n_primes * 16 / 1e9  # int64 primes + float64 csum
    if need_gib < 0.55 * (ram / 1e9):
        print(
            f"loading primes+csum into RAM (~{need_gib:.1f} GiB) for residual…",
            flush=True,
        )
        WORKER_PRIMES = np.load(
            prime_dir / f"primes_le_{int(x_max)}.npy"
        ).astype(np.int64, copy=False)
        WORKER_CSUM = np.load(csum_path)
    else:
        print(
            f"RAM tight for full load (need~{need_gib:.1f}GiB); using memmap + 4 workers",
            flush=True,
        )
        WORKER_PRIMES = np.load(
            prime_dir / f"primes_le_{int(x_max)}.npy", mmap_mode="r"
        )
        WORKER_CSUM = np.load(csum_path, mmap_mode="r")
    WORKER_N_POINTS = int(args.n_points)

    T_max = float(np.log(x_max))
    if args.T_list.strip():
        T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    else:
        T_values = list(
            np.unique(
                np.round(
                    np.linspace(max(10.0, T_max - 8), T_max, 10),
                    3,
                )
            )
        )
    degrees = [int(x) for x in args.degrees.split(",") if x.strip()]
    detrends = [x.strip() for x in args.detrends.split(",") if x.strip()]
    WORKER_DEGREES = degrees

    payloads = []
    for T in T_values:
        if T > T_max + 1e-9:
            continue
        for det in detrends:
            payloads.append({"T": float(T), "detrend": det})

    # Grouped jobs; fewer workers if arrays memmapped
    if isinstance(WORKER_PRIMES, np.memmap) or (
        hasattr(WORKER_PRIMES, "base") and isinstance(getattr(WORKER_PRIMES, "base", None), np.memmap)
    ):
        n_job_workers = min(4, workers, max(1, len(payloads)))
    else:
        n_job_workers = min(16, workers, max(1, len(payloads)))
    print(
        f"multi-T residual groups={len(payloads)} workers={n_job_workers} "
        f"n_primes={n_primes} degrees={degrees}",
        flush=True,
    )
    ctx = mp.get_context("fork")
    rows = []
    with ProcessPoolExecutor(max_workers=n_job_workers, mp_context=ctx) as ex:
        for group in ex.map(_arith_group_job, payloads, chunksize=1):
            rows.extend(group)
    for row in sorted(rows, key=lambda z: (z["T"], z["detrend"], z["degree"])):
        print(
            f"T={row['T']:7.3f} d={row['degree']} {row['detrend']:5s} "
            f"R_d={row['R_d']:.4e} x_win={row['x_max_window']:.3e}",
            flush=True,
        )

    elapsed = time.time() - t0
    summary = {
        "status": "completed",
        "elapsed_s": elapsed,
        "sieve_seconds": sieve_s,
        "x_max": float(x_max),
        "x_max_gt_1e10": bool(x_max > 1e10),
        "n_primes": int(pmeta.get("n_primes", 0)),
        "choose_note": choose_note,
        "ram_bytes": ram,
        "disk_free_bytes": disk,
        "prime_meta": {
            k: pmeta.get(k)
            for k in ("x_max", "n_primes", "path", "method", "extended_from", "sieve_seconds")
        },
        "T_values": T_values,
        "degrees": degrees,
        "detrends": detrends,
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "interpretation": (
            f"Multi-T arithmetic residual through x_max={x_max:.3e} (>1e10). "
            "Not a proof of RH or full Theorem A."
        ),
    }
    (out_dir / "extend_x_scan.json").write_text(json.dumps(summary, indent=2))
    lines = [
        "PBSS extend-x multi-T scan",
        f"status=completed elapsed_s={elapsed:.1f} x_max={x_max:.3e} n_primes={summary['n_primes']}",
        f"chooser: {choose_note}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
        "",
        f"{'T':>8} {'d':>3} {'det':>6} {'R_d':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['T']:8.3f} {r['degree']:3d} {r['detrend']:>6} {r['R_d']:12.4e}"
        )
    text = "\n".join(lines) + "\n"
    (out_dir / "extend_x_scan.txt").write_text(text)
    print(text, flush=True)

    if not args.no_plot and rows:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            focus = [r for r in rows if r["degree"] == 4 and r["detrend"] == "deg1"]
            focus = sorted(focus, key=lambda z: z["T"])
            fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
            if focus:
                ax.plot([r["T"] for r in focus], [r["R_d"] for r in focus], "o-")
            ax.set_xlabel("T")
            ax.set_ylabel(r"$R_4$ arith deg1")
            ax.set_title(f"Extend-x scan x_max={x_max:.2e} (not RH proof)")
            ax.grid(True, alpha=0.35)
            fig.tight_layout()
            fig.savefig(out_dir / "extend_x_Rd_vs_T.png", bbox_inches="tight")
            plt.close()
        except Exception as exc:
            print(f"plot skipped: {exc}", flush=True)

    if args.scratch:
        import shutil

        sc = Path(args.scratch)
        sc.mkdir(parents=True, exist_ok=True)
        for name in ("extend_x_scan.json", "extend_x_scan.txt", "extend_x_Rd_vs_T.png"):
            p = out_dir / name
            if p.exists():
                shutil.copy(p, sc / name)

    if not summary["x_max_gt_1e10"]:
        raise SystemExit("x_max not extended beyond 1e10")
    print(f"EXTEND_X_COMPLETE elapsed_s={elapsed:.1f} x_max={x_max:.3e}", flush=True)


if __name__ == "__main__":
    main()
