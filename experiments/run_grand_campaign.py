#!/usr/bin/env python3
"""
PBSS grand multi-hour campaign with resume support.

Phases:
  A) Segmented sieve + on-disk prime checkpoint (target x_max>=1e10)
  B) Arithmetic residual multi-T × (d, detrend, smooth) ablations
  C) Controls per T: critical-line, off-critical, persistent defect
  D) MC >= 2000 spectral-defect trials per T (default much larger for wall time)
     each trial → R_d for every d in the ablation list

Checkpoints after each phase/T so a restart resumes.
Clean shutdown: workers exit; stamp file written.

No RH proof claim.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import signal
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_MP_CTX = mp.get_context("fork")
_SHUTDOWN = False
# Set in parent before forking workers so children inherit via COW (no multi-GB pickle).
WORKER_PRIMES = None
WORKER_CSUM = None


def _handle_sig(sig, frame):
    global _SHUTDOWN
    _SHUTDOWN = True
    print(f"\nSignal {sig} — will checkpoint and stop after current batch…", flush=True)


def _env1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _state_path(out_dir: Path) -> Path:
    return out_dir / "grand_campaign_state.json"


def _save_state(out_dir: Path, state: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = _state_path(out_dir)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2))
    tmp.replace(p)


def _load_state(out_dir: Path) -> Optional[dict]:
    p = _state_path(out_dir)
    if not p.exists():
        return None
    return json.loads(p.read_text())


def choose_xmax(candidates: List[float], avail_bytes: int) -> tuple[float, str]:
    """Pick largest x_max with estimated prime storage under ~40% of avail RAM."""
    for x in sorted(candidates, reverse=True):
        n = int(x)
        n_primes_est = n / max(math.log(n), 2.0)
        # int64 primes + memmap overhead
        need = n_primes_est * 8 * 1.3 + 2e9
        if need < 0.40 * avail_bytes:
            return float(n), (
                f"x_max={n:.3e} est_primes~{n_primes_est:.3e} "
                f"need~{need/1e9:.1f}GiB avail~{avail_bytes/1e9:.1f}GiB"
            )
    x = float(candidates[-1])
    return x, f"fallback smallest candidate x_max={x:.3e}"


def avail_ram() -> int:
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except Exception:
        pass
    return 40 * 1024**3


# ----- worker jobs (must be top-level for pickling under spawn; fork ok too) -----


def _arith_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection_backend import energy_ratio_auto

    u = sample_grid(int(payload["n_points"]))
    primes = WORKER_PRIMES
    csum = WORKER_CSUM
    if primes is None or csum is None:
        raise RuntimeError("WORKER_PRIMES/CSUM not set (fork inherit failed)")
    q, T, meta = arithmetic_residual(
        u,
        T=float(payload["T"]),
        primes=primes,
        csum=csum,
        detrend=str(payload["detrend"]),
        smooth=int(payload["smooth"]),
    )
    r, backend = energy_ratio_auto(
        q, u, int(payload["degree"]), prefer_gpu=bool(payload.get("prefer_gpu", False))
    )
    return {
        "T": T,
        "x_max": meta["x_max"],
        "degree": int(payload["degree"]),
        "detrend": payload["detrend"],
        "smooth": int(payload["smooth"]),
        "R_d": r,
        "backend": backend,
        "n_primes": meta["n_primes"],
    }


def _control_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import (
        probe_critical_line_mode,
        probe_off_critical_mode,
        probe_persistent_defect,
        sample_grid,
    )
    from pbss.projection_backend import energy_ratio_auto

    T = float(payload["T"])
    n = int(payload["n_points"])
    d = int(payload["degree"])
    u = sample_grid(n)
    gpu = bool(payload.get("prefer_gpu", False))
    r_cl, b1 = energy_ratio_auto(probe_critical_line_mode(u, T=T), u, d, prefer_gpu=gpu)
    r_off, b2 = energy_ratio_auto(
        probe_off_critical_mode(u, T=T, sigma=float(payload["sigma"])),
        u,
        d,
        prefer_gpu=gpu,
    )
    r_def, b3 = energy_ratio_auto(
        probe_persistent_defect(u, eps=float(payload["eps"]), waves=120),
        u,
        d,
        prefer_gpu=gpu,
    )
    return {
        "T": T,
        "degree": d,
        "critical_line_R_d": r_cl,
        "off_critical_R_d": r_off,
        "persistent_defect_R_d": r_def,
        "backends": [b1, b2, b3],
    }


def _mc_batch_job(payload: dict) -> dict:
    """
    Run a batch of MC spectral-defect trials; return summary stats per degree.
    Each trial: random orthogonal-ish defect mixture (Lemma M2 family).
    """
    _env1()
    from pbss.probes import probe_defective, sample_grid
    from pbss.projection_backend import energy_ratio_auto

    n = int(payload["n_points"])
    degrees = list(payload["degrees"])
    n_trials = int(payload["n_trials"])
    seed0 = int(payload["seed0"])
    T = float(payload["T"])  # recorded for bookkeeping; defect is T-independent
    gpu = bool(payload.get("prefer_gpu", False))

    u = sample_grid(n)
    rng = np.random.default_rng(seed0)
    # accumulate per degree
    sums = {d: 0.0 for d in degrees}
    sumsq = {d: 0.0 for d in degrees}
    backend_counts = {"numpy": 0, "cupy": 0}
    eps2_sum = 0.0

    for i in range(n_trials):
        # Controlled spectral defect: HF + random low-degree weight (fast, shipped probe)
        weight = float(rng.uniform(0.5, 3.0))
        waves = int(rng.integers(30, 120))
        deg_def = int(rng.integers(0, 3))
        q = probe_defective(u, waves=waves, defect_degree=deg_def, defect_weight=weight)
        # proxy "eps2" bookkeeping: not exact M2 eps
        eps2_sum += weight * weight / (1.0 + weight * weight)
        for d in degrees:
            r, backend = energy_ratio_auto(q, u, d, prefer_gpu=gpu)
            sums[d] += r
            sumsq[d] += r * r
            backend_counts[backend] = backend_counts.get(backend, 0) + 1
    out_stats = {}
    for d in degrees:
        mean = sums[d] / n_trials
        var = max(0.0, sumsq[d] / n_trials - mean * mean)
        out_stats[str(d)] = {
            "mean_R_d": mean,
            "std_R_d": float(math.sqrt(var)),
            "n": n_trials,
        }
    return {
        "T": T,
        "seed0": seed0,
        "n_trials": n_trials,
        "mean_eps2": float(eps2_sum / n_trials),
        "per_degree": out_stats,
        "backend_counts": backend_counts,
    }


def make_plots(state: dict, out_dir: Path) -> List[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = []
    arith = state.get("arithmetic_rows", [])
    controls = state.get("control_rows", [])
    if not arith:
        return paths

    # focus: median degree, deg1, smooth=1
    degrees = state.get("degrees") or sorted(
        {int(r["degree"]) for r in arith if "degree" in r}
    )
    if not degrees:
        degrees = [4]
    detrends = state.get("detrends") or sorted(
        {str(r.get("detrend", "deg1")) for r in arith}
    )
    d_focus = sorted(degrees)[len(degrees) // 2]
    focus = [
        r
        for r in arith
        if r["degree"] == d_focus and r["detrend"] == "deg1" and r["smooth"] == 1
    ]
    focus = sorted(focus, key=lambda z: z["T"])

    fig, ax = plt.subplots(figsize=(11, 6.5), dpi=140)
    # all detrend for d_focus smooth=1
    for det in detrends:
        rows = sorted(
            [
                r
                for r in arith
                if r["degree"] == d_focus and r["detrend"] == det and r["smooth"] == 1
            ],
            key=lambda z: z["T"],
        )
        if rows:
            ax.plot(
                [r["T"] for r in rows],
                [r["R_d"] for r in rows],
                "o-",
                label=f"arith d={d_focus} {det}/s1",
                ms=4,
            )
    if controls:
        # controls may have multiple degrees — take d_focus
        crows = sorted(
            [c for c in controls if c["degree"] == d_focus], key=lambda z: z["T"]
        )
        if crows:
            ax.plot(
                [c["T"] for c in crows],
                [c["critical_line_R_d"] for c in crows],
                "s--",
                color="green",
                label="critical-line",
            )
            ax.plot(
                [c["T"] for c in crows],
                [c["persistent_defect_R_d"] for c in crows],
                "^--",
                color="red",
                label="persistent defect",
            )
            ax.plot(
                [c["T"] for c in crows],
                [c["off_critical_R_d"] for c in crows],
                "d--",
                color="orange",
                label="off-critical",
            )

    # MC means if present
    mc = state.get("mc_by_T", {})
    if mc:
        Ts = sorted(float(t) for t in mc.keys())
        means = [mc[str(t)]["per_degree"][str(d_focus)]["mean_R_d"] for t in Ts]
        stds = [mc[str(t)]["per_degree"][str(d_focus)]["std_R_d"] for t in Ts]
        ax.errorbar(
            Ts,
            means,
            yerr=stds,
            fmt="x-",
            color="purple",
            label=f"MC defect mean±std d={d_focus}",
            capsize=2,
        )

    ax.set_xlabel("T = log(x_max)")
    ax.set_ylabel("R_d")
    ax.set_yscale("log")
    ax.set_title(
        f"PBSS grand campaign  x_max={state.get('x_max', '?')}  "
        f"MC/T={state.get('mc_per_t', '?')}"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    p1 = out_dir / "grand_Rd_vs_T.png"
    fig.savefig(p1)
    plt.close(fig)
    paths.append(str(p1))

    if focus:
        fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
        ax.plot([r["T"] for r in focus], [r["R_d"] for r in focus], "o-", color="C0")
        ax.set_xlabel("T")
        ax.set_ylabel(f"R_{d_focus} arithmetic deg1/s1")
        ax.set_title("Arithmetic residual (linear scale)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p2 = out_dir / "grand_arith_focus_linear.png"
        fig.savefig(p2)
        plt.close(fig)
        paths.append(str(p2))

    return paths


def main() -> None:
    _env1()
    signal.signal(signal.SIGTERM, _handle_sig)
    signal.signal(signal.SIGINT, _handle_sig)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "grand_campaign"))
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--prime-dir", type=str, default=str(ROOT / "results" / "prime_checkpoints"))
    ap.add_argument(
        "--xmax-candidates",
        type=str,
        default="1e10,5e9,2e9,1e9",
    )
    ap.add_argument("--mc-per-t", type=int, default=50000, help=">=2000 required; default sized for long wall time")
    ap.add_argument("--mc-batch", type=int, default=500, help="trials per worker batch")
    ap.add_argument("--n-points", type=int, default=65536)
    ap.add_argument("--degrees", type=str, default="2,4,6,8")
    ap.add_argument("--detrends", type=str, default="none,deg0,deg1")
    ap.add_argument("--smooths", type=str, default="1,5,15")
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument("--sigma", type=float, default=0.75)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--prefer-gpu", action="store_true")
    ap.add_argument("--force-resieve", action="store_true")
    ap.add_argument(
        "--T-list",
        type=str,
        default="",
        help="comma T values; default auto from x_max",
    )
    args = ap.parse_args()

    if args.mc_per_t < 2000:
        raise SystemExit("AC requires mc-per-t >= 2000")

    out_dir = Path(args.out_dir)
    scratch = Path(args.scratch) if args.scratch else None
    prime_dir = Path(args.prime_dir)
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)
    degrees = [int(x) for x in args.degrees.split(",") if x.strip()]
    detrends = [x.strip() for x in args.detrends.split(",") if x.strip()]
    smooths = [int(x) for x in args.smooths.split(",") if x.strip()]

    candidates = [float(x) for x in args.xmax_candidates.split(",") if x.strip()]
    x_max, xmax_note = choose_xmax(candidates, avail_ram())

    # resume state
    state = _load_state(out_dir) or {}
    if state.get("status") == "completed" and not args.force_resieve:
        print("Campaign already completed — nothing to do.", flush=True)
        print(state.get("narrative", ""), flush=True)
        return

    if not state:
        state = {
            "status": "running",
            "started": time.strftime("%Y-%m-%d %H:%M:%S"),
            "started_unix": time.time(),
            "x_max": x_max,
            "xmax_note": xmax_note,
            "mc_per_t": args.mc_per_t,
            "n_points": args.n_points,
            "degrees": degrees,
            "detrends": detrends,
            "smooths": smooths,
            "workers": workers,
            "prefer_gpu": bool(args.prefer_gpu),
            "arithmetic_rows": [],
            "control_rows": [],
            "mc_by_T": {},
            "mc_completed_trials": {},
            "phase": "init",
            "not_an_RH_proof": True,
        }
    else:
        # allow upgrading mc target on resume
        state["mc_per_t"] = max(int(state.get("mc_per_t", 0)), args.mc_per_t)
        x_max = float(state.get("x_max", x_max))
        print(f"Resuming campaign from phase={state.get('phase')}", flush=True)

    if args.T_list.strip():
        T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    else:
        T_max = math.log(x_max)
        raw = list(np.linspace(10.0, min(16.0, T_max), 5)) + list(
            np.linspace(min(16.0, T_max), T_max, 10)
        )
        T_values = sorted({float(t) for t in raw if math.exp(t) <= x_max * 1.001})
    state["T_values"] = T_values
    state["x_max"] = x_max
    state["xmax_note"] = xmax_note
    _save_state(out_dir, state)
    if scratch:
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "grand_campaign_state.json").write_text(json.dumps(state, indent=2))

    print(xmax_note, flush=True)
    print(
        f"T={T_values}\nmc_per_t={args.mc_per_t} n_points={args.n_points} "
        f"degrees={degrees} workers={workers} gpu={args.prefer_gpu}",
        flush=True,
    )

    # ----- Phase A: primes -----
    from pbss.primes_io import ensure_primes

    state["phase"] = "sieve"
    _save_state(out_dir, state)
    print(f"Ensuring primes ≤ {x_max:.3e} (checkpoint dir {prime_dir})…", flush=True)
    t0 = time.time()
    try:
        primes, pmeta = ensure_primes(
            int(x_max),
            prime_dir,
            force_resieve=args.force_resieve,
            segment_size=25_000_000,
        )
    except MemoryError:
        # step down
        for fb in [5e9, 2e9, 1e9, 5e8]:
            if fb >= x_max:
                continue
            print(f"MemoryError — falling back to x_max={fb:.3e}", flush=True)
            x_max = fb
            T_values = [t for t in T_values if math.exp(t) <= x_max * 1.001]
            state["T_values"] = T_values
            state["x_max"] = x_max
            state["xmax_note"] = f"OOM fallback to {x_max:.3e}"
            try:
                primes, pmeta = ensure_primes(int(x_max), prime_dir, segment_size=25_000_000)
                break
            except MemoryError:
                continue
        else:
            raise
    print(
        f"  primes n={len(primes)}  meta={pmeta.get('loaded_from_checkpoint')=} "
        f"sieve_s={pmeta.get('sieve_seconds', 'ckpt')}  elapsed={time.time()-t0:.1f}s",
        flush=True,
    )
    state["n_primes"] = int(len(primes))
    state["prime_meta"] = {k: pmeta[k] for k in pmeta if k != "path"}
    state["phase"] = "arithmetic"
    _save_state(out_dir, state)

    # convert memmap to ndarray for fork COW reliability across large jobs
    # (memmap is fine read-only; keep as-is)
    primes_arr = primes

    # ----- Phase B: arithmetic ablations -----
    done_keys = {
        (round(r["T"], 6), r["degree"], r["detrend"], r["smooth"])
        for r in state.get("arithmetic_rows", [])
    }
    arith_payloads = []
    for T in T_values:
        for d in degrees:
            for det in detrends:
                for sm in smooths:
                    key = (round(T, 6), d, det, sm)
                    if key in done_keys:
                        continue
                    arith_payloads.append(
                        {
                            "T": T,
                            "degree": d,
                            "detrend": det,
                            "smooth": sm,
                            "n_points": args.n_points,
                            "prefer_gpu": args.prefer_gpu,
                        }
                    )
    print(f"Arithmetic jobs remaining: {len(arith_payloads)}", flush=True)
    # Sequential in parent: multi-GB primes/csum must NOT be forked to 70 workers
    # (even COW fails if any write/copy triggers). MC/controls stay parallel (no primes).
    from pbss.probes import arithmetic_residual, prime_log_cumsum, sample_grid
    from pbss.projection_backend import energy_ratio_auto

    print("Building log-cumsum(θ prefix) once in parent…", flush=True)
    t_cs = time.time()
    primes_np = np.asarray(primes_arr)
    csum_np = prime_log_cumsum(primes_np)
    print(f"  csum done in {time.time()-t_cs:.1f}s  len={len(csum_np)}", flush=True)

    if arith_payloads and not _SHUTDOWN:
        u_cache: dict[int, np.ndarray] = {}
        for i, p in enumerate(arith_payloads, 1):
            if _SHUTDOWN:
                break
            n = int(p["n_points"])
            if n not in u_cache:
                u_cache[n] = sample_grid(n)
            u = u_cache[n]
            q, T, meta = arithmetic_residual(
                u,
                T=float(p["T"]),
                primes=primes_np,
                csum=csum_np,
                detrend=str(p["detrend"]),
                smooth=int(p["smooth"]),
            )
            r, backend = energy_ratio_auto(
                q, u, int(p["degree"]), prefer_gpu=bool(p.get("prefer_gpu", False))
            )
            row = {
                "T": T,
                "x_max": meta["x_max"],
                "degree": int(p["degree"]),
                "detrend": p["detrend"],
                "smooth": int(p["smooth"]),
                "R_d": r,
                "backend": backend,
                "n_primes": meta["n_primes"],
            }
            state["arithmetic_rows"].append(row)
            if i % 25 == 0 or i == len(arith_payloads):
                _save_state(out_dir, state)
                print(
                    f"  arith {i}/{len(arith_payloads)} T={row['T']:.2f} d={row['degree']} "
                    f"{row['detrend']}/s{row['smooth']} R={row['R_d']:.4e}",
                    flush=True,
                )
        _save_state(out_dir, state)

    # free multi-GB tables before forking MC/control workers
    del primes_np, csum_np, primes_arr
    try:
        del primes
    except Exception:
        pass
    import gc

    gc.collect()
    print("Freed primes/csum before parallel MC/controls", flush=True)

    # ----- Phase C: controls -----
    state["phase"] = "controls"
    _save_state(out_dir, state)
    done_c = {(round(c["T"], 6), c["degree"]) for c in state.get("control_rows", [])}
    ctrl_payloads = []
    for T in T_values:
        for d in degrees:
            if (round(T, 6), d) in done_c:
                continue
            ctrl_payloads.append(
                {
                    "T": T,
                    "degree": d,
                    "n_points": args.n_points,
                    "eps": args.eps,
                    "sigma": args.sigma,
                    "prefer_gpu": args.prefer_gpu,
                }
            )
    print(f"Control jobs remaining: {len(ctrl_payloads)}", flush=True)
    if ctrl_payloads and not _SHUTDOWN:
        with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CTX) as ex:
            for row in ex.map(_control_job, ctrl_payloads):
                state["control_rows"].append(row)
        _save_state(out_dir, state)

    # ----- Phase D: MC defects -----
    state["phase"] = "mc"
    _save_state(out_dir, state)
    mc_completed: Dict[str, int] = {
        str(k): int(v) for k, v in state.get("mc_completed_trials", {}).items()
    }
    # Separate accumulator (never alias state["mc_by_T"] — finalized overwrites that)
    mc_accum: Dict[str, Any] = {}

    for T in T_values:
        if _SHUTDOWN:
            break
        key = str(round(T, 6))
        have = mc_completed.get(key, 0)
        need = args.mc_per_t - have
        if need <= 0:
            print(f"MC T={T:.3f} already complete ({have})", flush=True)
            continue
        print(f"MC T={T:.3f}: need {need} more trials (have {have})", flush=True)
        batches = []
        left = need
        seed_base = int(1_000_003 * T + have)
        bi = 0
        while left > 0:
            nb = min(args.mc_batch, left)
            batches.append(
                {
                    "T": T,
                    "n_trials": nb,
                    "seed0": seed_base + bi * 10_007,
                    "n_points": args.n_points,
                    "degrees": degrees,
                    "prefer_gpu": args.prefer_gpu,
                }
            )
            left -= nb
            bi += 1

        mc_accum[key] = {
            "T": T,
            "n_trials": have,
            "mean_eps2_num": 0.0,
            "per_degree": {
                str(d): {"sum": 0.0, "sumsq": 0.0, "n": 0} for d in degrees
            },
            "backend_counts": {},
        }
        # seed accum from prior finalized means if resuming mid-T
        if have > 0 and key in state.get("mc_by_T", {}):
            prev = state["mc_by_T"][key]
            mc_accum[key]["mean_eps2_num"] = float(prev.get("mean_eps2", 0.0)) * have
            for d_str, st in prev.get("per_degree", {}).items():
                n = int(st.get("n", 0))
                m = float(st.get("mean_R_d", 0.0))
                s = float(st.get("std_R_d", 0.0))
                mc_accum[key]["per_degree"][d_str] = {
                    "sum": m * n,
                    "sumsq": (s * s + m * m) * n,
                    "n": n,
                }

        with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CTX) as ex:
            futs = [ex.submit(_mc_batch_job, b) for b in batches]
            for j, fut in enumerate(as_completed(futs), 1):
                if _SHUTDOWN:
                    break
                batch = fut.result()
                rec = mc_accum[key]
                nt = batch["n_trials"]
                rec["n_trials"] += nt
                rec["mean_eps2_num"] += batch["mean_eps2"] * nt
                for d_str, st in batch["per_degree"].items():
                    pd = rec["per_degree"][d_str]
                    n = st["n"]
                    m = st["mean_R_d"]
                    s = st["std_R_d"]
                    pd["sum"] += m * n
                    pd["sumsq"] += (s * s + m * m) * n
                    pd["n"] += n
                for bk, cnt in batch["backend_counts"].items():
                    rec["backend_counts"][bk] = rec["backend_counts"].get(bk, 0) + cnt
                mc_completed[key] = rec["n_trials"]
                if j % 5 == 0 or j == len(futs):
                    finalized = {
                        "T": T,
                        "n_trials": rec["n_trials"],
                        "mean_eps2": rec["mean_eps2_num"] / max(rec["n_trials"], 1),
                        "per_degree": {},
                        "backend_counts": dict(rec["backend_counts"]),
                    }
                    for d_str, pd in rec["per_degree"].items():
                        n = max(pd["n"], 1)
                        mean = pd["sum"] / n
                        var = max(0.0, pd["sumsq"] / n - mean * mean)
                        finalized["per_degree"][d_str] = {
                            "mean_R_d": mean,
                            "std_R_d": math.sqrt(var),
                            "n": pd["n"],
                        }
                    state["mc_by_T"][key] = finalized
                    state["mc_completed_trials"] = mc_completed
                    _save_state(out_dir, state)
                    if scratch:
                        (scratch / "grand_campaign_state.json").write_text(
                            json.dumps(state, indent=2, default=str)
                        )
                    print(
                        f"  MC T={T:.2f} batch {j}/{len(futs)}  "
                        f"trials={rec['n_trials']}/{args.mc_per_t}",
                        flush=True,
                    )

        # finalize this T from accum
        rec = mc_accum.get(key)
        if rec is not None:
            finalized = {
                "T": T,
                "n_trials": rec["n_trials"],
                "mean_eps2": rec["mean_eps2_num"] / max(rec["n_trials"], 1),
                "per_degree": {},
                "backend_counts": dict(rec["backend_counts"]),
            }
            for d_str, pd in rec["per_degree"].items():
                n = max(pd["n"], 1)
                mean = pd["sum"] / n
                var = max(0.0, pd["sumsq"] / n - mean * mean)
                finalized["per_degree"][d_str] = {
                    "mean_R_d": mean,
                    "std_R_d": math.sqrt(var),
                    "n": pd["n"],
                }
            state["mc_by_T"][key] = finalized
        state["mc_completed_trials"] = mc_completed
        _save_state(out_dir, state)
    # ----- narrative -----
    # ensure metadata present even when resuming old state blobs
    state["degrees"] = degrees
    state["detrends"] = detrends
    state["smooths"] = smooths
    state["n_points"] = args.n_points
    state["workers"] = workers
    d_focus = sorted(degrees)[len(degrees) // 2]
    focus = sorted(
        [
            r
            for r in state["arithmetic_rows"]
            if r["degree"] == d_focus and r["detrend"] == "deg1" and r["smooth"] == 1
        ],
        key=lambda z: z["T"],
    )
    if len(focus) >= 2:
        R = np.array([r["R_d"] for r in focus])
        reading = {
            "d_focus": d_focus,
            "R_first": float(R[0]),
            "R_last": float(R[-1]),
            "R_min": float(R.min()),
            "R_max": float(R.max()),
            "max_T": float(focus[-1]["T"]),
            "max_x_max": float(focus[-1]["x_max"]),
        }
        narrative = (
            f"Grand campaign x_max={state['x_max']:.3e}, focus d={d_focus} deg1/s1: "
            f"R_d {reading['R_first']:.4e} → {reading['R_last']:.4e} "
            f"(min {reading['R_min']:.4e}, max {reading['R_max']:.4e}) at T≤{reading['max_T']:.2f}. "
            f"MC defect trials target {args.mc_per_t}/T. "
            "Not a proof of RH or full Theorem A."
        )
    else:
        reading = {}
        narrative = "Insufficient arithmetic focus rows. Not a proof of RH."

    state["reading"] = reading
    state["narrative"] = narrative
    state["phase"] = "plotting"
    plot_paths = make_plots(state, out_dir)
    state["plot_paths"] = plot_paths

    complete = (not _SHUTDOWN) and all(
        mc_completed.get(str(round(T, 6)), 0) >= args.mc_per_t for T in T_values
    )
    state["status"] = "completed" if complete else "stopped_partial"
    state["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    state["elapsed_s"] = time.time() - float(state.get("started_unix", time.time()))
    state["mc_completed_trials"] = mc_completed
    _save_state(out_dir, state)

    # human summary
    lines = [
        "PBSS GRAND CAMPAIGN",
        f"status={state['status']}  x_max={state['x_max']:.3e}  n_primes={state.get('n_primes')}",
        f"mc_per_t={args.mc_per_t}  n_points={args.n_points}  workers={workers}",
        f"elapsed_s={state['elapsed_s']:.1f}",
        "",
        narrative,
        "",
        "MC trials completed per T:",
    ]
    for T in T_values:
        k = str(round(T, 6))
        lines.append(f"  T={T:.3f}  trials={mc_completed.get(k, 0)}")
    lines += ["", "Not an unconditional RH proof.", f"Plots: {plot_paths}"]
    text = "\n".join(lines) + "\n"
    (out_dir / "grand_campaign_summary.txt").write_text(text)
    print(text)

    if scratch:
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "grand_campaign_state.json").write_text(json.dumps(state, indent=2, default=str))
        (scratch / "grand_campaign_summary.txt").write_text(text)
        for p in plot_paths:
            src = Path(p)
            if src.exists():
                (scratch / src.name).write_bytes(src.read_bytes())
        stamp = "COMPLETE" if complete else "PARTIAL"
        (scratch / f"GRAND_CAMPAIGN_{stamp}").write_text(text)

    print("Workers exited. Shutdown clean.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
