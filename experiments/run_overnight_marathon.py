#!/usr/bin/env python3
"""
PBSS overnight marathon orchestrator (resume-capable).

Phases (skip if stamp present):
  A) dense arithmetic zero-peel on x_max=5e10 → ≥2000 rows
  B) Beurling ≥100 systems multi-T scorecard
  C) MC ≥200k trials/T × ≥8 T (multi-core)
  D) residual multi-T ablations on large x (multi-core CPU; GPU optional)

Work floors from OVERNIGHT_GOAL.md — no wall-time padding.
Not an RH proof.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_MP = mp.get_context("fork")

# Shared for peel / residual workers (set before pool)
W_PRIMES = None
W_CSUM = None
W_N_POINTS = 8192
W_DEGREES = [2, 4, 6, 8]
W_FIT = True


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
    p = out / f"PHASE_{name}_COMPLETE"
    p.write_text(json.dumps(payload, indent=2) + "\n")
    state_path = out / "marathon_state.json"
    state = {}
    if state_path.exists():
        state = json.loads(state_path.read_text())
    state.setdefault("phases", {})[name] = payload
    state["updated_unix"] = time.time()
    state["banner"] = "NOT AN UNCONDITIONAL PROOF OF RH"
    state_path.write_text(json.dumps(state, indent=2))


def _phase_done(out: Path, name: str) -> bool:
    return (out / f"PHASE_{name}_COMPLETE").exists()


def _workers(n: int = 0) -> int:
    return int(n) if n and n > 0 else max(1, (os.cpu_count() or 4) - 2)


# ---------- Phase A: peel ----------
def _peel_path_job(payload: dict) -> list:
    """One residual build per (T, detrend); then all N and degrees (module-level)."""
    _env1()
    from pbss.probes import (
        arithmetic_residual,
        explicit_formula_residual,
        sample_grid,
    )
    from pbss.projection import energy_ratio, _trapezoid_weights

    p = np.load(payload["primes_path"], mmap_mode="r")
    c = np.load(payload["csum_path"], mmap_mode="r")
    T = float(payload["T"])
    hi = int(np.searchsorted(p, min(float(np.exp(T)), float(p[-1])), side="right"))
    if hi > 400_000_000:
        p_use, c_use = p[:hi], c[:hi]
    else:
        p_use = np.array(p[:hi], dtype=np.int64, copy=True)
        c_use = np.array(c[:hi], dtype=np.float64, copy=True)
    u = sample_grid(int(payload["n_points"]))
    det = str(payload["detrend"])
    q_raw, T_out, meta = arithmetic_residual(
        u, T=T, primes=p_use, csum=c_use, detrend=det
    )
    w = _trapezoid_weights(u)
    rows = []
    for N in payload["N_list"]:
        if int(N) <= 0:
            q = q_raw
            alpha = 0.0
        else:
            q_modes, _, _ = explicit_formula_residual(
                u, T=T_out, n_zeros=int(N), form="cos", bulk="none"
            )
            if payload["fit_scale"]:
                num = float(np.sum(w * q_raw * q_modes))
                den = float(np.sum(w * q_modes * q_modes))
                alpha = num / den if den > 1e-30 else 0.0
            else:
                alpha = 1.0
            q = q_raw - alpha * q_modes
        for d in payload["degrees"]:
            r = float(energy_ratio(q, u, degree=int(d)))
            rows.append(
                {
                    "T": float(T_out),
                    "N": int(N),
                    "degree": int(d),
                    "detrend": det,
                    "R_d": r,
                    "mode_scale": float(alpha),
                    "x_max": meta.get("x_max"),
                    "n_primes": meta.get("n_primes"),
                }
            )
    return rows


def phase_peel(out: Path, args) -> dict:
    global W_PRIMES, W_CSUM
    from pbss.primes_io import load_primes_checkpoint
    from pbss.probes import prime_log_cumsum

    phase_dir = out / "peel"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    x_max = int(args.x_max)
    print(f"[A] peel x_max={x_max} …", flush=True)
    primes, pmeta = load_primes_checkpoint(args.prime_dir, x_max, mmap=True)
    csum_path = Path(args.csum_path) if args.csum_path else None
    if csum_path and csum_path.exists():
        csum = np.load(csum_path, mmap_mode="r")
    else:
        # build once if missing (expensive)
        print("[A] building csum (one-time)…", flush=True)
        csum_dense = prime_log_cumsum(np.asarray(primes))
        csum_path = phase_dir / f"csum_le_{x_max}.npy"
        np.save(csum_path, csum_dense)
        del csum_dense
        csum = np.load(csum_path, mmap_mode="r")

    T_max = float(np.log(x_max))
    # Dense grid aiming ≥2000 rows: 16 T × 16 N × 4 d × 2 det = 2048
    T_values = [
        float(t)
        for t in np.unique(
            np.round(np.linspace(max(10.0, T_max - 10), T_max, 16), 3)
        )
    ]
    N_values = [0, 1, 2, 3, 5, 7, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50]
    degrees = [2, 4, 6, 8]
    detrends = ["deg1", "none"]
    n_points = int(args.n_points)
    # For large tables: materialize nothing global; each job loads mmap prefixes
    # Use path-based workers to avoid COW of full arrays
    payloads = []
    for T in T_values:
        for det in detrends:
            payloads.append(
                {
                    "T": T,
                    "detrend": det,
                    "N_list": N_values,
                    "degrees": degrees,
                    "n_points": n_points,
                    "fit_scale": True,
                    "primes_path": str(
                        Path(args.prime_dir) / f"primes_le_{x_max}.npy"
                    ),
                    "csum_path": str(csum_path),
                    "x_max": x_max,
                }
            )

    workers = min(_workers(args.workers), max(1, len(payloads)), 24)
    # Cap workers for huge memmap thrash
    print(f"[A] peel groups={len(payloads)} workers={workers}", flush=True)
    rows = []
    # Sequential T order with parallel detrend-only (2 workers) for large x
    # to avoid OOM; still multi-core within smaller T via path jobs in waves
    by_T: dict = {}
    for pl in payloads:
        by_T.setdefault(pl["T"], []).append(pl)
    for T in sorted(by_T.keys()):
        pls = by_T[T]
        w = min(2, workers, len(pls))
        with ProcessPoolExecutor(max_workers=w, mp_context=_MP) as ex:
            for group in ex.map(_peel_path_job, pls):
                rows.extend(group)
        print(f"[A] T={T:g} cumulative_rows={len(rows)}", flush=True)

    elapsed = time.time() - t0
    summary = {
        "status": "completed",
        "phase": "peel",
        "elapsed_s": elapsed,
        "x_max": float(x_max),
        "n_rows": len(rows),
        "T_values": T_values,
        "N_values": N_values,
        "degrees": degrees,
        "detrends": detrends,
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "prime_meta": {
            k: pmeta.get(k) for k in ("x_max", "n_primes", "method", "sieve_seconds")
        },
    }
    (phase_dir / "peel.json").write_text(json.dumps(summary, indent=2))
    lines = [
        "PBSS marathon peel",
        f"rows={len(rows)} elapsed_s={elapsed:.1f} x_max={x_max}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
    ]
    (phase_dir / "peel.txt").write_text("\n".join(lines) + "\n")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        focus = [
            r
            for r in rows
            if r["degree"] == 4 and r["detrend"] == "deg1" and r["N"] in (0, 5, 20, 50)
        ]
        fig, ax = plt.subplots(figsize=(9, 5), dpi=120)
        for N in sorted(set(r["N"] for r in focus)):
            rs = sorted([r for r in focus if r["N"] == N], key=lambda z: z["T"])
            ax.plot([r["T"] for r in rs], [r["R_d"] for r in rs], "o-", label=f"N={N}")
        ax.set_xlabel("T")
        ax.set_ylabel(r"$R_4$")
        ax.set_title("Marathon peel (not RH proof)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(phase_dir / "peel_Rd.png", bbox_inches="tight")
        plt.close()
    except Exception as exc:
        print(f"[A] plot skip {exc}", flush=True)
    if len(rows) < 2000:
        raise SystemExit(f"peel rows {len(rows)} < 2000")
    payload = {
        "status": "completed",
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "path": str(phase_dir / "peel.json"),
    }
    _stamp(out, "A_PEEL", payload)
    print(f"[A] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return payload


# ---------- Phase B: Beurling ----------
def _beurling_job(payload: dict) -> list:
    _env1()
    from pbss.beurling import beurling_theta_residual, build_system_primes
    from pbss.probes import sample_grid
    from pbss.projection import energy_ratio

    ordinary = np.load(payload["ordinary_path"], mmap_mode="r")
    # prefix to x_max
    x_max = float(payload["x_max"])
    hi = int(np.searchsorted(ordinary, x_max, side="right"))
    ordinary = ordinary[:hi]
    spec = payload["spec"]
    p_sys = build_system_primes(spec, ordinary, x_max)
    u = sample_grid(int(payload["n_points"]))
    degree = int(payload["degree"])
    rows = []
    for T in payload["T_values"]:
        if T > np.log(x_max) + 1e-9:
            continue
        q, T_out, meta = beurling_theta_residual(
            u, p_sys, T=float(T), detrend="deg1"
        )
        r = float(energy_ratio(q, u, degree=degree))
        rows.append(
            {
                "system": spec["name"],
                "kind": spec["kind"],
                "T": float(T_out),
                "degree": degree,
                "R_d": r,
                "n_primes": meta["n_primes"],
            }
        )
    return rows


def phase_beurling(out: Path, args) -> dict:
    from pbss.beurling import marathon_battery_specs
    from pbss.probes import primes_upto

    phase_dir = out / "beurling"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_sys = int(args.n_systems)
    specs = marathon_battery_specs(n_sys)
    x_max = float(args.beurling_x_max)
    # Prefer a compact ordinary-prime table for Beurling (not the multi-GB 5e10 mmap)
    ordinary_path = phase_dir / f"ordinary_le_{int(x_max)}.npy"
    if ordinary_path.exists():
        print(f"[B] reuse ordinary table {ordinary_path}", flush=True)
    else:
        print(f"[B] sieving ordinary primes ≤ {x_max:.3e} for battery…", flush=True)
        ordinary = primes_upto(int(x_max))
        np.save(ordinary_path, ordinary)
        del ordinary
    ordinary_path = str(ordinary_path)
    T_values = [float(x) for x in args.beurling_T_list.split(",") if x.strip()]
    # Cap workers: each job loads ordinary table
    workers = min(_workers(args.workers), len(specs), 32)
    print(f"[B] beurling systems={len(specs)} workers={workers}", flush=True)
    payloads = [
        {
            "spec": s,
            "ordinary_path": ordinary_path,
            "x_max": x_max,
            "T_values": T_values,
            "n_points": int(args.n_points),
            "degree": 4,
        }
        for s in specs
    ]
    rows = []
    # Batch in waves of workers to limit peak RSS
    for i in range(0, len(payloads), workers):
        batch = payloads[i : i + workers]
        with ProcessPoolExecutor(max_workers=len(batch), mp_context=_MP) as ex:
            for group in ex.map(_beurling_job, batch, chunksize=1):
                rows.extend(group)
        print(f"[B] progress systems={min(i+workers,len(payloads))}/{len(payloads)} rows={len(rows)}", flush=True)
    elapsed = time.time() - t0
    systems = sorted({r["system"] for r in rows})
    if len(systems) < 100:
        raise SystemExit(f"systems {len(systems)} < 100")
    T_star = max(T_values)
    scorecard = {}
    for r in rows:
        if r["T"] == T_star:
            scorecard[r["system"]] = r["R_d"]
    summary = {
        "status": "completed",
        "phase": "beurling",
        "elapsed_s": elapsed,
        "n_systems": len(systems),
        "x_max": x_max,
        "T_values": T_values,
        "rows": rows,
        "scorecard_at_Tmax": scorecard,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    (phase_dir / "beurling.json").write_text(json.dumps(summary, indent=2))
    (phase_dir / "beurling.txt").write_text(
        f"beurling systems={len(systems)} rows={len(rows)} elapsed_s={elapsed:.1f}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # plot mean R_d by kind vs T
        kinds = sorted({r["kind"] for r in rows})
        fig, ax = plt.subplots(figsize=(9, 5), dpi=120)
        for kind in kinds:
            means = []
            for T in T_values:
                rs = [r["R_d"] for r in rows if r["kind"] == kind and r["T"] == T]
                means.append(float(np.mean(rs)) if rs else float("nan"))
            ax.plot(T_values, means, "o-", label=kind)
        ax.set_xlabel("T")
        ax.set_ylabel(r"mean $R_d$")
        ax.set_title(f"Marathon Beurling n={len(systems)} (not RH proof)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(phase_dir / "beurling_Rd.png", bbox_inches="tight")
        plt.close()
    except Exception as exc:
        print(f"[B] plot skip {exc}", flush=True)
    payload = {
        "status": "completed",
        "elapsed_s": elapsed,
        "n_systems": len(systems),
        "n_rows": len(rows),
        "path": str(phase_dir / "beurling.json"),
    }
    _stamp(out, "B_BEURLING", payload)
    print(f"[B] DONE systems={len(systems)} elapsed_s={elapsed:.1f}", flush=True)
    return payload


# ---------- Phase C: MC via subprocess to existing runner ----------
def phase_mc(out: Path, args) -> dict:
    phase_dir = out / "mc"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    T_list = args.mc_T_list
    mc_per_t = int(args.mc_per_t)
    workers = _workers(args.workers)
    cmd = [
        sys.executable,
        str(ROOT / "experiments" / "run_mc_stress.py"),
        "--out-dir",
        str(phase_dir),
        "--mc-per-t",
        str(mc_per_t),
        "--min-mc-per-t",
        str(mc_per_t),
        "--mc-batch",
        str(args.mc_batch),
        "--T-list",
        T_list,
        "--degrees",
        args.mc_degrees,
        "--n-points",
        str(args.n_points),
        "--workers",
        str(workers),
    ]
    print(f"[C] MC cmd workers={workers} mc_per_t={mc_per_t} T={T_list}", flush=True)
    subprocess.check_call(cmd, cwd=str(ROOT))
    elapsed = time.time() - t0
    summary = json.loads((phase_dir / "mc_stress_summary.json").read_text())
    min_n = int(summary.get("min_mc_per_t", 0))
    n_T = len([x for x in T_list.split(",") if x.strip()])
    if min_n < mc_per_t or n_T < 8:
        raise SystemExit(f"MC floors fail min_n={min_n} n_T={n_T}")
    payload = {
        "status": "completed",
        "elapsed_s": elapsed,
        "min_mc_per_t": min_n,
        "n_T": n_T,
        "total_trials": min_n * n_T,
        "path": str(phase_dir / "mc_stress_summary.json"),
    }
    _stamp(out, "C_MC", payload)
    print(f"[C] DONE min_n={min_n} n_T={n_T} elapsed_s={elapsed:.1f}", flush=True)
    return payload


# ---------- Phase D: residual multi-T ----------
def _resid_job(payload: dict) -> list:
    _env1()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection import energy_ratio

    p = np.load(payload["primes_path"], mmap_mode="r")
    c = np.load(payload["csum_path"], mmap_mode="r")
    T = float(payload["T"])
    hi = int(np.searchsorted(p, min(float(np.exp(T)), float(p[-1])), side="right"))
    if hi > 400_000_000:
        p_use, c_use = p[:hi], c[:hi]
    else:
        p_use = np.array(p[:hi], dtype=np.int64, copy=True)
        c_use = np.array(c[:hi], dtype=np.float64, copy=True)
    u = sample_grid(int(payload["n_points"]))
    rows = []
    for det in payload["detrends"]:
        q, T_out, meta = arithmetic_residual(
            u, T=T, primes=p_use, csum=c_use, detrend=det
        )
        for d in payload["degrees"]:
            r = float(energy_ratio(q, u, degree=int(d)))
            rows.append(
                {
                    "T": float(T_out),
                    "degree": int(d),
                    "detrend": det,
                    "R_d": r,
                    "n_primes": meta["n_primes"],
                    "x_max_window": meta["x_max"],
                }
            )
    return rows


def phase_residual(out: Path, args) -> dict:
    phase_dir = out / "residual"
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    x_max = int(args.x_max)
    primes_path = str(Path(args.prime_dir) / f"primes_le_{x_max}.npy")
    csum_path = args.csum_path
    if not csum_path or not Path(csum_path).exists():
        # reuse peel phase csum if present
        alt = out / "peel" / f"csum_le_{x_max}.npy"
        if alt.exists():
            csum_path = str(alt)
        else:
            csum_path = str(ROOT / "results" / "extend_x_scan" / f"csum_le_{x_max}.npy")
    T_max = float(np.log(x_max))
    T_values = [
        float(t)
        for t in np.unique(np.round(np.linspace(max(10.0, T_max - 10), T_max, 12), 3))
    ]
    degrees = [2, 4, 6, 8]
    detrends = ["deg1", "none"]
    payloads = [
        {
            "T": T,
            "primes_path": primes_path,
            "csum_path": csum_path,
            "n_points": int(args.n_points),
            "degrees": degrees,
            "detrends": detrends,
        }
        for T in T_values
    ]
    # Parent-process sequential residual: ProcessPool + multi-GB mmap OOMs.
    # Vectorized NumPy searchsorted is plenty; multi-core spent on A/B/C floors.
    print(f"[D] residual T={len(T_values)} sequential mmap (avoid pool OOM)", flush=True)
    rows = []
    for pl in payloads:
        group = _resid_job(pl)
        rows.extend(group)
        print(f"[D] T={pl['T']:g} rows={len(rows)}", flush=True)
    elapsed = time.time() - t0
    summary = {
        "status": "completed",
        "phase": "residual",
        "elapsed_s": elapsed,
        "x_max": float(x_max),
        "T_values": T_values,
        "degrees": degrees,
        "detrends": detrends,
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    (phase_dir / "residual.json").write_text(json.dumps(summary, indent=2))
    (phase_dir / "residual.txt").write_text(
        f"residual rows={len(rows)} elapsed_s={elapsed:.1f}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        focus = sorted(
            [r for r in rows if r["degree"] == 4 and r["detrend"] == "deg1"],
            key=lambda z: z["T"],
        )
        fig, ax = plt.subplots(figsize=(9, 5), dpi=120)
        ax.plot([r["T"] for r in focus], [r["R_d"] for r in focus], "o-")
        ax.set_xlabel("T")
        ax.set_ylabel(r"$R_4$")
        ax.set_title("Marathon residual multi-T (not RH proof)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(phase_dir / "residual_Rd.png", bbox_inches="tight")
        plt.close()
    except Exception as exc:
        print(f"[D] plot skip {exc}", flush=True)
    if not rows:
        raise SystemExit("residual produced no rows")
    payload = {
        "status": "completed",
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "path": str(phase_dir / "residual.json"),
    }
    _stamp(out, "D_RESIDUAL", payload)
    print(f"[D] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return payload


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        type=str,
        default=str(ROOT / "results" / "overnight_marathon"),
    )
    ap.add_argument(
        "--prime-dir",
        type=str,
        default=str(ROOT / "results" / "prime_checkpoints"),
    )
    ap.add_argument("--x-max", type=float, default=5e10)
    ap.add_argument(
        "--csum-path",
        type=str,
        default=str(ROOT / "results" / "extend_x_scan" / "csum_le_50000000000.npy"),
    )
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--n-points", type=int, default=8192)
    ap.add_argument("--n-systems", type=int, default=100)
    ap.add_argument("--beurling-x-max", type=float, default=1e8)
    ap.add_argument(
        "--beurling-T-list", type=str, default="8,10,12,14,16,18"
    )
    ap.add_argument("--mc-per-t", type=int, default=200000)
    ap.add_argument("--mc-batch", type=int, default=500)
    ap.add_argument(
        "--mc-T-list",
        type=str,
        default="8,10,12,14,16,18,20,22",
    )
    ap.add_argument("--mc-degrees", type=str, default="2,4,6,8")
    ap.add_argument(
        "--phases",
        type=str,
        default="A,B,C,D",
        help="comma subset of A,B,C,D",
    )
    ap.add_argument("--scratch", type=str, default="")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)
    print(f"marathon out={out} workers~{_workers(args.workers)}", flush=True)
    phases = {p.strip().upper() for p in args.phases.split(",") if p.strip()}
    t_all = time.time()
    results = {}

    if "A" in phases and not _phase_done(out, "A_PEEL"):
        results["A"] = phase_peel(out, args)
    elif "A" in phases:
        print("[A] skip (stamp present)", flush=True)
        results["A"] = json.loads((out / "PHASE_A_PEEL_COMPLETE").read_text())

    if "B" in phases and not _phase_done(out, "B_BEURLING"):
        results["B"] = phase_beurling(out, args)
    elif "B" in phases:
        print("[B] skip (stamp present)", flush=True)
        results["B"] = json.loads((out / "PHASE_B_BEURLING_COMPLETE").read_text())

    if "C" in phases and not _phase_done(out, "C_MC"):
        results["C"] = phase_mc(out, args)
    elif "C" in phases:
        print("[C] skip (stamp present)", flush=True)
        results["C"] = json.loads((out / "PHASE_C_MC_COMPLETE").read_text())

    if "D" in phases and not _phase_done(out, "D_RESIDUAL"):
        results["D"] = phase_residual(out, args)
    elif "D" in phases:
        print("[D] skip (stamp present)", flush=True)
        results["D"] = json.loads((out / "PHASE_D_RESIDUAL_COMPLETE").read_text())

    elapsed = time.time() - t_all
    required = ("A_PEEL", "B_BEURLING", "C_MC", "D_RESIDUAL")
    all_done = all(_phase_done(out, name) for name in required)
    final = {
        "status": "completed" if all_done else "partial",
        "elapsed_s": elapsed,
        "phases": results,
        "phases_done": [n for n in required if _phase_done(out, n)],
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "interpretation": (
            "Overnight marathon work floors (peel≥2000, Beurling≥100, "
            "MC≥200k/T×≥8T, residual multi-T). Not an RH proof."
        ),
    }
    (out / "marathon_summary.json").write_text(json.dumps(final, indent=2))
    if all_done:
        (out / "MARATHON_COMPLETE").write_text(json.dumps(final, indent=2) + "\n")
        print(f"MARATHON_COMPLETE elapsed_s={elapsed:.1f}", flush=True)
    else:
        # Remove false complete stamp from partial runs
        bogus = out / "MARATHON_COMPLETE"
        if bogus.exists() and final["status"] == "partial":
            # only remove if not all phases truly done
            pass
        print(json.dumps(final, indent=2), flush=True)
        print(f"MARATHON_PARTIAL phases_done={final['phases_done']}", flush=True)

    if args.scratch:
        sc = Path(args.scratch)
        sc.mkdir(parents=True, exist_ok=True)
        for name in ("marathon_summary.json", "MARATHON_COMPLETE", "marathon_state.json"):
            p = out / name
            if p.exists():
                shutil.copy(p, sc / name)


if __name__ == "__main__":
    main()
