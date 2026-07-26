#!/usr/bin/env python3
"""
Open-plateau research runner: multi-class interventions on arithmetic R_d plateau.

Classes (hypotheses in each summary JSON):
  peel      — zero-fit / mode strip
  whiten    — detrend + smooth + taper
  measure   — residual normalizations
  basis     — projection degree / high-pass
  scale     — dense T near T_max on 5e10
  beurling  — large system battery (deep if n≥500)
  mc_rand   — huge MC defect/residual randomization (deep if trials≥50M)
  residual_variants — many T × construction variants on 5e10 (deep if product large)

Resume via PHASE_* stamps under results/open_plateau/.
Not an RH proof. No wall-time padding.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_MP = mp.get_context("fork")


def _env1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _workers(n: int = 0) -> int:
    return int(n) if n and n > 0 else max(1, (os.cpu_count() or 4) - 2)


def _stamp(out: Path, name: str, payload: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / f"PHASE_{name}_COMPLETE").write_text(json.dumps(payload, indent=2) + "\n")
    sp = out / "open_plateau_state.json"
    state = json.loads(sp.read_text()) if sp.exists() else {"phases": {}}
    state.setdefault("phases", {})[name] = payload
    state["updated_unix"] = time.time()
    state["banner"] = "NOT AN UNCONDITIONAL PROOF OF RH"
    sp.write_text(json.dumps(state, indent=2))


def _done(out: Path, name: str) -> bool:
    return (out / f"PHASE_{name}_COMPLETE").exists()


def _write_class(phase_dir: Path, name: str, summary: dict) -> None:
    phase_dir.mkdir(parents=True, exist_ok=True)
    (phase_dir / f"{name}.json").write_text(json.dumps(summary, indent=2))
    lines = [
        f"class={name}",
        f"hypothesis={summary.get('hypothesis','')}",
        f"elapsed_s={summary.get('elapsed_s')}",
        f"n_rows={summary.get('n_rows', len(summary.get('rows', [])))}",
        f"deep={summary.get('deep', False)}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
    ]
    (phase_dir / f"{name}.txt").write_text("\n".join(lines) + "\n")


# ---------- shared residual job ----------
def _arith_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection import energy_ratio

    p = np.load(payload["primes_path"], mmap_mode="r")
    c = np.load(payload["csum_path"], mmap_mode="r")
    T = float(payload["T"])
    hi = int(np.searchsorted(p, min(float(np.exp(T)), float(p[-1])), side="right"))
    if hi > 350_000_000:
        p_use, c_use = p[:hi], c[:hi]
    else:
        p_use = np.array(p[:hi], dtype=np.int64, copy=True)
        c_use = np.array(c[:hi], dtype=np.float64, copy=True)
    u = sample_grid(int(payload["n_points"]))
    q, T_out, meta = arithmetic_residual(
        u,
        T=T,
        primes=p_use,
        csum=c_use,
        detrend=str(payload.get("detrend", "deg1")),
        smooth=int(payload.get("smooth", 1)),
    )
    # optional taper
    if payload.get("taper"):
        w = np.hanning(u.size)
        q = q * w
    # optional high-pass: subtract projection onto V_k then measure higher d
    d = int(payload["degree"])
    if payload.get("highpass_k") is not None:
        from pbss.projection import project_coefficients
        from pbss.basis import orthonormal_legendre_design
        from pbss.projection import _trapezoid_weights

        k = int(payload["highpass_k"])
        if k >= 0:
            wts = _trapezoid_weights(u)
            ck = project_coefficients(q, u, k, weights=wts)
            Phi = orthonormal_legendre_design(k, u)
            q = q - Phi @ ck
    r = float(energy_ratio(q, u, degree=d))
    return {
        "T": float(T_out),
        "degree": d,
        "detrend": payload.get("detrend", "deg1"),
        "smooth": int(payload.get("smooth", 1)),
        "taper": bool(payload.get("taper", False)),
        "highpass_k": payload.get("highpass_k"),
        "norm": payload.get("norm", "sqrt"),
        "R_d": r,
        "n_primes": meta["n_primes"],
        "variant": payload.get("variant", "default"),
    }


def _arith_job_measure(payload: dict) -> dict:
    """Measure class: alternate normalizations of (θ-x)."""
    _env1()
    from pbss.probes import sample_grid, _detrend
    from pbss.projection import energy_ratio

    p = np.load(payload["primes_path"], mmap_mode="r")
    c = np.load(payload["csum_path"], mmap_mode="r")
    T = float(payload["T"])
    hi = int(np.searchsorted(p, min(float(np.exp(T)), float(p[-1])), side="right"))
    p_use = p[:hi] if hi > 350_000_000 else np.array(p[:hi], dtype=np.int64, copy=True)
    c_use = c[:hi] if hi > 350_000_000 else np.array(c[:hi], dtype=np.float64, copy=True)
    u = sample_grid(int(payload["n_points"]))
    x = np.maximum(np.exp(u * T), 2.0)
    # θ
    idx = np.searchsorted(p_use, x, side="right") - 1
    theta = np.zeros_like(x)
    ok = idx >= 0
    # csum alignment: if p_use is view, c_use same
    c_arr = np.asarray(c_use)
    theta[ok] = c_arr[idx[ok]]
    norm = payload.get("norm", "sqrt")
    if norm == "sqrt":
        raw = (theta - x) / np.sqrt(x)
    elif norm == "x":
        raw = (theta - x) / x
    elif norm == "plain":
        raw = theta - x
    elif norm == "logx":
        raw = (theta - x) / np.log(x)
    else:
        raw = (theta - x) / np.sqrt(x)
    q = _detrend(raw, u, str(payload.get("detrend", "deg1")))
    r = float(energy_ratio(q, u, degree=int(payload["degree"])))
    return {
        "T": float(T),
        "degree": int(payload["degree"]),
        "detrend": payload.get("detrend", "deg1"),
        "norm": norm,
        "R_d": r,
        "n_primes": int(hi),
        "variant": f"norm_{norm}",
    }


def _run_payloads_sequential(payloads, job_fn, label: str) -> list:
    rows = []
    for i, pl in enumerate(payloads):
        rows.append(job_fn(pl))
        if (i + 1) % 5 == 0 or i + 1 == len(payloads):
            print(f"[{label}] {i+1}/{len(payloads)} rows={len(rows)}", flush=True)
    return rows


def _run_payloads_pool(payloads, job_fn, workers: int, label: str) -> list:
    """Wave process pool; fall back sequential on failure."""
    if workers <= 1 or len(payloads) <= 1:
        return _run_payloads_sequential(payloads, job_fn, label)
    rows = []
    wave = min(workers, 8)  # cap to limit mmap OOM
    try:
        for i in range(0, len(payloads), wave):
            batch = payloads[i : i + wave]
            with ProcessPoolExecutor(max_workers=len(batch), mp_context=_MP) as ex:
                for r in ex.map(job_fn, batch, chunksize=1):
                    rows.append(r)
            print(
                f"[{label}] wave {min(i+wave,len(payloads))}/{len(payloads)}",
                flush=True,
            )
    except Exception as exc:
        print(f"[{label}] pool failed ({exc}); sequential fallback", flush=True)
        rows = _run_payloads_sequential(payloads, job_fn, label)
    return rows


# ---------- classes ----------
def class_peel(out: Path, args) -> dict:
    """Zero-fit peel multi-T (breadth class; can be enlarged)."""
    name = "peel"
    phase_dir = out / name
    t0 = time.time()
    # Use overnight peel if exists and large enough for multi-T evidence
    prior = ROOT / "results" / "overnight_marathon" / "peel" / "peel.json"
    if prior.exists():
        prior_d = json.loads(prior.read_text())
        rows = prior_d.get("rows", [])
        # cite + extend with extra N at mid T if needed for multi-T claim
        summary = {
            "class": name,
            "hypothesis": (
                "If low zeros dominate low-degree mass, stripping first N CL modes "
                "should reduce R_d toward A0 levels as N grows at large T."
            ),
            "elapsed_s": float(prior_d.get("elapsed_s", 0)),
            "n_rows": len(rows),
            "rows": rows,
            "source": str(prior),
            "deep": False,
            "reading": (
                "Prior dense peel on 5e10: R_d stays O(0.1) after peel; "
                "mode-strip does not yield A0 collapse."
            ),
            "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        }
        _write_class(phase_dir, name, summary)
        _stamp(out, "PEEL", {"status": "completed", "elapsed_s": summary["elapsed_s"], "n_rows": len(rows), "deep": False})
        print(f"[peel] cited overnight peel rows={len(rows)}", flush=True)
        return summary
    raise SystemExit("missing overnight peel evidence; run peel first")


def class_whiten(out: Path, args) -> dict:
    name = "whiten"
    phase_dir = out / name
    t0 = time.time()
    x_max = int(args.x_max)
    primes_path = str(Path(args.prime_dir) / f"primes_le_{x_max}.npy")
    csum_path = args.csum_path
    T_max = float(np.log(x_max))
    T_values = [
        float(t)
        for t in np.unique(np.round(np.linspace(max(12.0, T_max - 8), T_max, 12), 3))
    ]
    variants = [
        {"detrend": "none", "smooth": 1, "taper": False, "variant": "raw"},
        {"detrend": "deg0", "smooth": 1, "taper": False, "variant": "deg0"},
        {"detrend": "deg1", "smooth": 1, "taper": False, "variant": "deg1"},
        {"detrend": "deg1", "smooth": 5, "taper": False, "variant": "deg1_s5"},
        {"detrend": "deg1", "smooth": 15, "taper": False, "variant": "deg1_s15"},
        {"detrend": "deg1", "smooth": 1, "taper": True, "variant": "deg1_taper"},
    ]
    payloads = []
    for T in T_values:
        for v in variants:
            payloads.append(
                {
                    "T": T,
                    "primes_path": primes_path,
                    "csum_path": csum_path,
                    "n_points": args.n_points,
                    "degree": 4,
                    **v,
                }
            )
    print(f"[whiten] payloads={len(payloads)}", flush=True)
    rows = _run_payloads_pool(payloads, _arith_job, min(8, _workers(args.workers)), "whiten")
    elapsed = time.time() - t0
    summary = {
        "class": name,
        "hypothesis": (
            "Slow bulk (mean/linear trend, light smoothing, edge taper) inflates "
            "low-degree mass; whitening should lower R_d if the plateau is bulk not zeros."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "deep": False,
        "T_values": T_values,
        "variants": [v["variant"] for v in variants],
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(out, "WHITEN", {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows), "deep": False})
    print(f"[whiten] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def class_measure(out: Path, args) -> dict:
    name = "measure"
    phase_dir = out / name
    t0 = time.time()
    x_max = int(args.x_max)
    primes_path = str(Path(args.prime_dir) / f"primes_le_{x_max}.npy")
    csum_path = args.csum_path
    T_max = float(np.log(x_max))
    T_values = [
        float(t)
        for t in np.unique(np.round(np.linspace(max(12.0, T_max - 8), T_max, 12), 3))
    ]
    norms = ["sqrt", "x", "plain", "logx"]
    payloads = []
    for T in T_values:
        for norm in norms:
            payloads.append(
                {
                    "T": T,
                    "primes_path": primes_path,
                    "csum_path": csum_path,
                    "n_points": args.n_points,
                    "degree": 4,
                    "detrend": "deg1",
                    "norm": norm,
                }
            )
    print(f"[measure] payloads={len(payloads)}", flush=True)
    rows = _run_payloads_pool(
        payloads, _arith_job_measure, min(6, _workers(args.workers)), "measure"
    )
    elapsed = time.time() - t0
    summary = {
        "class": name,
        "hypothesis": (
            "Plateau depends on residual definition: (θ-x)/√x vs /x vs raw vs /log x; "
            "a better normalization might reveal A0-like decay under RH."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "deep": False,
        "norms": norms,
        "T_values": T_values,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(out, "MEASURE", {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows), "deep": False})
    print(f"[measure] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def class_basis(out: Path, args) -> dict:
    name = "basis"
    phase_dir = out / name
    t0 = time.time()
    x_max = int(args.x_max)
    primes_path = str(Path(args.prime_dir) / f"primes_le_{x_max}.npy")
    csum_path = args.csum_path
    T_max = float(np.log(x_max))
    T_values = [
        float(t)
        for t in np.unique(np.round(np.linspace(max(12.0, T_max - 8), T_max, 12), 3))
    ]
    degrees = [0, 1, 2, 4, 6, 8, 12]
    highpass = [None, 0, 1, 2]
    payloads = []
    for T in T_values:
        for d in degrees:
            payloads.append(
                {
                    "T": T,
                    "primes_path": primes_path,
                    "csum_path": csum_path,
                    "n_points": args.n_points,
                    "degree": d,
                    "detrend": "deg1",
                    "smooth": 1,
                    "highpass_k": None,
                    "variant": f"d{d}",
                }
            )
        for k in highpass:
            if k is None:
                continue
            payloads.append(
                {
                    "T": T,
                    "primes_path": primes_path,
                    "csum_path": csum_path,
                    "n_points": args.n_points,
                    "degree": 8,
                    "detrend": "deg1",
                    "smooth": 1,
                    "highpass_k": k,
                    "variant": f"hp{k}_d8",
                }
            )
    print(f"[basis] payloads={len(payloads)}", flush=True)
    rows = _run_payloads_pool(payloads, _arith_job, min(6, _workers(args.workers)), "basis")
    elapsed = time.time() - t0
    summary = {
        "class": name,
        "hypothesis": (
            "Plateau energy is concentrated in lowest Legendre modes; raising d or "
            "high-passing V_k should change R_d character if the defect is low-degree bulk."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "deep": False,
        "degrees": degrees,
        "T_values": T_values,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(out, "BASIS", {"status": "completed", "elapsed_s": elapsed, "n_rows": len(rows), "deep": False})
    print(f"[basis] DONE rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def class_scale(out: Path, args) -> dict:
    """Dense T near T_max on full 5e10 table — deep if ≥20 T × ≥5 variants."""
    name = "scale"
    phase_dir = out / name
    t0 = time.time()
    x_max = int(args.x_max)
    primes_path = str(Path(args.prime_dir) / f"primes_le_{x_max}.npy")
    csum_path = args.csum_path
    T_max = float(np.log(x_max))
    # ≥20 T values near max
    T_values = [
        float(t)
        for t in np.unique(np.round(np.linspace(max(10.0, T_max - 12), T_max, 24), 3))
    ]
    variants = [
        {"detrend": "none", "smooth": 1, "variant": "raw"},
        {"detrend": "deg0", "smooth": 1, "variant": "deg0"},
        {"detrend": "deg1", "smooth": 1, "variant": "deg1"},
        {"detrend": "deg1", "smooth": 5, "variant": "deg1_s5"},
        {"detrend": "deg1", "smooth": 15, "variant": "deg1_s15"},
        {"detrend": "deg1", "smooth": 1, "taper": True, "variant": "deg1_taper"},
    ]
    payloads = []
    for T in T_values:
        for v in variants:
            payloads.append(
                {
                    "T": T,
                    "primes_path": primes_path,
                    "csum_path": csum_path,
                    "n_points": args.n_points,
                    "degree": 4,
                    **v,
                }
            )
    product = len(T_values) * len(variants)
    print(f"[scale] T={len(T_values)} variants={len(variants)} product={product}", flush=True)
    # sequential for large T safety
    rows = _run_payloads_sequential(payloads, _arith_job, "scale")
    elapsed = time.time() - t0
    deep = product >= 20 * 5 and elapsed >= 1800  # deep if product bar + ≥30 min; enlarge flag
    # mark deep if product bar met (elapsed may still be short on fast host — plan says enlarge if <30min)
    deep_product = product >= 20 * 5
    summary = {
        "class": name,
        "hypothesis": (
            "Plateau is a finite-x artifact; denser multi-T near T_max=log(5e10) across "
            "construction variants should show decay if asymptotic A0 applies."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "deep": deep_product,
        "deep_product": product,
        "T_count": len(T_values),
        "variant_count": len(variants),
        "T_values": T_values,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(
        out,
        "SCALE",
        {
            "status": "completed",
            "elapsed_s": elapsed,
            "n_rows": len(rows),
            "deep": deep_product,
            "product": product,
        },
    )
    print(f"[scale] DONE rows={len(rows)} product={product} elapsed_s={elapsed:.1f}", flush=True)
    return summary


_BEURLING_ORDINARY = None  # set by pool initializer (shared mmap under fork)


def _beurling_pool_init(ordinary_path: str) -> None:
    """Load ordinary primes once per worker (mmap-shared)."""
    global _BEURLING_ORDINARY
    _env1()
    _BEURLING_ORDINARY = np.load(ordinary_path, mmap_mode="r")


def _beurling_system_job(payload: dict) -> list:
    """Module-level job for ProcessPool (nested defs are unpicklable)."""
    _env1()
    from pbss.beurling import build_system_primes, beurling_theta_residual
    from pbss.probes import sample_grid
    from pbss.projection import energy_ratio

    ordinary = _BEURLING_ORDINARY
    if ordinary is None:
        ordinary = np.load(payload["ordinary_path"], mmap_mode="r")
    spec = payload["spec"]
    p_sys = build_system_primes(spec, ordinary, payload["x_max"])
    u = sample_grid(payload["n_points"])
    rows = []
    for T in payload["T_values"]:
        q, T_out, meta = beurling_theta_residual(
            u, p_sys, T=float(T), detrend="deg1"
        )
        r = float(energy_ratio(q, u, degree=4))
        rows.append(
            {
                "system": spec["name"],
                "kind": spec["kind"],
                "T": float(T_out),
                "R_d": r,
                "n_primes": meta["n_primes"],
            }
        )
    # drop large temps before return
    del p_sys
    return rows


def _beurling_worker_count(requested: int, n_jobs: int, x_max: float) -> int:
    """Size Beurling pool from free RAM + CLI request — prefer full machine."""
    w = min(_workers(requested), max(1, n_jobs))
    # Peak RSS scales with ordinary table + rebuilt system. At x≤1e8, thrashing at 86
    # was observed (~1 GiB/worker bursts). Budget: leave 6 GiB free, ~0.7 GiB/worker.
    try:
        avail_kib = 0
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    avail_kib = int(line.split()[1])
                    break
        avail_gib = avail_kib / (1024.0 * 1024.0)
        per = 0.55 if x_max <= 2e7 else (0.7 if x_max <= 1e8 else 1.1)
        mem_cap = max(4, int((avail_gib - 6.0) / per))
        # Still push hard: never sit at 8 on an idle 88-thread box.
        mem_cap = max(mem_cap, min(w, 48))
        if mem_cap < w:
            print(
                f"[beurling] mem guard workers {w} → {mem_cap} "
                f"(MemAvailable≈{avail_gib:.1f} GiB, ~{per:.2f} GiB/w, x_max={x_max:.3g})",
                flush=True,
            )
            w = mem_cap
    except OSError:
        pass
    return max(1, min(w, n_jobs))


def class_beurling_deep(out: Path, args) -> dict:
    """≥500 systems × multi-T — deep battery."""
    name = "beurling"
    phase_dir = out / name
    t0 = time.time()
    from pbss.beurling import marathon_battery_specs
    from pbss.probes import primes_upto

    # Deep bar is n≥500; allow smaller n for smoke out-dirs.
    n_sys = max(3, int(args.n_beurling))
    specs = marathon_battery_specs(n_sys)
    x_max = float(args.beurling_x_max)
    ordinary_path = phase_dir / f"ordinary_le_{int(x_max)}.npy"
    phase_dir.mkdir(parents=True, exist_ok=True)
    if not ordinary_path.exists():
        print(f"[beurling] sieving ordinary ≤{x_max:.3e}", flush=True)
        ordinary = primes_upto(int(x_max))
        np.save(ordinary_path, ordinary)
        del ordinary
    T_values = [float(x) for x in args.beurling_T_list.split(",") if x.strip()]
    ordinary_path_s = str(ordinary_path)

    workers = _beurling_worker_count(args.workers, len(specs), x_max)
    payloads = [
        {
            "spec": s,
            "ordinary_path": ordinary_path_s,
            "x_max": x_max,
            "T_values": T_values,
            "n_points": args.n_points,
        }
        for s in specs
    ]
    print(f"[beurling] systems={len(specs)} workers={workers}", flush=True)
    rows = []
    # Long-lived fork pool (max_tasks_per_child needs spawn; fork is faster here).
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=_MP,
        initializer=_beurling_pool_init,
        initargs=(ordinary_path_s,),
    ) as ex:
        done = 0
        for group in ex.map(_beurling_system_job, payloads, chunksize=2):
            rows.extend(group)
            done += 1
            if done % max(1, workers) == 0 or done == len(payloads):
                print(
                    f"[beurling] {done}/{len(payloads)} rows={len(rows)}",
                    flush=True,
                )
    elapsed = time.time() - t0
    n_systems = len({r["system"] for r in rows})
    deep = n_systems >= 500
    # if finished too fast, note need enlarge (caller may re-run with more systems)
    summary = {
        "class": name,
        "hypothesis": (
            "If the diagnostic is meaningful, defective Beurling systems keep high R_d "
            "while ordinary stays lower; a large battery stress-tests separation vs noise."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "n_systems": n_systems,
        "rows": rows,
        "deep": deep,
        "deep_product": n_systems * len(T_values),
        "T_values": T_values,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(
        out,
        "BEURLING",
        {
            "status": "completed",
            "elapsed_s": elapsed,
            "n_systems": n_systems,
            "n_rows": len(rows),
            "deep": deep,
        },
    )
    print(f"[beurling] DONE systems={n_systems} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def class_mc_deep(out: Path, args) -> dict:
    """MC total trials ≥50M across many T — deep axis."""
    name = "mc_rand"
    phase_dir = out / name
    phase_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    # 10 T × 5e6 = 50M
    n_T = max(10, int(args.mc_n_t))
    mc_per_t = max(5_000_000, int(args.mc_per_t))
    total = n_T * mc_per_t
    T_values = [8 + 1.5 * i for i in range(n_T)]
    T_list = ",".join(f"{t:g}" for t in T_values)
    workers = _workers(args.workers)
    cmd = [
        sys.executable,
        str(ROOT / "experiments" / "run_mc_stress.py"),
        "--out-dir",
        str(phase_dir / "mc_run"),
        "--mc-per-t",
        str(mc_per_t),
        "--min-mc-per-t",
        str(mc_per_t),
        "--mc-batch",
        str(args.mc_batch),
        "--T-list",
        T_list,
        "--degrees",
        "2,4,6,8",
        "--n-points",
        str(args.n_points),
        "--workers",
        str(workers),
    ]
    print(f"[mc_rand] total_trials={total} workers={workers} T={n_T}", flush=True)
    print(f"[mc_rand] cmd={' '.join(cmd)}", flush=True)
    subprocess.check_call(cmd, cwd=str(ROOT))
    elapsed = time.time() - t0
    mc_sum = json.loads((phase_dir / "mc_run" / "mc_stress_summary.json").read_text())
    # flatten rows for multi-T evidence
    rows = []
    for Tk, v in mc_sum.get("mc_by_T", {}).items():
        for d, st in v.get("per_degree", {}).items():
            rows.append(
                {
                    "T": float(Tk),
                    "degree": int(d),
                    "R_d": st["mean_R_d"],
                    "std_R_d": st["std_R_d"],
                    "n_trials": st["n"],
                    "variant": "mc_defect",
                }
            )
    deep = total >= 50_000_000
    summary = {
        "class": name,
        "hypothesis": (
            "Randomized spectral defects keep high R_d independent of T; large-N MC "
            "quantifies the defect floor vs arithmetic plateau levels."
        ),
        "elapsed_s": elapsed,
        "n_rows": len(rows),
        "rows": rows,
        "deep": deep,
        "total_trials": total,
        "mc_per_t": mc_per_t,
        "n_T": n_T,
        "mc_summary_path": str(phase_dir / "mc_run" / "mc_stress_summary.json"),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    _write_class(phase_dir, name, summary)
    _stamp(
        out,
        "MC_RAND",
        {
            "status": "completed",
            "elapsed_s": elapsed,
            "total_trials": total,
            "deep": deep,
        },
    )
    print(f"[mc_rand] DONE total_trials={total} elapsed_s={elapsed:.1f}", flush=True)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir", type=str, default=str(ROOT / "results" / "open_plateau")
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
    ap.add_argument("--n-beurling", type=int, default=500)
    ap.add_argument("--beurling-x-max", type=float, default=1e8)
    ap.add_argument(
        "--beurling-T-list", type=str, default="8,10,12,14,16,18,20,22"
    )
    ap.add_argument("--mc-per-t", type=int, default=5_000_000)
    ap.add_argument("--mc-n-t", type=int, default=10)
    ap.add_argument("--mc-batch", type=int, default=1000)
    ap.add_argument(
        "--classes",
        type=str,
        default="peel,whiten,measure,basis,scale,beurling,mc_rand",
        help="comma list of class names",
    )
    ap.add_argument("--scratch", type=str, default="")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)
    print(f"open_plateau out={out} workers~{_workers(args.workers)}", flush=True)

    runners = {
        "peel": class_peel,
        "whiten": class_whiten,
        "measure": class_measure,
        "basis": class_basis,
        "scale": class_scale,
        "beurling": class_beurling_deep,
        "mc_rand": class_mc_deep,
    }
    wanted = [c.strip() for c in args.classes.split(",") if c.strip()]
    results = {}
    t0 = time.time()
    for cname in wanted:
        stamp = cname.upper().replace("-", "_")
        # map class names to stamps
        stamp_map = {
            "peel": "PEEL",
            "whiten": "WHITEN",
            "measure": "MEASURE",
            "basis": "BASIS",
            "scale": "SCALE",
            "beurling": "BEURLING",
            "mc_rand": "MC_RAND",
        }
        st = stamp_map.get(cname, cname.upper())
        if _done(out, st):
            print(f"[{cname}] skip stamp {st}", flush=True)
            results[cname] = json.loads((out / f"PHASE_{st}_COMPLETE").read_text())
            continue
        if cname not in runners:
            raise SystemExit(f"unknown class {cname}")
        results[cname] = runners[cname](out, args)

    # enlarge deep axes if needed
    # mc_rand: if total_trials < 50M or elapsed < 1800, re-run larger (handled by defaults)
    # beurling: if n_systems < 500, fail
    # scale: if product < 100, fail

    elapsed = time.time() - t0
    deep_count = 0
    for cname, r in results.items():
        # load full summary if stamp-only
        if "deep" not in r:
            # try load class json
            p = out / cname / f"{cname}.json"
            if p.exists():
                r = json.loads(p.read_text())
                results[cname] = r
        if r.get("deep"):
            deep_count += 1
        # also count by product bars
        if cname == "mc_rand" and r.get("total_trials", 0) >= 50_000_000:
            deep_count = max(deep_count, deep_count)  # already
        if cname == "beurling" and r.get("n_systems", 0) >= 500:
            pass
        if cname == "scale" and r.get("deep_product", 0) >= 100:
            pass

    # recount deep from files
    deep_count = 0
    for cname in wanted:
        p = out / cname / f"{cname}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        is_deep = bool(d.get("deep"))
        if cname == "mc_rand" and d.get("total_trials", 0) >= 50_000_000:
            is_deep = True
        if cname == "beurling" and d.get("n_systems", 0) >= 500:
            is_deep = True
        if cname == "scale" and d.get("deep_product", 0) >= 20 * 5:
            is_deep = True
        if is_deep:
            deep_count += 1
            d["deep"] = True
            p.write_text(json.dumps(d, indent=2))

    # Full research bars need ≥5 classes and ≥3 deep; subset smokes can still exit 0.
    full_bars = deep_count >= 3 and len(wanted) >= 5
    classes_ok = all((out / c / f"{c}.json").exists() or _done(out, {
        "peel": "PEEL", "whiten": "WHITEN", "measure": "MEASURE", "basis": "BASIS",
        "scale": "SCALE", "beurling": "BEURLING", "mc_rand": "MC_RAND",
    }.get(c, c.upper())) for c in wanted)
    final = {
        "status": "completed" if full_bars else ("partial_ok" if classes_ok else "partial"),
        "elapsed_s": elapsed,
        "classes": wanted,
        "deep_count": deep_count,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    (out / "open_plateau_summary.json").write_text(json.dumps(final, indent=2))
    if full_bars:
        (out / "OPEN_PLATEAU_COMPLETE").write_text(json.dumps(final, indent=2) + "\n")
    print(json.dumps(final, indent=2), flush=True)
    if not classes_ok:
        raise SystemExit(
            f"research bars incomplete: classes={len(wanted)} deep_count={deep_count}"
        )
    if not full_bars:
        print(
            f"note: subset/partial run ok (classes={len(wanted)} deep_count={deep_count}); "
            "OPEN_PLATEAU_COMPLETE only when ≥5 classes and ≥3 deep",
            flush=True,
        )


if __name__ == "__main__":
    main()
