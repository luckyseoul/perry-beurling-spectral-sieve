#!/usr/bin/env python3
"""
Unattended large-T PBSS campaign.

1) Choose max x_max (try 1e9, fall back 1e8 / lower if memory fails)
2) Segmented sieve once
3) Multi-T × detrend × smooth arithmetic R_d (multi-core)
4) Controls at each T: critical-line, persistent defect, off-critical
5) Write JSON/TXT + plots; completion stamp; clean exit

Uses shipped arithmetic_residual + project only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

# Prefer fork so large primes arrays are COW-shared, not pickled per worker
_MP_CTX = mp.get_context("fork")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _env_single_thread() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _avail_ram_bytes() -> int:
    try:
        import psutil

        return int(psutil.virtual_memory().available)
    except Exception:
        pass
    # /proc/meminfo MemAvailable
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except Exception:
        pass
    return 32 * 1024**3  # assume 32 GiB


def choose_xmax(candidates: list[float]) -> tuple[float, str]:
    """Pick largest x_max that fits a comfortable RAM budget for primes storage."""
    avail = _avail_ram_bytes()
    # rough: primes ≤ n use ~ (n/log n)*8 bytes; segmented sieve peak ~ segment + √n
    for x in sorted(candidates, reverse=True):
        n = int(x)
        n_primes_est = n / max(np.log(n), 2.0)
        need = n_primes_est * 8 * 2.5 + 500e6  # safety + scratch
        if need < 0.45 * avail:
            return (
                float(n),
                f"selected x_max={n:.3e}; est_primes~{n_primes_est:.3e}; "
                f"avail_RAM={avail/1e9:.1f}GiB",
            )
    return float(candidates[-1]), "fallback to smallest candidate"


def _arith_job(payload: dict) -> dict:
    _env_single_thread()
    from pbss.probes import arithmetic_residual, sample_grid
    from pbss.projection import project

    T = float(payload["T"])
    degree = int(payload["degree"])
    n_points = int(payload["n_points"])
    detrend = str(payload["detrend"])
    smooth = int(payload["smooth"])
    primes = np.asarray(payload["primes"], dtype=np.int64)

    u = sample_grid(n_points)
    q, T_ar, meta = arithmetic_residual(
        u, T=T, primes=primes, detrend=detrend, smooth=smooth
    )
    r = project(q, u, degree=degree, T=T_ar)
    return {
        "kind": "arithmetic",
        "T": T_ar,
        "x_max": meta["x_max"],
        "detrend": detrend,
        "smooth": smooth,
        "n_primes": meta["n_primes"],
        "R_d": r.energy_ratio,
        "S_d": r.scaled_strength,
        "energy": r.energy,
        "l2_norm_sq": r.l2_norm_sq,
    }


def _control_job(payload: dict) -> dict:
    _env_single_thread()
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
        probe_persistent_defect(u, eps=eps, j=0, waves=80), u, degree=degree, T=T
    )
    r_off = project(
        probe_off_critical_mode(u, T=T, sigma=sigma), u, degree=degree, T=T
    )
    return {
        "kind": "controls",
        "T": T,
        "critical_line": {"R_d": r_cl.energy_ratio, "S_d": r_cl.scaled_strength},
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
    }


def make_plots(summary: dict, out_dir: Path) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = []
    arith = summary["arithmetic_rows"]
    controls = {c["T"]: c for c in summary["control_rows"]}

    # group by (detrend, smooth)
    configs = sorted({(r["detrend"], r["smooth"]) for r in arith})
    fig, ax = plt.subplots(figsize=(10, 6), dpi=140)
    for det, sm in configs:
        rows = [r for r in arith if r["detrend"] == det and r["smooth"] == sm]
        rows = sorted(rows, key=lambda z: z["T"])
        Ts = [r["T"] for r in rows]
        Rs = [r["R_d"] for r in rows]
        ax.plot(Ts, Rs, "o-", label=f"arith detrend={det} smooth={sm}", linewidth=1.5)

    if controls:
        Ts_c = sorted(controls.keys())
        ax.plot(
            Ts_c,
            [controls[t]["critical_line"]["R_d"] for t in Ts_c],
            "s--",
            label="critical-line (A0)",
            color="green",
        )
        ax.plot(
            Ts_c,
            [controls[t]["persistent_defect"]["R_d"] for t in Ts_c],
            "^--",
            label="persistent defect (B0)",
            color="red",
        )
        ax.plot(
            Ts_c,
            [controls[t]["off_critical"]["R_d"] for t in Ts_c],
            "d--",
            label=f"off-critical σ={summary['sigma']}",
            color="orange",
        )

    ax.set_xlabel("T = log(x_max)")
    ax.set_ylabel(f"R_d  (d={summary['degree']})")
    ax.set_title("PBSS overnight: arithmetic residual vs controls")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    p1 = out_dir / "overnight_Rd_vs_T.png"
    fig.savefig(p1)
    plt.close(fig)
    paths.append(str(p1))

    # linear-scale focus on arithmetic deg1 smooth=1 if present
    focus = [r for r in arith if r["detrend"] == "deg1" and r["smooth"] == 1]
    if focus:
        focus = sorted(focus, key=lambda z: z["T"])
        fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
        ax.plot([r["T"] for r in focus], [r["R_d"] for r in focus], "o-", color="C0")
        if controls:
            Ts_c = sorted(controls.keys())
            ax.axhline(
                controls[Ts_c[0]]["persistent_defect"]["R_d"],
                color="red",
                ls="--",
                label="defect floor",
            )
        ax.set_xlabel("T")
        ax.set_ylabel("R_d arithmetic (deg1, smooth=1)")
        ax.set_title("Arithmetic residual plateau / trend (linear scale)")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        p2 = out_dir / "overnight_arith_deg1_linear.png"
        fig.savefig(p2)
        plt.close(fig)
        paths.append(str(p2))

    return paths


def main() -> None:
    _env_single_thread()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument("--sigma", type=float, default=0.75)
    ap.add_argument(
        "--xmax-candidates",
        type=str,
        default="1e9,5e8,1e8,5e7",
        help="try largest first",
    )
    ap.add_argument(
        "--T-list",
        type=str,
        default="",
        help="override T grid; default auto from x_max",
    )
    ap.add_argument(
        "--detrends",
        type=str,
        default="none,deg0,deg1",
    )
    ap.add_argument(
        "--smooths",
        type=str,
        default="1,5,15",
    )
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--scratch", type=str, default="")
    args = ap.parse_args()

    scratch = Path(args.scratch) if args.scratch else None
    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)

    candidates = [float(x) for x in args.xmax_candidates.split(",") if x.strip()]
    x_max, xmax_note = choose_xmax(candidates)
    print(xmax_note, flush=True)

    # auto T grid: denser at large end
    if args.T_list.strip():
        T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    else:
        T_max = float(np.log(x_max))
        # from ~8 to T_max
        T_values = sorted(
            set(
                list(np.linspace(8.0, min(16.0, T_max), 6))
                + list(np.linspace(min(16.0, T_max), T_max, 8))
            )
        )
        T_values = [float(t) for t in T_values if np.exp(t) <= x_max * 1.001]

    detrends = [s.strip() for s in args.detrends.split(",") if s.strip()]
    smooths = [int(s) for s in args.smooths.split(",") if s.strip()]

    stamp = {
        "started_unix": time.time(),
        "started": time.strftime("%Y-%m-%d %H:%M:%S"),
        "x_max": x_max,
        "xmax_note": xmax_note,
        "T_values": T_values,
        "detrends": detrends,
        "smooths": smooths,
        "workers": workers,
        "degree": args.degree,
        "n_points": args.n_points,
        "status": "running",
    }
    if scratch:
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "campaign_status.json").write_text(json.dumps(stamp, indent=2))

    # --- sieve ---
    from pbss.probes import primes_upto

    print(f"Sieving primes ≤ {x_max:.3e} (segmented if large)…", flush=True)
    t0 = time.time()
    try:
        primes = primes_upto(int(x_max))
    except MemoryError:
        # step down
        for fb in [1e8, 5e7, 2e7]:
            if fb >= x_max:
                continue
            print(f"OOM at {x_max:.3e}, falling back to {fb:.3e}", flush=True)
            x_max = fb
            T_values = [t for t in T_values if np.exp(t) <= x_max * 1.001]
            try:
                primes = primes_upto(int(x_max))
                xmax_note += f" | OOM fallback to {x_max:.3e}"
                break
            except MemoryError:
                continue
        else:
            raise
    sieve_s = time.time() - t0
    print(f"  n_primes={primes.size} in {sieve_s:.1f}s", flush=True)

    # --- arithmetic jobs ---
    arith_payloads = []
    for T in T_values:
        for det in detrends:
            for sm in smooths:
                arith_payloads.append(
                    {
                        "T": T,
                        "degree": args.degree,
                        "n_points": args.n_points,
                        "detrend": det,
                        "smooth": sm,
                        "primes": primes,
                    }
                )

    print(
        f"Arithmetic jobs: {len(arith_payloads)} "
        f"(T×detrend×smooth); workers={workers}",
        flush=True,
    )
    arith_rows = []
    t1 = time.time()
    # ProcessPool: on Linux fork, primes array is COW-friendly if not written
    with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CTX) as ex:
        futs = [ex.submit(_arith_job, p) for p in arith_payloads]
        done = 0
        for f in as_completed(futs):
            row = f.result()
            arith_rows.append(row)
            done += 1
            if done % max(1, len(arith_payloads) // 20) == 0 or done == len(arith_payloads):
                print(
                    f"  arith {done}/{len(arith_payloads)}  "
                    f"last T={row['T']:.2f} {row['detrend']}/s{row['smooth']} "
                    f"R_d={row['R_d']:.4e}",
                    flush=True,
                )
    print(f"Arithmetic done in {time.time()-t1:.1f}s", flush=True)

    # --- controls ---
    ctrl_payloads = [
        {
            "T": T,
            "degree": args.degree,
            "n_points": args.n_points,
            "eps": args.eps,
            "sigma": args.sigma,
        }
        for T in T_values
    ]
    print(f"Control jobs: {len(ctrl_payloads)}", flush=True)
    control_rows = []
    with ProcessPoolExecutor(
        max_workers=min(workers, max(4, len(ctrl_payloads))),
        mp_context=_MP_CTX,
    ) as ex:
        for row in ex.map(_control_job, ctrl_payloads):
            control_rows.append(row)
            print(
                f"  controls T={row['T']:.2f}  "
                f"R_cl={row['critical_line']['R_d']:.3e}  "
                f"R_def={row['persistent_defect']['R_d']:.3e}  "
                f"R_off={row['off_critical']['R_d']:.3e}",
                flush=True,
            )

    # --- analysis ---
    focus = [
        r
        for r in arith_rows
        if r["detrend"] == "deg1" and r["smooth"] == 1
    ]
    focus = sorted(focus, key=lambda z: z["T"])
    if len(focus) >= 2:
        R = np.array([r["R_d"] for r in focus])
        T_arr = np.array([r["T"] for r in focus])
        slope = float(np.polyfit(np.log(T_arr), np.log(np.maximum(R, 1e-30)), 1)[0])
        reading = {
            "focus": "detrend=deg1 smooth=1",
            "R_first": float(R[0]),
            "R_last": float(R[-1]),
            "R_min": float(R.min()),
            "R_max": float(R.max()),
            "loglog_slope": slope,
            "decreases": bool(R[-1] < R[0]),
            "plateau_std_last_half": float(np.std(R[len(R) // 2 :])),
        }
        if reading["decreases"] and slope < -0.3:
            narrative = (
                f"Focus arithmetic R_d decreases over T "
                f"({reading['R_first']:.4e} → {reading['R_last']:.4e}, "
                f"log-log slope {slope:.2f}). Not a proof of Theorem A / RH."
            )
        elif reading["plateau_std_last_half"] < 0.02:
            narrative = (
                f"Focus arithmetic R_d plateaus near {reading['R_last']:.3f} "
                f"at large T (max T={focus[-1]['T']:.2f}, "
                f"x_max={focus[-1]['x_max']:.3e}). "
                "No A0-style decay observed for arithmetic residual at this scale. "
                "Not a proof of RH."
            )
        else:
            narrative = (
                f"Focus arithmetic R_d drifts without clean decay "
                f"({reading['R_first']:.4e} → {reading['R_last']:.4e}). "
                "Not a proof of RH or full Theorem A."
            )
    else:
        reading = {}
        narrative = "Insufficient focus rows for trend analysis."

    summary = {
        "campaign": "overnight_large_T",
        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": time.time() - stamp["started_unix"],
        "x_max": x_max,
        "xmax_note": xmax_note,
        "n_primes": int(primes.size),
        "sieve_seconds": sieve_s,
        "degree": args.degree,
        "n_points": args.n_points,
        "eps": args.eps,
        "sigma": args.sigma,
        "workers": workers,
        "T_values": T_values,
        "detrends": detrends,
        "smooths": smooths,
        "arithmetic_rows": arith_rows,
        "control_rows": control_rows,
        "reading": reading,
        "narrative": narrative,
        "not_an_RH_proof": True,
        "not_full_theorem_A": True,
    }

    json_path = results_dir / "overnight_campaign.json"
    txt_path = results_dir / "overnight_campaign.txt"
    json_path.write_text(json.dumps(summary, indent=2))

    lines = [
        "PBSS overnight large-T campaign",
        f"x_max={x_max:.3e}  n_primes={primes.size}  sieve_s={sieve_s:.1f}",
        f"d={args.degree} n_points={args.n_points} workers={workers}",
        f"T grid: {T_values}",
        f"detrends={detrends}  smooths={smooths}",
        "",
        narrative,
        "",
        "Focus (deg1, smooth=1):",
    ]
    for r in focus:
        lines.append(f"  T={r['T']:7.3f}  x={r['x_max']:.3e}  R_d={r['R_d']:.6e}")
    lines += ["", "Controls (sample):"]
    for c in sorted(control_rows, key=lambda z: z["T"])[:: max(1, len(control_rows)//5)]:
        lines.append(
            f"  T={c['T']:7.3f}  R_cl={c['critical_line']['R_d']:.3e}  "
            f"R_def={c['persistent_defect']['R_d']:.3e}  "
            f"R_off={c['off_critical']['R_d']:.3e}"
        )
    lines += [
        "",
        "Not an unconditional RH proof; not a complete proof of Theorem A.",
        f"Elapsed {summary['elapsed_s']:.1f}s",
    ]
    text = "\n".join(lines) + "\n"
    txt_path.write_text(text)
    print(text)

    plot_paths = make_plots(summary, results_dir)
    summary["plot_paths"] = plot_paths
    json_path.write_text(json.dumps(summary, indent=2))
    print("Plots:", plot_paths, flush=True)

    if scratch:
        (scratch / "overnight_campaign.json").write_text(json.dumps(summary, indent=2))
        (scratch / "overnight_campaign.txt").write_text(text)
        for p in plot_paths:
            # copy plot basenames
            src = Path(p)
            if src.exists():
                (scratch / src.name).write_bytes(src.read_bytes())
        stamp["status"] = "completed"
        stamp["finished"] = summary["finished"]
        stamp["elapsed_s"] = summary["elapsed_s"]
        stamp["x_max"] = x_max
        stamp["narrative"] = narrative
        (scratch / "campaign_status.json").write_text(json.dumps(stamp, indent=2))
        (scratch / "CAMPAIGN_COMPLETE").write_text(
            f"OK {summary['finished']}\n{narrative}\n"
        )

    print("CAMPAIGN COMPLETE — workers exited", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        # still write failure stamp if possible
        scratch = None
        for i, a in enumerate(sys.argv):
            if a == "--scratch" and i + 1 < len(sys.argv):
                scratch = Path(sys.argv[i + 1])
        if scratch:
            scratch.mkdir(parents=True, exist_ok=True)
            (scratch / "CAMPAIGN_FAILED").write_text(traceback.format_exc())
        raise
