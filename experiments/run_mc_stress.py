#!/usr/bin/env python3
"""
MC / ablation stress campaign: ≥50_000 spectral-defect trials per T.

Resume-capable state under results/mc_stress/. Not an RH proof.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_MP_CTX = mp.get_context("fork")


def _env1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _mc_batch_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import probe_defective, sample_grid
    from pbss.projection import energy_ratio

    n = int(payload["n_points"])
    degrees = list(payload["degrees"])
    n_trials = int(payload["n_trials"])
    seed0 = int(payload["seed0"])
    T = float(payload["T"])
    u = sample_grid(n)
    rng = np.random.default_rng(seed0)
    sums = {d: 0.0 for d in degrees}
    sumsq = {d: 0.0 for d in degrees}
    for _ in range(n_trials):
        weight = float(rng.uniform(0.5, 3.0))
        waves = int(rng.integers(30, 120))
        deg_def = int(rng.integers(0, 3))
        q = probe_defective(u, waves=waves, defect_degree=deg_def, defect_weight=weight)
        for d in degrees:
            r = float(energy_ratio(q, u, degree=d))
            sums[d] += r
            sumsq[d] += r * r
    out = {}
    for d in degrees:
        mean = sums[d] / n_trials
        var = max(0.0, sumsq[d] / n_trials - mean * mean)
        out[str(d)] = {"mean_R_d": mean, "std_R_d": float(math.sqrt(var)), "n": n_trials}
    return {"T": T, "seed0": seed0, "n_trials": n_trials, "per_degree": out}


def _control_job(payload: dict) -> dict:
    _env1()
    from pbss.probes import (
        probe_critical_line_mode,
        probe_persistent_defect,
        sample_grid,
    )
    from pbss.projection import energy_ratio

    u = sample_grid(int(payload["n_points"]))
    T = float(payload["T"])
    d = int(payload["degree"])
    return {
        "T": T,
        "degree": d,
        "critical_line_R_d": float(
            energy_ratio(probe_critical_line_mode(u, T=T), u, degree=d)
        ),
        "persistent_defect_R_d": float(
            energy_ratio(probe_persistent_defect(u, eps=0.5), u, degree=d)
        ),
    }


def main() -> None:
    _env1()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "mc_stress"))
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--mc-per-t", type=int, default=50000)
    ap.add_argument("--mc-batch", type=int, default=500)
    ap.add_argument("--n-points", type=int, default=8192)
    ap.add_argument("--degrees", type=str, default="2,4,6,8")
    ap.add_argument("--T-list", type=str, default="10,14,18,22")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument(
        "--min-mc-per-t",
        type=int,
        default=50000,
        help="refuse to complete if achieved min trials/T is below this",
    )
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if args.mc_per_t < int(args.min_mc_per_t):
        raise SystemExit(
            f"mc-per-t={args.mc_per_t} must be >= min-mc-per-t={args.min_mc_per_t}"
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    state_path = out_dir / "mc_stress_state.json"
    degrees = [int(x) for x in args.degrees.split(",") if x.strip()]
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)
    batch = max(50, int(args.mc_batch))

    if state_path.exists():
        state = json.loads(state_path.read_text())
    else:
        state = {
            "status": "running",
            "mc_per_t": args.mc_per_t,
            "degrees": degrees,
            "T_values": T_values,
            "n_points": args.n_points,
            "mc_by_T": {},
            "control_rows": [],
            "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
            "started_unix": time.time(),
        }

    t0 = time.time()
    print(
        f"MC stress mc_per_t={args.mc_per_t} T={T_values} d={degrees} workers={workers}",
        flush=True,
    )
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    # controls (cheap)
    if not state.get("control_rows"):
        ctrl_payloads = [
            {"T": T, "degree": d, "n_points": args.n_points}
            for T in T_values
            for d in degrees
        ]
        with ProcessPoolExecutor(max_workers=min(workers, len(ctrl_payloads) or 1), mp_context=_MP_CTX) as ex:
            state["control_rows"] = list(ex.map(_control_job, ctrl_payloads))
        state_path.write_text(json.dumps(state, indent=2))

    for T in T_values:
        key = str(T)
        done = int(state.get("mc_by_T", {}).get(key, {}).get("n_trials_done", 0))
        if done >= args.mc_per_t:
            print(f"T={T} MC already done ({done})", flush=True)
            continue
        remaining = args.mc_per_t - done
        # accumulators
        acc = state["mc_by_T"].get(key, {})
        sums = {str(d): float(acc.get("sum", {}).get(str(d), 0.0)) for d in degrees}
        sumsq = {str(d): float(acc.get("sumsq", {}).get(str(d), 0.0)) for d in degrees}
        n_done = done
        seed_base = int(acc.get("next_seed", 1000 + int(T * 100)))

        batches = []
        left = remaining
        seed = seed_base
        while left > 0:
            nt = min(batch, left)
            batches.append(
                {
                    "T": T,
                    "n_points": args.n_points,
                    "degrees": degrees,
                    "n_trials": nt,
                    "seed0": seed,
                }
            )
            seed += 1
            left -= nt

        print(f"T={T} launching {len(batches)} batches ({remaining} trials)", flush=True)
        with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CTX) as ex:
            for res in ex.map(_mc_batch_job, batches):
                for d in degrees:
                    ds = str(d)
                    st = res["per_degree"][ds]
                    # reconstruct sum from mean*n
                    sums[ds] += st["mean_R_d"] * st["n"]
                    # var = E[x^2]-mean^2 => E[x^2]=var+mean^2
                    ex2 = st["std_R_d"] ** 2 + st["mean_R_d"] ** 2
                    sumsq[ds] += ex2 * st["n"]
                n_done += res["n_trials"]
        per_degree = {}
        for d in degrees:
            ds = str(d)
            mean = sums[ds] / n_done
            var = max(0.0, sumsq[ds] / n_done - mean * mean)
            per_degree[ds] = {
                "mean_R_d": mean,
                "std_R_d": float(math.sqrt(var)),
                "n": n_done,
            }
        state["mc_by_T"][key] = {
            "T": T,
            "n_trials_done": n_done,
            "per_degree": per_degree,
            "sum": sums,
            "sumsq": sumsq,
            "next_seed": seed,
        }
        state["status"] = "running"
        state_path.write_text(json.dumps(state, indent=2))
        print(
            f"T={T} done n={n_done} mean_R4={per_degree.get('4', per_degree[str(degrees[0])])['mean_R_d']:.4e}",
            flush=True,
        )

    elapsed = time.time() - t0
    # verify min trials
    min_n = min(
        int(state["mc_by_T"][str(T)]["n_trials_done"]) for T in T_values
    )
    state["status"] = "completed"
    state["elapsed_s"] = elapsed
    state["min_mc_per_t"] = min_n
    state["mc_per_t_target"] = args.mc_per_t
    state["workers"] = workers
    state["completed_unix"] = time.time()
    state["interpretation"] = (
        f"MC stress ≥{args.mc_per_t}/T (achieved min {min_n}/T). "
        "Instrument stress only; not an RH proof."
    )
    state_path.write_text(json.dumps(state, indent=2))

    # summary txt
    lines = [
        "PBSS MC stress campaign",
        f"status=completed elapsed_s={elapsed:.1f} mc_per_t={args.mc_per_t} min_achieved={min_n}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
        "",
    ]
    for T in T_values:
        pd = state["mc_by_T"][str(T)]["per_degree"]
        d0 = str(degrees[len(degrees) // 2])
        lines.append(
            f"T={T:g} n={state['mc_by_T'][str(T)]['n_trials_done']} "
            f"mean_R_d[{d0}]={pd[d0]['mean_R_d']:.6e} std={pd[d0]['std_R_d']:.6e}"
        )
    text = "\n".join(lines) + "\n"
    (out_dir / "mc_stress_summary.txt").write_text(text)
    (out_dir / "mc_stress_summary.json").write_text(json.dumps(state, indent=2))
    print(text, flush=True)

    if not args.no_plot:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            d_focus = str(degrees[len(degrees) // 2])
            Ts = T_values
            means = [state["mc_by_T"][str(T)]["per_degree"][d_focus]["mean_R_d"] for T in Ts]
            stds = [state["mc_by_T"][str(T)]["per_degree"][d_focus]["std_R_d"] for T in Ts]
            fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
            ax.errorbar(Ts, means, yerr=stds, fmt="o-", capsize=3, label=f"MC d={d_focus}")
            # controls
            for d in degrees:
                crows = [c for c in state["control_rows"] if c["degree"] == d]
                if d == int(d_focus) and crows:
                    crows = sorted(crows, key=lambda z: z["T"])
                    ax.plot(
                        [c["T"] for c in crows],
                        [c["critical_line_R_d"] for c in crows],
                        "s--",
                        label="CL",
                    )
                    ax.plot(
                        [c["T"] for c in crows],
                        [c["persistent_defect_R_d"] for c in crows],
                        "^--",
                        label="defect",
                    )
            ax.set_xlabel("T")
            ax.set_ylabel(r"$R_d$")
            ax.set_title(f"MC stress ≥{args.mc_per_t}/T (not RH proof)")
            ax.grid(True, alpha=0.35)
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_dir / "mc_stress_Rd_vs_T.png", bbox_inches="tight")
            plt.close()
        except Exception as exc:
            print(f"plot skipped: {exc}", flush=True)

    if args.scratch:
        import shutil

        sc = Path(args.scratch)
        sc.mkdir(parents=True, exist_ok=True)
        for name in (
            "mc_stress_state.json",
            "mc_stress_summary.json",
            "mc_stress_summary.txt",
            "mc_stress_Rd_vs_T.png",
        ):
            p = out_dir / name
            if p.exists():
                shutil.copy(p, sc / name)

    min_req = int(args.min_mc_per_t)
    if min_n < min_req:
        raise SystemExit(f"min mc per T {min_n} < required {min_req}")
    print(f"MC_STRESS_COMPLETE elapsed_s={elapsed:.1f} min_n={min_n}", flush=True)
    (out_dir / "MC_STRESS_COMPLETE").write_text(
        f"elapsed_s={elapsed:.1f}\nmin_mc_per_t={min_n}\nrequired={min_req}\n"
    )


if __name__ == "__main__":
    main()
