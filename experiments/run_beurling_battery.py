#!/usr/bin/env python3
"""
Beurling / generalized-prime PBSS battery: multi-T R_d scorecard.

Systems: ordinary primes (RH-like ambient), gapped, thinned (defective).
Not an RH proof.
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prime-dir", type=str, default=str(ROOT / "results" / "prime_checkpoints"))
    ap.add_argument("--x-max", type=float, default=1e8, help="window max for battery (uses primes ≤ this)")
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "beurling_battery"))
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--n-points", type=int, default=8192)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--T-list", type=str, default="8,10,12,14,16,18")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from pbss.beurling import (
        beurling_theta_residual,
        build_system_primes,
        default_battery_specs,
    )
    from pbss.primes_io import ensure_primes, find_largest_checkpoint, load_primes_checkpoint
    from pbss.probes import sample_grid
    from pbss.projection import energy_ratio

    t0 = time.time()
    print(f"beurling battery x_max={args.x_max:.3e} d={args.degree}", flush=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    # load ordinary primes: prefer existing largest checkpoint that covers x_max
    x_need = int(args.x_max)
    found = find_largest_checkpoint(args.prime_dir)
    if found and found[0] >= x_need:
        ordinary, pmeta = load_primes_checkpoint(args.prime_dir, found[0], mmap=True)
        # use prefix only
        ordinary = ordinary[: int(np.searchsorted(ordinary, x_need, side="right"))]
    else:
        ordinary, pmeta = ensure_primes(x_need, args.prime_dir)

    ordinary = np.asarray(ordinary)
    u = sample_grid(args.n_points)
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    T_cap = float(np.log(args.x_max))
    T_values = [T for T in T_values if T <= T_cap + 1e-9]
    specs = default_battery_specs()

    rows = []
    systems_meta = []
    for spec in specs:
        p_sys = build_system_primes(spec, ordinary, args.x_max)
        systems_meta.append(
            {
                "name": spec["name"],
                "kind": spec["kind"],
                "n_primes": int(np.asarray(p_sys).size),
                "p_max": float(np.asarray(p_sys)[-1]) if len(p_sys) else 0.0,
                "builder": spec.get("builder"),
            }
        )
        for T in T_values:
            q, T_out, meta = beurling_theta_residual(
                u, p_sys, T=T, detrend="deg1"
            )
            r = float(energy_ratio(q, u, degree=args.degree))
            row = {
                "system": spec["name"],
                "kind": spec["kind"],
                "T": float(T_out),
                "degree": int(args.degree),
                "R_d": r,
                "n_primes": meta["n_primes"],
            }
            rows.append(row)
            print(
                f"{spec['name']:16s} kind={spec['kind']:10s} "
                f"T={T_out:6.2f} R_d={r:.4e}",
                flush=True,
            )

    elapsed = time.time() - t0
    # separation: at largest T, mean R_d defective vs ordinary
    T_star = max(T_values) if T_values else 0.0
    by_sys = {}
    for r in rows:
        if r["T"] == T_star:
            by_sys[r["system"]] = r["R_d"]
    ordinary_R = by_sys.get("ordinary_primes")
    def_Rs = [by_sys[k] for k in by_sys if k != "ordinary_primes"]
    sep_ok = (
        ordinary_R is not None
        and def_Rs
        and all(d > ordinary_R for d in def_Rs)
    )

    summary = {
        "status": "completed",
        "elapsed_s": elapsed,
        "x_max": float(args.x_max),
        "degree": args.degree,
        "T_values": T_values,
        "systems": systems_meta,
        "rows": rows,
        "scorecard_at_Tmax": by_sys,
        "qualitative": {
            "T_star": T_star,
            "defective_above_ordinary": bool(sep_ok),
            "ordinary_R_d": ordinary_R,
            "defective_R_d": def_Rs,
        },
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "interpretation": (
            "Beurling battery: ordinary vs gapped/thinned systems. "
            "Defective systems should keep higher R_d. Not an RH proof."
        ),
        "n_systems": len(specs),
    }
    (out_dir / "beurling_battery.json").write_text(json.dumps(summary, indent=2))
    lines = [
        "PBSS Beurling battery scorecard",
        f"status=completed elapsed_s={elapsed:.1f} systems={len(specs)}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
        "",
        f"{'system':16s} {'kind':10s} {'T':>8} {'R_d':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['system']:16s} {r['kind']:10s} {r['T']:8.2f} {r['R_d']:12.4e}"
        )
    lines += [
        "",
        f"At T={T_star}: scorecard={by_sys}",
        f"defective_above_ordinary={sep_ok}",
        summary["interpretation"],
    ]
    text = "\n".join(lines) + "\n"
    (out_dir / "beurling_battery.txt").write_text(text)
    print(text, flush=True)

    if not args.no_plot and rows:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(9, 5), dpi=140)
            for name in sorted(set(r["system"] for r in rows)):
                rs = sorted([r for r in rows if r["system"] == name], key=lambda z: z["T"])
                ax.plot(
                    [r["T"] for r in rs],
                    [r["R_d"] for r in rs],
                    "o-",
                    label=name,
                )
            ax.set_xlabel("T")
            ax.set_ylabel(r"$R_d$")
            ax.set_title("Beurling battery multi-T (not RH proof)")
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.35)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(out_dir / "beurling_battery_Rd_vs_T.png", bbox_inches="tight")
            plt.close()
        except Exception as exc:
            print(f"plot skipped: {exc}", flush=True)

    if args.scratch:
        import shutil

        sc = Path(args.scratch)
        sc.mkdir(parents=True, exist_ok=True)
        for name in (
            "beurling_battery.json",
            "beurling_battery.txt",
            "beurling_battery_Rd_vs_T.png",
        ):
            p = out_dir / name
            if p.exists():
                shutil.copy(p, sc / name)

    if len(specs) < 2:
        raise SystemExit("need ≥2 systems")
    print(f"BEURLING_BATTERY_COMPLETE elapsed_s={elapsed:.1f}", flush=True)


if __name__ == "__main__":
    main()
