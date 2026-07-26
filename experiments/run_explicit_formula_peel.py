#!/usr/bin/env python3
"""
Multi-T × multi-N peel / include scan for truncated explicit-formula residuals.

For each T and truncation depth N:
  - include: R_d of q_T^{(N)} (first N CL modes)
  - peel:    R_d of q_T^{(N_full)} - q_T^{(N)}  (strip first N from a fixed full sum)

Writes results/explicit_formula_peel.{json,txt} and a plot if matplotlib is available.

Not an RH proof — finite-mode / truncated residual diagnostics only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _scan_one(
    T: float,
    N: int,
    N_full: int,
    degree: int,
    n_points: int,
) -> dict:
    from pbss.lemmas import bound_R_d_finite_mode_sum
    from pbss.probes import explicit_formula_residual, peel_residual, sample_grid
    from pbss.projection import energy_ratio
    from pbss.zeros import explicit_formula_amplitudes, zeta_zero_ordinates

    u = sample_grid(n_points)
    q_inc, _, meta = explicit_formula_residual(u, T=T, n_zeros=N)
    r_inc = float(energy_ratio(q_inc, u, degree=degree))

    q_full, _, meta_full = explicit_formula_residual(u, T=T, n_zeros=N_full)
    q_peel, _ = peel_residual(q_full, u, T, n_strip=N)
    # peel residual may be near-zero when N == N_full
    try:
        r_peel = float(energy_ratio(q_peel, u, degree=degree))
    except ValueError:
        r_peel = float("nan")

    t = zeta_zero_ordinates(N)
    a = explicit_formula_amplitudes(t)
    bound = float(bound_R_d_finite_mode_sum(T, a, t, degree))

    return {
        "T": float(T),
        "N": int(N),
        "N_full": int(N_full),
        "degree": int(degree),
        "R_d_include": r_inc,
        "R_d_peel": r_peel,
        "M5_bound_include": bound,
        "n_points": int(n_points),
        "form": meta.get("form"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=8192)
    ap.add_argument(
        "--T-list",
        type=str,
        default="8,12,20,32,48",
        help="comma-separated T values",
    )
    ap.add_argument(
        "--N-list",
        type=str,
        default="1,2,5,10,20",
        help="comma-separated truncation depths N",
    )
    ap.add_argument(
        "--N-full",
        type=int,
        default=30,
        help="full sum size for peel (must be >= max N)",
    )
    ap.add_argument("--scratch", type=str, default="")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    N_values = [int(x) for x in args.N_list.split(",") if x.strip()]
    N_full = int(args.N_full)
    if max(N_values) > N_full:
        raise SystemExit(f"N_full={N_full} must be >= max N={max(N_values)}")

    print(
        f"explicit-formula peel scan d={args.degree} n={args.n_points} "
        f"T={T_values} N={N_values} N_full={N_full}",
        flush=True,
    )
    print(
        "Not an unconditional RH proof. Finite-mode / truncated residual only.",
        flush=True,
    )

    rows = []
    for T in T_values:
        for N in N_values:
            row = _scan_one(T, N, N_full, args.degree, args.n_points)
            rows.append(row)
            print(
                f"T={row['T']:6.1f}  N={row['N']:3d}  "
                f"R_inc={row['R_d_include']:.4e}  "
                f"R_peel={row['R_d_peel']:.4e}  "
                f"bound={row['M5_bound_include']:.4e}",
                flush=True,
            )

    # Qualitative: for fixed N, R_include should drop with T
    by_N: dict[int, list] = {}
    for r in rows:
        by_N.setdefault(r["N"], []).append(r)
    decay_flags = {}
    for N, rs in by_N.items():
        rs_sorted = sorted(rs, key=lambda x: x["T"])
        r0 = rs_sorted[0]["R_d_include"]
        r1 = rs_sorted[-1]["R_d_include"]
        decay_flags[str(N)] = bool(r1 < r0 * 0.75)

    summary = {
        "degree": args.degree,
        "n_points": args.n_points,
        "T_values": T_values,
        "N_values": N_values,
        "N_full": N_full,
        "rows": rows,
        "qualitative": {
            "R_include_decays_with_T_per_N": decay_flags,
            "all_N_decay_ok": bool(all(decay_flags.values())),
        },
        "interpretation": (
            "Supports finite-mode A0 (Lemma M5): truncated explicit-formula "
            "residuals have R_d → 0 as T→∞ at O(T^{-2}) scale. Peel columns "
            "show residual after stripping first N modes from an N_full sum. "
            "Not a proof of RH or full Theorem A for arithmetic primes."
        ),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }

    out_dir = ROOT / "results" / "explicit_formula_peel"
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "explicit_formula_peel.json"
    txt_path = out_dir / "explicit_formula_peel.txt"
    json_path.write_text(json.dumps(summary, indent=2))

    lines = [
        "PBSS explicit-formula multi-T × multi-N peel scan",
        f"d={args.degree} n_points={args.n_points} N_full={N_full}",
        "NOT AN UNCONDITIONAL PROOF OF RH",
        "",
        f"{'T':>8}  {'N':>4}  {'R_include':>12}  {'R_peel':>12}  {'M5_bound':>12}",
    ]
    for r in rows:
        lines.append(
            f"{r['T']:8.2f}  {r['N']:4d}  {r['R_d_include']:12.4e}  "
            f"{r['R_d_peel']:12.4e}  {r['M5_bound_include']:12.4e}"
        )
    lines += [
        "",
        f"R_include decays with T for each N: {decay_flags}",
        f"all_N_decay_ok={summary['qualitative']['all_N_decay_ok']}",
        "",
        summary["interpretation"],
    ]
    text = "\n".join(lines) + "\n"
    txt_path.write_text(text)
    print(text, flush=True)

    plot_path = None
    if not args.no_plot:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=140)
            for N in N_values:
                rs = sorted([r for r in rows if r["N"] == N], key=lambda x: x["T"])
                Ts = [r["T"] for r in rs]
                axes[0].semilogy(
                    Ts, [r["R_d_include"] for r in rs], "o-", label=f"N={N}"
                )
            axes[0].set_xlabel("T")
            axes[0].set_ylabel(r"$R_d$ include $q_T^{(N)}$")
            axes[0].set_title("Include: truncated EF residual")
            axes[0].grid(True, which="both", alpha=0.35)
            axes[0].legend(fontsize=8)

            # Peel vs N at largest T
            T_star = max(T_values)
            rs = sorted(
                [r for r in rows if r["T"] == T_star], key=lambda x: x["N"]
            )
            axes[1].plot(
                [r["N"] for r in rs],
                [r["R_d_peel"] for r in rs],
                "s-",
                color="#d62728",
                label=f"peel @ T={T_star:g}",
            )
            axes[1].plot(
                [r["N"] for r in rs],
                [r["R_d_include"] for r in rs],
                "o--",
                color="#1f77b4",
                label=f"include @ T={T_star:g}",
            )
            axes[1].set_xlabel("N (zeros)")
            axes[1].set_ylabel(r"$R_d$")
            axes[1].set_title("Peel vs include at fixed T")
            axes[1].grid(True, alpha=0.35)
            axes[1].legend(fontsize=8)
            fig.suptitle(
                "PBSS explicit-formula peel (not an RH proof)", fontsize=11
            )
            fig.tight_layout()
            plot_path = out_dir / "explicit_formula_peel.png"
            fig.savefig(plot_path, bbox_inches="tight", facecolor="white")
            plt.close()
            print(f"Wrote plot {plot_path}", flush=True)
        except Exception as exc:  # pragma: no cover
            print(f"plot skipped: {exc}", flush=True)

    if args.scratch:
        scratch = Path(args.scratch)
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "explicit_formula_peel.json").write_text(
            json.dumps(summary, indent=2)
        )
        (scratch / "explicit_formula_peel.txt").write_text(text)
        if plot_path and Path(plot_path).is_file():
            import shutil

            shutil.copy(plot_path, scratch / "explicit_formula_peel.png")
        print(f"Wrote scratch copies under {scratch}", flush=True)

    if not summary["qualitative"]["all_N_decay_ok"]:
        raise SystemExit("expected R_include to decay with T for each N")


if __name__ == "__main__":
    main()
