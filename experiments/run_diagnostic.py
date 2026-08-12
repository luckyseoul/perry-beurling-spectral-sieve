#!/usr/bin/env python3
"""
End-to-end PBSS diagnostic experiment.

Compares:
  (A) RH-consistent probes: high-frequency synthetic + prime residual
  (B) Defective control: low-degree contaminated oscillation

Writes JSON + text under results/ and optional extra output path.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.probes import (  # noqa: E402
    probe_critical_line_mode,
    probe_defective,
    probe_high_frequency,
    probe_prime_residual,
    sample_grid,
)
from pbss.projection import project  # noqa: E402


def run(degree: int = 4, n_points: int = 4096, x_max: float = 1e5) -> dict:
    u = sample_grid(n_points)
    T_syn = 20.0  # synthetic log-window scale for S_d

    # (A1) RH-like synthetic: pure high frequency
    q_hf = probe_high_frequency(u, waves=48)
    r_hf = project(q_hf, u, degree=degree, T=T_syn)

    # (A2) RH-like: critical-line oscillatory mode (first zeta zero height)
    q_cl = probe_critical_line_mode(u, T=T_syn)
    r_cl = project(q_cl, u, degree=degree, T=T_syn)

    # (A3) prime residual (demeaned θ(x)-x); modest x_max still carries
    # staircase / slow components — reported honestly
    q_pr, T_pr = probe_prime_residual(u, x_max=x_max)
    r_pr = project(q_pr, u, degree=degree, T=T_pr)

    # (B) defective control
    q_def = probe_defective(u, waves=48, defect_degree=1, defect_weight=2.5)
    r_def = project(q_def, u, degree=degree, T=T_syn)

    out = {
        "degree": degree,
        "n_points": n_points,
        "definitions": {
            "energy_ratio": "R_d = ||P_d q||^2 / ||q||^2",
            "scaled_strength_P": "P := S_d = T^{2(d+1)} R_d  (working projection strength)",
            "basis": "orthonormal shifted Legendre on [0,1]",
            "note": (
                "Legacy archive quoted P≈3.92 for zeta with lost normalization; "
                "this reconstruction reports honest R_d and S_d from the shipped code."
            ),
        },
        "rh_like_high_frequency": {
            "label": "synthetic high-frequency (RH-like)",
            "T": T_syn,
            "energy": r_hf.energy,
            "energy_ratio": r_hf.energy_ratio,
            "P_scaled_strength": r_hf.P,
        },
        "rh_like_critical_line_mode": {
            "label": "sin(t T u) with t=Im(first zeta zero) (RH-consistent form)",
            "T": T_syn,
            "energy": r_cl.energy,
            "energy_ratio": r_cl.energy_ratio,
            "P_scaled_strength": r_cl.P,
        },
        "rh_like_prime_residual": {
            "label": f"demeaned prime Chebyshev residual x_max={x_max:g}",
            "T": T_pr,
            "x_max": x_max,
            "energy": r_pr.energy,
            "energy_ratio": r_pr.energy_ratio,
            "P_scaled_strength": r_pr.P,
            "note": (
                "At modest x_max the residual still has substantial low-degree "
                "mass; not expected to match pure high-frequency R_d."
            ),
        },
        "defective_control": {
            "label": "high-frequency + strong degree-1 defect (non-RH control)",
            "T": T_syn,
            "energy": r_def.energy,
            "energy_ratio": r_def.energy_ratio,
            "P_scaled_strength": r_def.P,
        },
        "separation": {
            "energy_ratio_defective_minus_hf": r_def.energy_ratio - r_hf.energy_ratio,
            "energy_ratio_defective_minus_critical": r_def.energy_ratio
            - r_cl.energy_ratio,
            "classifier_ok": bool(
                r_def.energy_ratio > r_hf.energy_ratio
                and r_def.energy_ratio > r_cl.energy_ratio
            ),
        },
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--x-max", type=float, default=1e5)
    ap.add_argument(
        "--scratch",
        type=str,
        default="",
        help="optional directory for extra copies of results",
    )
    args = ap.parse_args()

    result = run(degree=args.degree, n_points=args.n_points, x_max=args.x_max)

    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / "diagnostic_run.json"
    txt_path = results_dir / "diagnostic_run.txt"
    json_path.write_text(json.dumps(result, indent=2))

    lines = [
        "Perry–Beurling Spectral Sieve — diagnostic run",
        f"degree d={result['degree']}, n_points={result['n_points']}",
        "",
        "RH-like high-frequency:",
        f"  R_d={result['rh_like_high_frequency']['energy_ratio']:.6e}",
        f"  P=S_d={result['rh_like_high_frequency']['P_scaled_strength']:.6e}",
        "",
        "RH-like critical-line mode:",
        f"  R_d={result['rh_like_critical_line_mode']['energy_ratio']:.6e}",
        f"  P=S_d={result['rh_like_critical_line_mode']['P_scaled_strength']:.6e}",
        "",
        "Prime residual (demeaned, modest x_max — see note in JSON):",
        f"  R_d={result['rh_like_prime_residual']['energy_ratio']:.6e}",
        f"  P=S_d={result['rh_like_prime_residual']['P_scaled_strength']:.6e}",
        f"  T=log(x_max)={result['rh_like_prime_residual']['T']:.4f}",
        "",
        "Defective control:",
        f"  R_d={result['defective_control']['energy_ratio']:.6e}",
        f"  P=S_d={result['defective_control']['P_scaled_strength']:.6e}",
        "",
        f"Separation (R_def - R_hf)="
        f"{result['separation']['energy_ratio_defective_minus_hf']:.6e}",
        f"Separation (R_def - R_critical)="
        f"{result['separation']['energy_ratio_defective_minus_critical']:.6e}",
        f"classifier_ok={result['separation']['classifier_ok']}",
        "",
        "Not a proof of RH. Conditional diagnostic only.",
    ]
    text = "\n".join(lines) + "\n"
    txt_path.write_text(text)
    print(text)

    if args.scratch:
        scratch = Path(args.scratch)
        scratch.mkdir(parents=True, exist_ok=True)
        (scratch / "diagnostic_run.json").write_text(json.dumps(result, indent=2))
        (scratch / "diagnostic_run.txt").write_text(text)
        print(f"Also wrote to {scratch}", flush=True)


if __name__ == "__main__":
    main()
