#!/usr/bin/env python3
"""Reproducible PBSS control experiment with an independent finer-grid check.

Measures low-degree energy, not the truth of RH. Four independent probes run
in a process pool; each is evaluated on n and 2n-1 samples of its defining
function (never by interpolating the coarse samples).
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.diagnostic import assess_demo
from pbss.probes import (
    probe_critical_line_mode, probe_defective, probe_high_frequency,
    probe_prime_residual, sample_grid,
)
from pbss.projection import project


PROBES = ("high_frequency", "critical_line", "prime_residual", "defective")
KEYS = ("rh_like_high_frequency", "rh_like_critical_line_mode",
        "rh_like_prime_residual", "defective_control")  # historical JSON keys


def _evaluate(job):
    name, degree, n_points, x_max = job
    results = []
    for n in (n_points, 2 * n_points - 1):
        u = sample_grid(n)
        T = 20.0
        if name == "high_frequency":
            q = probe_high_frequency(u, waves=48)
        elif name == "critical_line":
            q = probe_critical_line_mode(u, T=T)
        elif name == "defective":
            q = probe_defective(u, waves=48, defect_degree=1, defect_weight=2.5)
        else:
            q, T = probe_prime_residual(u, x_max=x_max, smooth=1)
        r = project(q, u, degree=degree, T=T)
        results.append(r)
    coarse, fine = results
    delta = abs(coarse.energy_ratio - fine.energy_ratio)
    tolerance = 1e-3 + 1e-2 * abs(fine.energy_ratio)
    row = {
        "name": name, "T": coarse.T, "energy": coarse.energy,
        "l2_norm_sq": coarse.l2_norm_sq, "energy_ratio": coarse.energy_ratio,
        "P_scaled_strength": coarse.P, "coeffs": coarse.coeffs.tolist(),
        "refinement": {
            "n_points": fine.n_points, "energy_ratio": fine.energy_ratio,
            "absolute_change": delta, "tolerance": tolerance,
            "stable": bool(delta <= tolerance),
        },
    }
    if name == "prime_residual":
        row["x_max"] = x_max
        row["note"] = "Prime Chebyshev residual with degree-1 detrending and no smoothing; staircase and slow components remain."
    return row


def run(degree: int = 4, n_points: int = 4096, x_max: float = 1e5,
        workers: int | None = None) -> dict:
    if isinstance(degree, bool) or not isinstance(degree, (int, np.integer)) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    if isinstance(n_points, bool) or not isinstance(n_points, (int, np.integer)) or n_points < 2:
        raise ValueError("n_points must be an integer >= 2")
    if degree >= n_points:
        raise ValueError("degree must be smaller than n_points")
    if not np.isfinite(x_max) or x_max < 3:
        raise ValueError("x_max must be finite and >= 3")
    if workers is None:
        workers = max(1, (os.cpu_count() or 1) - 2)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    # Four independent jobs; additional workers would have no work to perform.
    workers = min(workers, len(PROBES))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_evaluate, [(name, degree, n_points, x_max) for name in PROBES]))
    assessment = assess_demo(
        [{"name": row["name"], "Rd": row["energy_ratio"]} for row in rows],
        degree=degree, n_points=n_points, T=20.0,
    )
    fine_assessment = assess_demo(
        [{"name": row["name"], "Rd": row["refinement"]["energy_ratio"]} for row in rows],
        degree=degree, n_points=2 * n_points - 1, T=20.0,
    )
    stable = all(row["refinement"]["stable"] for row in rows)
    consistent = assessment["status"] == fine_assessment["status"]
    out = {
        "schema_version": 2,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "degree": int(degree), "n_points": int(n_points), "x_max": float(x_max),
        "compute": {"backend": "numpy", "workers": workers, "jobs": len(PROBES)},
        "definitions": {
            "energy_ratio": "R_d = ||P_d q||_W^2 / ||q||_W^2",
            "scaled_strength_P": "P := S_d = T^{2(d+1)} R_d (scaling convention only)",
            "basis": "shifted Legendre on [0,1]; weighted discrete least-squares fit",
            "coeffs": "Legendre expansion of the fitted polynomial; energy is c^T G c",
            "quadrature": "composite trapezoid",
            "note": "Generic model-mode decay is O(T^-2), not O(T^-2(d+1)); no universal RH threshold.",
        },
        "probe_parameters": {"waves": 48, "synthetic_T": 20.0,
                             "critical_ordinate": 14.134725,
                             "defect_degree": 1, "defect_weight": 2.5,
                             "prime_detrend": "deg1", "prime_smooth_samples": 1},
        "separation": assessment,
        "refined_separation": fine_assessment,
        "grid_refinement": {
            "stable": stable, "absolute_tolerance": 1e-3, "relative_tolerance": 1e-2,
            "note": "One refinement is a sensitivity check, not a certified continuum error bound.",
        },
        "status": assessment["status"] if stable and consistent else "inconclusive",
    }
    out.update(zip(KEYS, rows))
    return out


def render_text(result: dict) -> str:
    lines = ["Perry–Beurling Spectral Sieve — diagnostic run", result["banner"],
             f"degree d={result['degree']}, n_points={result['n_points']}",
             "Energy measures polynomial content; it does not classify RH.", ""]
    for key in KEYS:
        row = result[key]
        lines.append(f"{row['name']}: R_d={row['energy_ratio']:.6e}  S_d={row['P_scaled_strength']:.6e}")
        ref = row["refinement"]
        lines.append(f"  finer n={ref['n_points']}: R_d={ref['energy_ratio']:.6e}, change={ref['absolute_change']:.3e}, stable={ref['stable']}")
    separation = result["separation"]
    lines.extend(["", f"Synthetic control check: {separation['status']}",
                  f"Finer-grid control check: {result['refined_separation']['status']}",
                  f"Both separations must reach {separation['minimum_separation']:.1f}.",
                  f"Grid refinement stable: {result['grid_refinement']['stable']}",
                  f"Overall diagnostic status: {result['status']}"])
    lines.extend(separation["warnings"])
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--x-max", type=float, default=1e5)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--output-dir", type=Path, default=ROOT / "results")
    ap.add_argument("--scratch", type=Path, default=None, help="optional extra report directory")
    args = ap.parse_args()
    try:
        result = run(args.degree, args.n_points, args.x_max, args.workers)
        payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
        report = render_text(result)
        for directory in dict.fromkeys(p for p in (args.output_dir, args.scratch) if p is not None):
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "diagnostic_run.json").write_text(payload, encoding="utf-8")
            (directory / "diagnostic_run.txt").write_text(report, encoding="utf-8")
    except (ValueError, OSError) as exc:
        ap.error(str(exc))
    print(report, end="")
    return 0 if result["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
