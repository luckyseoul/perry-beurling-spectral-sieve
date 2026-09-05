#!/usr/bin/env python3
"""Reproduce the README's finite ordinary-prime moment figure.

Evaluates 36 equally spaced log-window lengths from log(2) to log(100000)
using analytic prime-interval integrals and continuous affine detrending.
The accompanying JSON retains the complete decimal-valued API reports.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

for _thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
WINDOW_COUNT = 36
DEGREE = 4
PRIME_LIMIT = 100_000


def _window(job):
    from pbss.arithmetic_moments import arithmetic_moments

    index, T = job
    report = arithmetic_moments(T, degree=DEGREE, max_prime_limit=PRIME_LIMIT).to_dict()
    ratio = report["affine_detrended_energy_ratio"]
    if ratio is None or not 0 < ratio <= 1:
        raise ValueError(f"window {index}: affine ratio cannot be plotted on a log scale")
    return {"window_index": index, **report}


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _relative(path):
    return os.path.relpath(path, ROOT)


def _plot(rows, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    T = np.array([float(row["T_decimal"]) for row in rows])
    c2 = np.array([float(row["coefficients_decimal"][2]) for row in rows])
    ratio = np.array([row["affine_detrended_energy_ratio"] for row in rows])
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 13, "axes.labelsize": 12,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#6b7280", "axes.labelcolor": "#172033",
        "xtick.color": "#344054", "ytick.color": "#344054",
        "text.color": "#172033", "svg.hashsalt": "pbss-arithmetic-moments-v1",
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4))
    fig.subplots_adjust(left=0.075, right=0.975, bottom=0.20, top=0.74, wspace=0.30)
    fig.suptitle("Continuous moments of the ordinary-prime residual", x=0.075,
                 y=0.975, ha="left", fontsize=18, fontweight="bold")
    fig.text(0.075, 0.90,
             r"$Q(y)=e^{-y/2}[\theta(e^y)-e^y]$   ·   36 finite windows, $2\leq e^T\leq10^5$",
             fontsize=12, color="#475467")

    axes[0].axhline(0, color="#667085", linewidth=1.0, zorder=1)
    axes[0].plot(T, c2, color="#087e8b", linewidth=1.6, marker="o", markersize=4.0,
                 markeredgewidth=0, zorder=2)
    axes[0].set_title("(a) Signed degree-two moment", loc="left", pad=13)
    axes[0].set_ylabel(r"Signed absolute moment $c_2(T)$")

    axes[1].plot(T, ratio, color="#7652a3", linewidth=1.6, marker="o", markersize=4.0,
                 markeredgewidth=0)
    axes[1].set_yscale("log")
    axes[1].set_title("(b) Continuous affine projection ratio", loc="left", pad=13)
    axes[1].set_ylabel(r"Projection energy fraction $R_4(q_T)$")
    axes[1].set_ylim(min(ratio) / 1.6, 1.25)
    for ax in axes:
        ax.set_xlabel(r"Log-window length $T$")
        ax.set_xlim(T[0] - 0.15, T[-1] + 0.15)
        ax.grid(axis="y", color="#e4e7ec", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(length=4)

    fig.text(0.075, 0.075,
             r"Analytic prime-interval integrals · $q_T=(I-P_1)Q(T\,\cdot)$ · continuous $L^2([0,1],du)$",
             fontsize=10, color="#475467")
    fig.text(0.075, 0.032,
             "Each marker is one evaluated finite window; line segments connect adjacent markers.",
             fontsize=10, color="#667085")
    stem.parent.mkdir(parents=True, exist_ok=True)
    png, svg = stem.with_suffix(".png"), stem.with_suffix(".svg")
    fig.savefig(png, dpi=220, metadata={
        "Title": "Continuous finite ordinary-prime arithmetic moments",
        "Description": "36 finite log-windows; exact prime-interval identities in longdouble.",
    })
    fig.savefig(svg, metadata={"Date": None, "Title": "Continuous finite ordinary-prime arithmetic moments"})
    plt.close(fig)
    return png, svg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 2),
                        help="CPU worker budget; capped at the 36 independent windows")
    parser.add_argument("--figure-stem", type=Path,
                        default=ROOT / "docs/figures/continuous_arithmetic_moments")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "results/readme_figures_20260905/arithmetic_moments.json")
    parser.add_argument("--reuse-data", type=Path,
                        help="Re-render retained rows after checking all windows and evaluator source hashes")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")

    import matplotlib
    import numpy as np
    import scipy

    windows = np.linspace(np.log(np.longdouble(2)), np.log(np.longdouble(PRIME_LIMIT)),
                          WINDOW_COUNT, dtype=np.longdouble)
    jobs = [(index, np.format_float_scientific(T, unique=True, trim="-"))
            for index, T in enumerate(windows)]
    workers = min(args.workers, len(jobs))
    sources = ["experiments/plot_arithmetic_moments.py", "src/pbss/arithmetic_moments.py",
               "src/pbss/basis.py", "src/pbss/probes.py"]
    previous = None
    if args.reuse_data:
        previous = json.loads(args.reuse_data.read_text())
        for name in sources[1:]:
            if previous["source_sha256"][name] != _sha256(ROOT / name):
                raise ValueError(f"retained rows use a different evaluator dependency: {name}")
        rows = previous["rows"]
        if (previous["window_count"], previous["degree"], previous["max_prime_limit"]) != (WINDOW_COUNT, DEGREE, PRIME_LIMIT):
            raise ValueError("retained campaign parameters do not match")
        if [(row["window_index"], row["T_decimal"]) for row in rows] != jobs:
            raise ValueError("retained windows do not match the exact campaign grid")
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_window, jobs))
    figures = _plot(rows, args.figure_stem)
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              text=True, capture_output=True, check=False)
    report = {
        "status": "finite_ordinary_prime_windows_evaluated",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(), "python": sys.version,
        "source_git_head": revision.stdout.strip() if revision.returncode == 0 else None,
        "source_sha256": {name: _sha256(ROOT / name) for name in sources},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "backend": "ProcessPool_analytic_prime_intervals_longdouble",
        "requested_worker_budget": args.workers, "workers": workers,
        "worker_selection": "minimum_of_worker_budget_and_36_independent_windows",
        "gpu_used": False,
        "cpu_reason": "The frozen evaluator uses extended-precision longdouble analytic identities; each window is independent.",
        "window_count": WINDOW_COUNT, "degree": DEGREE,
        "window_spacing": "uniform_in_T_using_numpy_longdouble_linspace",
        "endpoint_T_decimals": [jobs[0][1], jobs[-1][1]],
        "max_prime_limit": PRIME_LIMIT,
        "residual": "Q(y)=exp(-y/2)*(theta(exp(y))-exp(y))",
        "integration": "analytic_prime_interval_identities_without_sampled_residual_interpolation",
        "affine_detrending": "continuous_orthogonal_projection_onto_span(1,u)_removed_exactly",
        "figure_interpretation": "Markers are finite-window evaluations; connecting segments only guide the eye.",
        "certified_interval_result": False,
        "figures": [{"path": _relative(path), "sha256": _sha256(path)} for path in figures],
        "rows": rows,
    }
    if previous is not None:
        # Preserve the original numerical evaluation provenance when only the
        # rendering changes. The new script and figure hashes describe this run.
        for key in ("timestamp_utc", "host", "python", "source_git_head", "versions",
                    "requested_worker_budget", "workers", "worker_selection"):
            report[key] = previous[key]
    report["figure_render"] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(), "matplotlib_version": matplotlib.__version__,
        "reused_window_reports": previous is not None,
        "window_reports_sha256": hashlib.sha256(
            json.dumps(rows, sort_keys=True, allow_nan=False).encode()).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"artifact": _relative(args.out), "sha256": _sha256(args.out),
                      "windows": len(rows), "workers": workers,
                      "figures": report["figures"]}, indent=2))


if __name__ == "__main__":
    main()
