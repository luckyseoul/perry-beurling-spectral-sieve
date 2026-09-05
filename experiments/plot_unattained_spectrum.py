#!/usr/bin/env python3
"""Plot the explicit analytic spectrum and its proved normalized rate factor.

These are formula illustrations, not measured R_d values or zeta zeros.
Only small NumPy arrays are evaluated; no numerical campaign is involved.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 13, "axes.labelsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#8892a2", "text.color": "#18263b",
        "axes.labelcolor": "#18263b", "xtick.color": "#46556a",
        "ytick.color": "#46556a", "savefig.facecolor": "white",
        "svg.hashsalt": "pbss-unattained-spectrum-20260905",
    })
    n = np.arange(8, 81)
    real_parts = 0.25 - 1 / n
    windows = np.geomspace(64, 10000, 300)
    log10_rate = (np.log(windows) - 4 * (np.sqrt(2) - 1) * np.sqrt(windows)) / np.log(10)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(left=.08, right=.97, bottom=.24, top=.74, wspace=.31)
    fig.suptitle("An unattained boundary can evade polynomial-rate tests",
                 x=.06, y=.96, ha="left", fontsize=17, fontweight="bold")
    fig.text(.06, .885, "Explicit infinite-series theorem · arithmetic residue coefficients $-1/\\rho_n$",
             fontsize=11, color="#526277")

    ax = axes[0]
    ax.scatter(real_parts, n, s=18, color="#156f91", alpha=.88, zorder=3)
    ax.axvline(.25, color="#c77135", linestyle="--", linewidth=1.6)
    ax.text(.253, 45, "unattained boundary 1/4", ha="left", va="center",
            rotation=90, color="#a75222", fontsize=10)
    ax.set(xlim=(.115, .261), ylim=(4, 85), xlabel="Real part of Laplace pole, $a_n$",
           ylabel=r"Logarithmic ordinate, $\log \gamma_n=n$",
           title=r"Genuine poles: $a_n=1/4-1/n,\ \gamma_n=e^n$")
    ax.grid(alpha=.16)
    ax.text(.03, .89, "First 73 positive-ordinate poles shown", transform=ax.transAxes,
            fontsize=9, color="#526277")

    ax = axes[1]
    ax.plot(windows, log10_rate, color="#6c55a2", linewidth=2.5)
    ax.set_xscale("log")
    ax.set(xlim=(64, 10000), xlabel="Log-window length, $T$",
           ylabel=r"$\log_{10} f(T)$", title="Proved rate factor for every fixed degree")
    ax.grid(alpha=.16)
    ax.text(.04, .12, r"$R_d=O_d(f(T))$" + "\n" + r"$f(T)=T\,e^{-4(\sqrt{2}-1)\sqrt{T}}$",
            transform=ax.transAxes, fontsize=13, color="#554184")
    fig.text(.06, .105, "Right: formula illustration; the degree-dependent constant is omitted. This is not a measured ratio.",
             fontsize=10, color="#526277")
    fig.text(.06, .05, "Analytic spectrum example. No prime realization, Euler product, or functional equation is asserted.",
             fontsize=10, color="#526277")

    figure_dir = ROOT / "docs/figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    outputs = [figure_dir / f"unattained_spectrum.{ext}" for ext in ("png", "svg")]
    fig.savefig(outputs[0], dpi=180)
    fig.savefig(outputs[1], metadata={"Date": None})
    plt.close(fig)
    sources = ["experiments/plot_unattained_spectrum.py", "docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md"]
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "kind": "analytic_formula_illustration_not_measured_ratio",
        "proof": "docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md",
        "prime_realization_asserted": False, "RH_proved": False,
        "formulas": {"a_n": "1/4-1/n", "gamma_n": "exp(n)",
                     "residue": "-1/(1/2+a_n+i*gamma_n)",
                     "log_f": "log(T)-4*(sqrt(2)-1)*sqrt(T)",
                     "meaning": "R_d=O_d(f(T)) for every fixed d>=2"},
        "n": n.tolist(), "pole_real_parts": real_parts.tolist(),
        "T": windows.tolist(), "log10_rate_factor": log10_rate.tolist(),
        "source_sha256": {name: sha(ROOT / name) for name in sources},
        "figure_sha256": {str(path.relative_to(ROOT)): sha(path) for path in outputs},
        "versions": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
    }
    destination = ROOT / "results/readme_figures_20260905/unattained_spectrum.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"figures": list(report["figure_sha256"]), "data": str(destination.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
