"""
pbss — command-line tool for spectral residual diagnostics.

General mathematics usage (not RH-only):
  - project an arbitrary 1D residual q(u) on [0,1]
  - scorecard RH-like vs defective synthetic controls
  - Gamma-weight sensitivity (offline vs online discriminability)
  - Beurling-style ordinary vs defective separation (optional)

Examples:
  pbss project --input q.npy --degree 4 --T 20
  pbss diagnose --demo
  pbss sensitivity --confirm-53
  pbss scorecard --x-max 1e6
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"


def _load_array(path: str) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    if p.suffix.lower() == ".npy":
        return np.load(p)
    if p.suffix.lower() in (".csv", ".txt"):
        return np.loadtxt(p, delimiter="," if p.suffix.lower() == ".csv" else None)
    # fallback
    return np.loadtxt(p)


def _unit_grid(n: int) -> np.ndarray:
    from .probes import sample_grid

    return sample_grid(int(n))


def cmd_project(args: argparse.Namespace) -> int:
    from .projection import project

    q = np.asarray(_load_array(args.input), dtype=np.float64).ravel()
    if args.u:
        u = np.asarray(_load_array(args.u), dtype=np.float64).ravel()
    else:
        u = _unit_grid(q.size)
    if u.size != q.size:
        raise SystemExit("u and q must have the same length")
    r = project(q, u, degree=int(args.degree), T=float(args.T))
    out = {
        "banner": BANNER,
        "degree": r.degree,
        "T": r.T,
        "n_points": r.n_points,
        "energy": r.energy,
        "energy_ratio_Rd": r.energy_ratio,
        "scaled_strength_Sd": r.scaled_strength,
        "coeffs": r.coeffs.tolist(),
    }
    text = (
        f"PBSS project\n"
        f"  R_d = {r.energy_ratio:.6e}\n"
        f"  S_d = {r.scaled_strength:.6e}\n"
        f"  d={r.degree}  T={r.T}  n={r.n_points}\n"
        f"{BANNER}\n"
    )
    print(text, end="")
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(out, indent=2))
        print(f"Wrote {args.json_out}", file=sys.stderr)
    return 0


def cmd_diagnose(args: argparse.Namespace) -> int:
    from .probes import (
        probe_critical_line_mode,
        probe_defective,
        probe_high_frequency,
        probe_prime_residual,
        sample_grid,
    )
    from .projection import project

    u = sample_grid(int(args.n_points))
    degree = int(args.degree)
    T = float(args.T)
    rows = []
    if args.demo or not args.input:
        specs = [
            ("high_frequency", probe_high_frequency(u, waves=48), T),
            ("critical_line", probe_critical_line_mode(u, T=T), T),
            ("defective", probe_defective(u, waves=48, defect_degree=1, defect_weight=2.5), T),
        ]
        if args.with_primes:
            q_pr, T_pr = probe_prime_residual(u, x_max=float(args.x_max))
            specs.append(("prime_residual", q_pr, T_pr))
    else:
        q = np.asarray(_load_array(args.input), dtype=np.float64).ravel()
        specs = [("input", q, T)]

    print("PBSS diagnose")
    print(f"degree={degree}  n_points={u.size}")
    for name, q, Tv in specs:
        r = project(q, u, degree=degree, T=float(Tv))
        rows.append(
            {
                "name": name,
                "T": float(Tv),
                "Rd": r.energy_ratio,
                "Sd": r.scaled_strength,
            }
        )
        print(f"  {name:16s}  R_d={r.energy_ratio:.6e}  S_d={r.scaled_strength:.6e}  T={Tv:.4g}")
    # separation if demo
    names = {r["name"]: r["Rd"] for r in rows}
    if "defective" in names and "high_frequency" in names:
        sep = names["defective"] - names["high_frequency"]
        print(f"  separation (defective - HF) = {sep:.6e}")
        print(f"  classifier_ok = {sep > 0.5}")
    print(BANNER)
    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps({"banner": BANNER, "rows": rows}, indent=2)
        )
    return 0


def cmd_sensitivity(args: argparse.Namespace) -> int:
    from .measure_sensitivity import confirm_sensitivity_claim, sensitivity_experiment

    if args.confirm_53:
        rep = confirm_sensitivity_claim(
            min_gain=0.53,
            seed=int(args.seed),
            noisy_noise=float(args.noise),
        )
        print("PBSS Gamma-weight sensitivity confirmation")
        print(f"  verdict: {rep['verdict']}")
        print(
            f"  noisy:  flat d′={rep['noisy']['flat_dprime']:.4f}  "
            f"gamma d′={rep['noisy']['gamma_dprime']:.4f}  "
            f"gain={rep['noisy']['relative_gain_percent']:.1f}%  "
            f"meets≥53%={rep['noisy_meets_53pct']}"
        )
        print(
            f"  clean:  flat d′={rep['clean']['flat_dprime']:.4f}  "
            f"gamma d′={rep['clean']['gamma_dprime']:.4f}  "
            f"gain={rep['clean']['relative_gain_percent']:.1f}%"
        )
        print(f"  confirmed: {rep['confirmed']}")
        print(BANNER)
        if args.json_out:
            Path(args.json_out).write_text(json.dumps(rep, indent=2, default=str))
        return 0 if rep["confirmed"] else 2

    rep = sensitivity_experiment(
        n_per_class=int(args.n_per_class),
        noise=float(args.noise),
        seed=int(args.seed),
        k=float(args.k),
        sigma=float(args.sigma),
    )
    print("PBSS Gamma-weight sensitivity")
    print(f"  flat d′  = {rep['flat_dprime']:.6f}")
    print(f"  gamma d′ = {rep['gamma_dprime']:.6f}")
    print(f"  relative gain = {rep['relative_gain_percent']:.2f}%")
    print(f"  peak u⋆ = {rep['peak_u_star']:.4f}  (k={rep['k']}, σ={rep['sigma']})")
    print(BANNER)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(rep, indent=2, default=str))
    return 0


def cmd_scorecard(args: argparse.Namespace) -> int:
    from .beurling import beurling_theta_residual, build_system_primes
    from .probes import sample_grid
    from .projection import project

    u = sample_grid(int(args.n_points))
    x_max = float(args.x_max)
    degree = int(args.degree)
    systems = [
        ("ordinary_primes", "ordinary_primes"),
        ("gapped_gap3", "gapped_gap3"),
        ("thinned_every3", "thinned_every3"),
    ]
    print(f"PBSS Beurling scorecard  x_max={x_max:g}  d={degree}")
    rows = []
    for label, name in systems:
        try:
            primes = build_system_primes(name, x_max=x_max)
            q, Tv = beurling_theta_residual(u, primes, x_max=x_max)
            r = project(q, u, degree=degree, T=Tv)
            rows.append({"system": label, "Rd": r.energy_ratio, "T": Tv})
            print(f"  {label:20s}  R_d={r.energy_ratio:.6f}  T={Tv:.4g}")
        except Exception as e:
            print(f"  {label:20s}  ERROR: {e}")
    print(BANNER)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps({"rows": rows, "banner": BANNER}, indent=2))
    return 0


def cmd_version(_: argparse.Namespace) -> int:
    from . import __version__

    print(f"pbss {__version__}")
    print(BANNER)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pbss",
        description=(
            "Perry–Beurling Spectral Sieve — spectral residual diagnostic tool "
            "for general mathematics use (not a proof of RH)."
        ),
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("version", help="print version")
    v.set_defaults(func=cmd_version)

    pr = sub.add_parser("project", help="project q(u) onto shifted Legendre; report R_d, S_d")
    pr.add_argument("--input", "-i", required=True, help="q array (.npy/.csv/.txt)")
    pr.add_argument("--u", default="", help="optional u grid (default linspace [0,1])")
    pr.add_argument("--degree", "-d", type=int, default=4)
    pr.add_argument("--T", type=float, default=20.0, help="log-window length for S_d")
    pr.add_argument("--json-out", default="", help="write JSON report")
    pr.set_defaults(func=cmd_project)

    di = sub.add_parser("diagnose", help="demo scorecard or diagnose one residual")
    di.add_argument("--demo", action="store_true", help="built-in HF / CL / defective probes")
    di.add_argument("--input", "-i", default="", help="optional q array")
    di.add_argument("--degree", "-d", type=int, default=4)
    di.add_argument("--T", type=float, default=20.0)
    di.add_argument("--n-points", type=int, default=2048)
    di.add_argument("--with-primes", action="store_true")
    di.add_argument("--x-max", type=float, default=1e5)
    di.add_argument("--json-out", default="")
    di.set_defaults(func=cmd_diagnose)

    se = sub.add_parser(
        "sensitivity",
        help="Gamma vs flat weight discriminability (offline vs online)",
    )
    se.add_argument("--confirm-53", action="store_true", help="confirm ≥53%% relative gain claim")
    se.add_argument("--noise", type=float, default=0.55)
    se.add_argument("--seed", type=int, default=20260522)
    se.add_argument("--n-per-class", type=int, default=200)
    se.add_argument("--k", type=float, default=4.0)
    se.add_argument("--sigma", type=float, default=6.0)
    se.add_argument("--json-out", default="")
    se.set_defaults(func=cmd_sensitivity)

    sc = sub.add_parser("scorecard", help="ordinary vs gapped/thinned Beurling R_d")
    sc.add_argument("--x-max", type=float, default=1e6)
    sc.add_argument("--degree", "-d", type=int, default=4)
    sc.add_argument("--n-points", type=int, default=2048)
    sc.add_argument("--json-out", default="")
    sc.set_defaults(func=cmd_scorecard)

    return p


def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
