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
  pbss moments --T 10 --degree 4
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

BANNER = "NOT AN UNCONDITIONAL PROOF OF RH"


def _load_array(path: str) -> np.ndarray:
    """Load one real sample vector without silently flattening or discarding data."""
    p = Path(path).expanduser()
    try:
        if p.suffix.lower() == ".npy":
            values = np.load(p, allow_pickle=False)
        else:
            values = np.loadtxt(
                p, delimiter="," if p.suffix.lower() == ".csv" else None, ndmin=1
            )
    except (OSError, ValueError, UnicodeError, EOFError) as exc:
        raise ValueError(f"cannot load {p}: {exc}") from exc
    if not isinstance(values, np.ndarray):
        values.close()  # np.load can return an archive even with a .npy suffix.
        raise ValueError(f"{p}: expected one sample array, not an archive")
    if values.ndim != 1:
        raise ValueError(f"{p}: expected a one-dimensional sample array, got {values.shape}")
    if np.iscomplexobj(values) or values.dtype.kind not in "biuf":
        raise ValueError(f"{p}: samples must be real numbers")
    values = np.asarray(values, dtype=np.float64)
    if values.size < 2:
        raise ValueError(f"{p}: need at least 2 sample points")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{p}: samples must all be finite")
    return values


def _write_json(path: str, report: Dict[str, Any]) -> None:
    """Reject nonfinite results before creating or replacing a report."""
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(payload, encoding="utf-8")


def _report(command: str) -> Dict[str, Any]:
    from . import __version__

    return {
        "schema_version": 1,
        "pbss_version": __version__,
        "command": command,
        "banner": BANNER,
    }


def _array_metadata(values: np.ndarray) -> Dict[str, Any]:
    canonical = np.asarray(values, dtype="<f8")
    return {
        "n_points": int(canonical.size),
        "sha256_float64_le": hashlib.sha256(canonical.tobytes(order="C")).hexdigest(),
    }


def _grid_metadata(u: np.ndarray, source: str = "") -> Dict[str, Any]:
    return {
        **_array_metadata(u),
        "source": source or "uniform_linspace",
        "minimum": float(u[0]),
        "maximum": float(u[-1]),
        "quadrature": "composite_trapezoid",
    }


def _integer_at_least(minimum: int):
    def parse(value: str) -> int:
        try:
            result = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError("must be an integer") from exc
        if result < minimum:
            raise argparse.ArgumentTypeError(f"must be at least {minimum}")
        return result

    return parse


def _finite_float(value: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a finite number") from exc
    if not np.isfinite(result):
        raise argparse.ArgumentTypeError("must be a finite number")
    return result


def _positive_float(value: str) -> float:
    result = _finite_float(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def _positive_decimal(value: str) -> str:
    """Keep a decimal window input for the longdouble arithmetic evaluator."""
    try:
        parsed = np.longdouble(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise argparse.ArgumentTypeError("must be a finite positive number") from exc
    if not np.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return value


def _nonnegative_float(value: str) -> float:
    result = _finite_float(value)
    if result < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return result


def _x_max(value: str) -> float:
    result = _finite_float(value)
    if result < 3:
        raise argparse.ArgumentTypeError("must be at least 3")
    return result


def _unit_grid(n: int) -> np.ndarray:
    from .probes import sample_grid

    return sample_grid(int(n))


def cmd_project(args: argparse.Namespace) -> int:
    from .projection import project

    q = _load_array(args.input)
    if args.u:
        u = _load_array(args.u)
    else:
        u = _unit_grid(q.size)
    if u.size != q.size:
        raise ValueError("u and q must have the same length")
    r = project(q, u, degree=int(args.degree), T=float(args.T))
    out = {
        **_report("project"),
        "degree": r.degree,
        "T": r.T,
        "n_points": r.n_points,
        "energy": r.energy,
        "l2_norm_sq": r.l2_norm_sq,
        "energy_ratio_Rd": r.energy_ratio,
        "scaled_strength_Sd": r.scaled_strength,
        "coeffs": r.coeffs.tolist(),
        "coefficient_basis": "orthonormal_shifted_legendre",
        "coefficient_meaning": "weighted_least_squares_polynomial_coefficients",
        "projection_method": "weighted_qr",
        "input": {"source": args.input, **_array_metadata(q)},
        "grid": _grid_metadata(u, args.u),
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
        _write_json(args.json_out, out)
        print(f"Wrote {args.json_out}", file=sys.stderr)
    return 0


def cmd_diagnose(args: argparse.Namespace) -> int:
    from .diagnostic import assess_demo
    from .probes import (
        probe_critical_line_mode,
        probe_defective,
        probe_high_frequency,
        probe_prime_residual,
        sample_grid,
    )
    from .projection import project

    degree = int(args.degree)
    T = float(args.T)
    rows = []
    is_demo = not bool(args.input)
    if is_demo:
        u = sample_grid(args.n_points if args.n_points is not None else 2048)
        specs = [
            ("high_frequency", probe_high_frequency(u, waves=48), T),
            ("critical_line", probe_critical_line_mode(u, T=T), T),
            ("defective", probe_defective(u, waves=48, defect_degree=1, defect_weight=2.5), T),
        ]
        if args.with_primes:
            q_pr, T_pr = probe_prime_residual(
                u, x_max=args.x_max if args.x_max is not None else 1e5
            )
            specs.append(("prime_residual", q_pr, T_pr))
    else:
        q = _load_array(args.input)
        u = _load_array(args.u) if args.u else _unit_grid(q.size)
        specs = [("input", q, T)]

    print("PBSS diagnose")
    print(f"degree={degree}  n_points={u.size}")
    for name, q, Tv in specs:
        r = project(q, u, degree=degree, T=float(Tv))
        rows.append(
            {
                "name": name,
                "T": float(Tv),
                "degree": r.degree,
                "n_points": r.n_points,
                "energy": r.energy,
                "l2_norm_sq": r.l2_norm_sq,
                "Rd": r.energy_ratio,
                "Sd": r.scaled_strength,
                "samples": _array_metadata(q),
            }
        )
        print(f"  {name:16s}  R_d={r.energy_ratio:.6e}  S_d={r.scaled_strength:.6e}  T={Tv:.4g}")
    assessment = None
    if is_demo:
        assessment = assess_demo(rows, degree=degree, n_points=u.size, T=T)
        sep = assessment["energy_ratio_defective_minus_hf"]
        print(f"  separation (defective - HF) = {sep:.6e}")
        print(f"  classifier_ok = {assessment['classifier_ok']} ({assessment['status']})")
        for warning in assessment["warnings"]:
            print(f"  warning: {warning}")
    print(BANNER)
    if args.json_out:
        report = {
            **_report("diagnose"),
            "mode": "demo" if is_demo else "input",
            "degree": degree,
            "T": T,
            "n_points": int(u.size),
            "projection_method": "weighted_qr",
            "grid": _grid_metadata(u, args.u),
            "rows": rows,
        }
        if assessment is not None:
            report["assessment"] = assessment
            report["demo_parameters"] = {
                "waves": 48,
                "critical_ordinate": 14.134725,
                "defect_degree": 1,
                "defect_weight": 2.5,
                "with_primes": bool(args.with_primes),
            }
            if args.with_primes:
                report["prime_parameters"] = {
                    "x_max": args.x_max if args.x_max is not None else 1e5,
                    "detrend": "deg1",
                    "smooth": 3,
                }
        else:
            report["input"] = {"source": args.input, **_array_metadata(specs[0][1])}
        _write_json(args.json_out, report)
    return 0 if assessment is None or assessment["classifier_ok"] else 2


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
            _write_json(args.json_out, {**_report("sensitivity"), **rep})
        return 0 if rep["confirmed"] else 2

    rep = sensitivity_experiment(
        n_per_class=args.n_per_class if args.n_per_class is not None else 200,
        noise=float(args.noise),
        seed=int(args.seed),
        k=args.k if args.k is not None else 4.0,
        sigma=args.sigma if args.sigma is not None else 6.0,
    )
    print("PBSS Gamma-weight sensitivity")
    print(f"  flat d′  = {rep['flat_dprime']:.6f}")
    print(f"  gamma d′ = {rep['gamma_dprime']:.6f}")
    print(f"  relative gain = {rep['relative_gain_percent']:.2f}%")
    print(f"  peak u⋆ = {rep['peak_u_star']:.4f}  (k={rep['k']}, σ={rep['sigma']})")
    print(BANNER)
    if args.json_out:
        _write_json(args.json_out, {**_report("sensitivity"), **rep})
    return 0


def cmd_scorecard(args: argparse.Namespace) -> int:
    from .beurling import beurling_theta_residual, build_system_primes, default_battery_specs
    from .probes import primes_upto, sample_grid
    from .projection import project

    u = sample_grid(int(args.n_points))
    x_max = float(args.x_max)
    degree = int(args.degree)
    systems = default_battery_specs()
    ordinary = primes_upto(int(np.floor(x_max)))
    print(f"PBSS Beurling scorecard  x_max={x_max:g}  d={degree}")
    rows = []
    failed = False
    for spec in systems:
        label = spec["name"]
        try:
            primes = build_system_primes(spec, ordinary, x_max=x_max)
            q, Tv, residual_metadata = beurling_theta_residual(
                u, primes, x_max=x_max, detrend="deg1"
            )
            r = project(q, u, degree=degree, T=Tv)
            rows.append({
                "system": label,
                "status": "ok",
                "Rd": r.energy_ratio,
                "Sd": r.scaled_strength,
                "T": Tv,
                "energy": r.energy,
                "l2_norm_sq": r.l2_norm_sq,
                "n_primes": int(primes.size),
                "residual": residual_metadata,
            })
            print(f"  {label:20s}  R_d={r.energy_ratio:.6f}  T={Tv:.4g}")
        except Exception as e:
            failed = True
            rows.append({"system": label, "status": "error", "error": str(e)})
            print(f"  {label:20s}  ERROR: {e}")
    print(BANNER)
    if args.json_out:
        _write_json(args.json_out, {
            **_report("scorecard"),
            "status": "failed" if failed else "ok",
            "x_max": x_max,
            "degree": degree,
            "n_points": int(u.size),
            "projection_method": "weighted_qr",
            "system_specs": systems,
            "grid": _grid_metadata(u),
            "rows": rows,
        })
    return 2 if failed else 0


def cmd_version(_: argparse.Namespace) -> int:
    from . import __version__

    print(f"pbss {__version__}")
    print(BANNER)
    return 0


def cmd_hstar(args: argparse.Namespace) -> int:
    from .hstar_anticancellation import hstar_anticancellation_report

    rep = hstar_anticancellation_report(
        write_results=bool(args.write_results),
        results_dir=args.results_dir or None,
        heavy=bool(args.heavy) or bool(args.write_results),
    )
    print("PBSS hstar (model cancellation audit)")
    print(f"  lemma={rep['lemma']['id']}  status={rep['candidate_status']}")
    print(f"  axis={rep['axis']}")
    print(f"  grid_n_rows={rep.get('grid_n_rows')}  constraint_frac={rep.get('grid_constraint_fraction')}")
    print(f"  all_constraint_signals={rep['all_constraint_signals']}")
    print(f"  partial_cancel_monotone={rep.get('partial_cancel_monotone')}")
    print(f"  rh_claimed={rep['rh_claimed']}  b_res_solved={rep['b_res_solved']}")
    eps = rep.get("eps_star") or {}
    print(f"  eps_star_half_min={eps.get('eps_star_global_half_min')}")
    for r in rep["multi_T_proxy"]:
        print(
            f"  T={r['T']:.1f}  R_pure={r['R_d_pure']:.4e}  "
            f"R_cancel={r['R_d_cancelled']:.4e}  C_pure={r['C_d_pure']:.4e}"
        )
    print(BANNER)
    if args.json_out:
        _write_json(args.json_out, {**_report("hstar"), **rep})
        print(f"Wrote {args.json_out}", file=sys.stderr)
    if rep.get("artifacts"):
        print(f"Artifacts: {rep['artifacts']}", file=sys.stderr)
    return 0


def cmd_moments(args: argparse.Namespace) -> int:
    from .arithmetic_moments import arithmetic_moments

    result = arithmetic_moments(
        args.T, degree=args.degree, max_prime_limit=args.max_prime_limit,
    )
    report = {**_report("moments"), **result.to_dict()}
    metadata = report["metadata"]
    print("PBSS moments (continuous arithmetic integrals)")
    print(f"  T={report['T_decimal']}  degree={result.degree}  primes={metadata['n_primes']}")
    print(f"  working_precision_bits={metadata['working_precision_bits']}")
    for k, coefficient in enumerate(report["coefficients_decimal"]):
        print(f"  c_{k}={coefficient}")
    print(f"  l2_norm_sq={report['l2_norm_sq_decimal']}")
    print(f"  affine_detrended_l2_norm_sq={report['affine_detrended_l2_norm_sq_decimal']}")
    print(f"  affine_detrended_Rd={report['affine_detrended_energy_ratio']}")
    print("  Finite-window evaluation; no all-window growth bound established.")
    print(BANNER)
    if args.json_out:
        _write_json(args.json_out, report)
        print(f"Wrote {args.json_out}", file=sys.stderr)
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
    pr.add_argument("--degree", "-d", type=_integer_at_least(0), default=4)
    pr.add_argument("--T", type=_positive_float, default=20.0, help="log-window length for S_d")
    pr.add_argument("--json-out", default="", help="write JSON report")
    pr.set_defaults(func=cmd_project)

    di = sub.add_parser("diagnose", help="demo scorecard or diagnose one residual")
    mode = di.add_mutually_exclusive_group()
    mode.add_argument("--demo", action="store_true", help="built-in HF / CL / defective probes (default)")
    mode.add_argument("--input", "-i", default="", help="q array; grid length inferred from samples")
    di.add_argument("--degree", "-d", type=_integer_at_least(0), default=4)
    di.add_argument("--T", type=_positive_float, default=20.0)
    grid = di.add_mutually_exclusive_group()
    grid.add_argument("--u", default="", help="u grid for --input (default linspace [0,1])")
    grid.add_argument("--n-points", type=_integer_at_least(2), default=None, help="demo grid size (default 2048)")
    di.add_argument("--with-primes", action="store_true", help="add a prime residual to the demo")
    di.add_argument("--x-max", type=_x_max, default=None, help="prime range with --with-primes (default 1e5)")
    di.add_argument("--json-out", default="")
    di.set_defaults(func=cmd_diagnose)

    se = sub.add_parser(
        "sensitivity",
        help="Gamma vs flat weight discriminability (offline vs online)",
    )
    se.add_argument("--confirm-53", action="store_true", help="confirm ≥53%% relative gain claim")
    se.add_argument("--noise", type=_nonnegative_float, default=0.55)
    se.add_argument("--seed", type=_integer_at_least(0), default=20260522)
    se.add_argument("--n-per-class", type=_integer_at_least(2), default=None, help="ensemble size per class (default 200)")
    se.add_argument("--k", type=_positive_float, default=None, help="Gamma shape (default 4)")
    se.add_argument("--sigma", type=_positive_float, default=None, help="Gamma rate (default 6)")
    se.add_argument("--json-out", default="")
    se.set_defaults(func=cmd_sensitivity)

    sc = sub.add_parser("scorecard", help="ordinary vs gapped/thinned Beurling R_d")
    sc.add_argument("--x-max", type=_x_max, default=1e6)
    sc.add_argument("--degree", "-d", type=_integer_at_least(0), default=4)
    sc.add_argument("--n-points", type=_integer_at_least(2), default=2048)
    sc.add_argument("--json-out", default="")
    sc.set_defaults(func=cmd_scorecard)

    hs = sub.add_parser(
        "hstar",
        help="corrected model cancellation audit; no arithmetic coercivity claim",
    )
    hs.add_argument("--write-results", action="store_true", help="write corrected results/hstar_corrected/")
    hs.add_argument("--heavy", action="store_true", help="full (T,σ,d) ProcessPool grid")
    hs.add_argument("--results-dir", default="", help="override results directory")
    hs.add_argument("--json-out", default="")
    hs.set_defaults(func=cmd_hstar)

    mo = sub.add_parser("moments", help="signed absolute theta moments from complete prime intervals")
    mo.add_argument("--T", type=_positive_decimal, default="10", help="positive logarithmic window length")
    mo.add_argument("--degree", "-d", type=_integer_at_least(0), default=4, help="maximum degree (0..32)")
    mo.add_argument("--max-prime-limit", type=_integer_at_least(2), default=1_000_000,
                    help="explicit bound on the complete prime prefix to generate")
    mo.add_argument("--json-out", default="", help="write decimal coefficients, norms and precision metadata")
    mo.set_defaults(func=cmd_moments)

    return p


def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cmd == "diagnose":
        if args.input:
            if args.n_points is not None:
                parser.error("--n-points applies only to the demo; input grid length comes from the samples")
            if args.with_primes:
                parser.error("--with-primes applies only to the demo")
        elif args.u:
            parser.error("--u requires --input")
        if args.x_max is not None and not args.with_primes:
            parser.error("--x-max requires --with-primes")
    if args.cmd == "sensitivity" and args.confirm_53:
        if any(getattr(args, key) is not None for key in ("n_per_class", "k", "sigma")):
            parser.error("--confirm-53 uses a fixed ensemble; --n-per-class, --k and --sigma cannot override it")
    if args.cmd in ("project", "diagnose") and args.json_out:
        destination = Path(args.json_out).expanduser().resolve()
        for source in (args.input, args.u):
            if source and Path(source).expanduser().resolve() == destination:
                parser.error("--json-out must not overwrite an input array")
    try:
        return int(args.func(args))
    except (OSError, ValueError, OverflowError) as exc:
        parser.exit(2, f"pbss: error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
