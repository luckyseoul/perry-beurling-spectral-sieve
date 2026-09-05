#!/usr/bin/env python3
"""Corroborate the continuous phase obstruction with independent mpmath sums.

Optional dependency: mpmath. Mathematical statements are proved in
OFF_CRITICAL_PHASE_OBSTRUCTION.md; finite checks here are not certificates of
an asymptotic theorem or claims about arithmetic RH. Independent cases run in
a ProcessPool; use --workers according to the budget on the executing node.
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
import sys

# Pin numerical-library inner threads; process workers own independent cases.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _exact_coefficients(mp, T, A, t, degree):
    z = T * mp.mpc(A, t)
    theta = t * T
    coefficients = []
    for k in range(degree + 1):
        upper = mp.mpc(0)
        lower = mp.mpc(0)
        for j in range(k + 1):
            derivative = mp.sqrt(2 * k + 1) * mp.factorial(k + j) / (mp.factorial(j) * mp.factorial(k - j))
            term = derivative / z**(j + 1)
            upper += (-1)**j * term
            lower += term
        coefficients.append(mp.im(mp.exp(1j * theta) * upper - mp.exp(-A * T) * (-1)**k * lower))
    lam = mp.mpc(A, t)
    H = 1 / A - mp.re(mp.exp(2j * theta) / lam)
    H0 = 1 / A - mp.re(1 / lam)
    norm = (H - mp.exp(-2 * A * T) * H0) / (4 * T)
    return coefficients, norm


def _case(job):
    import mpmath as mp
    import numpy as np
    from pbss.off_critical_exact import off_critical_projection

    kind, n, degree, digits, A_text, t_text = job
    with mp.workdps(digits):
        A, t = mp.mpf(A_text), mp.mpf(t_text)
        radius_sq = A * A + t * t
        beta = mp.atan2(t, A)
        if kind == "maximal":
            T = n * mp.pi / t
        elif kind == "leading_cancel":
            T = (beta + n * mp.pi) / t
        else:
            theta = beta + n * mp.pi
            mu = mp.mpf(degree * (degree + 2)) / 2
            T = (theta + mp.sqrt(theta**2 - 4 * t * t * mu / radius_sq)) / (2 * t)
        coefficients, norm = _exact_coefficients(mp, T, A, t, degree)
        ratio = mp.fsum(c * c for c in coefficients) / norm
        limit_TR = 4 * A * (degree + 1)**2 / radius_sq
        limit_T3 = A * (degree * (degree + 1) * (degree + 2))**2 / (3 * radius_sq**2)
        expected = limit_TR if kind == "maximal" else limit_T3 * (4 if kind == "leading_cancel" else 1)
        scaled = T * ratio if kind == "maximal" else T**3 * ratio

        # Compare at the same binary64 input, independently of phase rounding.
        float_T, float_A, float_t = float(T), float(A), float(t)
        reference_coefficients, reference_norm = _exact_coefficients(
            mp, mp.mpf(float_T), mp.mpf(float_A), mp.mpf(float_t), degree)
        reference_ratio = mp.fsum(c * c for c in reference_coefficients) / reference_norm
        library = off_critical_projection(float_T, A=float_A, t=float_t, degree=degree)
        max_coefficient_error = max(abs(float(c) - actual) for c, actual in zip(reference_coefficients, library.scaled_coeffs))
        ratio_error = abs(float(reference_ratio) - library.energy_ratio)
        assert np.allclose([float(c) for c in reference_coefficients], library.scaled_coeffs, rtol=1e-8, atol=2e-13)
        assert np.isclose(float(reference_ratio), library.energy_ratio, rtol=1e-8, atol=2e-13)
        return {
            "kind": kind, "n": n, "degree": degree,
            "T_decimal": mp.nstr(T, digits), "R_decimal": mp.nstr(ratio, digits),
            "scaled_statistic": mp.nstr(scaled, digits), "proved_limit": mp.nstr(expected, digits),
            "relative_distance_to_limit": float(abs(scaled / expected - 1)),
            "binary64_T": float_T, "float64_R": library.energy_ratio,
            "float64_R_absolute_error": ratio_error,
            "max_float64_coefficient_absolute_error": max_coefficient_error,
        }


def _quadrature(job):
    import mpmath as mp
    digits, A_text, t_text = job
    with mp.workdps(digits):
        A, t, T, degree = mp.mpf(A_text), mp.mpf(t_text), mp.mpf(12), 4
        coefficients, norm = _exact_coefficients(mp, T, A, t, degree)
        points = [mp.mpf(k) / 32 for k in range(33)]
        q = lambda u: mp.exp(A * T * (u - 1)) * mp.sin(t * T * u)
        errors = []
        for k, coefficient in enumerate(coefficients):
            integral = mp.quad(lambda u: q(u) * mp.sqrt(2 * k + 1) * mp.legendre(k, 2 * u - 1), points)
            errors.append(abs(integral - coefficient))
        norm_error = abs(mp.quad(lambda u: q(u)**2, points) - norm)
        threshold = mp.power(10, -digits + 15)
        assert max(errors) < threshold and norm_error < threshold
        return {
            "T": 12.0, "degree": degree, "subintervals": 32,
            "max_coefficient_absolute_error": mp.nstr(max(errors), digits),
            "norm_absolute_error": mp.nstr(norm_error, digits),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 2))
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--A", default="0.4")
    parser.add_argument("--t", default="14.134725")
    parser.add_argument("--out", type=Path, default=ROOT / "results/off_critical_phase_obstruction_20260905/exact_checks.json")
    args = parser.parse_args()
    if args.workers < 1 or args.digits < 40:
        parser.error("workers must be positive and digits must be at least 40")
    import mpmath
    import numpy
    import scipy
    jobs = [(kind, n, d, args.digits, args.A, args.t)
            for d in (1, 4)
            for kind in ("maximal", "leading_cancel", "minimum_T3")
            for n in (100, 400, 1600, 6400)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        quadrature = pool.submit(_quadrature, (args.digits, args.A, args.t))
        rows = list(pool.map(_case, jobs))
        integral_check = quadrature.result()
    inputs = ["src/pbss/off_critical_exact.py", "experiments/check_off_critical_exact.py", "docs/OFF_CRITICAL_PHASE_OBSTRUCTION.md"]
    report = {
        "status": "analytic_model_theorem_with_numerical_corroboration",
        "arithmetic_B_RES_solved": False, "RH_proved": False,
        "certification": "80-digit numerical checks are not interval certificates; proof is in the document",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(), "python": sys.version, "workers": args.workers,
        "backend": "ProcessPool_mpmath", "digits": args.digits,
        "versions": {"mpmath": mpmath.__version__, "numpy": numpy.__version__, "scipy": scipy.__version__},
        "A_decimal": args.A, "t_decimal": args.t,
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in inputs},
        "independent_quadrature": integral_check, "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"artifact": str(args.out), "cases": len(rows), "host": report["host"],
                      "workers": args.workers, "max_float64_R_error": max(row["float64_R_absolute_error"] for row in rows)}, indent=2))


if __name__ == "__main__":
    main()
