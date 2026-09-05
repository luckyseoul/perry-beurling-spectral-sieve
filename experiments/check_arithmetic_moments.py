#!/usr/bin/env python3
"""Independent arbitrary-precision corroboration of finite arithmetic moments.

Requires optional mpmath. Evaluates global finite prime antiderivatives,
compares the interval implementation, and checks a small case by direct
integration. None of these finite checks proves the B-ABS growth condition.
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

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _primitive_polynomial(mp, coefficients, exponent):
    """Solve G'+sG=P downward; exp(s*y)*G(y) is an exact antiderivative."""
    result = [mp.mpf(0)] * len(coefficients)
    for j in range(len(coefficients) - 1, -1, -1):
        derivative = (j + 1) * result[j + 1] if j + 1 < len(result) else 0
        result[j] = (coefficients[j] - derivative) / exponent
    return result


def _reference(mp, T, degree, primes):
    prime_logs = [(int(p), mp.log(int(p))) for p in primes]
    prime_logs = [(p, lp) for p, lp in prime_logs if lp <= T]
    theta_T = mp.fsum(lp for p, lp in prime_logs)
    coefficients = []
    for k in range(max(degree, 1) + 1):
        polynomial = [mp.sqrt(2 * k + 1) * (-1)**(k + j) * mp.binomial(k, j)
                      * mp.binomial(k + j, j) / T**j for j in range(k + 1)]
        minus = _primitive_polynomial(mp, polynomial, -mp.mpf("0.5"))
        plus = _primitive_polynomial(mp, polynomial, mp.mpf("0.5"))
        upper_prime_term = theta_T * mp.exp(-T / 2) * mp.polyval(minus[::-1], T)
        lower_prime_terms = mp.fsum(lp / mp.sqrt(p) * mp.polyval(minus[::-1], lp) for p, lp in prime_logs)
        main_term = mp.exp(T / 2) * mp.polyval(plus[::-1], T) - plus[0]
        coefficients.append((upper_prime_term - lower_prime_terms - main_term) / T)
    # Uncompleted, globally cancelling norm antiderivative, at high precision.
    edges = [mp.mpf(0)] + [lp for p, lp in prime_logs] + [T]
    cumulative = mp.mpf(0)
    intervals = []
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        if i:
            cumulative += prime_logs[i - 1][1]
        intervals.append(cumulative**2 * (mp.exp(-a) - mp.exp(-b))
                         - 2 * cumulative * (b - a) + mp.exp(b) - mp.exp(a))
    norm = mp.fsum(intervals) / T
    affine_norm = norm - coefficients[0]**2 - coefficients[1]**2
    return coefficients[:degree + 1], norm, affine_norm


def _case(job):
    import mpmath as mp
    import numpy as np
    from pbss.arithmetic_moments import arithmetic_moments
    from pbss.probes import primes_upto

    x, degree, digits = job
    T = np.log(np.longdouble(x))
    result = arithmetic_moments(T, degree=degree)
    report = result.to_dict()
    primes = primes_upto(x + 1)
    with mp.workdps(digits):
        precise_T = mp.mpf(report["T_decimal"])
        reference, norm, affine_norm = _reference(mp, precise_T, degree, primes)
        errors = [abs(mp.mpf(actual) - expected) for actual, expected in zip(report["coefficients_decimal"], reference)]
        norm_error = abs(mp.mpf(report["l2_norm_sq_decimal"]) - norm)
        affine_error = abs(mp.mpf(report["affine_detrended_l2_norm_sq_decimal"]) - affine_norm)
        assert max(errors) < mp.mpf("5e-12")
        assert norm_error < mp.mpf("5e-12") and affine_error < mp.mpf("5e-12")
        with mp.workdps(max(40, digits // 2)):
            lower_precision, _, _ = _reference(mp, mp.mpf(report["T_decimal"]), degree, primes)
        precision_difference = max(abs(a - b) for a, b in zip(reference, lower_precision))
        return {
            "x_nominal": x, "degree": degree, "interval_result": report,
            "reference_coefficients_decimal": [mp.nstr(c, digits) for c in reference],
            "reference_l2_norm_sq_decimal": mp.nstr(norm, digits),
            "reference_affine_l2_norm_sq_decimal": mp.nstr(affine_norm, digits),
            "coefficient_absolute_errors_decimal": [mp.nstr(e, digits) for e in errors],
            "max_coefficient_absolute_error": float(max(errors)),
            "norm_absolute_error": float(norm_error),
            "affine_norm_absolute_error": float(affine_error),
            "reference_precision_comparison_decimal": mp.nstr(precision_difference, digits),
        }


def _quadrature(digits):
    import mpmath as mp
    from pbss.probes import primes_upto

    with mp.workdps(digits):
        T, degree = mp.log(30), 4
        primes = primes_upto(30)
        reference, norm, _ = _reference(mp, T, degree, primes)
        logs = [mp.log(int(p)) for p in primes]
        edges = [mp.mpf(0)] + logs + [T]
        cumulative = mp.mpf(0)
        coefficients = [mp.mpf(0)] * (degree + 1)
        integrated_norm = mp.mpf(0)
        for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            if i:
                cumulative += logs[i - 1]
            q = lambda y: cumulative * mp.exp(-y / 2) - mp.exp(y / 2)
            for k in range(degree + 1):
                coefficients[k] += mp.quad(lambda y: q(y) * mp.sqrt(2 * k + 1) * mp.legendre(k, 2 * y / T - 1), [a, b]) / T
            integrated_norm += mp.quad(lambda y: q(y)**2, [a, b]) / T
        coefficient_error = max(abs(a - b) for a, b in zip(coefficients, reference))
        norm_error = abs(integrated_norm - norm)
        assert max(coefficient_error, norm_error) < mp.power(10, -digits + 15)
        return {
            "x": 30, "degree": degree, "n_intervals": len(edges) - 1,
            "max_coefficient_error_decimal": mp.nstr(coefficient_error, digits),
            "norm_error_decimal": mp.nstr(norm_error, digits),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 2))
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--out", type=Path, default=ROOT / "results/arithmetic_moments_20260905/checks.json")
    args = parser.parse_args()
    if args.workers < 1 or args.digits < 60:
        parser.error("workers must be positive and digits must be at least 60")
    import mpmath
    import numpy
    import scipy
    jobs = [(x, degree, args.digits) for x in (2, 10, 100, 1000, 10000) for degree in (2, 4, 6)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        quadrature = pool.submit(_quadrature, args.digits)
        rows = list(pool.map(_case, jobs))
        direct_check = quadrature.result()
    sources = ["src/pbss/arithmetic_moments.py", "experiments/check_arithmetic_moments.py", "docs/EXACT_ARITHMETIC_MOMENTS.md"]
    report = {
        "status": "finite_prime_integral_identities_numerically_corroborated",
        "certified_interval_result": False, "RH_proved": False, "B_ABS_growth_bound_proved": False,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "python": sys.version, "workers": args.workers, "backend": "ProcessPool_mpmath",
        "digits": args.digits,
        "versions": {"mpmath": mpmath.__version__, "numpy": numpy.__version__, "scipy": scipy.__version__},
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources},
        "independent_direct_quadrature": direct_check, "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"artifact": str(args.out), "cases": len(rows), "host": report["host"],
                      "workers": args.workers, "max_coefficient_error": max(r["max_coefficient_absolute_error"] for r in rows),
                      "max_norm_error": max(r["norm_absolute_error"] for r in rows)}, indent=2))


if __name__ == "__main__":
    main()
