"""Independent checks of the psi - theta correction (docs/PSI_THETA_CORRECTION.md).

Each test compares the experiment's code path with a reference that does not share it:
brute-force factorization, scipy quadrature over every jump, prime data for the zeta-side
transform, a published constant, the shipped theta evaluator, or the other implementation.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from pbss.probes import primes_upto

mpmath = pytest.importorskip("mpmath")

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def exp():
    path = ROOT / "experiments" / "check_psi_theta_correction.py"
    spec = importlib.util.spec_from_file_location("check_psi_theta_correction", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _factor(n: int) -> dict[int, int]:
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def _higher_prime_powers(x: int) -> list[tuple[int, int]]:
    """(n, p) for every prime power n = p^j <= x with j >= 2, by trial factorization."""
    rows = []
    for n in range(4, x + 1):
        f = _factor(n)
        if len(f) == 1:
            ((p, j),) = f.items()
            if j >= 2:
                rows.append((n, p))
    return rows


def _delta_by_quadrature(x: int, k: int) -> float:
    """Delta_k = (1/T) int_0^T e^{-y/2} (psi - theta)(e^y) phi_k(y/T) dy, split at every jump."""
    T = math.log(x)
    jumps = _higher_prime_powers(x)
    edges = [math.log(n) for n, _ in jumps] + [T]
    phi = lambda y: math.sqrt(2 * k + 1) * eval_legendre(k, 2 * y / T - 1)
    level, total = 0.0, 0.0
    for (_, p), a, b in zip(jumps, edges[:-1], edges[1:]):
        level += math.log(p)
        total += quad(lambda y: level * math.exp(-y / 2) * phi(y), a, b, epsabs=1e-14, epsrel=1e-13)[0]
    return total / T


def test_mobius_matches_factorization(exp):
    for n in range(1, 301):
        f = _factor(n)
        expected = 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)
        assert exp.mobius(n) == expected, n


def test_jumps_are_exactly_the_higher_prime_powers(exp):
    x = 20_000
    with mpmath.workdps(30):
        events = exp.prime_power_events(mpmath.mpf(x))
        got = [(int(mpmath.nint(mpmath.exp(pos))), int(mpmath.nint(mpmath.exp(lp)))) for pos, lp in events]
    assert got == _higher_prime_powers(x)


def test_psi_steps_reproduce_the_von_mangoldt_sum(exp):
    x = 5_000
    values, weights = exp.step_events(x, "psi")
    prime_powers = [n for n in range(2, x + 1) if len(_factor(n)) == 1]
    assert values.tolist() == prime_powers
    psi_x = sum(math.log(next(iter(_factor(n)))) for n in prime_powers)
    assert float(np.sum(weights)) == pytest.approx(psi_x, rel=1e-14)
    theta_values, _ = exp.step_events(x, "theta")
    assert theta_values.tolist() == primes_upto(x).tolist()


@pytest.mark.parametrize("x", [1_000, 30_000])
def test_exact_correction_matches_direct_quadrature(exp, x):
    result = exp.exact_corrections(x, dps=30)
    for k in range(1, 5):
        assert float(result["delta"][k]) == pytest.approx(_delta_by_quadrature(x, k), abs=1e-12)


def test_three_term_identity_for_delta_2(exp):
    result = exp.exact_corrections(10**6, dps=40)
    with mpmath.workdps(40):
        assert abs(sum(result["split"]) - result["delta"][2]) < mpmath.mpf(10) ** -30


def test_campaign_endpoint_value(exp):
    result = exp.exact_corrections(5 * 10**10, dps=30)
    assert result["events"] == 20_639
    assert float(result["delta"][2]) == pytest.approx(-0.201121262828, abs=5e-13)
    expected_terms = [0.017070, -0.498214, 0.280023]
    assert [float(t) for t in result["split"]] == pytest.approx(expected_terms, abs=5e-7)


def test_closed_form_constants(exp):
    c = exp.closed_form_constants(dps=25)
    # Rosser-Schoenfeld: sum_{p<=x} log p / p = log x + E + o(1), with -E = gamma + B = 1.3325822757332208...
    with mpmath.workdps(25):
        assert float(mpmath.euler + c["B"]) == pytest.approx(1.33258227573322088, abs=1e-15)
    assert float(c["C3"]) == pytest.approx(2.476836385022314, abs=1e-14)
    assert float(c["K0"]) == pytest.approx(0.288508218578187, abs=1e-14)
    assert float(c["sqrt5_K0"]) == pytest.approx(0.645123988808194, abs=1e-14)


@pytest.mark.parametrize("eps,x,tol", [(1.0, 10**8, 1e-8), (0.6, 10**10, 1e-6)])
def test_zeta_side_transform_matches_prime_data(exp, eps, x, tol):
    """G(eps), the formula whose eps -> 0 limit gives K0, against int_0^T (g - 1) e^{-eps y} dy from primes.
    The omitted tail beyond T is below tol for these (eps, x)."""
    with mpmath.workdps(25):
        closed = exp.G_eps(mpmath.mpf(eps))
        X = mpmath.mpf(x)
        T = mpmath.log(X)
        s = mpmath.mpf(1) / 2 + eps
        events = exp.prime_power_events(X)
        level, total = mpmath.mpf(0), mpmath.mpf(0)
        for i, (y0, lp) in enumerate(events):
            level += lp
            y1 = events[i + 1][0] if i + 1 < len(events) else T
            total += level * (mpmath.exp(-s * y0) - mpmath.exp(-s * y1)) / s
        direct = total - (1 - mpmath.exp(-eps * T)) / eps
        assert float(abs(closed - direct)) < tol


def test_theta_path_matches_shipped_evaluator(exp):
    check = exp.self_check_theta(10**4)
    assert check["status"] == "ok", check


def test_two_implementations_of_the_correction_agree(exp):
    """Prime-interval longdouble evaluator (theta and psi run separately) versus the mpmath jump formula."""
    x = 10**5
    theta, psi = exp.continuous_moments(x, "theta"), exp.continuous_moments(x, "psi")
    result = exp.exact_corrections(x, dps=30)
    for k in range(1, 5):
        assert psi["c"][k] - theta["c"][k] == pytest.approx(float(result["delta"][k]), abs=1e-12)


def test_psi_residual_has_larger_R4_than_theta(exp):
    """Documented plateau-test values: adding the prime powers back raises R_4."""
    expected = {10**3: (0.029509, 0.252845), 10**4: (0.073942, 0.355687), 10**5: (0.102701, 0.422086)}
    for x, (r4_theta, r4_psi) in expected.items():
        theta, psi = exp.continuous_moments(x, "theta"), exp.continuous_moments(x, "psi")
        assert theta["R4"] == pytest.approx(r4_theta, abs=5e-7)
        assert psi["R4"] == pytest.approx(r4_psi, abs=5e-7)
        assert psi["R4"] > 4 * theta["R4"]
