"""
Tests for src/pbss/rh_moment_asymptotics.py (docs/RH_PROOF_SIMPLIFICATIONS.md).

References are independent of the module under test: published values from Rubinstein's
exposition, the shipped interval evaluator pbss.arithmetic_moments, the PR #2 artifact
(results/psi_theta_correction_20261120), scipy quadrature and mpmath.
RH-conditional checks (envelopes) are labelled; they test consistency, not RH.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_legendre

from pbss import rh_moment_asymptotics as R
from pbss.arithmetic_moments import arithmetic_moments

ROOT = Path(__file__).resolve().parents[1]
PR2 = ROOT / "results" / "psi_theta_correction_20261120" / "psi_theta_correction.json"


@pytest.fixture(scope="module")
def pr2():
    return json.loads(PR2.read_text())


def test_phi_coefficients_match_shifted_legendre():
    u = np.linspace(0, 1, 7)
    for k in range(7):
        a = R.phi_coefficients(k)
        mine = sum(a[j] * u ** j for j in range(k + 1))
        assert np.allclose(mine, math.sqrt(2 * k + 1) * eval_legendre(k, 2 * u - 1), atol=1e-11)
        d = R.phi_derivatives(k, 1.0)
        assert d[0] == pytest.approx(math.sqrt(2 * k + 1))


def test_prime_power_events_brute_force():
    n, w = R.prime_power_events(1000, "psi")
    brute = []
    for m in range(2, 1001):
        p = next(q for q in range(2, m + 1) if m % q == 0)
        x = m
        while x % p == 0:
            x //= p
        if x == 1:
            brute.append((m, math.log(p)))
    assert [int(v) for v in n] == [b[0] for b in brute]
    assert np.allclose(w, [b[1] for b in brute])
    nd, _ = R.prime_power_events(1000, "psi_minus_theta")
    nt, _ = R.prime_power_events(1000, "theta")
    assert sorted(np.concatenate([nd, nt]).tolist()) == n.tolist()


def test_theta_prime_sum_matches_interval_library():
    X = 10 ** 5
    ref = arithmetic_moments(T=math.log(X), degree=4, max_prime_limit=X).to_dict()
    mine = R.prime_sum_moments(X, 4, "theta")
    for a, b in zip(mine["c"], ref["coefficients_decimal"]):
        assert a == pytest.approx(float(b), abs=1e-12)


def test_psi_prime_sum_matches_pr2_artifact(pr2):
    for row in pr2["part_c_continuous_R4_theta_vs_psi"][:5]:
        mine = R.prime_sum_moments(row["psi"]["X"], 4, "psi")
        assert np.allclose(mine["c"], row["psi"]["c"], atol=1e-12)


def test_delta_matches_pr2_at_campaign_endpoint(pr2):
    for row in pr2["part_a_exact_corrections"]:
        if float(row["X"]) in (5e10, 1e12):
            d = R.prime_sum_moments(float(row["X"]), 4, "psi_minus_theta")
            for k in range(1, 5):
                assert d["c"][k] == pytest.approx(float(row["delta_k"][str(k)]), abs=1e-13)


def test_psi_minus_theta_identity():
    X = 3 * 10 ** 4
    p, t, d = (R.prime_sum_moments(X, 4, kind) for kind in ("psi", "theta", "psi_minus_theta"))
    assert np.allclose(np.array(p["c"]) - np.array(t["c"]), d["c"], atol=1e-13)


def test_prime_sum_matches_quadrature():
    X = 2000
    n, w = R.prime_power_events(X, "psi")
    T = math.log(X)
    knots = [0.0] + [math.log(v) for v in n] + [T]
    levels = np.concatenate([[0.0], np.cumsum(w)])  # psi on each interval
    mine = R.prime_sum_moments(X, 3, "psi")
    for k in range(4):
        def f(y, s):
            return math.exp(-y / 2) * (s - math.exp(y)) * math.sqrt(2 * k + 1) * eval_legendre(k, 2 * y / T - 1)
        tot = sum(quad(f, a, b, args=(s,))[0] for a, b, s in zip(knots[:-1], knots[1:], levels) if b > a)
        assert mine["Tc"][k] == pytest.approx(tot, abs=1e-9)


def test_ramp_statistic_reproduces_published_values():
    """Rubinstein's exposition, section I.1: G(1e4)=1.749, G(1e5)=3.671, G(1e6)=11.040, G(2e6)=-12.268."""
    G = R.ramp_statistic([1e4, 1e5, 1e6, 2e6])
    assert np.allclose(G, [1.749, 3.671, 11.040, -12.268], atol=6e-4)


def test_ramp_constant_and_explicit_formula():
    assert R.RH_A2 == pytest.approx(0.0461914179322, abs=1e-12)
    gam = R.zeta_zero_ordinates()
    partial = float(np.sum(2.0 / (0.25 + gam ** 2)))
    assert 0.03 < partial < R.RH_A2
    for x in (1e3, 1e4, 1e5, 1e6, 2e6):  # RH-conditional envelope on the zero tail
        diff = float(R.ramp_statistic([x])[0]) - R.ramp_explicit_formula(x)
        assert abs(diff) <= R.ramp_zero_tail_envelope(x)


def test_ramp_identity_links_to_psi_residual():
    """e^{-y} G(e^y) - (1 + gamma) = int_0^y F(u) e^{-u/2} du, checked by quadrature on prime intervals."""
    X = 5000
    n, w = R.prime_power_events(X, "psi")
    knots = [1.0] + [float(v) for v in n] + [float(X)]
    psi = np.concatenate([[0.0], np.cumsum(w)])
    integral = sum(quad(lambda t, s=s: (s - t) / t ** 2, a, b)[0] for a, b, s in zip(knots[:-1], knots[1:], psi))
    G = float(R.ramp_statistic([X])[0])
    assert G / X - (1 + R.EULER_GAMMA) == pytest.approx(integral, abs=1e-10)


def test_ramp_windows_within_rh_envelope():
    for Y in (1e2, 1e3, 1e4):
        wd = R.ramp_window(Y, 1.5)
        assert wd["within_rh_envelope"], wd
        assert wd["points"] > 10


def test_mellin_constants():
    mu = R.psi_mellin_constants()
    with mp.workdps(30):
        half = mp.mpf(1) / 2
        assert mu[0] == pytest.approx(float(2 - 2 * mp.zeta(half, 1, 1) / mp.zeta(half)), abs=1e-13)
        tay = mp.taylor(R._Phi, half, 4)
        for j in range(5):
            assert mu[j] == pytest.approx(float((-1) ** j * mp.factorial(j) * tay[j]), rel=1e-9)
    # the s = 0 pole of Phi (the -log(2pi) e^{-y/2} term) dominates the high moments
    assert mu[8] / (-2 * R.LOG_2PI * 2 ** 8 * math.factorial(8)) == pytest.approx(1.0, abs=1e-5)


def test_mellin_constants_against_damped_prime_sums():
    """Unconditional Mellin identity behind mu_j, at eps = 2 (s = 5/2):
    int_0^inf y^j F(y) e^{-2y} dy = (-1)^j Phi^{(j)}(5/2).  The prime sum stops at 1e5; the part beyond is
    replaced by its density-one value, which leaves an error below 1e-8."""
    n, w = R.prime_power_events(10 ** 5, "psi")
    y = np.log(n.astype(float))
    a = 2.5
    for j in range(3):
        z = a * y
        upper = math.factorial(j) * np.exp(-z) * sum(z ** i / math.factorial(i) for i in range(j + 1))
        lhs = float(np.sum(w * upper)) / a ** (j + 1) - math.factorial(j) / 1.5 ** (j + 1)
        with mp.workdps(30):
            T = math.log(10 ** 5)
            lhs += float(mp.quad(lambda u: mp.exp(u) * mp.gammainc(j + 1, a * u), [T, mp.inf]) / a ** (j + 1))
            rhs = float((-1) ** j * mp.diff(R._Phi, mp.mpf(5) / 2, j))
        assert lhs == pytest.approx(rhs, abs=1e-8)


@pytest.mark.parametrize("X", [1e4, 1e5, 1e6, 1e7])
def test_rh_expansion_against_exact_moments(X):
    """RH-conditional: |T c_k - S_k| <= envelope, and the first 50 zeros leave less than the tail bound."""
    m = R.prime_sum_moments(X, 4, "psi")
    T = m["T"]
    for k in range(5):
        s = R.psi_smooth_part(T, k)
        res = m["Tc"][k] - s
        assert abs(res) <= R.rh_oscillation_envelope(T, k)
        assert abs(res - R.zero_sum_oscillation(T, k)) <= R.zero_tail_envelope(T, k)


def test_limits_and_t2rd_constant():
    mu0 = R.psi_mellin_constants()[0]
    for k in range(5):
        assert R.psi_smooth_part(1e5, k) == pytest.approx(R.psi_limit(k), rel=2e-3 * (k + 1) ** 2)
    assert R.t2_rd_limit_psi(4) == pytest.approx(21 * mu0 ** 2 / R.RH_A2)
    assert 5100 < R.t2_rd_limit_psi(4) < 5250


def test_predicted_r4_matches_measured(pr2):
    for row in pr2["part_c_continuous_R4_theta_vs_psi"][2:6]:
        X = row["psi"]["X"]
        pred = R.predict_moments(X, 4)
        T = pred["T"]
        for kind, key in (("psi", "Tc_psi_pred"), ("theta", "Tc_theta_pred")):
            num = sum(r[key] ** 2 for r in pred["rows"][2:5]) / T ** 2
            assert num / row[kind]["N_aff"] == pytest.approx(row[kind]["R4"], rel=0.02)
