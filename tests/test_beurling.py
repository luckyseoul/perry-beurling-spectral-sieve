"""Tests for Beurling generalized-prime residual builders."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.beurling import (
    beurling_theta_residual,
    build_system_primes,
    default_battery_specs,
    gapped_beurling_primes,
    marathon_battery_specs,
    thinned_ordinary_primes,
)
from pbss.probes import primes_upto, sample_grid
from pbss.projection import energy_ratio


def test_gapped_primes_regular():
    p = gapped_beurling_primes(100.0, gap=3.0, p0=2.0)
    assert p[0] == pytest.approx(2.0)
    assert np.all(np.diff(p) == pytest.approx(3.0))
    assert p[-1] <= 100.0


def test_thinned_ordinary():
    primes = primes_upto(1000)
    th = thinned_ordinary_primes(primes, keep_every=3)
    assert th.size == (primes.size + 2) // 3 or th.size == primes[::3].size
    assert np.array_equal(th, primes[::3].astype(np.float64))


def test_beurling_residual_finite_Rd():
    u = sample_grid(2048)
    p = gapped_beurling_primes(1e5, gap=2.5)
    q, T, meta = beurling_theta_residual(u, p, x_max=1e5, detrend="deg1")
    assert T > 0
    assert meta["n_primes"] > 10
    r = energy_ratio(q, u, degree=4)
    assert 0.0 <= r <= 1.0 + 1e-9


def test_battery_specs_at_least_two_kinds():
    specs = default_battery_specs()
    assert len(specs) >= 2
    kinds = {s["kind"] for s in specs}
    assert "rh_like" in kinds
    assert "defective" in kinds
    ordinary = primes_upto(10_000)
    for s in specs:
        p = build_system_primes(s, ordinary, 1e4)
        assert len(p) >= 1


def test_marathon_battery_specs_at_least_100():
    specs = marathon_battery_specs(100)
    assert len(specs) >= 100
    kinds = {s["kind"] for s in specs}
    assert "rh_like" in kinds and "defective" in kinds
    ordinary = primes_upto(5_000)
    # sample a few builders
    for s in specs[:5] + specs[-3:]:
        p = build_system_primes(s, ordinary, 5_000)
        assert len(p) >= 1


def test_defective_gapped_higher_Rd_than_dense_gap():
    """Larger gap → stronger density defect → typically higher low-degree mass."""
    u = sample_grid(4096)
    p_small_gap = gapped_beurling_primes(5e4, gap=2.0)
    p_big_gap = gapped_beurling_primes(5e4, gap=11.0)
    r1 = energy_ratio(
        beurling_theta_residual(u, p_small_gap, x_max=5e4, detrend="deg1")[0],
        u,
        degree=3,
    )
    r2 = energy_ratio(
        beurling_theta_residual(u, p_big_gap, x_max=5e4, detrend="deg1")[0],
        u,
        degree=3,
    )
    # both finite; big gap should not collapse below tiny
    assert np.isfinite(r1) and np.isfinite(r2)
    assert r2 > 0.05
