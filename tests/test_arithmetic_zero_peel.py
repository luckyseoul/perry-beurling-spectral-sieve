"""Tests for arithmetic residual zero-peel (shipped path)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.probes import (
    arithmetic_residual,
    arithmetic_zero_peel,
    peel_residual,
    primes_upto,
    sample_grid,
)
from pbss.projection import energy_ratio


def test_arithmetic_zero_peel_n0_matches_raw():
    u = sample_grid(2048)
    primes = primes_upto(50_000)
    q0, T0, m0 = arithmetic_zero_peel(
        u, x_max=5e4, primes=primes, n_strip=0, detrend="deg1"
    )
    q1, T1, _ = arithmetic_residual(u, x_max=5e4, primes=primes, detrend="deg1")
    assert T0 == pytest.approx(T1)
    assert np.allclose(q0, q1)
    assert m0["n_strip"] == 0
    assert "not an rh" in m0["note"].lower() or "not an unconditional" in m0["note"].lower()


def test_arithmetic_zero_peel_strips_and_finite_Rd():
    u = sample_grid(4096)
    primes = primes_upto(200_000)
    q, T, meta = arithmetic_zero_peel(
        u,
        x_max=2e5,
        primes=primes,
        n_strip=5,
        detrend="deg1",
        fit_scale=True,
    )
    assert meta["n_strip"] == 5
    assert np.all(np.isfinite(q))
    r = energy_ratio(q, u, degree=4)
    assert 0.0 <= r <= 1.0 + 1e-9
    assert np.isfinite(meta["mode_scale"])


def test_arithmetic_zero_peel_uses_peel_residual_semantics():
    """Stripping modes from raw arithmetic equals peel_residual helper path."""
    u = sample_grid(2048)
    primes = primes_upto(80_000)
    q_raw, T, _ = arithmetic_residual(u, x_max=8e4, primes=primes, detrend="none")
    q_peel, pmeta = peel_residual(q_raw, u, T, n_strip=3, mode_scale=1.0)
    q2, _, meta = arithmetic_zero_peel(
        u,
        T=T,
        primes=primes,
        n_strip=3,
        detrend="none",
        fit_scale=False,
        mode_scale=1.0,
    )
    assert np.allclose(q_peel, q2, rtol=1e-10, atol=1e-10)
    assert meta["mode_scale"] == pytest.approx(1.0)
