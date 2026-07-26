"""Tests for GPU-assisted residual path (CuPy when available, NumPy fallback)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.gpu_residual import arithmetic_residual_fast, energy_ratios_multi_degree
from pbss.probes import arithmetic_residual, primes_upto, sample_grid
from pbss.projection_backend import cupy_available


def test_arithmetic_residual_fast_matches_cpu_on_small():
    u = sample_grid(2048)
    primes = primes_upto(100_000)
    csum = np.cumsum(np.log(primes.astype(np.float64)))
    q_fast, T1, meta = arithmetic_residual_fast(
        u, T=float(np.log(1e5)), primes=primes, csum=csum, detrend="deg1", prefer_gpu=True
    )
    q_cpu, T2, _ = arithmetic_residual(
        u, T=float(np.log(1e5)), primes=primes, csum=csum, detrend="deg1"
    )
    assert T1 == pytest.approx(T2)
    assert np.allclose(q_fast, q_cpu, rtol=1e-10, atol=1e-9)
    assert meta["kind"] == "arithmetic_residual_fast"
    assert "not an rh" in meta["note"].lower() or "not an rh proof" in meta["note"].lower()


def test_energy_ratios_multi_degree_in_unit_interval():
    u = sample_grid(4096)
    primes = primes_upto(50_000)
    csum = np.cumsum(np.log(primes.astype(np.float64)))
    q, _, _ = arithmetic_residual_fast(
        u, T=10.0, primes=primes, csum=csum, detrend="deg1", prefer_gpu=True
    )
    out = energy_ratios_multi_degree(q, u, [2, 4], prefer_gpu=True)
    assert set(out.keys()) == {2, 4}
    for d, info in out.items():
        assert 0.0 <= info["R_d"] <= 1.0 + 1e-9
        assert info["backend"] in ("numpy", "cupy")


def test_cupy_available_bool():
    assert isinstance(cupy_available(), bool)
