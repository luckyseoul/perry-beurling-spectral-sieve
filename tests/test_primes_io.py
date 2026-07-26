"""Tests for prime checkpoint I/O and GPU-optional projection backend."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.primes_io import ensure_primes, load_primes_checkpoint, save_primes_checkpoint
from pbss.probes import primes_upto, sample_grid
from pbss.projection_backend import cupy_available, energy_ratio_auto


def test_ensure_primes_checkpoint_roundtrip(tmp_path):
    x_max = 5000
    primes, meta = ensure_primes(x_max, tmp_path, force_resieve=True)
    assert meta["n_primes"] == len(primes)
    assert primes[0] == 2
    # second call loads checkpoint
    primes2, meta2 = ensure_primes(x_max, tmp_path, force_resieve=False)
    assert meta2.get("loaded_from_checkpoint") is True
    assert np.array_equal(np.asarray(primes), np.asarray(primes2))


def test_energy_ratio_auto_numpy_path():
    u = sample_grid(512)
    q = np.sin(20 * np.pi * u)
    r, backend = energy_ratio_auto(q, u, degree=3, prefer_gpu=False)
    assert backend == "numpy"
    assert 0.0 <= r <= 1.0


def test_cupy_available_is_bool():
    assert isinstance(cupy_available(), bool)
