"""Tests for Gamma-weight sensitivity + ≥53% claim confirmation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.measure_sensitivity import (  # noqa: E402
    confirm_sensitivity_claim,
    fisher_dprime,
    gamma_weight,
    sensitivity_experiment,
    weighted_energy,
)
from pbss.probes import sample_grid  # noqa: E402


def test_gamma_weight_peaks_at_half_for_k4_sigma6():
    u = sample_grid(2048)
    w = gamma_weight(u, k=4.0, sigma=6.0, normalize=True)
    assert abs(float(u[np.argmax(w)]) - 0.5) < 0.02
    assert float(np.sum(w) * (u[1] - u[0])) == pytest.approx(1.0, rel=1e-3)


def test_fisher_dprime_separated_vs_identical():
    rng = np.random.default_rng(0)
    a = rng.normal(0.0, 1.0, size=500)
    b = rng.normal(3.0, 1.0, size=500)
    assert fisher_dprime(a, b) > 2.0
    assert fisher_dprime(a, a) < 0.2


def test_sensitivity_gamma_beats_flat():
    rep = sensitivity_experiment(
        n_per_class=120, noise=0.4, seed=20260522, k=4.0, sigma=6.0
    )
    assert rep["gamma_dprime"] > rep["flat_dprime"]
    assert rep["relative_gain"] > 0.0
    assert rep["rh_claimed"] is False


def test_confirm_53_percent_lower_bound_achieved():
    """Project-record claim: ≥53% relative gain on a noisy ensemble."""
    rep = confirm_sensitivity_claim(min_gain=0.53, seed=20260522)
    assert rep["noisy_meets_53pct"] is True
    assert rep["noisy"]["relative_gain"] >= 0.53
    assert rep["clean_improves"] is True
    assert rep["confirmed"] is True
    assert rep["verdict"] == "CONFIRMED"


def test_weighted_energy_positive():
    u = sample_grid(512)
    q = np.sin(40 * np.pi * u)
    w = gamma_weight(u, k=4, sigma=6)
    assert weighted_energy(q, u, w) > 0.0
