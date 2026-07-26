"""Tests for admissible weight class and endpoint contribution (shipped APIs)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.probes import (  # noqa: E402
    finite_cl_superposition,
    probe_critical_line_mode,
    sample_grid,
)
from pbss.projection import energy_ratio  # noqa: E402
from pbss.weights import (  # noqa: E402
    admissible_weight,
    apply_weight,
    bulk_vs_weighted_report,
    endpoint_contribution,
    hanning_weight,
    is_admissible_weight,
    tukey_weight,
    weighted_energy_ratio,
)
from pbss.zeros import ZETA_ZERO_ORDINATES_50  # noqa: E402


def test_tukey_vanishes_at_endpoints():
    u = sample_grid(2048)
    w = tukey_weight(u, alpha=0.1)
    assert w[0] == pytest.approx(0.0, abs=1e-12)
    assert w[-1] == pytest.approx(0.0, abs=1e-12)
    # interior near 0.5 should be ~1
    mid = w[len(w) // 2]
    assert mid == pytest.approx(1.0, abs=1e-9)


def test_hanning_is_admissible_mass():
    u = sample_grid(1024)
    w = hanning_weight(u)
    adm = is_admissible_weight(w, u, alpha=0.5)
    assert adm["ok"] is True
    assert adm["l2_mass"] > 0.1


def test_endpoint_contribution_positive_for_bulk_trend():
    """Linear bulk has mass at ends; E_end should be positive."""
    u = sample_grid(4096)
    q = (u - 0.5) + 0.01 * np.sin(40 * np.pi * u)
    stats = endpoint_contribution(q, u, degree=2, alpha=0.1)
    assert stats["E_end"] >= 0.0
    assert stats["R_d"] > 0.0
    assert np.isfinite(stats["E_bulk"])


def test_weight_reduces_endpoint_effect_on_linear_plus_oscillation():
    """
    Admissible weight should lower R_d for a residual dominated by slow bulk
    relative to raw (open-plateau taper phenomenon on a synthetic).
    """
    u = sample_grid(8192)
    # slow bulk + high-frequency
    q = 3.0 * (u - 0.5) + 0.2 * np.sin(80 * np.pi * u)
    r_raw = energy_ratio(q, u, degree=4)
    r_w, _, end = weighted_energy_ratio(
        q, u, degree=4, weight_name="tukey", alpha=0.15
    )
    assert r_raw > 0.5  # bulk-dominated
    assert r_w < r_raw
    assert end["E_end"] > 0.0


def test_bulk_vs_weighted_report_multi_T_cl_modes():
    """CL modes: weighted R_d stays small for several T (M3-compatible)."""
    u = sample_grid(8192)
    t0 = ZETA_ZERO_ORDINATES_50[0]
    rows = []
    for T in (15.0, 25.0, 40.0):
        q = probe_critical_line_mode(u, T=T, t=t0)
        rep = bulk_vs_weighted_report(q, u, degree=4, alpha=0.1, weight_name="hanning")
        rows.append(rep)
        assert 0.0 <= rep["R_d_raw"] <= 1.0 + 1e-9
        assert 0.0 <= rep["R_d_weighted"] <= 1.0 + 1e-9
        assert "NOT AN UNCONDITIONAL PROOF OF RH" in rep["banner"]
    assert len(rows) == 3
    # higher T → smaller raw R_d for pure CL mode
    assert rows[-1]["R_d_raw"] < rows[0]["R_d_raw"]


def test_admissible_weight_names():
    u = sample_grid(512)
    for name in ("tukey", "hanning", "flat"):
        w = admissible_weight(u, name=name, alpha=0.1)
        assert w.shape == u.shape
        assert np.all(np.isfinite(w))
    with pytest.raises(ValueError):
        admissible_weight(u, name="not_a_window")


def test_apply_weight_matches_multiply():
    u = sample_grid(256)
    q = probe_critical_line_mode(u, T=20.0)
    w = tukey_weight(u, 0.2)
    assert np.allclose(apply_weight(q, w), q * w)
