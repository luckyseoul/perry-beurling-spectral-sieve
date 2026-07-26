"""
Tests for explicit-formula residual builder and peel helpers.

Drive shipped ``explicit_formula_residual`` / ``peel_residual`` + ``energy_ratio``
only — no reimplementation of the residual formula inside the tests.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.lemmas import bound_R_d_finite_mode_sum
from pbss.probes import (
    explicit_formula_residual,
    finite_cl_superposition,
    peel_residual,
    probe_critical_line_mode,
    sample_grid,
)
from pbss.projection import energy_ratio, project
from pbss.zeros import (
    ZETA_ZERO_ORDINATES_50,
    explicit_formula_amplitudes,
    zeta_zero_ordinates,
)


def test_zeta_zero_table_first_entries():
    t = zeta_zero_ordinates(3)
    assert t.shape == (3,)
    assert t[0] == pytest.approx(14.134725141734693, rel=0, abs=1e-9)
    assert t[1] == pytest.approx(21.022039638771554, rel=0, abs=1e-9)
    assert len(ZETA_ZERO_ORDINATES_50) == 50


def test_explicit_formula_residual_finite_nonzero():
    u = sample_grid(4096)
    q, T, meta = explicit_formula_residual(u, T=20.0, n_zeros=8)
    assert T == 20.0
    assert meta["n_zeros"] == 8
    assert np.all(np.isfinite(q))
    assert float(np.max(np.abs(q))) > 0.0
    assert "not an RH proof" in meta["note"].lower() or "not an RH" in meta["note"]


def test_explicit_formula_residual_Rd_in_unit_interval():
    u = sample_grid(8192)
    q, T, meta = explicit_formula_residual(u, T=25.0, n_zeros=10)
    r = energy_ratio(q, u, degree=4)
    assert 0.0 < r <= 1.0 + 1e-9
    res = project(q, u, degree=4, T=T)
    assert res.energy_ratio == pytest.approx(r, rel=0, abs=1e-12)
    assert np.isfinite(res.scaled_strength)


def test_explicit_formula_single_zero_matches_scaled_cl_mode():
    """N=1 cos form matches amplitude * cos(t T u - α); energy ratio equals pure mode."""
    u = sample_grid(10000)
    t0 = ZETA_ZERO_ORDINATES_50[0]
    q, _, meta = explicit_formula_residual(u, T=30.0, n_zeros=1, form="cos")
    a0 = meta["amplitudes"][0]
    # pure cos mode (phase from meta) has same R_d as unit sinusoid of same freq
    q_unit = finite_cl_superposition(
        u,
        T=30.0,
        amplitudes=np.array([1.0]),
        ordinates=np.array([t0]),
        phases=np.array([meta["phases"][0]]),
        form="cos",
    )
    r = energy_ratio(q, u, degree=3)
    r_unit = energy_ratio(q_unit, u, degree=3)
    assert r == pytest.approx(r_unit, rel=1e-6, abs=1e-8)
    assert abs(a0) > 0


def test_explicit_formula_Rd_decays_with_T():
    u = sample_grid(10000)
    d = 4
    r_lo = energy_ratio(
        explicit_formula_residual(u, T=10.0, n_zeros=6)[0], u, degree=d
    )
    r_hi = energy_ratio(
        explicit_formula_residual(u, T=50.0, n_zeros=6)[0], u, degree=d
    )
    assert r_hi < r_lo
    assert r_hi < 0.1
    t = zeta_zero_ordinates(6)
    a = explicit_formula_amplitudes(t)
    assert r_hi <= bound_R_d_finite_mode_sum(50.0, a, t, d) + 0.05


def test_peel_residual_strips_modes():
    u = sample_grid(4096)
    T = 22.0
    q_full, _, meta = explicit_formula_residual(u, T=T, n_zeros=12)
    q_peel, pmeta = peel_residual(q_full, u, T, n_strip=5)
    assert pmeta["n_strip"] == 5
    # after stripping first 5 of a 12-mode sum, residual ≈ modes 6..12
    q_tail, _, _ = explicit_formula_residual(
        u,
        T=T,
        n_zeros=7,
        ordinates=np.asarray(meta["ordinates"][5:]),
        amplitudes=np.asarray(meta["amplitudes"][5:]),
        phases=np.asarray(meta["phases"][5:]),
        form="cos",
    )
    assert np.allclose(q_peel, q_tail, rtol=1e-10, atol=1e-10)


def test_peel_n_strip_zero_is_identity():
    u = sample_grid(512)
    q, _, _ = explicit_formula_residual(u, T=15.0, n_zeros=4)
    q2, meta = peel_residual(q, u, 15.0, n_strip=0)
    assert meta["stripped"] is False
    assert np.array_equal(q, q2)


def test_bulk_term_shifts_residual():
    u = sample_grid(2048)
    q0, _, _ = explicit_formula_residual(
        u, T=18.0, n_zeros=3, bulk="none", bulk_scale=0.0
    )
    q_bulk, _, meta = explicit_formula_residual(
        u, T=18.0, n_zeros=3, bulk="deg0", bulk_scale=0.5
    )
    assert meta["bulk"] == "deg0"
    assert np.allclose(q_bulk, q0 + 0.5)
