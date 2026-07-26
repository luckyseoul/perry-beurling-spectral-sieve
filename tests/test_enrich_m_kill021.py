"""Tests: enrich m only, multi-N Ed(r) on H_theta_sqrt (shipped APIs)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.ef_identify import (  # noqa: E402
    M_ENRICHMENTS,
    build_m_columns,
    identify_ef,
    multi_N_enrich_scan,
    summarize_enrich_kill021,
)
from pbss.probes import explicit_formula_residual, primes_upto, sample_grid  # noqa: E402


def test_build_m_columns_zeros_and_enrich():
    u = sample_grid(256)
    cols, names, meta = build_m_columns(u, T=20.0, n_zeros=5, enrich="zeros")
    assert len(cols) == 1 and names == ["zeros"]
    cols2, names2, _ = build_m_columns(
        u, T=20.0, n_zeros=5, enrich="zeros_Vd", degree=4
    )
    assert len(cols2) == 1 + 5  # zeros + φ0..φ4
    assert "phi0" in names2


def test_model_zeros_identity_still_holds():
    u = sample_grid(1024)
    q, _, _ = explicit_formula_residual(u, T=25.0, n_zeros=6)
    idn = identify_ef(q, u, T=25.0, n_zeros=6, m_enrich="zeros")
    assert idn["E_d_remainder_over_l2q"] < 1e-10
    assert idn["frac_l2_remainder"] < 1e-10


def test_zeros_Vd_kills_Ed_on_arithmetic_small():
    """Spanning V_d inside m must drive Ed rem ~0 (construction)."""
    primes = primes_upto(50_000)
    rows = multi_N_enrich_scan(
        T=10.0,
        n_zeros_list=[5, 10],
        primes=primes,
        enrichments=("zeros", "zeros_Vd"),
        n_points=512,
        degree=4,
        detrend="deg1",
    )
    by = {}
    for r in rows:
        by.setdefault(r["m_enrich"], []).append(r["Ed_r_over_l2q"])
    assert np.mean(by["zeros"]) > 0.05
    assert np.mean(by["zeros_Vd"]) < 1e-6


def test_summarize_detects_structure():
    primes = primes_upto(40_000)
    rows = multi_N_enrich_scan(
        T=10.0,
        n_zeros_list=[5, 10, 15],
        primes=primes,
        enrichments=("zeros", "zeros_smooth", "zeros_Vd"),
        n_points=512,
        degree=4,
    )
    s = summarize_enrich_kill021(rows)
    assert "means_Ed_r_over_l2q" in s
    assert "zeros" in s["means_Ed_r_over_l2q"]
    assert s["n_rows"] == len(rows)
    # Vd should be near zero mean
    assert s["means_Ed_r_over_l2q"]["zeros_Vd"] < 0.05
