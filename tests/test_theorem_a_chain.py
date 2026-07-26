"""Tests for Conditional Theorem A model chain (shipped theorem_a_chain API)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.theorem_a_chain import (  # noqa: E402
    model_chain_report,
    multi_T_model_chain,
    package_status,
)


def test_package_status_does_not_claim_rh_or_full_a():
    s = package_status()
    assert s["conditional_theorem_a_package"] == "complete"
    assert s["full_arithmetic_A"] == "open"
    assert s["rh"] == "open"
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in s["banner"]


def test_model_chain_report_fields_and_labels():
    row = model_chain_report(20.0, degree=4, n_zeros=8, n_points=2048)
    assert row["T"] == 20.0
    assert row["full_arithmetic_A_status"] == "open"
    assert row["rh_status"] == "open"
    emp = row["empirical"]
    assert 0.0 <= emp["R_d_cl"] <= 1.0 + 1e-9
    assert 0.0 <= emp["R_d_ef"] <= 1.0 + 1e-9
    maj = row["majorants"]
    assert "proved_style" in maj["M5_ef_flat"]["label"]
    assert "scaffolding" in maj["infinite_tail_scaffold"]["label"]
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in row["banner"]


def test_model_chain_multi_T_decay_cl():
    """Pure CL empirical R_d decreases with T; majorants shrink."""
    rows = multi_T_model_chain(
        [15.0, 30.0, 60.0], degree=4, n_zeros=6, n_points=4096
    )
    assert len(rows) == 3
    r_cl = [r["empirical"]["R_d_cl"] for r in rows]
    assert r_cl[2] < r_cl[0]
    m5 = [r["majorants"]["M5_ef_flat"]["value"] for r in rows]
    assert m5[2] < m5[0]
    # weighted CL also decays
    r_w = [r["empirical"]["R_d_cl_weighted"] for r in rows]
    assert r_w[2] < r_w[0]
    assert all(r["proved_model_decay_ok"] for r in rows)


def test_model_chain_rejects_nonpositive_T():
    with pytest.raises(ValueError):
        model_chain_report(0.0)
