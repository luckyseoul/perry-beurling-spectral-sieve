"""Tests for Full A/B package closure surface (real shipped ab_closure API)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.ab_closure import (  # noqa: E402
    B_RESIDUAL_STEP_ID,
    ant_citations,
    conditional_full_a_report,
    energy_ratio_perturbation_bound,
    full_a_gap_table,
    full_b_gap_table,
    off_critical_model_obstruction,
    package_status,
    verify_m7_on_grid,
)
from pbss.theorem_a_chain import model_chain_report, package_status as chain_status


def test_package_status_full_a_closed_conditional_not_rh():
    s = package_status()
    assert s["full_arithmetic_A"] == "closed_conditional"
    assert s["full_B"] == "package_complete_single_residual"
    assert s["full_B_residual_step_id"] == B_RESIDUAL_STEP_ID
    assert s["rh"] == "open"
    assert s["model_A0"] == "proved"
    assert s["model_B0"] == "proved"
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in s["banner"]
    assert "RH remains open" in s["note"] or s["rh"] == "open"


def test_chain_package_status_matches_ab_closure():
    s = chain_status()
    assert s["full_arithmetic_A"] == "closed_conditional"
    assert s["rh"] == "open"
    assert s["conditional_theorem_a_package"] == "complete"
    assert s["full_B"] == "package_complete_single_residual"


def test_full_a_gap_table_no_open_required_steps():
    rows = full_a_gap_table()
    assert any(r["step"].startswith("Full A") for r in rows)
    # Required ANT/model steps must not be bare "open" or "scaffolding"
    for r in rows:
        if r["step"] == "Unconditional RH":
            assert r["disposition"] == "open"
            continue
        if "Scaffold" in r["step"] or "scaffold" in r["step"].lower():
            continue
        assert r["disposition"] not in ("open", "scaffolding", "scaffolding only")


def test_full_b_exactly_one_open_residual():
    rows = full_b_gap_table()
    openish = [
        r
        for r in rows
        if r["disposition"] in ("open_single_residual", "open")
        and "RH" not in r["step"]
        and "via Full B" not in r["step"]
    ]
    # Exactly one residual gap named B-RES
    assert any("B-RES" in r["step"] for r in rows)
    residual = [r for r in rows if "B-RES" in r["step"]]
    assert len(residual) == 1
    assert residual[0]["disposition"] == "open_single_residual"


def test_ant_citations_have_hypotheses_and_refs():
    cites = ant_citations()
    ids = {c["id"] for c in cites}
    assert {"ANT-1", "ANT-2", "ANT-3"}.issubset(ids)
    for c in cites:
        assert c["status"] in ("cited", "cited_optional")
        assert len(c["classical_refs"]) >= 1
        assert len(c["hypotheses"]) >= 1
        assert "adapted_conclusion" in c
        assert c["not_proved_in_repo"] is True


def test_m7_majorant_dominates_empirical_on_real_path():
    out = verify_m7_on_grid(24.0, degree=4, n_points=2048, defect_weight=0.12)
    assert out["holds"] == 1.0
    assert out["R_empirical"] <= out["majorant_M7"] + 1e-9
    assert 0.0 <= out["R0"] <= 1.0


def test_m7_formula_limit_small_pieces():
    # R0=0, delta=0 → majorant 0
    assert energy_ratio_perturbation_bound(0.0, q0_norm=1.0, r_norm=0.0) == 0.0
    # small delta, small R0 → small majorant
    maj = energy_ratio_perturbation_bound(0.01, q0_norm=1.0, r_norm=0.05)
    assert maj < 0.1
    with pytest.raises(ValueError):
        energy_ratio_perturbation_bound(0.5, q0_norm=0.0, r_norm=0.1)


def test_off_critical_ratio_grows_with_T():
    a = off_critical_model_obstruction(8.0, sigma=0.9, n_points=2048)
    b = off_critical_model_obstruction(32.0, sigma=0.9, n_points=2048)
    assert b["ratio_off_over_cl"] > a["ratio_off_over_cl"]
    assert a["R_d_off"] > a["R_d_cl"]


def test_conditional_full_a_report_and_model_chain_labels():
    rep = conditional_full_a_report(18.0, degree=4, n_zeros=6, n_points=1024)
    assert rep["full_a_status"] == "closed_conditional"
    assert rep["full_b_status"] == "package_complete_single_residual"
    assert rep["rh_status"] == "open"
    assert rep["m7_grid_check"]["holds"] == 1.0
    row = rep["model_chain"]
    assert row["full_arithmetic_A_status"] == "closed_conditional"
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in row["banner"]
    # model CL still decays vs majorants field present
    assert 0.0 <= row["empirical"]["R_d_cl"] <= 1.0


def test_model_chain_legacy_fields_updated():
    row = model_chain_report(20.0, degree=4, n_zeros=5, n_points=1024)
    assert row["full_arithmetic_A_status"] == "closed_conditional"
    assert row["rh_status"] == "open"
