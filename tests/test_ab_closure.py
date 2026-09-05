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


def test_package_status_distinguishes_conditional_equivalence_and_unconditional_claims():
    s = package_status()
    assert s["full_arithmetic_A"] == "proved_conditional_continuous"
    assert s["full_B"] == "open_normalized_converse"
    assert s["full_B_residual_step_id"] == B_RESIDUAL_STEP_ID
    assert s["rh"] == "open"
    assert s["model_A0"] == "proved"
    assert s["model_B0"] == "proved"
    assert s["rh_claimed"] is False
    assert s["absolute_moment_converse"] == "proved_rh_equivalence"
    assert s["absolute_moment_growth_bound"] == "open_unconditional"
    assert s["normalized_recurrence"] == "proved_unconditional_continuous"
    assert s["normalized_full_limit"] == "open_unconditional"
    assert s["normalized_rate"] == "proved_only_under_rh"
    assert s["normalized_b_res_floor"] == "rh_equivalent_open"
    scope = s["continuous_scope"]
    assert scope["degree"] == "fixed integer d>=2"
    assert scope["detrend"] == "continuous orthogonal affine projection P_1"
    assert scope["sampled_transfer_proved"] is False
    assert scope["weighted_transfer_proved"] is False
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in s["banner"]
    assert "RH remains open" in s["note"] or s["rh"] == "open"


def test_chain_package_status_matches_ab_closure():
    s = chain_status()
    assert s["full_arithmetic_A"] == "proved_conditional_continuous"
    assert s["rh"] == "open"
    assert s["conditional_theorem_a_package"] == "complete"
    assert s["full_B"] == "open_normalized_converse"


def test_full_a_gap_table_separates_proved_continuous_result_from_open_transfers():
    rows = full_a_gap_table()
    assert any(r["step"].startswith("Full A") for r in rows)
    required = [r for r in rows if r["required_for_continuous_A"]]
    assert required
    assert all(r["disposition"] in {
        "cited_exact_statement", "proved_continuous_adaptation", "proved_conditional_continuous"
    } for r in required)
    transfers = [r for r in rows if r["disposition"] == "open_transfer"]
    assert len(transfers) == 2
    assert all(not r["required_for_continuous_A"] for r in transfers)


def test_full_b_does_not_equate_absolute_converse_with_normalized_converse():
    rows = full_b_gap_table()
    residual = [r for r in rows if "B-RES" in r["step"]]
    assert len(residual) == 1
    assert residual[0]["disposition"] == "rh_equivalent_open"
    assert next(r for r in rows if "B-ABS" in r["step"])["disposition"] == "proved_rh_equivalence"
    assert next(r for r in rows if "Full B:" in r["step"])["disposition"] == "open_normalized_converse"
    assert next(r for r in rows if "recurrence" in r["step"])["disposition"] == "proved_unconditional_continuous"


def test_ant_citations_have_hypotheses_and_refs():
    cites = ant_citations()
    ids = {c["id"] for c in cites}
    assert {"ANT-1", "ANT-2", "ANT-3", "ANT-MS"}.issubset(ids)
    for c in cites:
        assert len(c["classical_refs"]) >= 1
        assert len(c["hypotheses"]) >= 1
        assert "adapted_conclusion" in c
        assert (ROOT / c["proof"].split(",")[0]).exists()
    external = [c for c in cites if c["status"] == "cited_exact_statement"]
    assert {c["id"] for c in external} == {"ANT-3", "ANT-MS"}
    assert all(c["source_url"].endswith("Akbary-Ng-Shahabi.pdf") for c in external)
    assert all(c["not_proved_in_repo"] is True for c in external)
    assert next(c for c in cites if c["id"] == "ANT-4")["adaptation_proved_in_repo"] is False


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


def test_off_critical_comparison_is_only_sampled_model_evidence():
    a = off_critical_model_obstruction(8.0, sigma=0.9, n_points=2048)
    b = off_critical_model_obstruction(32.0, sigma=0.9, n_points=2048)
    assert b["ratio_off_over_cl"] > a["ratio_off_over_cl"]
    assert a["R_d_off"] > a["R_d_cl"]
    assert a["evidence_kind"] == "sampled_model_comparison"
    assert a["uniform_lower_bound_proved"] is False


def test_off_critical_comparison_large_T_uses_the_common_exponential_scale():
    from pbss.b_res_threshold import off_critical_rd_lower_model

    report = off_critical_model_obstruction(4000.0, degree=4, n_points=2048)
    other = off_critical_rd_lower_model(4000.0, degree=4, n_points=2048)
    assert 0.0 <= report["R_d_off"] <= 1.0
    assert report["R_d_off"] == pytest.approx(other["R_d_off"], abs=1e-12)
    assert report["uniform_lower_bound_proved"] is False


@pytest.mark.parametrize("T", [float("nan"), float("inf"), 0.0, -1.0])
def test_off_critical_comparison_rejects_nonfinite_or_nonpositive_T(T):
    with pytest.raises(ValueError, match="T"):
        off_critical_model_obstruction(T)


def test_conditional_full_a_report_and_model_chain_labels():
    rep = conditional_full_a_report(18.0, degree=4, n_zeros=6, n_points=1024)
    assert rep["full_a_status"] == "proved_conditional_continuous"
    assert rep["full_b_status"] == "open_normalized_converse"
    assert rep["rh_status"] == "open"
    assert rep["m7_grid_check"]["holds"] == 1.0
    assert rep["evidence_kind"] == "sampled_model_report"
    assert rep["continuous_limit_verified_by_this_run"] is False
    assert rep["rh_verified_by_this_run"] is False
    row = rep["model_chain"]
    assert row["full_arithmetic_A_status"] == "proved_conditional_continuous"
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in row["banner"]
    assert row["continuous_limit_verified_by_this_run"] is False
    # model CL still decays vs majorants field present
    assert 0.0 <= row["empirical"]["R_d_cl"] <= 1.0


def test_model_chain_legacy_fields_updated():
    row = model_chain_report(20.0, degree=4, n_zeros=5, n_points=1024)
    assert row["full_arithmetic_A_status"] == "proved_conditional_continuous"
    assert row["rh_status"] == "open"
