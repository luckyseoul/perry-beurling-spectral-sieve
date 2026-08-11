"""Tests for consensus ranks 2–5 (plateau, ANT audit, zero-proportion, B-RES threshold)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.ant_audit import ant_interface_audit  # noqa: E402
from pbss.b_res_threshold import (  # noqa: E402
    b_res_threshold_report,
    cancelled_off_critical_rd,
    off_critical_rd_lower_model,
)
from pbss.plateau_secondary import (  # noqa: E402
    enrich_scan_small,
    evaluate_predictions,
    plateau_secondary_report,
)
from pbss.zero_proportion_feasibility import (  # noqa: E402
    zero_proportion_feasibility_report,
)


def test_rank2_plateau_predictions_on_shipped_ef_path():
    rep = plateau_secondary_report(T=14.0, x_max=2e6, n_points=1024, degree=4)
    assert rep["rh_claimed"] is False
    assert "NOT AN UNCONDITIONAL" in rep["banner"]
    pred = rep["prediction_eval"]
    assert "zeros" in pred["means"]
    assert pred["predictions"]["P1_smooth_beats_zeros"]["pass"]
    assert pred["predictions"]["P3_Vd_is_oracle"]["pass"]
    assert pred["all_tested_pass"]


def test_rank2_enrich_scan_uses_real_ed_metric():
    rows = enrich_scan_small(
        T=14.0, x_max=1e5, n_points=512, n_zeros_list=(5, 10), enrichments=("zeros", "zeros_Vd")
    )
    assert len(rows) == 4
    for r in rows:
        assert "Ed_r_over_l2q" in r
        assert 0.0 <= r["Ed_r_over_l2q"] <= 1.0 + 1e-9


def test_rank3_ant_audit_freezes_full_a():
    aud = ant_interface_audit()
    assert aud["rh_claimed"] is False
    assert aud["unlabeled_count"] == 0
    assert aud["freeze_full_a_packaging"] is True
    assert aud["full_a_status"] == "closed_conditional"
    ids = {r["id"] for r in aud["checklist"]}
    assert {"ANT-1", "ANT-2", "ANT-3", "M7"}.issubset(ids)


def test_rank4_zero_proportion_stops():
    rep = zero_proportion_feasibility_report()
    assert rep["decision"] == "STOP"
    assert rep["rh_claimed"] is False
    assert len(rep["comparison"]) >= 3
    assert all(c["status"] != "ready_to_implement" for c in rep["incremental_candidates"])


def test_rank5_b_res_threshold_not_solved():
    rep = b_res_threshold_report(T=16.0, sigma=0.9, degree=4)
    assert rep["b_res_solved"] is False
    assert rep["rh_claimed"] is False
    assert rep["pure_above_cancel"] is True
    assert rep["model_off_critical"]["R_d_off"] > rep["model_cancellation_counterexample"][
        "R_d_after_Vd_kill"
    ]
    assert "H*" in rep["threshold_hypothesis"]["statement"] or "H*" in rep[
        "threshold_hypothesis"
    ]["id"]


def test_rank5_cancellation_near_zero_rd():
    c = cancelled_off_critical_rd(18.0, degree=4, n_points=2048)
    assert c["R_d_after_Vd_kill"] < 1e-6
    o = off_critical_rd_lower_model(18.0, degree=4, n_points=2048)
    assert o["R_d_off"] > 0.01
