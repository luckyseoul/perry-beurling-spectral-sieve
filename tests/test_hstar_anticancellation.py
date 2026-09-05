"""Regression checks for the retracted L* proxy and retained model identities."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.hstar_anticancellation import (  # noqa: E402
    CANDIDATE_STATUS,
    analytic_off_critical_Rd_leading,
    analytic_vs_numeric_row,
    arithmetic_injection_bridge_row,
    coercivity_Cd,
    hstar_anticancellation_report,
    lemma_statement,
    partial_cancel,
    partial_cancel_curve,
    weil_local_proxy_row,
    write_campaign_artifacts,
    write_multi_T_artifacts,
)
from pbss.probes import probe_off_critical_mode as _probe  # noqa: E402
from pbss.probes import sample_grid as _grid  # noqa: E402
from pbss.zero_proportion_feasibility import (  # noqa: E402
    zero_proportion_feasibility_report,
)


def test_lemma_named_and_coercivity_retracted():
    lem = lemma_statement()
    assert lem["name"] == "L*"
    assert lem["id"]
    assert lem["status"] == "retracted_coercivity_proxy"
    assert "No Guinand-Weil quadratic form" in lem["hypotheses"]
    assert "retracted" in lem["conclusion"]
    assert "proportion" in lem["not_a_proportion"].lower() or "Re=1/2" in lem["not_a_proportion"]
    assert lem["gaps_to_full_H_star"]


def test_hstar_report_status_and_rh_nonclaim():
    rep = hstar_anticancellation_report(
        T_values=(8.0, 16.0, 24.0, 32.0),
        sigma=0.9,
        degree=4,
    )
    assert rep["rh_claimed"] is False
    assert rep["b_res_solved"] is False
    assert "NOT AN UNCONDITIONAL" in rep["banner"]
    assert rep["candidate_status"] == "retracted_coercivity_proxy"
    assert rep["weil_form_evaluated"] is False
    assert rep["positive_limiting_model_coercivity"] is False
    assert rep["candidate_status"] == CANDIDATE_STATUS
    assert rep["axis"] == "anti_cancellation_toward_H_star"
    assert all(c["status"] == "retracted_coercivity_proxy" for c in rep["candidates"])
    assert rep["all_constraint_signals"] is False
    assert rep["all_sample_separation_signals"] is True
    assert len(rep["multi_T_proxy"]) >= 3
    assert rep["rank5_crosscheck"]["pure_above_cancel"] is True
    assert rep["partial_cancel_monotone"] is True


def test_rank4_proportion_still_stop():
    """L* must not flip Rank-4 proportion candidates to ready_to_implement."""
    rep = zero_proportion_feasibility_report()
    assert rep["decision"] == "STOP"
    assert rep["rh_claimed"] is False
    assert all(c["status"] != "ready_to_implement" for c in rep["incremental_candidates"])
    assert "anti-cancellation" in rep["related_but_distinct"].lower() or "H*" in rep["related_but_distinct"]


def test_full_cancel_drives_Cd_to_zero():
    """Shipped projector: α=1 ⇒ C_d ≈ 0 (drives real project/energy_ratio path)."""
    u = _grid(2048)
    q = _probe(u, T=20.0, sigma=0.9, t=14.134725)
    q1 = partial_cancel(q, u, degree=4, T=20.0, alpha=1.0)
    assert coercivity_Cd(q1, u, degree=4, T=20.0) < 1e-5
    assert coercivity_Cd(q, u, degree=4, T=20.0) > 0.05


def test_partial_cancel_monotone_Cd():
    curve = partial_cancel_curve(24.0, sigma=0.85, degree=4, alphas=(0.0, 0.25, 0.5, 0.75, 1.0))
    cds = [r["C_d"] for r in curve]
    for i in range(len(cds) - 1):
        assert cds[i] >= cds[i + 1] - 1e-9
    assert cds[0] > 0.05
    assert cds[-1] < 1e-5


def test_pure_beats_cl_on_row():
    row = weil_local_proxy_row(16.0, sigma=0.9, degree=4)
    assert row["C_d_pure"] > row["C_d_cl"]
    assert row["constraint_signal"] == 0.0
    assert row["sample_separation_signal"] >= 0.5
    assert row["weil_form_evaluated"] is False
    assert row["R_d_pure"] > row["R_d_cancelled"] + 1e-9


@pytest.fixture
def bounded_artifact_campaign(monkeypatch):
    """Exercise artifact plumbing with small real calibrations, not default sweeps."""
    import pbss.hstar_anticancellation as hstar
    monkeypatch.setattr(hstar, "analytic_envelope_campaign", lambda **kwargs: [
        hstar.analytic_vs_numeric_row(16.0, n_points=1024)])
    monkeypatch.setattr(hstar, "arithmetic_bridge_campaign", lambda **kwargs: [
        hstar.arithmetic_injection_bridge_row(1e4, n_strip=3, n_points=512)])


def test_multi_T_artifacts_durable(tmp_path, bounded_artifact_campaign):
    path = write_multi_T_artifacts(tmp_path, T_values=(8.0, 16.0, 32.0))
    assert path.exists()
    data = json.loads(path.read_text())
    assert data["rh_claimed"] is False
    assert len(data["rows"]) >= 3
    assert data["all_constraint_signals"] is False
    assert data["all_sample_separation_signals"] is True
    assert data["candidate_status"] == "retracted_coercivity_proxy"
    grid = tmp_path / "grid.json"
    assert grid.exists()
    g = json.loads(grid.read_text())
    assert g["n_rows"] >= 3
    assert g["eps_star"]["eps_star_global_half_min"] > 0
    assert g["fraction_meets_eps_star"] == 0.0
    assert g["fraction_meets_sample_reference"] > 0.9
    assert g["eps_star"]["is_uniform_lower_bound"] is False


def test_campaign_artifacts_include_partial(tmp_path, bounded_artifact_campaign):
    paths = write_campaign_artifacts(
        tmp_path,
        T_values=(8.0, 16.0),
        sigma_values=(0.8, 0.9),
        degrees=(4,),
        max_workers=4,
        make_plot=False,
    )
    partial = json.loads(Path(paths["partial"]).read_text())
    assert len(partial["rows"]) >= 5
    assert partial["rh_claimed"] is False


def test_phase_aware_continuous_identity_matches_resolved_grid():
    lead = analytic_off_critical_Rd_leading(32.0, sigma=0.9, degree=4)
    row = analytic_vs_numeric_row(32.0, sigma=0.9, degree=4, n_points=8192)
    assert lead["R_d_leading"] == lead["R_d_phase_leading"]
    assert lead["R_d_phase_leading"] <= lead["R_d_limsup_envelope"] * (1 + 1e-12)
    assert row["relative_discretization_error"] < 0.002
    assert lead["eps_star_analytic_is_lower_bound"] is False
    assert lead["continuous_Rd_limit"] == 0.0


def test_phase_trough_rejects_legacy_half_envelope_lower_bound():
    from pbss.off_critical_exact import phase_subsequence_T
    T = phase_subsequence_T(100, A=0.4, t=14.134725, degree=4)
    row = analytic_vs_numeric_row(T, sigma=0.9, degree=4, n_points=8192)
    assert row["R_d_continuous_exact"] < row["R_d_limsup_envelope"] / 4
    assert row["numeric_above_legacy_half_envelope"] == 0.0
    assert row["coercivity_status"] == "retracted_coercivity_proxy"


def test_partial_cancel_curve_matches_normalized_projector_identity():
    curve = partial_cancel_curve(24.0, sigma=0.85, degree=4)
    assert max(row["identity_absolute_error"] for row in curve) < 2e-14
    assert curve[2]["C_d_predicted"] != pytest.approx((1 - curve[2]["alpha"]) * curve[0]["C_d"], rel=1e-5)


def test_arithmetic_injection_bridge_restores_Cd():
    """Real arith leftover + off-mode bump: C_d matches prediction (shipped path)."""
    row = arithmetic_injection_bridge_row(1e6, amp=0.1, sigma=0.9, degree=4, n_strip=20)
    assert row["C_d_orth"] < 1e-4
    assert row["pred_ratio"] == pytest.approx(1.0, rel=2e-10)
    assert row["bridge_ok"] >= 0.5
    assert row["injected_model_component"] is True
    assert row["arithmetic_off_zero_forced"] is False
    assert row["bridge_role"] == "projection_identity_calibration_only"


def test_doc_hstar_records_retraction_and_valid_identity():
    doc = ROOT / "docs" / "HSTAR_ANTICANCELLATION.md"
    text = doc.read_text()
    assert "retracted_coercivity_proxy" in text
    assert "OFF_CRITICAL_PHASE_OBSTRUCTION.md" in text
    assert "not a uniform lower bound" in text
    assert "arithmetic injection bridge" in text.lower()
    assert "B-RES" in text and "remain open" in text
