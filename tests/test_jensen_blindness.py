"""Tests for Jensen/moment hierarchy blindness helpers (shipped path only)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.jensen_blindness import (  # noqa: E402
    ACCESSIBLE_MOMENT_ORDER,
    ACCESSIBLE_SHIFT_N_MAX,
    LEHMER_ZERO_INDEX,
    approx_height_of_nth_zero,
    central_shift_moment_order,
    index_argument_report,
    jensen_blindness_report,
    laguerre_L1_ratio,
    max_zero_ordinal_probed_by_even_moment_order,
    min_even_moment_order_for_zero_ordinal,
    model_false_hyperbolicity_demo,
    riemann_von_mangoldt_N,
)


def test_max_zero_ordinal_is_floor_half_order():
    # Drive the shipped function — do not hardcode only the answer without the call
    assert max_zero_ordinal_probed_by_even_moment_order(0) == 0
    assert max_zero_ordinal_probed_by_even_moment_order(1) == 0
    assert max_zero_ordinal_probed_by_even_moment_order(34) == 17
    assert max_zero_ordinal_probed_by_even_moment_order(13418) == 6709
    with pytest.raises(ValueError):
        max_zero_ordinal_probed_by_even_moment_order(-1)


def test_min_moment_order_for_lehmer_is_twice_index():
    m = min_even_moment_order_for_zero_ordinal(LEHMER_ZERO_INDEX)
    assert m == 2 * LEHMER_ZERO_INDEX
    # Accessible order cannot reach Lehmer under the index inequality
    assert max_zero_ordinal_probed_by_even_moment_order(ACCESSIBLE_MOMENT_ORDER) < LEHMER_ZERO_INDEX
    assert m > ACCESSIBLE_MOMENT_ORDER


def test_central_shift_maps_30_to_34():
    assert central_shift_moment_order(ACCESSIBLE_SHIFT_N_MAX, pad=4) == ACCESSIBLE_MOMENT_ORDER


def test_index_argument_report_blinds_lehmer_via_shipped_path():
    rep = index_argument_report()
    assert rep["banner"].startswith("NOT AN UNCONDITIONAL")
    assert rep["rh_claimed"] is False
    assert rep["lambda_upper_bound_claimed"] is False
    assert rep["index_blinds_lehmer"] is True
    # Consistency: probed ordinal from accessible order
    assert rep["max_zero_ordinal_probed"] == max_zero_ordinal_probed_by_even_moment_order(
        rep["moment_order_accessible"]
    )
    assert rep["min_moment_order_for_lehmer"] == min_even_moment_order_for_zero_ordinal(
        rep["lehmer_zero_index"]
    )
    assert rep["moment_order_gap"] == (
        rep["min_moment_order_for_lehmer"] - rep["moment_order_accessible"]
    )
    # Historical numerics labeled not re-run
    hist = rep["historical_false_hyperbolicity"]
    assert hist["rerun_here"] is False
    assert hist["t_min_certified"] == -0.7


def test_approx_height_monotone_in_n():
    h30 = approx_height_of_nth_zero(30)
    h6709 = approx_height_of_nth_zero(6709)
    assert h30 < 200.0  # asymptotic ballpark; ≲101 is project record with tables
    assert h6709 > 1000.0
    assert h6709 > h30
    # N(T) roughly recovers n
    assert riemann_von_mangoldt_N(h6709) > 6000.0


def test_laguerre_ratio_finite_on_real_gammas():
    # first few ordinates as pure real model
    from pbss.zeros import zeta_zero_ordinates

    g = zeta_zero_ordinates(10)
    tau = laguerre_L1_ratio(g)
    assert tau == tau  # not NaN
    assert tau > 0.0


def test_model_distant_sensitivity_scales_like_one_over_h2():
    demo = model_false_hyperbolicity_demo(n_low=20, complex_height=7005.0)
    assert demo["sensitivity_abs"] >= 0.0
    # Adding a height-h feature changes the certificate by O(1/h^2) scale
    assert demo["scale_1_over_h2"] == pytest.approx(1.0 / 7005.0**2, rel=1e-12)
    # Sensitivity should be on the order of a few / h^2 (not O(1))
    assert demo["sensitivity_abs"] < 50.0 * demo["scale_1_over_h2"]


def test_jensen_blindness_report_entry_point():
    rep = jensen_blindness_report(n_shift=30, include_model_demo=True)
    assert rep["rh_claimed"] is False
    assert "NOT AN UNCONDITIONAL" in rep["banner"]
    idx = rep["index_argument"]
    assert idx["index_blinds_lehmer"] is True
    assert "model_distant_sensitivity" in rep
