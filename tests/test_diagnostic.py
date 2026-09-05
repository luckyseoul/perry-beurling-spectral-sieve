"""Control checks must distinguish an instrument failure from missing evidence."""
import pytest

from pbss.diagnostic import assess_demo


def controls(hf=0.01, cl=0.02, defect=0.9):
    return [
        {"name": "high_frequency", "Rd": hf},
        {"name": "critical_line", "Rd": cl},
        {"name": "defective", "Rd": defect},
    ]


def test_separation_checks_both_controls_and_explicit_margin():
    assert assess_demo(controls(), degree=4, n_points=2048, T=20)["classifier_ok"]
    result = assess_demo(controls(cl=0.7), degree=4, n_points=2048, T=20)
    assert result["status"] == "failed"
    assert not result["classifier_ok"]


@pytest.mark.parametrize("degree,n_points,T", [(0, 2048, 20), (4, 97, 20), (4, 2048, 1000)])
def test_unsupported_or_aliased_controls_never_pass(degree, n_points, T):
    result = assess_demo(controls(), degree=degree, n_points=n_points, T=T)
    assert result["status"] == "inconclusive"
    assert not result["classifier_ok"]
    assert result["warnings"]


def test_nonfinite_control_cannot_pass():
    with pytest.raises(ValueError, match="finite"):
        assess_demo(controls(defect=float("nan")), degree=4, n_points=2048, T=20)
