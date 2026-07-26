"""Tests for infinite-tail scaffolding majorant (shipped remainder API)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.remainder import bound_infinite_zero_tail_scaffold  # noqa: E402


def test_infinite_tail_scaffold_decays_in_T():
    lo = bound_infinite_zero_tail_scaffold(12.0, n_kept=15, N_eff=2000, d=4)
    hi = bound_infinite_zero_tail_scaffold(48.0, n_kept=15, N_eff=2000, d=4)
    assert hi["bound_R_d_model_tail"] < lo["bound_R_d_model_tail"]
    assert "scaffolding" in lo["label"]
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in lo["banner"]
    assert "arithmetic" in lo["does_not_control"].lower() or "ψ" in lo["does_not_control"] or "psi" in lo["does_not_control"].lower() or "prime" in lo["does_not_control"].lower()


def test_infinite_tail_requires_positive_T():
    with pytest.raises(ValueError):
        bound_infinite_zero_tail_scaffold(0.0, n_kept=5, N_eff=10)
