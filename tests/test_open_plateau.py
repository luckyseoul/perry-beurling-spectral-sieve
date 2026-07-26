"""Tests for open-plateau intervention entry points (shipped runner jobs)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.probes import primes_upto  # noqa: E402


def _load_runner():
    path = ROOT / "experiments" / "run_open_plateau.py"
    spec = importlib.util.spec_from_file_location("run_open_plateau", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def runner():
    return _load_runner()


@pytest.fixture
def tiny_tables(tmp_path):
    primes = primes_upto(50_000)
    csum = np.cumsum(np.log(primes.astype(np.float64)))
    p_path = tmp_path / "primes_le_50000.npy"
    c_path = tmp_path / "csum.npy"
    np.save(p_path, primes)
    np.save(c_path, csum)
    return str(p_path), str(c_path)


def test_arith_job_returns_multi_T_Rd(runner, tiny_tables):
    primes_path, csum_path = tiny_tables
    rows = []
    for T in (8.0, 9.5, 10.5):
        row = runner._arith_job(
            {
                "T": T,
                "primes_path": primes_path,
                "csum_path": csum_path,
                "n_points": 1024,
                "degree": 4,
                "detrend": "deg1",
                "smooth": 1,
                "variant": "deg1",
            }
        )
        assert "R_d" in row
        assert 0.0 <= float(row["R_d"]) <= 1.0 + 1e-9
        assert float(row["T"]) > 0
        assert row["n_primes"] > 10
        rows.append(row)
    assert len({r["T"] for r in rows}) == 3


def test_arith_job_measure_norms(runner, tiny_tables):
    primes_path, csum_path = tiny_tables
    norms = []
    for norm in ("sqrt", "x", "plain", "logx"):
        row = runner._arith_job_measure(
            {
                "T": 9.0,
                "primes_path": primes_path,
                "csum_path": csum_path,
                "n_points": 1024,
                "degree": 4,
                "detrend": "deg1",
                "norm": norm,
            }
        )
        assert row["variant"] == f"norm_{norm}"
        assert 0.0 <= float(row["R_d"]) <= 1.0 + 1e-9
        norms.append(row["R_d"])
    # Distinct normalizations should not all be bit-identical in general
    assert len(set(round(x, 12) for x in norms)) >= 2


def test_resume_stamp_and_done(runner, tmp_path):
    out = tmp_path / "open_plateau_smoke"
    assert runner._done(out, "PEEL") is False
    runner._stamp(out, "PEEL", {"status": "completed", "n_rows": 3, "deep": False})
    assert runner._done(out, "PEEL") is True
    payload = json.loads((out / "PHASE_PEEL_COMPLETE").read_text())
    assert payload["n_rows"] == 3
    state = json.loads((out / "open_plateau_state.json").read_text())
    assert "PEEL" in state["phases"]
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in state["banner"]


def test_write_class_summary(runner, tmp_path):
    phase = tmp_path / "whiten"
    summary = {
        "class": "whiten",
        "hypothesis": "bulk inflates low-degree mass",
        "elapsed_s": 1.23,
        "n_rows": 2,
        "rows": [{"T": 8.0, "R_d": 0.2}, {"T": 10.0, "R_d": 0.18}],
        "deep": False,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }
    runner._write_class(phase, "whiten", summary)
    data = json.loads((phase / "whiten.json").read_text())
    assert data["n_rows"] == 2
    assert len(data["rows"]) == 2
    assert all("R_d" in r for r in data["rows"])
    txt = (phase / "whiten.txt").read_text()
    assert "NOT AN UNCONDITIONAL PROOF OF RH" in txt


def test_marathon_specs_for_deep_beurling_bar():
    """Deep Beurling axis requires ≥500 systems from shipped marathon_battery_specs."""
    from pbss.beurling import marathon_battery_specs

    specs = marathon_battery_specs(500)
    assert len(specs) >= 500
    kinds = {s["kind"] for s in specs}
    assert "rh_like" in kinds and "defective" in kinds
