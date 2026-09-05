"""Interpretation shared by the CLI and the reproducible control experiment.

This checks a specified synthetic control pair, never an RH classifier.
Sampling checks apply to the known sinusoidal probes only. They cannot certify
the resolution of arbitrary input data or the prime-counting staircase.
"""
from __future__ import annotations

import math


def assess_demo(rows, *, degree: int, n_points: int, T: float) -> dict:
    """Compare the degree-1 contaminated control with both oscillatory controls.

    A 0.5 separation is an explicit instrument check for this fixed mixture,
    not a universal detection threshold. Insufficient sampling or a projection
    excluding the injected degree makes the check inconclusive.
    """
    values = {row["name"]: float(row["Rd"]) for row in rows}
    required = {"high_frequency", "critical_line", "defective"}
    if not required <= values.keys():
        raise ValueError("demo assessment requires all three synthetic controls")
    if any(not math.isfinite(values[name]) or not 0 <= values[name] <= 1 for name in required):
        raise ValueError("control energy ratios must be finite and in [0,1]")
    if not math.isfinite(T) or T <= 0 or n_points < 2 or degree < 0:
        raise ValueError("demo requires T > 0, n_points >= 2, degree >= 0")
    max_cycles = max(48.0, 14.134725 * T / (2 * math.pi))
    points_per_cycle = (n_points - 1) / max_cycles
    resolved = points_per_cycle > 2.0
    supported = degree >= 1
    sep_hf = values["defective"] - values["high_frequency"]
    sep_cl = values["defective"] - values["critical_line"]
    warnings = []
    if not supported:
        warnings.append("Degree 1 defect is outside the selected projection subspace.")
    if not resolved:
        warnings.append("Known synthetic frequencies reach or exceed the grid Nyquist limit.")
    elif points_per_cycle < 8:
        warnings.append("Fewer than 8 samples per shortest cycle; refine the grid before interpreting separation.")
    passed = supported and resolved and min(sep_hf, sep_cl) >= 0.5
    return {
        "status": "inconclusive" if not supported or not resolved else ("passed" if passed else "failed"),
        "classifier_ok": bool(passed),  # retained for report consumers
        "minimum_separation": 0.5,
        "energy_ratio_defective_minus_hf": sep_hf,
        "energy_ratio_defective_minus_critical": sep_cl,
        "defect_in_subspace": supported,
        "sampling": {
            "max_cycles": max_cycles,
            "points_per_cycle": points_per_cycle,
            "nyquist_resolved": resolved,
        },
        "warnings": warnings,
        "scope": "Separation of the specified synthetic controls only; no RH conclusion or guarantee for arbitrary defects.",
    }
