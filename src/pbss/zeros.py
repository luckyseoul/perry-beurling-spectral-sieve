"""
Known nontrivial zeta zero ordinates t_n = Im(ρ_n) for offline explicit-formula
residuals. Values are standard tabulated imaginary parts (first 50).

Not a claim about RH: under the model residual we *place* zeros on the
critical line by construction (σ = 1/2). Full Theorem A for the arithmetic
prime residual remains open.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

# First 50 ordinates (standard tables / Odlyzko-class listings).
ZETA_ZERO_ORDINATES_50: tuple[float, ...] = (
    14.134725141734693,
    21.022039638771554,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
    37.586178158825671,
    40.918719012147495,
    43.327073280914999,
    48.005150881167159,
    49.773832477672302,
    52.970321477714460,
    56.446247697063394,
    59.347044002602353,
    60.831778524609809,
    65.112544048081606,
    67.079810529494173,
    69.546401711173979,
    72.067157674481907,
    75.704690699083933,
    77.144840068874805,
    79.337375020249367,
    82.910380854086030,
    84.735492980517050,
    87.425274613125229,
    88.809111207634465,
    92.491899270558484,
    94.651344040519812,
    95.870634228245309,
    98.831194218193692,
    101.31785100573139,
    103.72553804047833,
    105.44662305232609,
    107.16861118427640,
    111.02953554316967,
    111.87465917699263,
    114.32022091545271,
    116.22668032085755,
    118.79078286597621,
    121.37012500242064,
    122.94682929355258,
    124.25681855434594,
    127.51668387959649,
    129.57870419995605,
    131.08768853093265,
    133.49773720299758,
    134.75650975337387,
    138.11604205453344,
    139.73620895212138,
    141.12370740489593,
    143.11184580762063,
)


def zeta_zero_ordinates(
    n: Optional[int] = None,
    ordinates: Optional[Sequence[float]] = None,
) -> np.ndarray:
    """
    Return the first ``n`` positive imaginary parts of nontrivial zeros.

    Parameters
    ----------
    n : if set, take first n from the built-in table (or from ``ordinates``).
    ordinates : optional override table (must be sorted ascending positive).
    """
    base = (
        np.asarray(ordinates, dtype=np.float64)
        if ordinates is not None
        else np.asarray(ZETA_ZERO_ORDINATES_50, dtype=np.float64)
    )
    if base.size == 0:
        raise ValueError("empty zero table")
    if np.any(base <= 0):
        raise ValueError("ordinates must be positive")
    if n is None:
        return base.copy()
    n = int(n)
    if n < 1:
        raise ValueError("n >= 1")
    if n > base.size:
        raise ValueError(
            f"requested n={n} zeros but table has only {base.size}; "
            "pass a longer ordinates array"
        )
    return base[:n].copy()


def explicit_formula_amplitudes(
    ordinates: Union[Sequence[float], np.ndarray],
) -> np.ndarray:
    """
    Standard |ρ|^{-1} scale for a zero at ρ = 1/2 + i t:

      a_n = 2 / |ρ_n| = 2 / sqrt(1/4 + t_n²)

    (pair contribution folded into the real cosine form).
    """
    t = np.asarray(ordinates, dtype=np.float64)
    return 2.0 / np.sqrt(0.25 + t * t)


def explicit_formula_phases(
    ordinates: Union[Sequence[float], np.ndarray],
) -> np.ndarray:
    """Phase α_n = arg(1/2 + i t_n) for the real-part cosine form."""
    t = np.asarray(ordinates, dtype=np.float64)
    return np.arctan2(t, 0.5)
