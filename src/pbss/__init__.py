"""
Perry–Beurling Spectral Sieve (PBSS)

Diagnostic framework: project a density perturbation q onto low-degree
orthonormal polynomials and measure projection strength / L² energy ratio.

This is a classifier / conditional diagnostic for RH-like spectral structure
on Beurling generalized prime systems — **not** a proof of RH.
"""

from .basis import shifted_legendre_values, orthonormal_legendre_design
from .projection import (
    project,
    project_coefficients,
    projection_energy,
    energy_ratio,
    scaled_projection_strength,
    ProjectionResult,
)
from .probes import (
    sample_grid,
    probe_low_degree,
    probe_high_frequency,
    probe_critical_line_mode,
    probe_defective,
    probe_prime_residual,
    normalize_l2,
)

__all__ = [
    "shifted_legendre_values",
    "orthonormal_legendre_design",
    "project",
    "project_coefficients",
    "projection_energy",
    "energy_ratio",
    "scaled_projection_strength",
    "ProjectionResult",
    "sample_grid",
    "probe_low_degree",
    "probe_high_frequency",
    "probe_critical_line_mode",
    "probe_defective",
    "probe_prime_residual",
    "normalize_l2",
]

__version__ = "0.1.0"
