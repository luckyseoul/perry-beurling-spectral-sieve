"""
Perry–Beurling Spectral Sieve (PBSS)

Diagnostic framework: project a density perturbation q onto low-degree
orthonormal polynomials and measure projection strength / L² energy ratio.

This is a classifier / conditional diagnostic for RH-like spectral structure
on Beurling generalized prime systems — **not** an unconditional proof of RH.
See docs/THEOREMS_AB.md for precise A/B statements and proved lemmas.
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
    probe_off_critical_mode,
    probe_defective,
    probe_persistent_defect,
    probe_prime_residual,
    normalize_l2,
)
from .lemmas import (
    continuous_R_d_pure_mode,
    continuous_R_d_orthogonal_defect,
    synthetic_orthogonal_defect,
    predicted_R_d_critical_scaling,
    critical_line_omega,
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
    "probe_off_critical_mode",
    "probe_defective",
    "probe_persistent_defect",
    "probe_prime_residual",
    "normalize_l2",
    "continuous_R_d_pure_mode",
    "continuous_R_d_orthogonal_defect",
    "synthetic_orthogonal_defect",
    "predicted_R_d_critical_scaling",
    "critical_line_omega",
]

__version__ = "0.2.0"
