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
    arithmetic_residual,
    primes_upto,
    normalize_l2,
    finite_cl_superposition,
    explicit_formula_residual,
    peel_residual,
    arithmetic_zero_peel,
)
from .beurling import (
    beurling_theta_residual,
    gapped_beurling_primes,
    thinned_ordinary_primes,
    default_battery_specs,
    build_system_primes,
)
from .lemmas import (
    continuous_R_d_pure_mode,
    continuous_R_d_orthogonal_defect,
    synthetic_orthogonal_defect,
    predicted_R_d_critical_scaling,
    critical_line_omega,
    bound_R_d_finite_mode_sum,
    finite_mode_R_d_order_T,
)
from .zeros import zeta_zero_ordinates, ZETA_ZERO_ORDINATES_50
from .gpu_residual import arithmetic_residual_fast, energy_ratios_multi_degree
from .weights import (
    admissible_weight,
    apply_weight,
    bulk_vs_weighted_report,
    endpoint_contribution,
    weighted_energy_ratio,
)
from .remainder import (
    bound_R_d_mode_tail,
    peel_via_remainder,
    remainder_diagnostic,
    truncated_mode_sum,
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
    "arithmetic_residual",
    "primes_upto",
    "normalize_l2",
    "finite_cl_superposition",
    "explicit_formula_residual",
    "peel_residual",
    "arithmetic_zero_peel",
    "beurling_theta_residual",
    "gapped_beurling_primes",
    "thinned_ordinary_primes",
    "default_battery_specs",
    "build_system_primes",
    "continuous_R_d_pure_mode",
    "continuous_R_d_orthogonal_defect",
    "synthetic_orthogonal_defect",
    "predicted_R_d_critical_scaling",
    "critical_line_omega",
    "bound_R_d_finite_mode_sum",
    "finite_mode_R_d_order_T",
    "zeta_zero_ordinates",
    "ZETA_ZERO_ORDINATES_50",
    "arithmetic_residual_fast",
    "energy_ratios_multi_degree",
    "admissible_weight",
    "apply_weight",
    "bulk_vs_weighted_report",
    "endpoint_contribution",
    "weighted_energy_ratio",
    "bound_R_d_mode_tail",
    "peel_via_remainder",
    "remainder_diagnostic",
    "truncated_mode_sum",
]

__version__ = "0.6.0"
