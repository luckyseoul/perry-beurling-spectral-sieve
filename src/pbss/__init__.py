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
    bound_R_d_weighted_sine_order,
    bound_R_d_weighted_finite_mode_sum,
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
    bound_infinite_zero_tail_scaffold,
    peel_via_remainder,
    remainder_diagnostic,
    truncated_mode_sum,
)
from .theorem_a_chain import model_chain_report, multi_T_model_chain, package_status
from .ab_closure import (
    ant_citations,
    conditional_full_a_report,
    energy_ratio_perturbation_bound,
    full_a_gap_table,
    full_b_gap_table,
    off_critical_model_obstruction,
    package_status as full_ab_package_status,
    verify_m7_on_grid,
)
from .jensen_blindness import (
    index_argument_report,
    jensen_blindness_report,
    max_zero_ordinal_probed_by_even_moment_order,
    min_even_moment_order_for_zero_ordinal,
)
from .plateau_secondary import plateau_secondary_report
from .ant_audit import ant_interface_audit
from .zero_proportion_feasibility import zero_proportion_feasibility_report
from .b_res_threshold import b_res_threshold_report
from .ef_identify import (
    identify_ef,
    hypothesis_residual,
    multi_hypothesis_scan,
    model_sanity_identify,
    summarize_attack,
    multi_N_enrich_scan,
    summarize_enrich_kill021,
    build_m_columns,
    M_ENRICHMENTS,
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
    "bound_R_d_weighted_sine_order",
    "bound_R_d_weighted_finite_mode_sum",
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
    "bound_infinite_zero_tail_scaffold",
    "peel_via_remainder",
    "remainder_diagnostic",
    "truncated_mode_sum",
    "model_chain_report",
    "multi_T_model_chain",
    "package_status",
    "full_ab_package_status",
    "ant_citations",
    "conditional_full_a_report",
    "energy_ratio_perturbation_bound",
    "full_a_gap_table",
    "full_b_gap_table",
    "off_critical_model_obstruction",
    "verify_m7_on_grid",
    "index_argument_report",
    "jensen_blindness_report",
    "max_zero_ordinal_probed_by_even_moment_order",
    "min_even_moment_order_for_zero_ordinal",
    "plateau_secondary_report",
    "ant_interface_audit",
    "zero_proportion_feasibility_report",
    "b_res_threshold_report",
    "identify_ef",
    "hypothesis_residual",
    "multi_hypothesis_scan",
    "model_sanity_identify",
    "summarize_attack",
    "multi_N_enrich_scan",
    "summarize_enrich_kill021",
    "build_m_columns",
    "M_ENRICHMENTS",
]

__version__ = "0.9.1"
