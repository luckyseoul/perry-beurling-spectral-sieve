# Rank 3 — Continuous theorem and sampled-interface audit

**Date:** 2026-09-05 (correction of the 2026-08-11 registry)\
**Code:** `pbss.ant_audit.ant_interface_audit`\
**Related:** [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) ·
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md)

## Purpose

Record which analytic inputs establish the exact continuous theorem, and which
extensions still need their own estimates. Calling the API reads a theorem/source
registry; it does not execute a proof checker or verify an arithmetic limit.

The continuous theorem assumes RH and uses the exact θ residual, continuous
Lebesgue measure, orthogonal affine detrending, fixed degree at least two, and
no additional weight. It proves the normalized ratio is O(T^-2).

## Entry

```bash
PYTHONPATH=src python3 -c "from pbss.ant_audit import ant_interface_audit; import json; print(json.dumps(ant_interface_audit(), indent=2))"
```

## Disposition

| ID | Status | Evidence or remaining obligation |
|----|--------|----------------------------------|
| ANT-3 exact ψ formula | Precisely cited | Akbary–Ng–Shahabi equation (1.5), with sign and phase |
| ANT-1 uniform moments | Adaptation proved | Cutoff G=exp(T), absolutely summable coefficient bound |
| ANT-2 prime-square correction | Adaptation proved | Constant -1 removed; secondary error belongs to L1 and L2 |
| ANT-MS denominator | Precisely cited and adapted | Corollary 1.15(i), positive mean-square limit |
| ANT-4 weighted transfer | Named optional gap | Weighted numerator and positive denominator estimates |
| SAMPLED grid/regression transfer | Named separate gap | Quantitative quadrature and detrend errors as T grows |
| M7 perturbation | Model lemma proved | Its numerical checks do not prove vanishing remainders |

The exact source, hypotheses, and deductions appear in
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md).

## Status behavior

`freeze_continuous_theorem_packaging=true` records that the required continuous
proof obligations have been resolved. The legacy broad
`freeze_full_a_packaging` flag is **false** while the weighted and sampled
transfers remain open. A named required gap prevents even the continuous freeze;
zero unlabeled rows alone are not evidence of a complete proof.

The API reports `analytic_inputs_verified_by_this_call=false`,
`sampled_transfer_proved=false`, and `weighted_transfer_proved=false`.

## Correction record and non-claims

The August registry marked generic textbook interfaces as matched and used the
absence of unlabeled rows to freeze the package. That did not establish the
claimed sign/phase adaptation, truncation uniformity, secondary correction, or
norm lower bound. The September continuous proof supplies explicit replacements.

This audit does not prove RH, establish a sampled asymptotic rate, or assert
that a bounded endpoint taper preserves arithmetic norm comparability. Model
tapers and numerical reports remain available within their stated scope.
