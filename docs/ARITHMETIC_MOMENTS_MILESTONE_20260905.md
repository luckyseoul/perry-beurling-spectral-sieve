# Numerical milestone: signed continuous arithmetic moments

**Date:** 2026-09-05. **Base:** `c33105cae7d1c9d3791ac08707a18cf1bfe5f2f8`.

The new `pbss moments` command and `arithmetic_moments` library API evaluate
signed absolute Legendre moments of the exact ordinary-prime theta residual.
They integrate between consecutive prime jumps and retain the true continuous
norm. This makes the quantity in the [one-sided target](SIGNED_MOMENT_OSCILLATION.md)
directly inspectable without removing its amplitude by normalization.

The [derivation](EXACT_ARITHMETIC_MOMENTS.md) includes finite prime-sum
antiderivatives, stable local exponential integrals, a positive interval norm
identity, and continuous orthogonal affine detrending. Complete prime prefixes
are verified; there is no sample interpolation, smoothing, or clamp below 2.
Decimal JSON values retain platform `longdouble` precision. Numerical error
indicators are explicitly distinguished from certified interval bounds.

## Validation

- The full suite passed: **373 tests in 17.73 seconds**, using 86 workers.
  New checks compare the integral identities with independently split adaptive
  integration, verify continuity across a prime jump, and reject incomplete
  prime tables and malformed windows. A Rodrigues identity independently
  checks the sign and scale of the small-window second moment.
- Nuka ran **15 independent 80-digit cases** for nominal endpoints 2 through
  10,000 and degrees 2, 4, and 6. The checker uses globally cancelling finite
  antiderivatives, separately from the local implementation. Maximum absolute
  discrepancies were `2.015e-16` for a coefficient, `4.836e-18` for the full
  norm, and `6.430e-19` for the affine-detrended norm.
- An additional direct 80-digit integration at endpoint 30 checked the
  independent antiderivative evaluator itself; the largest coefficient
  discrepancy was below `1e-77`.
- The CLI replay at T=10 uses 2,466 primes and reports
  `c_2=0.05476066123266131636`, full squared norm
  `1.6097615914374053387`, and affine-detrended R_4 approximately
  `0.08137291741408735`.

Review caught and corrected two edge cases before this full run: complex or
non-scalar T inputs are rejected before conversion, and unresolved affine norm
subtraction retains valid absolute moments while marking the affine norm and
ratio unavailable. These reviews were internal; they are not formal
verification or external peer review.

Evidence and source hashes are in
[`results/arithmetic_moments_20260905/`](../results/arithmetic_moments_20260905/).
All 75 Python source, test, and experiment files were hashed before the full
suite and compared again before commit. The precision check records the
executing host, library versions, arithmetic precision, and its source hashes.
It exercises CPU arbitrary precision; a dense GPU projection is not involved.

## Replay

```bash
python3 -m pytest -n 86 tests -q
PYTHONPATH=src python3 -m pbss moments --T 10 --degree 4 \
  --json-out results/arithmetic_moments_20260905/cli_T10.json
python3 experiments/check_arithmetic_moments.py --workers 14 --digits 80 \
  --out results/arithmetic_moments_20260905/checks.json
```

The precision experiment requires optional `mpmath`. Choose the worker count
from the executing host's compute budget. The recorded precision run used Nuka
with 14 workers; the full test suite used the local host with 86 workers.

## Remaining implication

No finite window establishes the all-large-T one-sided subexponential estimate
for c_2. That unconditional bound remains RH-equivalent and open. The
normalized ordinary-prime converse also remains open when the supremum of zero
real parts is unattained. This milestone supplies a faithful finite evaluator,
not a proof of either asymptotic target.
