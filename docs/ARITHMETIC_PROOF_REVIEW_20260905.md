# Arithmetic proof audit and successor checkpoint

**Date:** 2026-09-05. **RH and the normalized converse remain open.**
This record follows the diagnostic repair (`cc3191c`) and the exact off-critical
model theorem (`9ec6f8a`), both pushed to `origin/main` at their milestones.

## Results established in this audit

| Result | Scope and proof |
|---|---|
| Conditional continuous A | Under RH, the exact θ residual after continuous orthogonal affine detrending has R_d=O(T^-2), fixed d>=2. [Proof](CONTINUOUS_THEOREM_A.md). |
| Absolute-moment converse B-ABS | For any one fixed Legendre degree, subexponential absolute moment growth for all sufficiently large windows is equivalent to RH. [Proof](ABSOLUTE_MOMENT_CONVERSE.md). |
| Unconditional normalized recurrence | For that continuously detrended arithmetic residual, liminf R_d=0 without RH. This does not establish a full limit or rate. [Proof](NORMALIZED_RECURRENCE_THEOREM.md). |
| Exact real off-mode obstruction | The intact continuous real exponential sine has R_d->0, a phase-dependent T^-1 envelope, and a sharp T^-3 trough constant for d>=1. [Proof](OFF_CRITICAL_PHASE_OBSTRUCTION.md). |
| Weil transfer obstruction | An injective restriction of an indefinite form does not supply positive projection coercivity. [Counterexample and replacement hypotheses](WEIL_TRANSFER_OBSTRUCTION.md). |

The arithmetic proofs retain the exact residual, all zeros, the prime-power
correction, and the absolute transform. They assume no zero simplicity or
linear independence. Continuous A uses the explicit formula and mean-square
statements in Akbary–Ng–Shahabi, with equation numbers and primary-source links
in its proof. B-ABS uses the Euler-product logarithmic derivative and its
meromorphic continuation. No historical novelty claim is made.

Three agents independently reviewed the load-bearing arithmetic arguments;
the exact model's sharp phase constant and the corrected package statements
also received an independent derivation review. These were internal reviews,
not external peer review or a machine-checked formal proof.

## Corrections to executable claims

The status APIs now separate continuous theorems from sampled model reports,
optional weighted extensions, and unresolved transfer estimates. A named
unresolved analytic input cannot be hidden by a zero count of unnamed gaps.

The former H* component bound omitted projected remainder cancellation. A
valid sufficient margin subtracts that remainder norm. The pure-mode L*
coercivity proxy is marked `retracted_coercivity_proxy`: its old expression was
an upper envelope, not a lower bound. Partial cancellation uses the changed
total norm. The arithmetic injection experiment is labeled as a projection
calibration and uses the same trapezoid inner product for every norm.

Fresh campaign output defaults to `results/hstar_corrected/`. Earlier H*
artifacts are retained byte-for-byte, with a dated superseding notice. Shared
scaled off-mode evaluation prevents avoidable exponential overflow in the
sampled comparison APIs. The CLI writes structured JSON using its common
validated writer.

## Validation and replay

The preceding exact-model milestone passed 26 regression tests. Nuka independently
checked 24 phase cases at 80 decimal digits with 14 CPU workers, plus direct
quadrature; maximum float64 ratio error was 8.674e-18. Source hashes, software
versions, and measured errors are in
[`results/off_critical_phase_obstruction_20260905/`](../results/off_critical_phase_obstruction_20260905/).
Those computations corroborate analytic identities; they are not interval
certificates or proofs of an all-window arithmetic estimate.

The integrated review records its test output, status reports, CLI JSON, and
source hashes in
[`results/arithmetic_proof_review_20260905/`](../results/arithmetic_proof_review_20260905/).
The full suite passed **332 tests in 17.45 seconds** with 86 workers. The CLI
JSON smoke check also passed and reports the retraction and open arithmetic
claims explicitly.
Replay the full regression suite from the repository root:

```bash
python3 -m pytest -n 86 tests -q
PYTHONPATH=src python3 -m pbss hstar --json-out /tmp/pbss-hstar-corrected.json
PYTHONPATH=src python3 experiments/check_off_critical_exact.py --workers 86 --out /tmp/pbss-phase-check.json
```

Use the local compute budget to adjust worker counts on other hosts. The
optional precision checker requires mpmath; the core continuous evaluator uses
NumPy and SciPy. Exact rational tests corroborate the complete finite Laurent
principal parts and the recurrence constants. They do not test RH itself.

## Remaining target

The proved B-ABS equivalence identifies a precise missing estimate:
for one fixed k>=2, establish |c_k(T)|=O_epsilon(exp(epsilon T)) for every
epsilon>0 and all sufficiently large T, **without assuming RH**. That estimate
is RH-equivalent. No argument in this repository currently establishes it
unconditionally.

For a normalized alternative, state one precise rate and control the absolute
norm or otherwise retain the arithmetic pole information discarded by division.
The positive-floor B-RES implication is itself RH-equivalent, since unconditional
recurrence already rules out such a floor for the exact continuous residual.
Repeating finite-grid separations, injecting a model bump, or extending the
same phase-blind sweep cannot supply the missing estimate.
