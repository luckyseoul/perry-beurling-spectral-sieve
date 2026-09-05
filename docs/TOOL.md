# PBSS as a general mathematics tool

**Not a proof of RH.** The same projection diagnostic is useful for any 1D residual
on a unit interval: measure how much energy sits in a low-degree polynomial subspace.

## Install

```bash
# Prefer if venv available:
# python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pbss version
```

Without a venv (editable path):

```bash
PYTHONPATH=src python3 -m pbss version
```

## Commands

| Command | Purpose |
|---------|---------|
| `pbss project -i q.npy -d 4 --T 20` | Project arbitrary \(q(u)\) → \(R_d\), \(S_d\), coeffs |
| `pbss diagnose --demo` | Built-in HF / critical-line / defective scorecard |
| `pbss diagnose -i q.csv --u u.csv` | Diagnose samples on a supplied grid |
| `pbss sensitivity` | Flat vs Gamma-weight Fisher \(d'\) |
| `pbss sensitivity --confirm-53` | Confirm ≥53% relative gain claim |
| `pbss scorecard --x-max 1e6` | Ordinary vs gapped/thinned Beurling systems |

## Library API (general use)

```python
import numpy as np
from pbss import project, sample_grid, confirm_sensitivity_claim

u = sample_grid(2048)
q = np.sin(80 * np.pi * u)          # any residual on [0,1]
r = project(q, u, degree=4, T=20.0)
print(r.energy_ratio, r.scaled_strength)

rep = confirm_sensitivity_claim()   # Gamma weight discriminability
print(rep["verdict"], rep["noisy"]["relative_gain_percent"])
```

## Numerical contract

Inputs must be finite, real, one-dimensional arrays with matching lengths.
The grid must be strictly increasing in \([0,1]\). With no grid file, the CLI
uses a uniform grid with the input file's length. It never resamples a file to
`--n-points`; that option sets the built-in controls' grid. Arrays with multiple
columns and complex arrays are rejected instead of flattened or truncated.

The default measure uses composite trapezoid weights. Library callers can pass
nonnegative quadrature weights, with at least \(d+1\) positive-weight nodes.
The requested degree must be identifiable on those nodes. For the default
trapezoid rule, include both endpoints to approximate the full unit interval;
interior-only quadrature nodes require suitable explicit weights.

The sampled Legendre basis is not generally orthonormal under these weights.
The code solves \(\min_c\|W^{1/2}(q-\Phi c)\|_2\) using QR, checks numerical
rank, and computes energy in the same weighted norm. `coeffs` reconstructs the
fitted polynomial as `Phi @ coeffs`; its squared Euclidean norm is **not** the
reported energy. With \(G=\Phi^T W\Phi\),

\[
 E_d=c^TGc,\qquad R_d=E_d/(q^TWq)\in[0,1].
\]

Zero-norm inputs have no energy ratio. Invalid grids, nonfinite data,
unidentifiable degrees, nonpositive/nonfinite `T`, and unrepresentable output
energies or strengths produce errors. JSON reports use finite numbers and
record the parameters and numerical method. CPU and optional CuPy projections
share the same fit and energy definition.

The fit is a finite-sample quantity. Its boundedness does not certify accuracy
of the continuum integral. This distinction follows from continuous Legendre
[orthogonality](https://dlmf.nist.gov/18.3) versus the sampled
[least-squares problem](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html).

## Interpreting a diagnostic

`R_d` measures energy in the selected degree range. It does not identify the
cause of that energy. Defects above degree `d`, cancellation, and unresolved
oscillations can all give small values. Detrending or tapering changes the
residual being measured and must be held fixed for comparisons.

`S_d = T^(2(d+1)) R_d` is retained as a scaling convention. The proved generic
pure/finite model-mode rate is `O(T^-2)`; the stronger archive heuristic is not
generally valid. There is no calibrated universal threshold for RH or arbitrary
input residuals.

The built-in demo checks a specific degree-1 contamination against **both**
oscillatory controls using a stated minimum separation of 0.5. An excluded
defect degree or frequencies at/above the grid's Nyquist limit yield
`inconclusive`; a resolved check missing the separation yields `failed`.
These statuses describe the synthetic instrument check only. The command
returns exit code 2 for either non-passing status. File inputs receive energy
measurements without a synthetic classification. Scorecard failures likewise
appear in JSON and return a nonzero exit code.

For a reproducible control run, including the prime residual and fresh function
evaluations on a nested `2n-1` grid:

```bash
PYTHONPATH=src python3 experiments/run_diagnostic.py \
  --output-dir results/diagnostic_review --workers 4
```

Four independent probes use up to four workers. The refinement check compares
`R_d` at tolerance `0.001 + 0.01 * abs(R_d_fine)` and reports every change.
Both grids must agree on control status; a threshold crossing is inconclusive.
The experiment disables prime-residual smoothing on both grids to hold that
preprocessing choice fixed (the CLI prime demo retains its 3-sample smoothing).
This is a sensitivity check, not a rigorous error bound; a staircase residual
can require further resolution studies. Existing experiment JSON keys are
retained for compatibility; the `rh_like_*` names are historical labels.

Historical campaign artifacts used raw quadrature moments as coefficients.
Fresh reports use the corrected weighted fit and identify that convention;
old numerical artifacts have not been silently recalculated.

## Gamma-weight sensitivity (≥53%)

Project record: Gamma bump \(w(s)=s^{k-1}e^{-\sigma s}\) improved offline/online
discriminability by **~53%** on a noisy Beurling-like ensemble (idealized run often
larger). Shipped confirmation:

```bash
pbss sensitivity --confirm-53 --json-out /tmp/sens.json
```

See [`MEASURE_SENSITIVITY.md`](MEASURE_SENSITIVITY.md).

## Typical non-RH uses

- Detect **low-mode contamination** in a high-frequency residual  
- Compare **weighted** vs flat energy for matched-filter style detection  
- Scorecard **ordinary vs defective** counting systems (Beurling-style)  
- Teaching / demos of orthonormal projection energy ratios  

## Non-claims

- Does not prove RH or close B-RES.  
- Sensitivity gain is **ensemble-dependent**; 53% is a **lower-bound claim** that the
  shipped noisy ensemble meets (often exceeds).
