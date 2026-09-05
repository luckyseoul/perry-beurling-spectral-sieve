# README figures, 2026-09-05

The arithmetic figure contains 36 finite-window evaluations of the exact
continuous ordinary-prime theta residual. Its JSON retains signed absolute
moments, norms, ratios, and precision metadata as decimal values. The
arithmetic evaluator was validated in the preceding 373-test milestone and
independently checked at 80 digits. The figure's connected points are finite
measurements, not a fit or an all-window bound.

The unattained-spectrum figure illustrates the explicit pole locations and
proved rate factor from
[`UNATTAINED_SPECTRUM_OBSTRUCTION.md`](../../docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md).
The displayed factor is f(T) in R_d=O_d(f(T)), with its degree-dependent
constant omitted. It is not a sampled R_d curve and does not represent zeros
of ordinary zeta. The example has no asserted prime realization.

Both figures are supplied as PNG and SVG, were visually inspected, and are
embedded in the repository README with their scope stated in the captions.
Each JSON records its generator and mathematical-source hashes and the image
hashes. The plotting generators were run successfully; the previously tested
library and test sources remain unchanged.

Replay from the repository root, with Matplotlib installed:

```bash
python3 experiments/plot_arithmetic_moments.py --workers 86
python3 experiments/plot_unattained_spectrum.py
```

The arithmetic run uses 36 independent workers, capped by the supplied budget.
The analytic figure uses only small vector calculations. Neither requires a
GPU projection or establishes an additional theorem through computation.
