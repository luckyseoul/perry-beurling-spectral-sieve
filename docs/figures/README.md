# Figures (canonical)

Single store for README and paper graphics. Do **not** duplicate under `docs/paper/figures/`.

| File | Description |
|------|-------------|
| [`continuous_arithmetic_moments.png`](continuous_arithmetic_moments.png) / [SVG](continuous_arithmetic_moments.svg) | Signed absolute c2 and continuous affine R4, 36 finite ordinary-prime windows through endpoint 100,000 |
| [`unattained_spectrum.png`](unattained_spectrum.png) / [SVG](unattained_spectrum.svg) | Explicit Laplace poles and the proved asymptotic rate factor; analytic example, no prime realization asserted |
| `pipeline.png` | Residual → Legendre → scorecard |
| `scorecard_Rd.png` | Diagnostic \(R_4\) scorecard |
| `residuals_waveforms.png` | Example residuals with \(R_4\) |
| `mode_spectrum.png` | Energy by Legendre mode |
| `grand_Rd_vs_T.png` | Grand-campaign multi-\(T\) |
| `grand_arith_focus_linear.png` | Arithmetic soft plateau |
| `beurling_battery_Rd_vs_T.png` | Ordinary vs defective Beurling |
| `status_map.png` | Theorem / lemma status map |
| `campaign_schematic.png` | Campaign schematic |

Paper TeX includes these as `../figures/...`.

The two 2026-09-05 figures have standalone generators and recorded inputs:

```bash
python3 experiments/plot_arithmetic_moments.py --workers 86
python3 experiments/plot_unattained_spectrum.py
```

The first command uses at most the 36 independent windows; choose the worker
budget for the host. The second evaluates small formula arrays. Both require
optional Matplotlib. Data, source hashes, figure hashes, and scope are recorded
in [`results/readme_figures_20260905/`](../../results/readme_figures_20260905/).
