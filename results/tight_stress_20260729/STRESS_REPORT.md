# PBSS Tight Stress Report — 2026-07-29

**NOT AN UNCONDITIONAL PROOF OF THE RIEMANN HYPOTHESIS.**

Session: Grok Build long-horizon formalization + stress.

## 1. Formalized diagnostics

| Symbol | Definition | Role |
|--------|------------|------|
| $R_d(q)$ | $\|P_d q\|^2/\|q\|^2$ | L² energy ratio in low-degree Legendre space |
| $S_d(q;T)$ | $T^{2(d+1)} R_d$ | scaled strength (working $P(q)$) |
| $P(q)$ | $:= S_d$ | reconstruction convention; **not** legacy 3.92 |
| A₀ / M3 | CL pure mode $R_d=O(T^{-2})$ | **proved** |
| finite A₀ / M5 | finite CL sum same order | **proved** |
| A arithmetic | residual under RH | **open** |
| B₀ / M2–M4 | persistent defect $R_d=\varepsilon^2$ | **proved** |
| B converse | fast decay $\Rightarrow$ RH | **open** |

Legacy archive numbers $P\approx 3.92$, threshold $\approx 29.5$ used lost
high-precision scripts and are **not** hard-coded here.

## 2. MC stress (tighter grids + ablations)

- Total trials: **192,000**
- Degrees: [1, 2, 4, 6, 8]
- T grid: [8.0, 14.0, 20.0, 28.0]
- Ablations: ['baseline', 'heavy_defect', 'light_defect', 'high_freq', 'low_freq', 'high_deg_defect']
- Elapsed: 127.8s

- Baseline flatness: mean $R_4$=0.793935, std across T=1.979008e-03, flat=True

Mid-T ablation $R_4$ means:

| Ablation | mean $R_4$ |
|----------|-----------:|
| baseline | 0.792001 |
| heavy_defect | 0.970319 |
| high_deg_defect | 0.495979 |
| high_freq | 0.791900 |
| light_defect | 0.284841 |
| low_freq | 0.793550 |

**Reading:** Defective MC mass stays high (~0.5–0.95 depending on
weight/wave settings) and **flat in T** — instrument is stable.
Heavy defect lifts $R_d$; light defect lowers it but remains ≫ CL controls.

## 3. Off-critical σ sweep (zero-free *diagnostic*)

- all_far_higher_than_cl: **True**

| T | d | $R_d(\sigma{=}1/2)$ | $R_d(\sigma{=}0.9)$ | ratio |
|--:|--:|-------------------:|-------------------:|------:|
| 8 | 4 | 6.222e-03 | 2.416e-02 | 3.88 |
| 12 | 4 | 2.774e-03 | 1.644e-02 | 5.93 |
| 16 | 4 | 1.562e-03 | 1.241e-02 | 7.94 |
| 20 | 4 | 9.994e-04 | 9.942e-03 | 9.95 |
| 24 | 4 | 6.936e-04 | 8.286e-03 | 11.95 |
| 32 | 4 | 3.894e-04 | 6.203e-03 | 15.93 |

Model off-critical modes vs CL: measured R_d ratios. This is a *diagnostic probe*, not a zero-free-region theorem. Full B (fast arithmetic R_d decay => RH) remains open.

## 4. Expanded Beurling constructions

- Systems: **35**
- $x_{\max}$=1.000e+06
- Defective all above ordinary @ $T_{\max},d4$: **True**
- Ordinary $R_4$: 0.18882323622251154
- Defective stats: {'min': 0.9930217632510576, 'max': 0.9963177510020347, 'mean': 0.9952161927001252, 'median': 0.9954733822254824}
- Failure modes (no separation): none
- Thin-margin systems: 0

## 5. CL rate check (M3 scale)

- m3_scale_ok: **True**
- d=0: R_d ~ numerical floor (max 3.600e-10); skip M3 scale check
- d=1: R_d 1.877e-03->2.874e-05; T2*R_d in [1.177e-01,1.201e-01] (M3-scale OK)
- d=2: R_d 1.877e-03->2.874e-05; T2*R_d in [1.177e-01,1.201e-01] (M3-scale OK)
- d=4: R_d 6.222e-03->9.591e-05; T2*R_d in [3.929e-01,3.998e-01] (M3-scale OK)
- d=6: R_d 1.269e-02->2.017e-04; T2*R_d in [8.123e-01,8.377e-01] (M3-scale OK)
- Stronger archive rate $O(T^{-2(d+1)})$: NOT proved; T^{2(d+1)} R_d often grows

## 6. Failure modes & strengthened bounds

### Holds under stress

1. **B₀ / instrument:** MC defective $R_d$ stays high and flat across wider T and ablations.
2. **A₀ / M3:** pure and finite-mode CL $R_d$ continue to track $O(T^{-2})$ on expanded T.
3. **Beurling separation:** ordinary primes remain well below gapped/thinned families at large T
   for the expanded construction set (when prime table covers $x_{\max}$).

### Degrades / remains weak

1. **Arithmetic soft plateau** (prior campaigns): $R_4\sim0.15$–$0.19$ through $10^{10}$–$5\times10^{10}$ —
   does **not** approach A₀ levels; zero-peel does not close the gap. Full A open.
2. **Theorem B:** off-critical model probes give a directional diagnostic only;
   no reduction from arithmetic residual to off-critical envelopes is proved.
3. **Sharp rate** $O(T^{-2(d+1)})$: not supported as a proved bound; $T^{2(d+1)}R_d$ often grows.
   Stick to M3 rate $O(T^{-2})$ for model modes.
4. **Legacy $P(q)\approx3.92$ / threshold 29.5:** still unrecovered; do not use as classifier thresholds.
5. **Thin-margin Beurling systems:** some mild thinnings / small gaps can approach ordinary $R_d$;
   separation is construction-dependent — battery must keep strongly defective controls.

### What this does *not* do

- Prove RH.
- Prove full Theorems A or B.
- Exclude zeros with $\mathrm{Re}\,\rho\neq 1/2$ from arithmetic data alone.

---

*Private research. Nicholas Perry / Perry–Beurling Spectral Sieve.*
