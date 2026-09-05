# Perry–Beurling Spectral Sieve

**PBSS** in this repository means **Perry–Beurling Spectral Sieve** — a projection
diagnostic for density residuals on logarithmic windows (not the pharmaceutical
society, IEEE WLAN “personal basic service set,” actuarial IAA section, or other
public expansions of the same acronym).

Open research archive and **runnable reconstruction** of the **Perry–Beurling Spectral Sieve** (spectral diagnostic / P(q) framework) for testing consistency with the Riemann Hypothesis on Beurling generalized prime systems.

**Author:** Nicholas Perry  
**Status:** Independent research. Reconstruction of core projection diagnostic (2026-07).  
**License:** MIT — see [`LICENSE`](LICENSE).  
**Not a proof of RH** — see [`docs/STATUS.md`](docs/STATUS.md).  
**Docs index:** [`docs/README.md`](docs/README.md).

## Overview

A spectral approach combining Beurling’s theory of generalized primes / Beurling zeta functions with a projection-based diagnostic. Analyze density perturbations \(q\) associated with prime systems and test whether their low-degree polynomial energy is consistent with all non-trivial zeros on \(\mathrm{Re}(s)=1/2\).

The framework is a **classifier / diagnostic**, not a full proof of RH.
Its finite control separation is not a universal zero-location criterion:
a [positive discrete Beurling system](docs/BEURLING_NORMALIZED_COUNTEREXAMPLE.md)
can have off-critical zeta zeros while all fixed-degree continuous detrended
ratios decay exponentially.

## At a glance

<p align="center">
  <img src="docs/figures/continuous_arithmetic_moments.png" alt="Signed absolute arithmetic moment c2 and continuous affine-detrended R4 across 36 finite logarithmic windows, with endpoints from 2 to 100000" width="1000"/>
</p>

The new continuous evaluator integrates the ordinary-prime residual between
prime jumps. These 36 finite windows show the signed absolute moment beside
its normalized projection diagnostic; the plotted segments connect evaluated
points. They establish no all-window growth law.
[Integral derivation](docs/EXACT_ARITHMETIC_MOMENTS.md) ·
[Decimal data and provenance](results/readme_figures_20260905/arithmetic_moments.json) ·
[SVG](docs/figures/continuous_arithmetic_moments.svg).

<p align="center">
  <img src="docs/figures/unattained_spectrum.png" alt="Genuine Laplace poles approach the unattained real-part boundary one quarter, while the proved normalized rate factor decays faster than any inverse power" width="1000"/>
</p>

The new analytic counterexample has genuine Laplace poles and arithmetic
residue coefficients, yet every fixed-degree normalized ratio decays faster
than every inverse power. The right panel illustrates the proved rate factor,
with its degree-dependent constant omitted; it is not a measured ratio.
No prime realization is asserted.
[Proof](docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md) ·
[Formula data](results/readme_figures_20260905/unattained_spectrum.json) ·
[SVG](docs/figures/unattained_spectrum.svg).

### Pipeline and earlier campaigns

<p align="center">
  <img src="docs/figures/pipeline.png" alt="PBSS pipeline: residual → Legendre projection → R_d scorecard" width="720"/>
</p>

<p align="center"><em>Pipeline — density residual on a log-window, projected onto low-degree shifted Legendre modes.</em></p>

<p align="center">
  <img src="docs/figures/scorecard_Rd.png" alt="R_d diagnostic scorecard: RH-like vs defective" width="640"/>
</p>

<p align="center"><em>Diagnostic scorecard — low <code>R<sub>d</sub></code> for high-frequency / critical-line-like probes; high <code>R<sub>d</sub></code> for defective controls.</em></p>

| Residual waveforms | Mode energy spectrum |
|:------------------:|:--------------------:|
| <img src="docs/figures/residuals_waveforms.png" alt="Four residual waveforms with R4 labels" width="420"/> | <img src="docs/figures/mode_spectrum.png" alt="Legendre mode energy: RH-like vs defective" width="420"/> |

| Grand campaign \(R_d\) vs \(T\) | Arithmetic soft plateau |
|:------------------------------:|:-----------------------:|
| <img src="docs/figures/grand_Rd_vs_T.png" alt="Grand campaign R_d versus T" width="420"/> | <img src="docs/figures/grand_arith_focus_linear.png" alt="Arithmetic residual R4 soft plateau" width="420"/> |

| Beurling battery | Historical status map |
|:----------------:|:----------:|
| <img src="docs/figures/beurling_battery_Rd_vs_T.png" alt="Beurling ordinary vs defective R_d vs T" width="420"/> | <img src="docs/figures/status_map.png" alt="Historical Theorem A/B status map; see current status table below" width="420"/> |

More figures: [`docs/figures/`](docs/figures/) (canonical; paper uses the same files).

## Quick start (tool)

```bash
pip install -e ".[dev]"
pbss version
pbss diagnose --demo
pbss project --input path/to/q.npy --degree 4 --T 20
pbss sensitivity --confirm-53
pbss scorecard --x-max 1e6
pbss moments --T 10 --degree 4 --json-out results/moments.json
```

Or without install: `PYTHONPATH=src python3 -m pbss …`

Tool guide: [`docs/TOOL.md`](docs/TOOL.md) · Sensitivity (≥53% gain): [`docs/MEASURE_SENSITIVITY.md`](docs/MEASURE_SENSITIVITY.md)

`pbss moments` evaluates signed absolute arithmetic moments and the continuous
norm from complete prime intervals. Its decimal output preserves the amplitude
needed for the [absolute-moment target](docs/ABSOLUTE_MOMENT_CONVERSE.md).
See the [integral identities and precision contract](docs/EXACT_ARITHMETIC_MOMENTS.md);
finite evaluations do not establish an asymptotic growth bound.

### Research campaigns (optional)

```bash
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_diagnostic.py
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
```

**Grand campaign** (shipped): `experiments/run_grand_campaign.py` — arithmetic residual to
\(x_{\max}=10^{10}\) (checkpointed primes), multi-\((d,\mathrm{detrend},\mathrm{smooth})\),
**≥2000 MC defect trials/T** (default 20k), CL/off-critical/defect controls, resume,
plots in `results/grand_campaign/`. Focus deg1: \(R_4\sim0.15\)–\(0.19\) through \(10^{10}\).
See `docs/STATUS.md`.

## Core math (shipped)

### Projection strength

1. Normalize a density perturbation \(q\) on the unit log-window \(u\in[0,1]\).
2. Project onto orthonormal **shifted Legendre** polynomials
   \(\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1)\).
3. Energy ratio
   \[
   R_d(q)=\frac{\|P_d q\|_{L^2}^2}{\|q\|_{L^2}^2}
   =\frac{\sum_{k=0}^d|\langle q,\varphi_k\rangle|^2}{\|q\|^2}.
   \]
4. Working projection strength (scaled)
   \[
   P(q)\;:=\;S_d(q)=T^{2(d+1)}\,R_d(q),
   \]
   with \(T\) the logarithmic window length. This is a scaling convention;
   generic model-mode decay is \(O(T^{-2})\), so \(S_d\) need not remain bounded.
   The archive's sharper \(O(T^{-2(d+1)})\) heuristic is not a general RH criterion.

- **Low** \(R_d\): little energy in the selected polynomial subspace.
- **High** \(R_d\): substantial energy in that subspace.

For sampled data, the shipped implementation fits the polynomial by weighted
least squares using trapezoid weights. Legendre polynomials are continuously
orthonormal, but their sampled Gram matrix need not be the identity. The energy
is therefore \(c^T G c\), with \(G=\Phi^T W\Phi\), rather than \(\sum c_k^2\).
This keeps the sampled energy fraction in \([0,1]\), including on nonuniform
grids. A low value can also hide an out-of-subspace defect or undersampling;
neither outcome classifies RH. See [`docs/TOOL.md`](docs/TOOL.md) for the input
contract and control/refinement checks.

### Legacy numbers

Earlier notes quoted **P(q)≈3.92** for zeta and threshold **≈29.5**. Those
used lost high-precision scripts. This reconstruction **does not hard-code
those values**; it reports \(R_d\) and \(S_d\) from the shipped path. See
`docs/STATUS.md`.

### Theorems A/B (precise status)

| Result | Status |
|--------|--------|
| **A₀** critical-line mode \(R_d\to0\) at \(O(T^{-2})\) | **Proved** (Lemma M3) |
| **Finite-mode A₀** finite CL sum \(R_d=O(T^{-2})\) | **Proved** (Lemma M5) |
| **Weighted model A₀** \(R_d(w q)=O(T^{-2})\) for \(w\in W_\alpha\) | **Proved** (Lemma M6) |
| **M7** \(R_d\) perturbation majorant | **Proved** |
| **A** exact continuous arithmetic residual under RH | **Proved conditionally**, with continuous affine detrending |
| **B₀** persistent defect \(\Rightarrow R_d=\varepsilon^2\not\to0\) | **Proved** (Lemmas M2–M4) |
| **Arithmetic recurrence** \(\liminf R_d=0\) | **Proved unconditionally**, for that continuous residual |
| **B-ABS** subexponential absolute moment \(\Leftrightarrow\) RH | **Proved equivalence**; unconditional growth bound remains open |
| **Signed moments** either one-sided subexponential bound \(\Leftrightarrow\) RH | **Proved equivalence**; bound remains open |
| **Normalized converse with an attained rightmost off-critical zero** | **Proved:** \(\limsup T R_d>0\) |
| **B** specified normalized decay \(\Rightarrow\) RH | **Open**; the positive-floor B-RES formulation is RH-equivalent |
| **RH** unconditional | **Open** |

The arithmetic statements use the exact continuous θ residual and orthogonal
degree-one detrending, with fixed \(d\ge2\). They do not identify a fixed sampled
array with its continuous limit. The proof package now includes
[continuous A](docs/CONTINUOUS_THEOREM_A.md),
[unconditional recurrence](docs/NORMALIZED_RECURRENCE_THEOREM.md), and the
[absolute-moment converse](docs/ABSOLUTE_MOMENT_CONVERSE.md).
The [exact off-critical model](docs/OFF_CRITICAL_PHASE_OBSTRUCTION.md) also has
\(R_d\to0\), with phase-dependent troughs. Its formerly claimed positive
floor was incorrect.
The [attained-boundary theorem](docs/ATTAINED_BOUNDARY_CONVERSE.md) now isolates
the remaining normalized case: a supremum of zero real parts that no zero
attains. The [signed criterion](docs/SIGNED_MOMENT_OSCILLATION.md) reduces the
absolute target to either one-sided bound, still unproved unconditionally.
An [explicit infinite spectrum](docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md) shows
why the remaining normalized case needs additional arithmetic input: its
genuine poles approach an unattained boundary while every fixed-degree ratio
decays faster than any inverse power, even with the arithmetic residue
coefficients. This is an analytic spectrum example with no prime realization
asserted.

Explicit-formula residual: `pbss.probes.explicit_formula_residual` · peel scan:
`experiments/run_explicit_formula_peel.py` → `results/explicit_formula_peel/`.

Open-plateau: [`docs/RESEARCH_PLATEAU.md`](docs/RESEARCH_PLATEAU.md) · Theorem-A scaffold:
[`docs/THEOREM_A_SCAFFOLD.md`](docs/THEOREM_A_SCAFFOLD.md) · Tail remainder note:
[`docs/INFINITE_TAIL_REMAINDER.md`](docs/INFINITE_TAIL_REMAINDER.md).

Details: [`docs/THEOREMS_AB.md`](docs/THEOREMS_AB.md) · Proofs: [`docs/PROOFS_LEMMAS.md`](docs/PROOFS_LEMMAS.md) · Status: [`docs/STATUS.md`](docs/STATUS.md) · Paper: [`docs/paper/`](docs/paper/).

**Not an unconditional RH proof.**

## Repository layout

```
src/pbss/           # installable library
tests/              # pytest (real shipped paths)
experiments/        # campaign entry scripts
results/            # campaign JSON/TXT/plots (large primes gitignored)
docs/               # index: docs/README.md
docs/figures/       # canonical images (README + paper)
docs/paper/         # status note TeX/PDF only
docs/goals/         # campaign briefs (plateau / marathon)
docs/related/       # third-party reference PDFs
```

## Limitations

- Diagnostic, not a decisive RH proof.  
- Finite windows cannot exclude extremely high zeros.  
- Not a practical local primality sieve.  
- Large \(d\) or \(T\) is expensive.  
- Converse (low energy ⇒ RH) is essentially as hard as RH.

## Related work

### External (zeta / RH — reference only)

- **Anthropic (Aug 2026):** improved lower bound on the fraction of \(\zeta\) zeros on the critical line (~41.6% → ~67.2%).  
  Blog: [anthropic.com/research/riemann-zeta](https://www.anthropic.com/research/riemann-zeta) ·  
  Local PDFs + notes: [`docs/related/anthropic-riemann-zeta/`](docs/related/anthropic-riemann-zeta/)  
  (Third-party; not part of PBSS; different theorem class — see that README.)

### Jensen / de Bruijn–Newman (same author, pre-PBSS survivor)

- **Jensen/moment hierarchy blindness for \(\Lambda\):** central certificates falsely look hyperbolic down to \(t=-0.7\) while the binding Lehmer pair needs moment order \(\sim 13400\).  
  Note: [`docs/JENSEN_MOMENT_HIERARCHY_BLINDNESS.md`](docs/JENSEN_MOMENT_HIERARCHY_BLINDNESS.md) · code: `pbss.jensen_blindness`  
  (Not RH; not a new sharp \(\Lambda\) bound beyond optional bookkeeping.)

### Post-Jensen roadmap ranks 2–5 (executed)

| Rank | Deliverable |
|-----:|-------------|
| 2 | Plateau / secondary EF: [`docs/ARITHMETIC_PLATEAU_SECONDARY.md`](docs/ARITHMETIC_PLATEAU_SECONDARY.md) · `plateau_secondary_report` |
| 3 | ANT interface audit (freeze Full A): [`docs/ANT_INTERFACE_AUDIT.md`](docs/ANT_INTERFACE_AUDIT.md) · `ant_interface_audit` |
| 4 | Zero-proportion feasibility **STOP**: [`docs/ZERO_PROPORTION_FEASIBILITY.md`](docs/ZERO_PROPORTION_FEASIBILITY.md) |
| 5 | B-RES threshold \(H^*\): [`docs/B_RES_THRESHOLD.md`](docs/B_RES_THRESHOLD.md) · `b_res_threshold_report` |

### By the same author

- Other number-theory / exploratory archives (e.g. `perry-spirals`, `wieferich-hunts`) may appear under the same GitHub account; they are separate projects and not required to run PBSS.

---

*Independent research. MIT licensed. Not a proof of the Riemann Hypothesis.*
