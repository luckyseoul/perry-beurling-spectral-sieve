# Perry–Beurling Spectral Sieve  
## and the RH Research Program That Preceded It

**Research history, mathematical development, computational campaigns, failures, and current status**

Nicholas D. Perry | Independent Researcher  
Compiled 11 August 2026

**Scope note.** This is a research-history document, not a claim of a proof of the Riemann Hypothesis. It deliberately preserves abandoned branches and negative results because they explain how the present PBSS formulation emerged.

**Related earlier record.** The full Nov 2025 – Jun 2026 project record (Spectral Sieve operator framework **killed** as μ-weighted \(L^2\); Wronskian / OT / Li kills; de Bruijn–Newman survivors including Jensen-hierarchy blindness and \(\Lambda\le 0.20\)) lives in [`PROJECT_RECORD_2025-11_to_2026-06.md`](PROJECT_RECORD_2025-11_to_2026-06.md). This file focuses on how the **projection-diagnostic** strand was reconstructed as PBSS in July 2026 and what that strand currently claims.

**External concurrent work (not PBSS).** Anthropic (Aug 2026) reported an improved unconditional lower bound on the fraction of zeros of \(\zeta\) on the critical line (~41.6% → ~67.2%). That is a different theorem class from PBSS’s residual projection diagnostic. Reference copies and links: [`related/anthropic-riemann-zeta/`](related/anthropic-riemann-zeta/).

---

## Executive Summary

The Perry–Beurling Spectral Sieve (PBSS) is the surviving, formalized branch of a longer sequence of experiments aimed at converting the Riemann Hypothesis into a measurable property of prime-density residuals. The central idea is not to locate zeros directly, but to ask whether the normalized prime-counting error has persistent low-complexity structure when viewed on expanding logarithmic windows.

The program went through several materially different phases. Early work explored apparent spectral and geometric regularities, including a Perry-Gamma / macro-resonance line of attack. Those claims were explicitly discarded in December 2025 after adversarial review exposed measurement dependence and numerical artifacts. A later spectral diagnostic based on the normalized Chebyshev error survived longer: project a residual onto a low-degree polynomial space and measure how much energy remains there. That diagnostic produced the legacy \(P(q)\approx 3.92\) observation and an intended forward/converse RH criterion, but its original high-precision implementation was lost and the strongest asymptotic rate was not justified.

PBSS, reconstructed in July 2026, is the disciplined version of that surviving idea. It fixes the window, basis, projection operator, energy ratio, model classes, theorem labels, computational controls, and non-claims. It then uses Beurling generalized-prime systems as adversarial synthetic universes: ordinary-prime-like systems should look spectrally high-frequency, while deliberately defective systems should retain low-degree mass. In current tests, that classifier separates the controls extremely well, but the true arithmetic residual does not yet exhibit the model decay required for a theorem.

The strongest current mathematical results are **model** results: critical-line modes and finite sums of them have low-degree projection ratio \(R_d=O(T^{-2})\); fixed low-degree defects prevent \(R_d\) from vanishing; admissible endpoint weights preserve the model \(O(T^{-2})\) behavior. The full arithmetic implication under RH and the converse implication back to RH remain open.

### At a glance

| Phase | Approx. date | Outcome |
|-------|--------------|---------|
| Early resonance / Perry-Gamma | 2025, discarded Dec. 8–9 | Abandoned as artifact / measurement-dependent |
| Spectral diagnostic \(P_d(q_T)\) | Late 2025 | Core projection idea survives; strongest claimed rate and legacy normalization not retained |
| Beurling generalization / PBSS reconstruction | July 2026 | Formal diagnostic, proofs M1–M6, reproducible code and controls |
| Large arithmetic + Beurling campaigns | July 26–29, 2026 | Strong control separation; arithmetic \(R_4\) soft plateau ~0.15–0.19 |
| Current status | Aug. 2026 | Useful classifier / research instrument; RH and full Theorems A/B open |

---

## 1. Research Lineage Before PBSS

### 1.1 The broader objective

The recurring objective was to find a representation in which the critical-line condition \(\mathrm{Re}(s)=1/2\) becomes a structural signature rather than a list of individually verified zeros. The hoped-for route was: transform prime-counting error into a signal, identify a property forced by critical-line oscillations, then show that an off-critical contribution necessarily violates that property.

This is closely aligned with the explicit-formula viewpoint: zeros contribute oscillatory terms to prime-counting discrepancies. The research question was whether a coarse but robust spectral statistic could distinguish a signal made only from critical-line modes from one containing an off-critical or otherwise defective component.

### 1.2 Early macro-resonance / Perry-Gamma branch

An earlier branch attempted to identify macro-scale resonances and a Perry-Gamma correspondence. It included proposed special constants and apparent isomorphisms between numerical structures. By December 8–9, 2025, adversarial review showed that the macro-resonance quantities depended on the measurement procedure and that the associated Perry-Gamma / \(\alpha=|a^*|\) claims were not reliable. The project was explicitly reset and those claims were discarded.

**Historical significance:** this failure changed the methodology. Later RH work increasingly required explicit controls, invariant definitions, synthetic defect systems, reproducible code paths, and a hard distinction between model lemmas, empirical observations, conditional arithmetic statements, and actual RH claims.

### 1.3 The late-2025 spectral diagnostic

The next important branch used a normalized Chebyshev-type residual on a logarithmic window and projected it onto low-degree polynomials. In the archived formulation, a representative residual was

\[
q_T(s) = \frac{\psi(x)-x}{\sqrt{x}},
\]

with the \(x\)-domain mapped to \(s\in[0,1]\).

The diagnostic \(P_d(q_T)\) measured low-degree polynomial content. The intended intuition was simple: critical-line zero contributions become increasingly oscillatory as the log-window grows, so a fixed low-dimensional polynomial space should capture progressively less of the signal. A persistent off-critical contribution was expected to leave detectable low-frequency or low-degree structure.

At this stage, two ambitious claims were entertained: a forward implication under RH with a rate quoted as \(O(T^{-2(d+1)})\), and a converse in which sufficiently fast decay of the diagnostic would imply RH. Numerical work also produced a legacy zeta value near \(P(q)\approx 3.92\) and a threshold near 29.5. These values came from high-precision scripts that were later lost; PBSS does not hard-code or treat them as established.

The late-2025 work also used Diamond/Beurling generalized-prime systems as stress tests. That was the conceptual bridge to PBSS: instead of testing only the ordinary primes, construct generalized-prime worlds whose analytic behavior can be intentionally made RH-like or defective and ask whether the diagnostic classifies them correctly.

---

## 2. Why Beurling Generalized Primes Matter

Beurling generalized-prime systems provide controlled counterfactual arithmetic. They allow the prime sequence or density law to be altered while retaining enough zeta-like structure to test whether a proposed RH diagnostic is responding to the intended analytic feature rather than to incidental properties of the ordinary primes.

For this research program, Beurling systems serve the role that adversarial and ablation datasets serve in machine learning. A useful instrument should respond correctly not only to the target data but also to deliberately malformed controls.

- **Ordinary / RH-like control:** the statistic should place little energy in the fixed low-degree subspace as oscillation increases.
- **Gapped or thinned generalized primes:** engineered density defects should create persistent low-degree mass.
- **Off-critical model modes:** exponential envelope \(e^{T(\sigma-1/2)u}\) should become increasingly distinguishable from \(\sigma=1/2\).
- **Orthogonal low-degree defects:** a known injected component provides an exactly solvable calibration target.

---

## 3. PBSS: Formal Reconstruction

The repository `luckyseoul/perry-beurling-spectral-sieve` describes PBSS as a projection diagnostic for density residuals on logarithmic windows and explicitly labels it a classifier/diagnostic rather than an RH proof. The reconstruction dates to July 2026.

### 3.1 Logarithmic window

\[
x = e^{uT},\qquad u\in[0,1].
\]

\(T\) is the logarithmic window length. Pulling the arithmetic residual back to the unit interval separates window growth from the fixed projection geometry.

### 3.2 Orthonormal shifted-Legendre basis

\[
\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1).
\]

Let \(V_d=\mathrm{span}\{\varphi_0,\ldots,\varphi_d\}\), and let \(P_d\) be orthogonal projection onto \(V_d\) in \(L^2([0,1])\). This makes “low-degree content” a precise, basis-controlled quantity.

### 3.3 Core energy statistic

\[
R_d(q)=\frac{\|P_d q\|^2}{\|q\|^2}
=\frac{\sum_{k=0}^d|\langle q,\varphi_k\rangle|^2}{\|q\|^2},\qquad 0\le R_d\le 1.
\]

\[
S_d(q;T)=T^{2(d+1)}R_d(q),\qquad\text{working notation }P(q):=S_d(q;T).
\]

\(R_d\) is the robust core observable. The scaled \(S_d/P(q)\) retains the historical scaling convention, but the repository explicitly separates it from the lost legacy \(P(q)\approx 3.92\) normalization. Later stress tests also showed that the stronger \(O(T^{-2(d+1)})\) rate should not be treated as proved or empirically sharp.

### 3.4 Model residuals

| Model | Form | Purpose |
|-------|------|---------|
| Critical-line mode | \(\sin(tTu)\) | Ideal contribution from \(\tfrac12+it\) |
| Off-critical mode | \(e^{T(\sigma-1/2)u}\sin(tTu)\) | Directional non-RH control |
| Fixed low-degree defect | \(\sqrt{1-\varepsilon^2}\,f+\varepsilon\varphi_j\) | Exactly calibrated persistent defect |
| Finite CL sum | \(\sum a_n\sin(t_nTu+\phi_n)\) | Truncated explicit-formula model |
| Arithmetic residual | \((\theta(x)-x)/\sqrt{x}\) | Actual prime-data target |

---

## 4. What Is Actually Proved in PBSS

### 4.1 M1: pure basis-mode energy

If \(q=\varphi_m\), then \(R_d(q)=1\) for \(m\le d\) and \(R_d(q)=0\) for \(m>d\).

This is the basic calibration of the projection instrument.

### 4.2 M2: exact orthogonal-defect formula

\[
q=\sqrt{1-\varepsilon^2}\,f+\varepsilon\varphi_j,\quad f\perp V_d,\ \|f\|=1,\ j\le d
\quad\Rightarrow\quad R_d(q)=\varepsilon^2.
\]

This gives an exact synthetic defect whose strength is known before the experiment.

### 4.3 M3: critical-line mode decay

\[
q_T(u)=\sin(tTu)\quad\Rightarrow\quad R_d(q_T)=O_d(T^{-2}).
\]

This is the rigorous model rate. It replaces the earlier stronger archive heuristic \(O(T^{-2(d+1)})\) as the rate that can presently be defended for the basic mode.

### 4.4 M4: persistent defects block vanishing

If \(\varepsilon\) remains bounded away from zero in the M2 construction, then \(R_d\) remains bounded away from zero regardless of how oscillatory the orthogonal part becomes. Thus \(R_d\to 0\) requires the disappearance of persistent low-degree mass.

### 4.5 M5: finite critical-line superpositions

\[
q_T^{(N)}(u)=\sum_{n=1}^N a_n\sin(t_nTu+\phi_n)
\quad\Rightarrow\quad R_d(q_T^{(N)})=O_d(T^{-2}),
\]

for fixed finite \(N\).

M5 is the key bridge from a single ideal zero contribution to a truncated explicit-formula signal. It does not control the infinite zero sum or the arithmetic remainder.

### 4.6 M6: admissible weights

Endpoint tapering became important because arithmetic residuals are sensitive to the ends of finite windows. PBSS formalized a weight class \(W_\alpha\) and proved, at the model level, that admissible weights preserve the \(O(T^{-2})\) decay for critical-line and finite explicit-formula residuals.

---

## 5. Theorem Program: A, B, and Their Model Versions

| Statement | Meaning | Status |
|-----------|---------|--------|
| **A₀** | Single critical-line mode has vanishing low-degree energy | **Proved** (M3) |
| **Finite A₀** | Finite CL sum has \(R_d=O(T^{-2})\) | **Proved** (M5) |
| **Weighted A₀** | Admissible weighted model retains \(O(T^{-2})\) | **Proved** (M6) |
| **A** | Under RH, arithmetic residual has \(R_d\to 0\) | **Open**; conditional package/scaffold complete |
| **B₀** | Persistent low-degree defect prevents \(R_d\to 0\) | **Proved** (M2–M4) |
| **B** | Sufficiently fast arithmetic decay implies RH | **Open**; essentially RH-hard |

The central unresolved bridge is **arithmetic identification**: proving that the true \(\theta\) or \(\psi\) residual, including infinitely many zeros, secondary main terms, truncation error, and endpoint behavior, falls into the model class strongly enough for A. The converse requires an equally rigorous reduction showing that any off-critical zero forces a nonvanishing or otherwise forbidden projection signature.

---

## 6. Computational Program

### 6.1 Arithmetic multi-T and grand campaign

PBSS constructed \(q_T=(\theta(x)-x)/\sqrt{x}\) over logarithmic windows and ran multi-T scans with different degrees, detrending choices, and smoothing. A grand campaign reached \(x_{\max}=10^{10}\) with approximately 455 million checkpointed primes and 20,000 Monte Carlo defect trials per \(T\) over 14 windows.

| \(T\) | \(x_{\max}\) | Arithmetic \(R_4\) (deg1) | CL \(R_4\) | Defect \(R_4\) | MC mean \(R_4\) |
|------:|-------------:|--------------------------:|-----------:|---------------:|----------------:|
| 10 | \(2\times10^4\) | 0.144 | ~0.006 | 0.250 | 0.793 |
| 16 | \(9\times10^6\) | 0.193 | ~0.002 | 0.250 | 0.793 |
| 23 | \(10^{10}\) | 0.155 | ~0.0005 | 0.250 | 0.795 |

The important **negative result** is the arithmetic soft plateau: \(R_4\) stayed roughly 0.15–0.19 instead of collapsing toward the ideal finite-mode levels. This does not refute RH. It shows that the current arithmetic observable contains substantial low-degree structure not explained by the simplest truncated zero model.

### 6.2 Extended-x campaign

The prime checkpoint was extended to \(x_{\max}=5\times10^{10}\), containing 2,119,654,578 primes. The parallel extension used 86 workers. At \(T\approx 24.6\) the focus \(R_4\) value remained about 0.146, so the soft plateau persisted. A \(10^{11}\) run was rejected on memory grounds and \(10^{12}\) was considered infeasible in the available environment.

### 6.3 Arithmetic zero peeling

The project explicitly subtracted the first \(N\) critical-line modes from the arithmetic residual, optionally fitting a scale \(\alpha\). If the plateau were mostly caused by a few known low zeros, peeling them should have driven \(R_d\) toward the finite-mode model. It did not. Through \(10^{10}\), arithmetic \(R_d\) remained far above A₀ levels.

This was a useful falsification of an easy explanation: the stubborn low-degree mass was not simply the visible contribution of the first handful of critical-line zeros.

### 6.4 Explicit-formula identification attack

A focused comparison found the best residual convention to be \(H_{\theta,\sqrt{}}=(\theta-x)/\sqrt{x}\). Correlation between the arithmetic residual \(q\) and a truncated zero model \(m\) reached mean \(|\mathrm{corr}(q,m)|\approx 0.65\) and about 0.72 at \(N=40\), while the \(L^2\) remainder fraction improved with \(N\).

However, a sharper invariant emerged: \(E_d(r)/\|q\|^2\approx 0.21\) stayed nearly flat in \(N\) even as total \(L^2\) capture improved. In other words, adding more zero modes explained more of the signal without explaining the stubborn low-degree component.

### 6.5 “Kill 0.21” enrichment

| Model enrichment | Mean low-degree remainder fraction |
|------------------|-----------------------------------:|
| Zeros only | ~0.207 |
| Zeros + smooth exponential/secondary mains | ~0.052 |
| Zeros + high-Legendre / free \(V_d\) span | ~0 |

This was genuine headway but also a warning. Smooth secondary terms cut the unexplained low-degree mass by roughly fourfold, showing that the explicit-formula model had omitted important nonzero secondary structure. Allowing the model to span \(V_d\) directly can trivially kill the diagnostic, so that route is not evidentially meaningful unless the added terms are independently justified.

---

## 7. Beurling Battery and Adversarial Validation

### 7.1 Initial battery

At \(x_{\max}=10^8\) and \(T\in\{8,10,12,14,16,18\}\), degree 4, the ordinary-prime control separated strongly from deliberately gapped and thinned systems.

| System | \(R_4\) at \(T=18\) |
|--------|--------------------:|
| ordinary_primes | ~0.181 |
| gapped_gap3 | ~0.983 |
| thinned_every3 | ~0.990 |

This established that the instrument can detect large engineered density defects.

### 7.2 Expanded tight-stress battery

On July 29, 2026 the battery was expanded to 35 systems at \(x_{\max}=10^6\). Ordinary \(R_4\) was about 0.189, while defective systems had min/median/max about 0.993/0.995/0.996. No failure systems or thin-margin cases appeared at that scale.

This is among the strongest empirical results in the project: PBSS is a very effective classifier for the defect families it was designed to distinguish. It is not evidence that the same separation theorem holds for all Beurling systems or for every possible off-critical zeta configuration.

---

## 8. Monte Carlo and Failure-Mode Stress

The project ran increasingly large Monte Carlo campaigns, including 50 million trials in the open-plateau work, 1.6 million trials in the overnight marathon, and a 192,000-trial focused ablation session. These tests were used to understand the statistic itself rather than to estimate a probability that RH is true.

| Ablation | Mid-T mean \(R_4\) |
|----------|-------------------:|
| baseline | 0.792 |
| heavy_defect | 0.970 |
| light_defect | 0.285 |
| high_freq | 0.792 |
| low_freq | 0.794 |
| high_deg_defect | 0.496 |

A crucial failure mode was discovered: a “high-degree defect” can lower \(R_d\) if its mass lies outside \(V_d\). PBSS therefore detects **low-degree defects relative to the chosen projection space**, not arbitrary defects. This limits any claim that a single fixed \(d\) is a universal classifier.

---

## 9. Off-Critical Model Tests

For model modes with \(\sigma>1/2\), the exponential envelope increasingly distinguishes them from critical-line modes as \(T\) grows. In the July 29 sweep, the \(R_4\) ratio between \(\sigma=0.9\) and \(\sigma=1/2\) grew from about 3.9 at \(T=8\) to about 15.9 at \(T=32\).

This is directionally consistent with the intended converse mechanism, but it is not a zero-free-region theorem. The missing step is a rigorous arithmetic reduction showing that an actual off-critical zero necessarily produces a surviving signature in the chosen normalized residual after all other terms and cancellations are accounted for.

---

## 10. Open-Plateau Campaign

The open-plateau campaign deliberately attacked the arithmetic plateau from multiple directions: zero peeling, whitening, measure choice, basis choice, scale, Monte Carlo randomization, and enlarged Beurling tests. Deep axes included 50 million MC trials, 24 \(T\) values across six scaling variants on the \(5\times10^{10}\) prime table, and at least 500 Beurling systems.

The synthesis was negative but clarifying: peeling did not collapse the arithmetic statistic; degree-1 detrending plus optional taper was the best residual recipe; and dense multi-T scans showed no A₀-like decay through \(x\approx 5\times10^{10}\). That result forced the project away from interpreting finite-window arithmetic behavior as a near-proof.

---

## 11. What Survived from the Pre-PBSS Work

- The core projection idea: convert prime-density error into a signal and quantify its low-complexity energy.
- Logarithmic windows: zero terms naturally oscillate in \(\log x\), making \(T\) the meaningful expansion parameter.
- Beurling systems as adversarial controls: generalized primes are central, not decorative.
- Forward/converse architecture: critical-line-only signals should lose low-degree mass; persistent off-critical structure should obstruct that loss.
- The value of synthetic defects: exact injected components make the instrument calibratable.
- The need for aggressive falsification: artifacts and easy explanations are treated as targets to kill, not results to preserve.

---

## 12. What Was Retracted, Weakened, or Reclassified

| Earlier item | Current treatment |
|--------------|-------------------|
| Perry-Gamma / macro-resonance claims | Discarded in Dec. 2025 as artifact / measurement-dependent |
| \(\alpha=|a^*|\) and associated resonance draft | Discarded with the above branch |
| Legacy \(P(q)\approx 3.92\) and threshold \(\approx 29.5\) | Historical only; original high-precision path lost |
| \(R_d=O(T^{-2(d+1)})\) as general rate | Not proved; model theorem supports \(O(T^{-2})\) |
| Fast decay \(\Rightarrow\) RH as established converse | Open |
| Arithmetic RH \(\Rightarrow\) observed finite-window A₀ decay | Open; current arithmetic data plateau |
| Beurling separation as proof of RH | Explicitly rejected; classifier validation only |

---

## 13. Current Research Boundary

PBSS currently has a clean internal boundary between what is known and what is hoped for.

- **Known in the model:** fixed-degree projection suppresses increasingly oscillatory critical-line modes at \(O(T^{-2})\).
- **Known in the model:** a fixed low-degree defect leaves a fixed nonzero \(R_d\).
- **Known for finite sums and admissible weights:** the \(O(T^{-2})\) model behavior survives.
- **Observed:** ordinary-prime and engineered Beurling defect families separate strongly in current batteries.
- **Observed:** the true arithmetic residual has a persistent low-degree plateau over all currently computed scales.
- **Open:** control of the infinite explicit-formula tail and the exact \(\theta/\psi\) remainder in the required norm.
- **Open:** a justified arithmetic decomposition that removes secondary low-degree terms without baking the answer into the model.
- **Open:** a converse showing that every off-critical zero forces a forbidden asymptotic PBSS signature.
- **Open:** RH itself.

---

## 14. Methodological Significance

Even without an RH proof, the project produced a useful research instrument and a much stronger methodology than the early resonance work. PBSS is falsifiable, parameterized, reproducible, and equipped with negative controls. Its failures are informative: the arithmetic plateau identifies precisely where the simple “critical-line modes become high-frequency” story stops being sufficient.

The Beurling component is especially important because it turns RH-inspired intuition into a family of classification experiments. A future theorem would need to explain why the empirical separation is structural, which generalized-prime pathologies evade it, and which norm/window/basis choices make the forward and converse statements invariant rather than instrument-specific.

---

## 15. Suggested Next Mathematical Bottlenecks

Based on the present repository, the work is bottlenecked less by additional brute-force computation than by analytic identification.

1. Close the arithmetic residual identity in the exact weighted \(L^2\) space used by PBSS, including secondary main terms and endpoint contributions.
2. Replace scaffold tail majorants with genuine analytic-number-theory bounds strong enough to pass from finite M5/M6 sums to the infinite explicit formula.
3. Characterize which off-critical terms necessarily inject energy into \(V_d\) after normalization, and whether \(d\) must grow with \(T\).
4. Construct adversarial Beurling systems specifically designed to fool fixed-\(d\) PBSS, rather than only obvious gap/thinning defects.
5. Determine whether a basis-independent low-frequency functional can preserve the useful separation while eliminating the high-degree-defect loophole.
6. Treat the arithmetic plateau as a phenomenon to explain analytically before interpreting larger-\(x\) runs as evidence for or against the theorem program.

---

## 16. Chronology

| Date / period | Development |
|---------------|-------------|
| 2025 | Exploratory resonance / Perry-Gamma and other spectral-geometric RH work |
| Dec. 8–9, 2025 | Macro-resonance, Perry-Gamma, \(\alpha=|a^*|\) branch explicitly discarded |
| Late 2025 | Normalized Chebyshev-error polynomial projection diagnostic; legacy \(P(q)\approx 3.92\); Beurling/Diamond controls |
| 2026-07 | PBSS reconstructed as reproducible open repository |
| 2026-07-26 | Marathon, extended-x, Beurling battery, MC stress, theorem scaffolding |
| 2026-07-29 | Tight stress: formalized status, 192k ablations, 35-system Beurling battery, off-critical sweep |
| 2026-08-11 | Current synthesis: model diagnostic established; arithmetic A, converse B, and RH remain open |

---

## 17. Source and Provenance Notes

Primary current source: GitHub repository `luckyseoul/perry-beurling-spectral-sieve`, especially `README.md`, `docs/STATUS.md`, and `docs/THEOREMS_AB.md` as inspected on 11 August 2026. The repository identifies itself as a July 2026 reconstruction of the earlier spectral diagnostic / \(P(q)\) framework.

Historical pre-PBSS material is reconstructed from prior research conversations and retained project history. In particular, the December 2025 reset explicitly discarded the Perry-Gamma / macro-resonance claims. Because the original high-precision scripts behind the legacy \(P(q)\approx 3.92\) and \(\approx 29.5\) threshold are not present in the current repository, those values are recorded only as historical observations, not reproduced results.

Repository status language is intentionally preserved: PBSS is not an unconditional proof of RH. Full arithmetic Theorem A, full converse Theorem B, infinite-zero control, and RH remain open.

---

## Appendix A. Compact Mathematical Map

| Layer | Object | Question |
|-------|--------|----------|
| Arithmetic | \(\theta(x)\), \(\psi(x)\), prime/generalized-prime systems | What discrepancy signal is generated? |
| Normalization | \(q_T\) on \(x=e^{uT}\) | How should growth and endpoints be removed? |
| Projection | \(P_d\) onto shifted Legendre \(V_d\) | How much low-degree energy remains? |
| Statistic | \(R_d=\|P_dq\|^2/\|q\|^2\) | Does low-complexity mass vanish? |
| Model forward | CL modes / finite sums | Proved \(O(T^{-2})\) |
| Arithmetic forward | True \(q_T\) under RH | Open |
| Model obstruction | Injected low-degree defect | Proved nonvanishing |
| Arithmetic converse | Off-critical zero \(\Rightarrow\) forbidden PBSS behavior | Open |
| Adversarial test | Beurling systems | Empirically strong separation for tested families |

---

## Appendix B. Repository Artifacts Mentioned in the Current Program

- `experiments/run_arithmetic_zero_peel.py`
- `experiments/run_beurling_battery.py`
- `experiments/run_mc_stress.py`
- `experiments/run_extend_x_scan.py`
- `experiments/run_explicit_formula_peel.py`
- `experiments/run_open_plateau.py`
- `experiments/run_theorem_a_scaffold.py`
- `experiments/run_tight_stress_20260729.py`
- `docs/PROOFS_LEMMAS.md`
- `docs/THEOREMS_AB.md`
- `docs/THEOREM_A_SCAFFOLD.md`
- `docs/THEOREM_A_PACKAGE.md`
- `docs/INFINITE_TAIL_REMAINDER.md`
- `docs/RH_CLOSEOUT_ROADMAP.md`
- `docs/RESEARCH_PLATEAU.md`
- `docs/EF_IDENTIFY_ATTACK.md`
- `docs/KILL021_ENRICH_M.md`
