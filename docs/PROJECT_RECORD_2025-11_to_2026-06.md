# The Spectral Sieve and Related Riemann Hypothesis Work

**Complete project record, November 2025 – June 2026**

Nicholas Perry · Perry Brands LLC · Council Bluffs, Iowa

> **Repo placement note (Aug 2026).** This document covers the **pre-July 2026** arc: the original Spectral Sieve operator framework (killed), Wronskian / OT / Li lines, and the de Bruijn–Newman survivors. The **July 2026 reconstruction** of a projection diagnostic under the name PBSS — model lemmas M1–M6, arithmetic plateau campaigns, open A/B — is a **later branch** documented in [`RESEARCH_HISTORY.md`](RESEARCH_HISTORY.md) and the live code in this repository. Read both for the full lineage; they are not the same phase.

---

## Executive summary

Five distinct lines of attack on the Riemann Hypothesis were built, tested, and resolved
over roughly seven months. Four were killed on rigorous mathematical grounds. One
produced a genuine, narrow, defensible result.

| Line of work | Period | Disposition |
|---|---|---|
| Spectral Sieve (Perry–Beurling) | Nov 2025 – Jun 2026 | **Killed.** Collapses to a μ-weighted L² norm |
| Wronskian positivity criterion | Jan 2026 | **Killed.** Logically equivalent to RH; no wedge |
| Wasserstein-2 / optimal transport framing | Jan 2026 | **Killed as a proof route.** Real negative result retained |
| Li's criterion / positivity cascade | Jan 2026 | **Killed.** Weaker than direct zero verification |
| de Bruijn–Newman constant Λ | Jun 2026 | **Two real outputs.** Jensen-hierarchy blindness; Λ ≤ 0.20 |

**What survived, in order of value:**

1. **Jensen/moment hierarchy blindness.** A precise, rigorous demonstration that the
   central moment/Jensen polynomial hierarchy is structurally blind to the constraint
   that actually binds Λ, with an explicit index argument. The cleanest publishable
   thing produced across the whole arc.
2. **Λ ≤ 0.20.** The published Polymath15 bound of 0.22 updated by feeding
   Platt–Trudgian's rigorous verification height 3×10¹² into the Ki–Kim–Lee dynamical
   inequality. Bookkeeping, not new mathematics, but correct and worth stating.
3. **Montgomery-blindness theorem.** No 1D spacing statistic can imply RH; the gap is
   exactly the off-line mass M_σ. A negative result with teeth.
4. **Gamma-shaped measure sensitivity.** μ(ds) = s^(k−1)e^(−σs) ds gives roughly 53%
   better sensitivity than unweighted L² for off-line zero detection in Beurling
   systems. Good engineering, not deep mathematics. Standalone note written.

**The meta-lesson, established across all five:** diagnostic and equivalence strategies
relocate the difficulty rather than reduce it. Every "new criterion" turned out to be RH
in disguise, or a detector for something already easy to detect. The only category that
moves a needle is a statement admitting partial, quantitative, *unconditional* progress —
which is why Λ was the one target where the work landed anywhere real.

---

# Part I — The Spectral Sieve

## 1. Core idea

A continuous analogue of a prime sieve. Instead of eliminating integers discretely, it
measures **spectral stiffness** — how resistant a system's harmonic structure is to
perturbation. The reframing:

> RH ⇔ vanishing of a measurable spectral defect 𝒟_k(t) for all t.

A generalized prime system's density deviation is encoded as a perturbation q(s) on
[0,1]. Two weighted operators are built over a μ-orthonormal basis — a baseline and a
perturbed version — and the operator norm of their difference is the diagnostic.

## 2. Construction

**Measure.**

$$\mu(ds) = s^{\,k-1} e^{-\sigma s}\, ds, \qquad s \in [0,1],\; k \ge 2,\; \sigma > 0$$

**Basis.** {ψᵢ} built by weighted Modified Gram–Schmidt against μ, so the Gram matrix is
the identity. Legendre polynomials were the working choice.

**Operators.**

$$J^0_{ij}(t) = \int_0^1 \psi_i(s)\,\psi_j(s)\cos(ts)\,\mu(ds)$$

$$J^\varepsilon_{ij}(t) = \int_0^1 \psi_i(s)\,\psi_j(s)\cos(ts)\,[1 + \varepsilon q(s)]\,\mu(ds)$$

**Defect.**

$$\mathcal{D}_k(t) = \bigl\| J^\varepsilon(t) - J^0(t) \bigr\|_{\mathrm{op}}$$

**Predictions.** Baseline 𝒟_k ≈ 0; perturbed 𝒟_k ~ ε·t² for small t.

**Classification rule.** Ratio of area under the defect curve to a baseline. Ratio below
1.5 → "RH-like".

**Cost.** O(d³) per t after 1D quadrature — trivially cheap, which was much of the appeal.

## 3. What differentiated it from a standard Beurling sieve

Beurling's 1937 generalized primes are a *construction* tool: pick any increasing
sequence 1 < p₁ ≤ p₂ ≤ …, build "integers" as products, study the associated zeta
function. Output is binary — the system satisfies RH or it doesn't.

The sieve's intended addition was a *measurement*: encode the density deviation as q(s),
build the two operators, take the norm of their difference, and get a continuous number
in [0,∞) instead of a yes/no. Beurling tells you *if*; the sieve was supposed to tell you
*how much*.

That distinction is the whole claim to novelty, and it is exactly what failed to survive
scrutiny — see §7.

## 4. Timeline

### Nov 2025 — Construction and validation

Three validation tests, all passed at machine precision:

| Test | Description | Result |
|---|---|---|
| 1 | Subspace invariance under basis rotation | Passed, 10⁻¹⁶ precision |
| 2 | Operator-orthogonal null suppression | Passed, ~100× suppression |
| 3 | Adversarial robustness / power-law scaling | Passed, α = 1.98 ± 0.01, R² > 0.9999 |

The quadratic law 𝒟_k ~ ε·t² held across bases, coordinate reparametrizations, and
perturbation families.

### Nov–Dec 2025 — Beurling and L-function classification

**Diamond's Beurling systems.** 5/5 correct on the initial set, 7/7 on the extended set.
Separation on the three-way synthetic comparison:

| System | Max 𝒟_k | AUC | Behavior |
|---|---|---|---|
| Regular (RH-like) | 0.12 | 1.8 | High-frequency modulation |
| Irregular (non-RH-like) | 0.31 | 5.2 | Low-frequency bias |
| Spike (strongly non-RH-like) | 0.38 | 7.1 | Localized excess |

Roughly 3× contrast by peak defect or AUC on this set; 3–17× across the broader sweep.

**Dirichlet L-functions.** ~1,150 characters tested, conductors from q = 5 up to ~10¹⁵.
99.75% classified RH-like. (Later verbal recaps quoted "up to conductor 10⁵" — the
primary session record is the larger range, and the anomaly investigated below sits at
q ≈ 1.7×10⁸, which confirms the sweep went well past 10⁵.)

**Elliptic curve L-functions.** ~130 curves via Frobenius trace encoding, ranks 0–3, CM
and non-CM. 95–100% RH-like.

**Riemann zeta zeros.** 95% detection accuracy on the first 20 zeros via spectral
resonances.

**Encoding artifact identified.** The L′/L encoding throws false positives at
approximately 74% of (q−1), with ratios in the 1.5–2.1× band. The specific case
q = 171,450,613 at index 127,348,284 (ratio 2.09×) was chased down by direct zero
computation and confirmed to be an encoding limitation, not a GRH violation. This was
honest work and remains the most useful negative detail from the classification phase.

**GPU implementation.** PyTorch with Intel Extension for PyTorch (IPEX) for an Intel Arc
A380, with CUDA and CPU backends. Batched eigensolves, projected 20–40 classifications/s
versus ~2.5 on CPU.

### Dec 2025 — Dimension reduction

- d = 3 is extremely stable (CoV < 3% across the (k,σ) grid) but retains only 64% of
  signal and drops contrast from 4× to 1.3×.
- d = 4–6 optimal — preserves full discrimination.
- Reformulated as a **two-channel classifier**: P(q), the polynomial projection strength,
  plus A(q), the integrated cosine defect.

Adversarial testing exposed boundary instability, disguise attacks (localized features
masked as oscillatory), and parameter sensitivity that flipped classifications.

### Dec 2025 — First deflation

AUC correlated 0.91–0.998 with P(q) alone. Grok's independent review confirmed: the
framework was functioning as a **polynomial rejection test**, not a spectral stiffness
metric. The project was killed here the first time.

### Jan 2026 — Failed rescue

Grok proposed a "positivity cascade" via the Weil–Guinand explicit formula. On inspection
it was **Li's criterion in disguise**, with sensitivity to near-critical zeros too small
to be practical. Dead again.

Around the same time, the one-pager was renamed from "Perry–Beurling Spectral Sieve" to
just "The Spectral Sieve" — the earlier name overclaimed for something still in
validation.

### Feb 2026 — Full autopsy

Fresh analysis found the December deflation was *partially premature*:

- The AUC–P(q) correlation was 0.91, not 0.998 — about 17% of variance unexplained.
- At matched P(q), AUC still varied 2×. The second channel was real.
- The residual norm ‖R_d‖ achieved **100%** classification accuracy against P(q)'s 92%.

But chasing that thread to its end closed the door properly: **R_d ≈ ‖q‖_μ**, because the
polynomial projection captures under 1% of total energy. The J-matrices, defect curves,
t-sweeps, and spectral norms were elaborate scaffolding around a μ-weighted L² norm of
the density error.

The deeper reason the whole enterprise was doomed: **the classification problem is
fundamentally easy.** Off-line zeros produce exponential growth; on-line zeros don't. Any
reasonable norm detects that. Choosing the optimal norm was the only non-trivial
contribution available.

### Jun 2026 — PSWF basis swap (final attempt)

Motivated by Connes–Moscovici's identification of prolate spheroidal wave functions as
the natural eigenbasis for zeta-zero spectral structure, the Legendre basis was replaced
with PSWFs. Rationale: the diagnosis said polynomial bases were the wrong choice.

Results, run in full:

- P–AUC correlation dropped 0.69 → 0.49, suggesting the PSWF basis carried
  non-polynomial information.
- **Zero detection collapsed from 95% to 5%.**

The second number is the only one that matters. The "new information" was noise, not
zero-location signal. The operator framework is trivial regardless of basis. Sieve stays
dead.

## 5. The one survivor

**The measure is not cosmetic.** The Gamma-shaped weight

$$w_\mu(s) = s^{\,k-1} e^{-\sigma s}$$

is a normalized bump peaking at s⋆ = (k−1)/σ. Choosing (k,σ) so s⋆ lands on the
intermediate-s band where the off-line-zero signature concentrates turns the quadratic
energy statistic

$$T_w = \int_0^1 w(s)\, q(s)^2\, ds$$

into an approximate matched filter for that excess.

**Mechanism.** The on-line and off-line populations differ *only* over an intermediate
band of s. A flat L² statistic integrates that discriminating band together with the
large flat regions on either side where the populations are identical — diluting signal
with pure variance. The s^(k−1) factor suppresses the s → 0 edge; e^(−σs) suppresses the
s → 1 tail. Matched filtering by another name.

**Measured gain.** Roughly **53%** on the original noisy Beurling ensemble. A clean
idealized reference run (seed 20260522, k = 4, σ = 6, so s⋆ = 0.5 sits at the band
center) gives flat d′ ≈ 2.37 against Gamma d′ ≈ 5.53 — more than double, because a
perfectly matched filter on a well-localized feature always does well. Both point the
same direction; the magnitude is ensemble-dependent and should be reported as such.

**Status.** A standalone technical note was written
(`weighted_measure_sensitivity_note.md`) with a self-contained reference implementation
in NumPy. Scope stated honestly in the note: a sensitivity result, not a proof tool. It
improves detector power; it says nothing about RH.

**Venue target.** 4–6 page computational note — *Integers*, a short *Journal of Number
Theory* communication, or *Experimental Mathematics*.

**Obvious next step if ever revisited.** Sweep (k, σ) and map realized gain against where
the weight peak s⋆ sits. That converts "approximately matched" into a quantitative tuning
curve — the figure that would anchor it as a standalone short paper.

---

# Part II — Related RH work

## 6. Wronskian positivity (Jan 2026)

**The criterion.**

$$(\Xi')^2 - \Xi\,\Xi'' > 0 \quad \text{for all } t$$

Rewritten, this says (Ξ′/Ξ)′ < 0 — the logarithmic derivative is strictly decreasing
wherever Ξ ≠ 0.

**Why it's equivalent to RH, not a route to it.** For an entire function of order 1 with
Hadamard factorization over its zeros, Ξ′/Ξ(t) = Σ_ρ 1/(t−ρ) + regular part. If all zeros
are real, each term is strictly decreasing away from its pole, and the sum is decreasing.
If some zero sits off the real line, the conjugate pair contributes local
non-monotonicity and the inequality fails somewhere.

So the criterion **is** RH in differential-inequality disguise. It characterizes the
Laguerre–Pólya class, and RH is precisely the statement that Ξ lies in that class.
Csordas–Norfolk–Varga (1986) had already established that the Turán inequalities and
generalizations are RH-equivalent. There is no known wedge between the two.

**What it later turned out to be.** The first Laguerre inequality — the d = 2 Jensen
condition — and the de Bruijn–Newman *lower*-bound mechanism. See §9.

## 7. Off-line mass and the optimal transport framing (Jan 2026)

**Off-line mass functional.**

$$M_\sigma = \sum_\rho \bigl(\mathrm{Re}(\rho) - \tfrac12\bigr)^2$$

Zero iff RH holds.

**Transport formulation.** Let μ be the empirical measure of zeros up to height T, and
R: s ↦ 1 − s̄ the reflection across the critical line. Then

$$T(\mu) = W_2^2(\mu,\, R_*\mu) = 4\sum_\rho \bigl(\mathrm{Re}(\rho) - \tfrac12\bigr)^2 = 4M_\sigma$$

giving the clean geometric statement:

> **RH ⇔ W₂(μ_zeros, R·μ_zeros) = 0** — the zeros are fixed points of reflection across
> Re(s) = 1/2.

Verified numerically: T = 4M_σ exactly. The small discrepancy in one test came from the
unconstrained OT solver matching zeros to non-partner reflections when close in
imaginary part; the functional equation constrains ρ to pair with 1 − ρ̄, and the exact
formula respects that.

Properties: deterministic (not statistical), exact (T = 4·n·σ² for n zeros at distance σ),
local (computable in windows), monotone (more contamination ⇒ strictly more cost).

**The theorem worth keeping.**

1. Montgomery's pair correlation statistic P(μ) depends only on |Im(ρᵢ) − Im(ρⱼ)|.
2. T(μ) = 4M_σ.
3. **P(μ) is blind to T(μ):** there exist measures with P ≈ GUE and T > 0. Explicit
   counterexamples constructed.

**Corollary.** No condition on P alone can imply T = 0. The minimal strengthening is
"P ≈ GUE **and** T = 0", but T = 0 is itself RH.

**What this kills.**

| Attack | Status |
|---|---|
| Montgomery ⟹ RH via Wronskian | Dead — explicit counterexamples |
| Any 1D spacing statistic ⟹ RH | Dead — 1D statistics are blind to 2D location |
| Statistical shortcuts to RH generally | Dead — any sufficient condition must encode Re = 1/2 directly |

Detector comparison established:

| Statistic | Detects off-line zeros? |
|---|---|
| 1D pair correlation (Montgomery) | No — completely blind |
| Off-line mass M_σ | Yes, exact |
| Transport cost T | Yes, T = 4M_σ |
| Wronskian W_min | Yes, goes negative |

**Assessment.** A negative result with teeth. It rules out a class of attacks and
identifies exactly why they fail. The transport formulation is novel framing but
equivalent content — it is M_σ = 0 restated, which is RH restated.

**Publishable as:** "Montgomery's Pair Correlation and the Wronskian Criterion: Why 1D
Statistics Cannot Reach RH" — 4–6 pages, *American Mathematical Monthly*, *Experimental
Mathematics*, or a *J. Number Theory* short communication.

**Code artifacts:** `wronskian_lab.py`, `detectors_2d.py`, `detectors_final.py`,
`transport_detector.py`.

## 8. Li's criterion (Jan 2026)

Reached twice — once through Grok's positivity cascade via Weil–Guinand, once directly.
Both times the conclusion was the same: computationally weaker than direct zero
verification, with sensitivity to near-critical zeros too small to be useful. Not a
route.

## 9. The de Bruijn–Newman constant Λ (Jun 2026)

The one target where the work interacted with genuinely open mathematics.

### Setup

$$H_t(z) = \int_0^\infty e^{tu^2}\,\Phi(u)\cos(zu)\,du, \qquad H_0(z) = \tfrac18\,\xi\!\left(\tfrac12 + \tfrac{iz}{2}\right)$$

H_t has all real zeros iff t ≥ Λ. **RH ⇔ Λ ≤ 0.**

Lower bound Λ ≥ 0 is Rodgers–Tao (*Forum of Mathematics Pi*, 2020). Upper bound engine is
Polymath15 (Tao et al., *Res. Math. Sci.* 2019): the Ki–Kim–Lee dynamical inequality

$$\Lambda \le t + \tfrac12\,\sigma_{\max}(t)^2$$

combined with a numerical "barrier" zero-free-region verification and rigorous
asymptotics.

### Thread 1 — H_t machinery and the Lehmer pair

Built H_t by direct quadrature, validated against the Polymath reference zero table to
10⁻⁷. Computed H₀ zeros at 28.27, 42.04, 50.02, 60.85, 65.87, 75.17 against the table's
28.25, 42.05, 50.05, 60.85, 65.85, 75.15, marching downward with increasing t exactly as
the table does.

Ran the Laguerre certificate L₁ = (H_t′)² − H_t H_t″ on the low window: **strictly
positive for every t from +0.4 down to −0.4, no dips.** That is the diagnosis, not a null
result — on the first eight zeros the all-real certificate survives deep into negative t.
The configuration does not start to fail at the bottom of the spectrum.

Then ran it on the actual **Lehmer pair at γ ≈ 7005.062866 and γ ≈ 7005.100565** (zero
#~6709). The pair leaves the real line at heat-time **t⋆ ≈ −7.13×10⁻⁴**. A generic
well-separated pair at the same height requires **>400× more backward flow** to go
complex.

Confirmed: the Wronskian condition is the first Laguerre inequality *and* the
de Bruijn–Newman lower-bound mechanism. Same object, arrived at from a different
direction eleven months earlier.

### Thread 2 — Exact heat smoothing

Built the exact representation of H_{−τ} as a convolution of H₀ with a Gaussian kernel,
proved via interchange of integration order, validated against direct quadrature to 10⁻⁷.
Confirmed exact at any height via Riemann–Siegel.

This mattered because direct quadrature dies at large height — the cos(zu) oscillation
and cancellation wreck float64, the same wall that forced Polymath onto an effective
A+B−C approximation. The convolution identity removes it without reconstructing
Polymath's effective sum from memory, where subtle constant errors hide.

### Thread 3 — Jensen hierarchy blindness (the best result of the arc)

Built the moment/Jensen polynomial hierarchy: moments

$$m_k(t) = \int u^{2k} e^{tu^2}\,\Phi(u)\,du$$

to order 34 at 40-digit precision. Validated m₀(0) = ξ(1/2)/8 to 10 digits. Confirmed the
Turán and higher-degree Jensen conditions pass at t = 0, consistent with
Csordas–Norfolk–Varga and Griffin–Ono–Rolen–Zagier (PNAS 2019).

Then demonstrated the sharp structural failure:

> The central certificate **falsely certifies hyperbolicity all the way down to
> t = −0.7**, with the minimum Turán ratio barely moving — 1.04343 → 1.04336.

**The index argument.** Accessible central Jensen shifts (n ≤ 30, moments to order ~34)
probe only the first ~30 zeros, heights ≲ 101. The Lehmer pair at zero #~6709 that
actually determines Λ would require moments to order **~13,400** to reach.

**Conclusion.** The upper bound on Λ is irreducibly a **local-at-height** phenomenon.
Bulk and central certificate reformulations cannot capture it — not as a matter of
computational budget, but structurally.

This is the cleanest publishable finding across the whole two-year arc: a precise,
rigorous, defensible statement about the *limits* of certificate methods for Λ. More
publishable as an honest note than anything pointed at the proof itself.

### The updated bound

Feeding Platt–Trudgian's rigorous verification of RH to height 3×10¹² (arXiv 2004.09765)
into the Polymath framework:

$$\boxed{\Lambda \le 0.20}$$

improving on the published Polymath value of 0.22. Bookkeeping rather than new
mathematics, but correct and worth stating explicitly.

### The ceiling, stated honestly

0.20 is gated entirely by rigorous verification height. Verify RH higher → the barrier
sits higher → σ_max(t) is squeezed smaller at the chosen t → tighter Λ. That is the
entire dependence.

Two hard facts bound the payoff:

- **Diminishing returns.** The jump to 3×10¹² — roughly 1.5–2 decades — bought 0.02.
  Order of magnitude: ~0.01 per decade of verified height. Λ ≤ 0.15 needs ~4 more
  decades; Λ ≤ 0.10 needs ~8 (height ~10²⁰). Infeasible fast.
- **It provably cannot reach Λ ≤ 0.** That requires verifying RH to infinite height —
  i.e. proving RH. The barrier route asymptotes strictly above zero.

Extending Platt–Trudgian past 3×10¹² is embarrassingly parallel and genuinely suited to
available hardware. It is a bounded compute grind that can chip but never close.

---

# Part III — Literature position

Scanned as of mid-2026, to make sure the work wasn't reinventing or ignoring the field.

**Guth–Maynard (May 2024, arXiv 2405.20552; accepted to *Annals*, 2025).** New bounds on
how often Dirichlet polynomials take large values, yielding the zero-density estimate
N(σ,T) ≤ T^{30(1−σ)/13+o(1)}. First improvement to Ingham's 1940 bound at σ = 3/4 in 84
years. Consequence: better control on zeros with real part 3/4, and correspondingly
shorter intervals for good prime-distribution estimates. **No known mechanism translates
an improved zero-density estimate into a better bound on σ_max(t) or on Λ.** That gap was
checked directly and is real.

**Connes, "Letter to Riemann" (Feb 2026).** 158-page survey with an original
contribution: extremizing a restriction of Weil's quadratic form using only primes below
13 yields approximations to the first 50 zeros with accuracies from 2.6×10⁻⁵⁵ to 10⁻³.
The prolate spheroidal operator plays a dual role — infrared, approximating the minimal
eigenvector of the Weil quadratic form; ultraviolet, modeling a self-adjoint operator
whose spectrum reflects the zeros. Proof strategy rests on convergence of zeros from
finite to infinite Euler products.

**Untested angle worth noting.** Nobody has published convergence curves for the finite
Euler product optimization. Connes showed P = 13; tracking eigenvalue convergence to
actual zeros as P grows (13 → 17 → 19 → … → 10⁶) is dense linear algebra on
GPU-friendly matrices, and the convergence *rate* is unknown. If exponential, that is
evidence for the proof strategy; if polynomial, it constrains which approaches can work.
This is the one clearly open computational question surfaced by the literature scan that
was never executed.

**Platt–Trudgian (arXiv 2004.09765).** Rigorous verification of RH to height 3×10¹² —
the current record, and the sole lever on the Λ upper bound.

---

# Part IV — Compute environment

| System | Configuration | Role |
|---|---|---|
| Intel Arc A380 | IPEX/PyTorch backend | Original sieve GPU implementation |
| NUKA | AMD RX 9070 XT | Effectively CPU-bound — ROCm/RDNA4 support unreliable |
| Luna | V100 16GB, 64GB DRAM, dual E5-2686 v4 | Intended host for heavy Λ / zero-verification work |

The sieve's GPU path targeted Intel Arc, CUDA, and CPU backends with batched eigensolves.
The Λ work is where serious compute would actually pay — zero verification past 3×10¹² is
embarrassingly parallel — but the payoff ceiling in §9 bounds how much it is worth.

---

# Part V — Publishable inventory

| Output | Length | Suggested venue | Status |
|---|---|---|---|
| Jensen/moment hierarchy blindness for Λ | Note | *Experimental Mathematics* / *J. Number Theory* | Result established, not written up |
| Montgomery blindness / 1D statistics cannot reach RH | 4–6 pp | *Amer. Math. Monthly*, *Exp. Math.*, *JNT* short comm. | Result established, outline drafted |
| Gamma-shaped measure sensitivity for Beurling detection | 4–6 pp | *Integers*, *JNT* short comm., *Exp. Math.* | **Note written**, with reference implementation |
| Λ ≤ 0.20 restatement | Short | Folds into the Jensen-blindness note | Stated, not written up |

The Jensen-blindness result is the strongest of the four and the only one that says
something new about the *limits* of an active research program rather than about an
in-house framework.

---

# Part VI — Post-mortem

**What went right.**

- Every kill was made on mathematical grounds, promptly, and stuck. The Dec 2025
  deflation, the Feb 2026 autopsy, and the Jun 2026 PSWF test all reached the same
  verdict independently.
- The encoding-artifact investigation at ~74% of (q−1) was chased to direct zero
  computation rather than being written off — that is the right instinct and it produced
  a correct answer.
- The Feb 2026 autopsy found a real second channel that had been killed prematurely, then
  followed it honestly to the conclusion that it didn't matter. Both halves of that are
  correct behavior.

**What went wrong, structurally.**

- **The classification problem was never hard.** Off-line zeros produce exponential
  growth; on-line zeros don't. Any reasonable norm detects that. Seven months of operator
  theory was built on top of a problem whose difficulty was already exhausted by the
  choice of norm. That should have been checked at the start with a single baseline: how
  well does a plain unweighted L² norm do?
- **Equivalence traps repeated four times.** The Wronskian, the transport framing, Li's
  criterion, and the positivity cascade are all RH restated. Each took real work to
  identify. The general test — *does this statement admit partial, quantitative,
  unconditional progress, or is it a biconditional?* — would have caught all four
  immediately.
- **Naming ran ahead of validation.** "Perry–Beurling Spectral Sieve" was attached to
  something still in validation and had to be walked back. Names last longer than
  results.

**The rule extracted.** Before building a framework: construct the simplest possible
baseline that solves the same problem, and measure against it. If the elaborate version
doesn't beat the baseline by a wide margin, the elaboration is scaffolding. That single
check would have collapsed the Spectral Sieve in a week instead of seven months — and it
is exactly the check that eventually killed it.

---

*Record compiled August 2026 from the full project history, November 2025 – June 2026.*
