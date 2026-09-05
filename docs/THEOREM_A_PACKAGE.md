# Full Theorem A — proved for the exact continuous residual under RH

**Date:** 2026-09-05 (corrected arithmetic proof); package history 2026-08-11\
**Status:** **Proved conditional continuous theorem**, using the exact explicit formula and a cited positive mean-square limit.\
**Unconditional RH:** **open** — this package does **not** prove RH.\
**Arithmetic numerics:** soft plateau \(R_d\sim0.15\)–\(0.19\) is **not** a counterexample to
conditional A (finite-\(T\) probe vs \(T\to\infty\) under RH + ANT).

Pointers: [`THEOREMS_AB.md`](THEOREMS_AB.md) · [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md) ·\
[`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md) · [`INFINITE_TAIL_REMAINDER.md`](INFINITE_TAIL_REMAINDER.md) ·\
[`RH_CLOSEOUT_ROADMAP.md`](RH_CLOSEOUT_ROADMAP.md) · Code: `pbss.ab_closure`, `pbss.theorem_a_chain`

---

## Explicit non-claim

**This document does not contain an unconditional proof of the Riemann Hypothesis.**

The proof now concerns the **exact continuous**, unweighted θ residual with
continuous orthogonal affine detrending and fixed degree \(d\ge2\). Under RH,
its ratio satisfies \(R_d=O_d(T^{-2})\).

**Proof:** [`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md), with exact primary
source statements and the complete adaptation. The sampled model API does not
verify a continuous limit theorem. Optional weighted and sampled transfers remain
separate obligations.

**Correction of the August package.** Generic references alone did not establish
its stated sign/phase identification, growing-mode uniformity, prime-square
remainder, or positive relative norm floor. The new proof supplies these steps
explicitly for the continuous unweighted object. The former broad claim that all
of the shipped residual/weight interfaces were already matched is superseded.
Scaffold-only tail majorants remain diagnostic and are not proof inputs.

---

## 1. Precise statement

### Definitions (as in THEOREMS_AB)

- Log window: \(x=e^{uT}\), \(u\in[0,1]\).
- Orthonormal shifted Legendre \(\{\varphi_k\}\); \(P_d\) onto \(V_d\); \(R_d=\|P_dq\|_2^2/\|q\|_2^2\).
- **Exact continuous arithmetic residual:**
  \(Q(y)=e^{-y/2}(\theta(e^y)-e^y)\),
  \(q_T=(I-P_1)Q(T\cdot)\), where \(P_1\) is orthogonal in continuous \(L^2(du)\).
- No additional weight or smoothing operator is included in this theorem.

### Theorem A (conditional continuous form)

**Assume RH.** For each fixed integer \(d\ge2\),
\[
R_d(q_T)=O_d(T^{-2}),\qquad
\|q_T\|_2^2\longrightarrow
2\sum_{\gamma>0}\frac{m_\gamma^2}{1/4+\gamma^2}>0.
\]
The external inputs are the classical explicit formula and the mean-square theorem
identified in §3. Simplicity and linear independence of ordinates are not assumed.
For \(d=0,1\), the numerator vanishes by detrending and supplies no converse information.

`probes.arithmetic_residual` computes samples and a sample least-squares detrend.
That implementation is a numerical diagnostic; transferring the theorem to a
refining grid requires quantitative quadrature and detrend error bounds.

---

## 2. Gap table (every step labeled)

| Step | Disposition | Support |
|------|-------------|---------|
| M1–M4 diagnostic lemmas | **Proved** | `PROOFS_LEMMAS.md`; `tests/test_lemmas.py` |
| M5 finite CL / truncated EF | **Proved** | `PROOFS_LEMMAS.md`; `lemmas.bound_R_d_finite_mode_sum` |
| M6 weighted model decay | **Proved** | `PROOFS_LEMMAS.md`; weighted majorants |
| M7 \(R_d\) perturbation majorant | **Proved** | `PROOFS_LEMMAS.md` (M7); `ab_closure.energy_ratio_perturbation_bound` |
| ANT-3 exact ψ formula | **Cited precisely** | Akbary–Ng–Shahabi (1.5) — §3 |
| ANT-1 uniform polynomial moments | **Proved adaptation** | Uniform coefficient series and cutoff \(G=e^T\) |
| ANT-2 prime-square correction | **Proved adaptation** | Explicit constant \(-1\) and \(L^1\cap L^2\) secondary error |
| ANT-MS positive denominator | **Cited + adapted** | Akbary–Ng–Shahabi Cor. 1.15(i) |
| ANT-4 optional weight transfer | **Not established here** | Separate extension; model M6 is insufficient alone |
| Sampled grid/detrend transfer | **Not established here** | Separate numerical error estimates required |
| Scaffold zero-tail model | **Diagnostic only** | Not a required Full-A step |
| Continuous Full A under RH | **Proved conditional** | `CONTINUOUS_THEOREM_A.md` |
| Unconditional RH | **Open** | The theorem assumes RH |

Machine-readable: `pbss.ab_closure.full_a_gap_table()`, `ant_citations()`.

---

## 3. Exact ANT inputs and continuous adaptations

The primary source is Akbary, Ng, and Shahabi, *Limiting distributions of the
classical error terms of prime number theory*, Q. J. Math. **65** (2014), 743–780.
Numbers below refer to [arXiv:1306.1657v1](https://arxiv.org/abs/1306.1657v1), also
available as the [author PDF](https://www.cs.uleth.ca/~akbary/Akbary-Ng-Shahabi.pdf).

### ANT-3 — Exact normalized ψ formula

Equation (1.5) gives, under RH,
\[
F_\psi(y)=\operatorname{Re}\sum_{0<\gamma\le G}
\frac{-2m_\gamma e^{i\gamma y}}{1/2+i\gamma}
+O\!\left(\frac{e^{y/2}(y+\log G)^2}{G}+ye^{-y/2}\right).
\]
The sign and phase are part of the identity. A generic positive-amplitude model
sum cannot be substituted with a vanishing relative error merely by matching
its ordinates.

### ANT-1 — Uniform polynomial moments

Choose \(G=e^T\). Integration by parts bounds each fixed polynomial moment of
a retained mode by \(C_g/(\gamma T)\). Since
\(\sum m_\gamma/(\gamma|1/2+i\gamma|)<\infty\), the summed bound is uniform in
\(G\); the integrated explicit-formula error is also \(O_g(T^{-1})\).
Thus \(\int_0^1F_\psi(Tu)g(u)du=O_g(T^{-1})\).
This replaces the unsupported passage from fixed-block M5 to an arbitrary
schedule \(N(T)\), and the former insufficient illustrative cutoff scales.

### ANT-2 — Explicit prime-square correction

The prime-power identity and RH estimates give
\[
Q(y)=F_\psi(y)-1+h(y),\qquad h\in L^1\cap L^2.
\]
The continuous affine detrend removes the constant \(-1\). The remaining term
has \(O(T^{-1})\) polynomial moments and \(O(T^{-1})\) squared norm.
The elementary bound \(\psi-\theta=O(\sqrt{x}\log x)\) alone would not prove this.

### ANT-MS — Positive mean-square limit

Corollary 1.15(i), specialized to \(L=\zeta\), gives
\[
T^{-1}\int_0^T F_\psi(y)^2dy\longrightarrow
V=2\sum_{\gamma>0}\frac{m_\gamma^2}{1/4+\gamma^2}\in(0,\infty).
\]
The explicit correction and detrend leave the same norm limit. This establishes
the denominator needed by the arithmetic ratio, with no guessed model norm floor.

### ANT-4 — Optional weight transfer

The continuous arithmetic theorem above is unweighted. M6 covers fixed regular
weights on model sums, but a bounded endpoint taper alone does not prove that
weighted arithmetic norms are comparable to unweighted norms. This extension is
not included in the proved statement.

## 4. Conditional deduction

Under RH, ANT-3 and the uniform calculation in ANT-1 give \(O_k(T^{-1})\)
coefficients. ANT-2 removes the prime-square constant explicitly. Therefore
\(\|P_dq_T\|_2^2=O_d(T^{-2})\). ANT-MS gives \(\|q_T\|_2^2\to V>0\), proving
\(R_d(q_T)=O_d(T^{-2})\).

Every estimate and its scope is detailed in
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md). The generic M7 perturbation
lemma remains valid, but finite sampled M7 checks are not the proof of these
arithmetic hypotheses.

---

## 5. Relation to numerics (do not misread)

| Campaign | Finding | Completes unconditional A? |
|----------|---------|----------------------------|
| Grand / extend-\(x\) plateau | Arith \(R_4\sim0.15\)–\(0.19\) | **No** — finite \(T\); conditional A is \(T\to\infty\) under RH+ANT |
| Open-plateau peel | Model zeros ≠ full arith identity | Consistent with needing ANT-2/3, not a refutation |
| Model CL/EF M5–M6 tests | Agreement with model majorants | **No** — numerical model checks |

---

## 6. Package status summary

| Item | Status |
|------|--------|
| Continuous Full A under RH | **Proved**, exact unweighted statement (§1) |
| Model lemmas M1–M6 + M7 | **Proved** |
| ANT-3 / ANT-MS | **Precisely cited**; ANT-1 / ANT-2 **adaptations proved** (§3) |
| Optional weighted and sampled transfers | **Separate, not established here** |
| Unconditional RH | **Open** |
| Machine status | `full_arithmetic_A = "proved_conditional_continuous"` |

**Continuous arithmetic Theorem A is proved under RH. RH remains open.**
