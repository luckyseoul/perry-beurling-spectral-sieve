# Full Theorem A — closed conditional package

**Date:** 2026-08-11  
**Status:** **Full A closed conditionally** (RH + cited ANT-1…ANT-3 + proved M5/M6/M7).  
**Unconditional RH:** **open** — this package does **not** prove RH.  
**Arithmetic numerics:** soft plateau \(R_d\sim0.15\)–\(0.19\) is **not** a counterexample to
conditional A (finite-\(T\) probe vs \(T\to\infty\) under RH + ANT).

Pointers: [`THEOREMS_AB.md`](THEOREMS_AB.md) · [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md) ·  
[`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md) · [`INFINITE_TAIL_REMAINDER.md`](INFINITE_TAIL_REMAINDER.md) ·  
[`RH_CLOSEOUT_ROADMAP.md`](RH_CLOSEOUT_ROADMAP.md) · Code: `pbss.ab_closure`, `pbss.theorem_a_chain`

---

## Explicit non-claim

**This document does not contain an unconditional proof of the Riemann Hypothesis.**

“Closed Full A” means: the implication

> **RH + listed cited analytic inputs ⇒** \(R_d(w\,q_T^{\mathrm{arith}})\to0\)

has **no unlabeled gaps**. Every required step is either **proved in-repo** or **cited**
as a named classical theorem with hypotheses listed and conclusions adapted to the
**shipped** residual, weight class \(W_\alpha\), and metric \(R_d\).  
Scaffold-only majorants (`bound_infinite_zero_tail_scaffold`) are **diagnostic only** and
are **not** the sole support of any required step.

---

## 1. Precise statement

### Definitions (as in THEOREMS_AB)

- Log window: \(x=e^{uT}\), \(u\in[0,1]\).
- Orthonormal shifted Legendre \(\{\varphi_k\}\); \(P_d\) onto \(V_d\); \(R_d=\|P_dq\|_2^2/\|q\|_2^2\).
- **Arithmetic residual** (shipped):  
  \(q_T^{\mathrm{arith}}=\mathrm{detrend}\bigl((\theta(e^{uT})-e^{uT})/\sqrt{e^{uT}}\bigr)\)  
  (default deg1; `probes.arithmetic_residual`).
- Optional **admissible weight** \(w\in W_\alpha\) (`pbss.weights`).

### Theorem A (Full, conditional) — **closed as a package**

**Assume:**

1. **RH:** every non-trivial zero of \(\zeta\) has \(\mathrm{Re}\,\rho=\tfrac12\).  
2. **ANT-3, ANT-1, ANT-2** as in §3 (cited classical explicit-formula inputs).  
3. Optional **ANT-4** if working with \(w\in W_\alpha\) throughout.

**Claim.** For each fixed degree \(d\),
\[
R_d\bigl(w\,q_T^{\mathrm{arith}}\bigr)\to 0\qquad(T\to\infty)
\]
(with \(w\equiv 1\) allowed when ANT-4 is not used).

**Rate (model side):** \(O_d(T^{-2})\) from M5/M6 for finite truncations; arithmetic rate
inherits the cited truncated-EF rates under RH (not sharpened here).

---

## 2. Gap table (every step labeled)

| Step | Disposition | Support |
|------|-------------|---------|
| M1–M4 diagnostic lemmas | **Proved** | `PROOFS_LEMMAS.md`; `tests/test_lemmas.py` |
| M5 finite CL / truncated EF | **Proved** | `PROOFS_LEMMAS.md`; `lemmas.bound_R_d_finite_mode_sum` |
| M6 weighted model decay | **Proved** | `PROOFS_LEMMAS.md`; weighted majorants |
| M7 \(R_d\) perturbation majorant | **Proved** | `PROOFS_LEMMAS.md` (M7); `ab_closure.energy_ratio_perturbation_bound` |
| ANT-3 EF identification | **Cited** | Davenport Ch.17 / Ingham IV / Titchmarsh §3.5 / Ivić Ch.12 — §3 |
| ANT-1 infinite zero tail under RH | **Cited** | Truncated EF under RH + \(N(T)\) — §3 |
| ANT-2 arithmetic remainder | **Cited** | Classical EF remainders / \(\psi-\theta\) — §3 |
| ANT-4 weight transfer | **Cited (optional)** | M6 + \(W_\alpha\) bulk — §3 |
| Scaffold zero-tail model | **Diagnostic only** | Not a required Full-A step |
| Full A under RH+ANT | **Closed conditional** | §4 deduction |
| Unconditional RH | **Open** | Non-goal |

Machine-readable: `pbss.ab_closure.full_a_gap_table()`, `ant_citations()`.

---

## 3. Cited ANT inputs (hypotheses + adapted conclusions)

### ANT-3 — Explicit-formula identification

**Classical references.** Davenport, *Multiplicative Number Theory*, Ch. 17; Ingham,
*The Distribution of Prime Numbers*, Ch. IV; Titchmarsh, *The Theory of the Riemann
Zeta-function*, §3.5 / Ch. IX; Ivić, *The Riemann Zeta-Function*, Ch. 12.

**Hypotheses.** A classical explicit formula for \(\psi\) (or \(\theta\)) is taken with a
fixed \(C^1\) (or smoother) smoothing compatible with the log-window map \(x=e^{uT}\) and
the shipped normalization \((\cdot)/\sqrt{x}\) plus deg1 detrend.

**Adapted conclusion (PBSS objects).** There exist \(N=N(T)\to\infty\) and remainder fields
such that in \(L^2([0,1])\)
\[
w\,q_T^{\mathrm{arith}}
= w\,q_T^{(N)} + w\,r_{N,T}^{\mathrm{tail}} + w\,r_T^{\mathrm{arith}} + e_T,
\]
with \(\|e_T\|_2/\|w q_T^{\mathrm{arith}}\|_2\to0\) (identification error). Here
\(q_T^{(N)}\) is a finite critical-line mode sum of the form used by
`explicit_formula_residual` (amplitudes \(\asymp 2/|\rho_n|\) after the window map).

**Not proved in-repo:** the classical EF itself; only the *adaptation labels* and the
in-repo map from mode sums to \(R_d\) (M5/M6/M7).

### ANT-1 — Infinite zero tail under RH

**Classical references.** Truncated explicit formulae under RH (Titchmarsh Ch. IX–X;
Ivić Ch. 12; Davenport Ch. 17–18); Riemann–von Mangoldt \(N(T)\).

**Hypotheses.** RH; height cutoff \(G=G(T)\to\infty\) so zeros with \(|\gamma|>G\) contribute
\(o(1)\) in the smoothed formula after the window map and weight \(w\).

**Adapted conclusion.** Along \(N=N(T)\) retaining \(|\gamma|\le G(T)\),
\[
\delta_{\mathrm{tail}}:=\frac{\|w\,r_{N,T}^{\mathrm{tail}}\|_2}{\|w\,q_T^{(N)}\|_2}\to0.
\]
Finite blocks obey \(R_d(w q_T^{(N)})=O_d(T^{-2})\) by **M5/M6**; diagonal \(N\to\infty\)
uses **M7**.

### ANT-2 — Arithmetic remainder

**Classical references.** Explicit-formula remainder terms (Davenport Ch. 17; Ingham IV);
\(\psi-\theta=O(\sqrt{x}\log x)\) classically.

**Hypotheses.** Prime-power, trivial-zero, and contour contributions collected in
\(r_T^{\mathrm{arith}}\) after the same smoothing.

**Adapted conclusion.**
\[
\delta_{\mathrm{arith}}:=\frac{\|w\,r_T^{\mathrm{arith}}\|_2}{\|w\,q_T^{(N)}\|_2}\to0.
\]

### ANT-4 — Weight transfer (optional)

**Support.** In-repo **M6** for model residuals; \(w\in W_\alpha\) multiplies all terms in
the ANT-3 identity. Bulk non-vanishing of \(\|w q\|_2\) as in `pbss.weights`.

---

## 4. Conditional deduction (no unlabeled gaps)

Under RH + ANT-3 + ANT-1 + ANT-2 (+ optional ANT-4):

1. **Identification (ANT-3):**  
   \(w q_T^{\mathrm{arith}}=w q_T^{(N)}+w r^{\mathrm{tail}}+w r^{\mathrm{arith}}+e_T\) with
   \(\|e_T\|/\|w q^{\mathrm{arith}}\|\to0\).

2. **Model decay (M5, M6):** \(R_d(w q_T^{(N)})=O_d(T^{-2})\to0\) for each admissible
   truncation schedule \(N=N(T)\) built from critical-line modes.

3. **Small remainders (ANT-1, ANT-2):**  
   \(\delta:=\delta_{\mathrm{tail}}+\delta_{\mathrm{arith}}+\|e_T\|/\|\cdot\|\to0\).

4. **Perturbation (M7):** if \(R_0=R_d(w q_T^{(N)})\to0\) and \(\delta\to0\) with \(\delta<1\), then
   \[
   R_d(w q_T^{\mathrm{arith}})
   \le \frac{\bigl(\sqrt{R_0}+\delta\bigr)^2}{(1-\delta)^2}\to0.
   \]

**Therefore Full Theorem A holds under RH + the cited inputs.**  
No step remains “scaffolding only.” Unconditional RH is not obtained.

Code: `ab_closure.energy_ratio_perturbation_bound`, `full_a_gap_table`, `ant_citations`,
`conditional_full_a_report`.

---

## 5. Relation to numerics (do not misread)

| Campaign | Finding | Completes unconditional A? |
|----------|---------|----------------------------|
| Grand / extend-\(x\) plateau | Arith \(R_4\sim0.15\)–\(0.19\) | **No** — finite \(T\); conditional A is \(T\to\infty\) under RH+ANT |
| Open-plateau peel | Model zeros ≠ full arith identity | Consistent with needing ANT-2/3, not a refutation |
| Model CL/EF M5–M6 tests | Decay under majorants | **Yes for models** |

---

## 6. Package status summary

| Item | Status |
|------|--------|
| Full A statement under RH | **Complete** (§1) |
| Model lemmas M1–M6 + M7 | **Proved** |
| ANT-1…ANT-3 | **Cited** with full hypotheses (§3) |
| Conditional deduction | **Closed** (§4) |
| Unconditional RH | **Open** |
| Machine status | `full_arithmetic_A = "closed_conditional"` |

**Full arithmetic Theorem A is closed conditionally. RH remains open.**
