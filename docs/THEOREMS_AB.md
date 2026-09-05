# Theorems A & B — precise statements (PBSS)

**Status date:** 2026-09-05 (continuous arithmetic proof and converse audit)\
**Not an unconditional proof of the Riemann Hypothesis.**

This document locks definitions so numerics, code, and math refer to the same objects.
Proved results in this repo are **about the diagnostic** (lemmas M1–M7, including
finite-mode A₀ and the \(R_d\) perturbation majorant M7). Full Theorem A is
**proved under RH for the exact continuous unweighted residual with orthogonal
affine detrending**. The absolute-moment converse is an RH equivalence; the
normalized fast-decay converse remains open. Unconditionally, the exact continuous
detrended ratio has lower limit zero. Unconditional RH remains **open**.

---

## 0. Definitions

### Log window

Fix \(T>0\) (logarithmic window length). Set
\[
x = e^{uT},\qquad u\in[0,1].
\]
Functions of \(x\) pull back to functions of \(u\) on the unit interval with
Lebesgue measure \(du\).

### Orthonormal basis

Let \(L_k\) be the classical Legendre polynomials on \([-1,1]\). Define
\[
\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1),\qquad u\in[0,1].
\]
Then \(\{\varphi_k\}_{k\ge0}\) is orthonormal in \(L^2([0,1],du)\):
\[
\langle f,g\rangle=\int_0^1 f(u)g(u)\,du.
\]
Write \(V_d=\mathrm{span}\{\varphi_0,\ldots,\varphi_d\}\) and \(P_d\) the
orthogonal projection onto \(V_d\).

### Density residual \(q_T\)

In the full arithmetic setting, \(q_T\) is a normalized residual built from a
Chebyshev / prime-counting discrepancy on \([1,e^T]\) (see `probes.probe_prime_residual`).
**Model residuals** used for proved lemmas:

| Name | Formula on \([0,1]\) | Interpretation |
|------|----------------------|----------------|
| Critical-line mode | \(q_T^{\mathrm{cl}}(u)=\sin(tTu)\) | zero at \(\tfrac12+it\) |
| Off-critical mode | \(q_T^{\sigma,t}(u)=e^{T(\sigma-1/2)u}\sin(tTu)\) | zero at \(\sigma+it\) |
| Orthogonal defect | \(q=\sqrt{1-\varepsilon^2}\,f+\varepsilon\varphi_j\) | fixed low-degree mass |

### Energy ratio and scaled strength

\[
R_d(q)=\frac{\|P_d q\|_{L^2}^2}{\|q\|_{L^2}^2}
=\frac{\sum_{k=0}^d|\langle q,\varphi_k\rangle|^2}{\|q\|^2}\in[0,1],
\]
\[
S_d(q;T)=T^{2(d+1)}R_d(q).
\]
Code: `pbss.projection.energy_ratio`, `scaled_projection_strength`.
Working name \(P(q):=S_d(q;T)\) is a scaling convention, not the lost legacy “3.92” normalization.

---

## 1. Proved in this repository (lemmas about the diagnostic)

Proofs: continuous \(L^2\) arguments in `src/pbss/lemmas.py` (module docstring)
and `docs/PROOFS_LEMMAS.md`. Discrete verification: `tests/test_lemmas.py`.

### Lemma M1 (pure mode energy)

If \(q=\varphi_m\), then
\[
R_d(q)=\begin{cases}1&m\le d,\\0&m>d.\end{cases}
\]

### Lemma M2 (orthogonal defect formula)

Let \(j\le d\), \(\varepsilon\in[0,1]\), \(f\perp V_d\), \(\|f\|_2=1\), and
\[
q=\sqrt{1-\varepsilon^2}\,f+\varepsilon\,\varphi_j.
\]
Then
\[
\boxed{R_d(q)=\varepsilon^2}.
\]

### Lemma M3 (critical-line pure mode decay)

For \(\omega>0\) and \(q_\omega(u)=\sin(\omega u)\),
\[
R_d(q_\omega)=O_d(\omega^{-2})\qquad(\omega\to\infty).
\]
Hence for \(q_T^{\mathrm{cl}}(u)=\sin(tTu)\) with fixed \(t>0\),
\[
\boxed{R_d(q_T^{\mathrm{cl}})=O(T^{-2})\qquad(T\to\infty)}.
\]

**Remark.** The archive heuristic \(R_d=O(T^{-2(d+1)})\) is **stronger** and is
**not** proved here. M3 is the rigorous model-mode rate we use.

### Lemma M4 (persistent defect blocks vanishing)

Under M2 with fixed \(\varepsilon>0\), \(R_d(q)=\varepsilon^2\not\to0\) no matter how
oscillatory the orthogonal part \(f\) becomes. Therefore, if a family satisfies
\(R_d(q_T)\to0\), its low-degree mass must tend to zero.

---

## 2. Theorem A (conditional on RH) — precise form

### Model theorem A₀ (**proved**)

For the critical-line pure mode \(q_T^{\mathrm{cl}}(u)=\sin(tTu)\),
\[
R_d\bigl(q_T^{\mathrm{cl}}\bigr)\to0\quad(T\to\infty)
\]
at rate \(O(T^{-2})\) (Lemma M3).

### Finite-mode A₀ (**proved**, Lemma M5)

For a **finite** superposition of critical-line modes
\[
q_T^{(N)}(u)=\sum_{n=1}^{N}a_n\sin(t_n T u+\phi_n)
\]
(\(N<\infty\), distinct \(t_n>0\), \(a\not\equiv0\)),
\[
R_d\bigl(q_T^{(N)}\bigr)=O_d(T^{-2})\qquad(T\to\infty)
\]
at the **same order** as pure-mode M3. Proof: `docs/PROOFS_LEMMAS.md` (M5);
code: `lemmas.bound_R_d_finite_mode_sum`, `probes.finite_cl_superposition`,
`probes.explicit_formula_residual`.

This is the bridge from pure A₀ to truncated explicit-formula residuals. Infinite
zero sums and the arithmetic residual are handled by the explicit uniform moment
bound and positive mean-square theorem in
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md), not by M5 alone.

### Scaffolding toward full A (2026-07-26)

Shipped (not full A):

1. **Weight class \(W_\alpha\)** — `pbss.weights`: Tukey/Hanning windows, endpoint
   contribution \(E_{\mathrm{end}}\), weighted \(R_d\). Documents the taper effect seen
   in open-plateau as a controlled object.
2. **Lemma M6 (proved model):** admissible weights preserve \(O(T^{-2})\) decay for pure
   CL and finite EF residuals — `bound_R_d_weighted_sine_order`,
   `docs/PROOFS_LEMMAS.md`.
3. **Truncation remainder path** — `pbss.remainder`: peel via \(q-\alpha q_T^{(N)}\),
   M5 majorant on finite stripped blocks, scaffolding infinite-tail majorant
   (`bound_infinite_zero_tail_scaffold`). Note: [`INFINITE_TAIL_REMAINDER.md`](INFINITE_TAIL_REMAINDER.md).
4. **Arithmetic weight multi-\(T\)** — `experiments/run_arithmetic_weights.py` →
   `results/arithmetic_weights/`.

Writeup: [`THEOREM_A_SCAFFOLD.md`](THEOREM_A_SCAFFOLD.md). Campaign:
`results/theorem_a_scaffold/`. Historical scaffolding path only — **Full A is closed
conditionally for the exact continuous unweighted residual** in
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md) (the scaffold is not a proof input).

### Full Theorem A (**proved under RH for the exact continuous residual**)

**Hypothesis (RH).** All non-trivial zeros of \(\zeta\) satisfy \(\mathrm{Re}\,\rho=\tfrac12\).

**Theorem under RH.** Let \(Q(y)=e^{-y/2}(\theta(e^y)-e^y)\) and
\(q_T=(I-P_1)Q(T\cdot)\), with continuous orthogonal affine detrend. For fixed
\(d\ge2\), with no added weight or smoothing,
\[
R_d(q_T)=O_d(T^{-2})\qquad(T\to\infty).
\]

**Package:** [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) · code `pbss.ab_closure`.

**Status:** **Proved conditional continuous theorem.** The proof explicitly handles
the formula's sign and phase, growing-zero-block uniformity, the prime-square
constant, and the positive denominator. See
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md) and its precise primary sources.
The earlier generic ANT adaptations are superseded. Optional weighted and sampled
transfers are separate obligations. Scaffold tail majorants remain diagnostic only.\
**Unconditional RH remains open.** Roadmap:
[`RH_CLOSEOUT_ROADMAP.md`](RH_CLOSEOUT_ROADMAP.md) (**not** a proof of RH).

### Arithmetic multi-\(T\) evidence (measured, not a proof)

Shipped builder: `pbss.probes.arithmetic_residual` — \(q_T=( \theta(x)-x)/\sqrt{x}\) on
\(x=e^{uT}\). Large campaign: `experiments/run_overnight_campaign.py`
→ `results/overnight_campaign.json` (also earlier mid-scale
`results/arithmetic_multi_T.json`).

**Grand campaign** (`experiments/run_grand_campaign.py`,
`results/grand_campaign/`): \(x_{\max}=10^{10}\), 455M primes checkpointed,
**20 000 MC defect trials per \(T\)** (14 windows → 280 000 trials), multi-\((d,\mathrm{detrend},\mathrm{smooth})\),
controls CL / off-critical / defect.

| \(T\) | \(x_{\max}\) | \(R_4\) arith (deg1) | \(R_4\) CL | \(R_4\) defect | MC mean \(R_4\) |
|------:|-------------:|---------------------:|-----------:|---------------:|----------------:|
| 10 | 2e4 | 0.144 | ~0.006 | 0.250 | 0.793 |
| 16 | 9e6 | 0.193 | ~0.002 | 0.250 | 0.793 |
| 23 | **1e10** | **0.155** | ~0.0005 | 0.250 | 0.795 |

**Reading:** arithmetic soft-plateau ~0.15–0.19 through \(10^{10}\); MC defects stay high (~0.79) and flat in \(T\). **No** full Theorem A / RH claim.

**Does not prove or disprove RH or full Theorem A.**

---

## 3. Theorem B (obstruction) — precise form

### Model theorem B₀ (**proved as obstruction for the diagnostic**)

If a residual family admits a uniform lower bound
\[
R_d(q_T)\ge\varepsilon_0^2>0
\]
for all large \(T\), then it **cannot** be of the form “high-frequency only”
in the sense of M2 with \(\varepsilon\to0\). Equivalently (M4): vanishing of
\(R_d\) is **necessary** for the absence of a persistent low-degree component.

### Full Theorem B (**normalized converse open**)

**Claim.** If \(R_d(q_T^{\mathrm{arith}})\) decays sufficiently rapidly as
\(T\to\infty\), then \(\zeta\) has no zero with \(\mathrm{Re}\,\rho\neq\tfrac12\).

**Package:** [`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md) · `pbss.ab_closure`.

**Status:** The normalized fast-decay converse remains **open**. The phrase
“sufficiently rapidly” requires a fixed rate and normalization; mere vanishing
and \(O(T^{-2})\) are different assertions. A pure off-critical real mode has
phase-dependent troughs and provides no uniform positive floor.

Two distinct arithmetic theorems are now available:

- **B-ABS:** one fixed absolute Legendre moment has subexponential growth if and
  only if RH holds; the equivalent growth bound remains unproved unconditionally.
  See [`ABSOLUTE_MOMENT_CONVERSE.md`](ABSOLUTE_MOMENT_CONVERSE.md).
- **Normalized recurrence:** \(\liminf R_d(q_T)=0\) unconditionally for the exact
  continuous affine-detrended θ residual. This does not prove a full limit or rate.
  The positive-floor B-RES formulation is consequently RH-equivalent and open.
  See [`NORMALIZED_RECURRENCE_THEOREM.md`](NORMALIZED_RECURRENCE_THEOREM.md).

Neither finite model checks nor an RH equivalence establishes its unproved side.
Unconditional RH remains **open**.

---

### Explicit-formula residual + peel scan (measured)

Builder: `pbss.probes.explicit_formula_residual` — truncated sum of first \(N\)
critical-line modes with amplitudes \(2/|\rho_n|\) on the log-window (offline
ordinates in `pbss.zeros`). Peel helper: `peel_residual`.

Campaign: `experiments/run_explicit_formula_peel.py` →
`results/explicit_formula_peel/` (JSON, TXT, plot). Multi-\(T\) × multi-\(N\):
include \(R_d(q_T^{(N)})\) decays with \(T\) for each fixed \(N\) (finite-mode A₀);
peel column records \(R_d\) after stripping the first \(N\) modes from a fixed
\(N_{\mathrm{full}}\) sum.

**Not a proof of RH or full Theorem A.**

---

## 4. What is open

1. **Unconditional RH**, including the equivalent absolute-moment growth bound.\
2. The **normalized fast-decay converse** and the RH-equivalent positive-floor B-RES formulation.\
3. A full normalized limit or rate **without RH**; the unconditional theorem proves only a lower limit of zero.\
4. Legacy normalization with \(P\approx3.92\).\
5. Quantitative transfer from the proved continuous theorem to sampled or optionally weighted residuals.\
6. Finite-\(T\) arithmetic plateau explanation (secondary terms) — consistent with
   conditional A, without constituting an asymptotic proof.

### Marathon campaigns (measured)

See `docs/STATUS.md` marathon section: arithmetic zero-peel, extend-\(x\) to
\(5\times10^{10}\) (parallel sieve), Beurling battery, MC ≥50k/T. **Not RH proofs.**


---

## 5. Numerical support (multi-\(T\))

See `results/multi_T_scan.json` .

Expected qualitative picture (must match the run):

- Critical-line \(R_d(T)\) **falls** roughly like \(T^{-2}\).\
- Persistent defect \(R_d\equiv\varepsilon^2\) **flat**.\
- Clear gap for large \(T\): defect stays high, RH-like mode drops.

This supports A₀/B₀. It does **not** prove RH.
