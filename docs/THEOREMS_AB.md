# Theorems A & B — precise statements (PBSS)

**Status date:** 2026-07-26  
**Not an unconditional proof of the Riemann Hypothesis.**

This document locks definitions so numerics, code, and math refer to the same objects.
Proved results in this repo are **about the diagnostic**. Claims that need RH or
that convert low \(R_d\) into zero-free regions are labeled **conditional / open**.

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

### Full Theorem A (**conditional / not proved in this repo**)

**Hypothesis (RH).** All non-trivial zeros of \(\zeta\) satisfy \(\mathrm{Re}\,\rho=\tfrac12\).

**Claim.** For the arithmetic residual \(q_T\) built from \(\theta(x)-x\) (or an
equivalent explicit-formula residual) on the window of length \(T\),
\[
R_d(q_T)\to0\qquad(T\to\infty).
\]

**Status:** Conditional on RH. Outline: under RH the residual is a superposition
of critical-line modes; each contributes \(O(T^{-2})\) (or better) to \(R_d\);
controlling the sum over zeros is analytic number theory **not completed here**.

### Arithmetic multi-\(T\) evidence (measured, not a proof)

Shipped builder: `pbss.probes.arithmetic_residual` — \(q_T=( \theta(x)-x)/\sqrt{x}\) on
\(x=e^{uT}\). Large campaign: `experiments/run_overnight_campaign.py`
→ `results/overnight_campaign.json` (also earlier mid-scale
`results/arithmetic_multi_T.json`).

**Overnight scan (d=4, n=4096, detrend=deg1, smooth=1, \(x_{\max}=10^9\)):**

| \(T\) | \(x_{\max}\) | \(R_4\) arith | \(R_4\) crit-line | \(R_4\) defect |
|------:|-------------:|-------------:|------------------:|---------------:|
| 8.0 | 3e3 | 0.108 | 0.0062 | 0.250 |
| 16.0 | 9e6 | 0.193 | 0.0016 | 0.250 |
| 18.7 | 1.3e8 | 0.180 | 0.0010 | 0.250 |
| 20.7 | **1e9** | **0.165** | 0.0010 | 0.250 |

Ablation at \(T=20.7\): raw residual \(R_4\sim 0.97\); deg0 ~0.20; deg1 ~0.17
(smooth=1). Heavy smooth increases \(R_d\).

**Reading (honest):**

- Arithmetic \(R_d\) (deg1) **rises then soft-plateaus** ~0.16–0.19 through \(x=10^9\).
- **No** A0-style \(O(T^{-2})\) collapse for the arithmetic residual at this scale.
- Always **below** defect floor 0.25, **far above** pure critical-line mode.
- A0 remains proved only for pure modes; full Theorem A stays **open**.

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

### Full Theorem B (**conditional / open**)

**Claim (archive intent).** If \(R_d(q_T)\) decays sufficiently rapidly as
\(T\to\infty\), then \(\zeta\) has no zero with \(\mathrm{Re}\,\rho\neq\tfrac12\).

**Status:** The implication
“fast decay of the prime residual’s \(R_d\) \(\Rightarrow\) RH”
is essentially as hard as RH. **Not proved here.** Off-critical model modes
\(e^{T(\sigma-1/2)u}\sin(tTu)\) are implemented for numerical comparison; a
complete reduction from arithmetic residuals to these modes is open.

---

## 4. What is open

1. Full Theorem A for the **arithmetic** residual under RH (zero-sum estimates).  
2. Full Theorem B (converse / obstruction ⇒ RH).  
3. Sharp rate \(O(T^{-2(d+1)})\) for model or arithmetic residuals.  
4. Legacy normalization with \(P\approx3.92\).  
5. Diamond-system battery and large-\(T\) prime scans beyond current \(x_{\max}\).

---

## 5. Numerical support (multi-\(T\))

See `results/multi_T_scan.json` and `{SCRATCH}/multi_T_scan.*`.

Expected qualitative picture (must match the run):

- Critical-line \(R_d(T)\) **falls** roughly like \(T^{-2}\).  
- Persistent defect \(R_d\equiv\varepsilon^2\) **flat**.  
- Clear gap for large \(T\): defect stays high, RH-like mode drops.

This supports A₀/B₀. It does **not** prove RH.
