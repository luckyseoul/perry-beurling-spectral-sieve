# Infinite-tail zeros and arithmetic remainder (scaffolding note)

**Date:** 2026-07-26  
**Related:** [`THEOREM_A_SCAFFOLD.md`](THEOREM_A_SCAFFOLD.md), [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md) (M5–M6)  
**Code:** `pbss.remainder.bound_R_d_mode_tail`, `bound_infinite_zero_tail_scaffold`

## Explicit non-claim

**Not a proof of RH.** Full Theorem A is **closed conditionally** under RH + cited
ANT-1…3 ([`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md)); this note documents the
**diagnostic scaffold** for infinite tails. The scaffold is **not** the sole support
for Full A — ANT-1 (cited truncated EF under RH) is.

---

## 1. Finite truncation (controlled — model)

Under the model residual
\[
q_T^{(N)}=\sum_{n=1}^{N}a_n\cos(t_n T u-\alpha_n),
\]
Lemma **M5** gives \(R_d(q_T^{(N)})=O_d(T^{-2})\). Peeling via
`peel_via_remainder` strips the first \(N\) modes. Full strip of a pure model
sum leaves \(\approx 0\).

**Proved-style majorant:** `lemmas.bound_R_d_finite_mode_sum`.

---

## 2. Zeros beyond \(N\) (scaffolding only)

### What is controlled

`bound_infinite_zero_tail_scaffold(T, n_kept=N, N_eff=…)` builds a **model**
tail of \(N_{\mathrm{eff}}\) fictitious zeros with
\[
t_n \ge t_{\mathrm{next}}\cdot\frac{n}{N+1},\qquad a_n=\frac{2}{t_n},
\]
and applies the same M5-style majorant. Conclusion under this **model**:
\[
R_d^{\mathrm{model\ tail}}=O_d(T^{-2})
\]
as \(T\to\infty\) for fixed \(N_{\mathrm{eff}}\).

Label in JSON: `scaffolding_model_zero_tail_not_arithmetic_remainder`.

### What is **not** controlled

| Gap | Why it matters |
|-----|----------------|
| True \(\zeta\) zero density | Model \(t_n\) growth is not a theorem about \(\zeta\) |
| Infinite sum \(N_{\mathrm{eff}}\to\infty\) | Need zero-density / large-value estimates to pass to the limit |
| RH / zero-free regions | Off-line zeros break the pure CL mode structure |
| Amplitude law \(a_n=2/t_n\) | Explicit-formula coefficients depend on smoothing and \(\rho\) |

---

## 3. Arithmetic \(\psi\)/\(\theta\) explicit-formula remainder (open)

The classical explicit formula schematically reads
\[
\psi(x)-x
=-\sum_{|\gamma|\le G}\frac{x^\rho}{\rho}
+R_{\mathrm{smooth}}(x;G)
+R_{\mathrm{arith}}(x),
\]
with \(R_{\mathrm{smooth}}\) depending on the smoothing / truncated height \(G\),
and \(R_{\mathrm{arith}}\) collecting prime-power and trivial zeros / contour
contributions (depending on formulation).

**Shipped PBSS residual** \(q_T=(\theta-x)/\sqrt{x}\) is **not** identical to a
fully expanded explicit formula with tracked \(R_{\mathrm{smooth}}\) and
\(R_{\mathrm{arith}}\). Open-plateau peel experiments show that subtracting a
**model** \(q_T^{(N)}\) from the arithmetic residual does **not** collapse \(R_d\).

### Package disposition (2026-08-11)

The checklist items below are **not left unlabeled**. In the Full A package they are
discharged as:

1. **RH** — **hypothesis** of conditional Full A (not proved).  
2. **Truncation / infinite tail** — **cited ANT-1** (truncated EF under RH); this
   module’s scaffold majorant is **diagnostic only**, not sole support.  
3. **Arithmetic remainder** — **cited ANT-2**.  
4. **Weight transfer** — **M6** (model, proved) + optional **cited ANT-4**.  
5. **Triangle / continuity** — **proved M7**.

See [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md). **Full A = closed conditional.**  
**Unconditional RH remains open.** The *in-repo re-derivation* of classical ANT
constants is optional research, not a package gap.

---

## 4. Constants in the scaffolding majorant

For the model tail with \(a_n=2/t_n\) and \(t_n\ge t_{\ast} n/(N+1)\),
\[
\sum_{n=1}^{N_{\mathrm{eff}}}\frac{|a_n|}{t_n T}
\le \frac{2}{T t_{\ast}^2}\sum_{n=1}^{N_{\mathrm{eff}}}\frac{(N+1)^2}{n^2}
\le \frac{C_N}{T},
\]
so the M5-style bound is \(O_d(T^{-2})\) with an explicit (large) constant
depending on \(d,N_{\mathrm{eff}},t_{\ast}\) — available numerically via
`bound_infinite_zero_tail_scaffold`.

**This constant is not claimed for \(\zeta\).**

---

## 5. How to run

```bash
PYTHONPATH=src python3 -c "
from pbss.remainder import bound_infinite_zero_tail_scaffold
print(bound_infinite_zero_tail_scaffold(30.0, n_kept=20, N_eff=5000, d=4))
"
```

**RH remains open.** Full Theorem A is **closed conditionally** (RH + cited ANT); this note is diagnostic scaffolding only.
