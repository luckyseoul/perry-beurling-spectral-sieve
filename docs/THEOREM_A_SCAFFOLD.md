# Theorem A scaffolding: weight class + truncation remainder

**Date:** 2026-07-26  
**Repo:** `perry-beurling-spectral-sieve`  
**Entry points:** `pbss.weights`, `pbss.remainder`, `experiments/run_theorem_a_scaffold.py`  
**Artifacts:** `results/theorem_a_scaffold/`

## Explicit non-claim

**This document does not prove the Riemann Hypothesis and does not claim full Theorem A
for the arithmetic residual.**  
What is shipped is *scaffolding*: admissible weights with checkable endpoint estimators,
and a zero-sum / truncation remainder path with proved-style majorants for **finite**
mode blocks (Lemma M5) plus **scaffolding** tail majorants.

---

## 1. Weight class / endpoint control (item 1)

### Problem (from open-plateau)

On \(q_T=(\theta-x)/\sqrt{x}\) with deg1 detrend, \(R_4\sim 0.15\)–\(0.19\) through
\(x\sim5\times10^{10}\). A Hanning taper dropped mean \(R_4\) to \(\sim0.08\) — evidence
that **endpoint / window pollution** feeds low-degree Legendre mass. Ad-hoc taper is not
a theorem; Theorem A needs a controlled weight class.

### Admissible class \(W_\alpha\)

A continuous weight \(w:[0,1]\to[0,1]\) is **admissible** (\(w\in W_\alpha\), \(\alpha\in(0,1/2]\)) if:

1. \(w\ge 0\) and \(\|w\|_{L^2}>0\)
2. \(w\) vanishes (or is strongly suppressed) near endpoints \(u\in[0,\alpha)\cup(1-\alpha,1]\)
3. \(w\) is sufficiently regular on the bulk \((\alpha,1-\alpha)\) for integration-by-parts
   against critical-line modes (C¹ / Lipschitz bulk)

**Shipped members** (`pbss.weights`):

| Name | Description |
|------|-------------|
| `tukey` | Cosine taper of half-width \(\alpha\) at each end (default \(\alpha=0.1\)) |
| `hanning` | Full Hann window (member of \(W_{1/2}\)) |
| `flat` / `none` | \(w\equiv 1\) (baseline, not endpoint-killing) |

### Endpoint contribution estimator

For residual sample \(q\) on grid \(u\),

\[
E_{\mathrm{end}}=\frac{\|P_d(q\cdot\mathbf{1}_{\mathrm{end}})\|_2^2}{\|q\|_2^2},\qquad
E_{\mathrm{bulk}}=\frac{\|P_d(q\cdot\mathbf{1}_{\mathrm{bulk}})\|_2^2}{\|q\|_2^2}.
\]

API: `endpoint_contribution`, `bulk_vs_weighted_report`, `weighted_energy_ratio`.

**What this buys for Theorem A:** a precise object (“endpoint mass feeding \(V_d\)”) that
must be controlled analytically (or by choosing \(w\in W_\alpha\)) before low \(R_d\) can be
attributed to zero structure alone.

**What it does *not* buy:** \(R_d(Wq)\to0\) for the arithmetic residual under RH. Weighted
\(R_d\) can drop while the bulk still plateaus (open-plateau).

### Multi-\(T\) evidence (`results/theorem_a_scaffold/weight/`)

Campaign: CL pure modes + truncated EF residuals × multi-\(T\) × weights
`flat` / `tukey` / `hanning`. Resume stamp `PHASE_WEIGHT_COMPLETE`.

Reading: pure CL modes keep small \(R_d\) under weights (M3-compatible). Bulk-sensitive
synthetics (unit tests) show \(R_d^{\mathrm{weighted}}<R_d^{\mathrm{raw}}\) when endpoints dominate.

---

## 2. Zero-sum / truncation remainder (item 2)

### Problem

Lemma **M5** (proved): any **finite** critical-line superposition \(q_T^{(N)}\) has
\(R_d=O_d(T^{-2})\). Full Theorem A needs the **infinite** zero sum plus arithmetic
remainders. Open-plateau peel numerics showed stripping more model zeros from the
*arithmetic* residual does **not** collapse \(R_d\) — so the gap is the tail + true
arithmetic error, not missing low zeros alone.

### Truncation path (shipped)

\[
q_T^{(N)}(u)=\sum_{n=1}^{N}a_n\cos(t_n T u-\alpha_n),\qquad
q_{\mathrm{rem}}=q-\alpha\,q_T^{(N)}.
\]

| API | Role |
|-----|------|
| `truncated_mode_sum` | \(q_T^{(N)}\) via `explicit_formula_residual` |
| `peel_via_remainder` | \(q-\alpha q_T^{(N)}\) (optional LS \(\alpha\)) |
| `bound_R_d_mode_tail` | M5-style majorant for a *model* tail of zeros beyond \(N\) |
| `remainder_diagnostic` | Multi-\((T,N)\) row: \(R_d\) full/rem, M5 bound, tail majorant |

**Proved-style:** M5 majorant on any finite stripped block
(`lemmas.bound_R_d_finite_mode_sum`).

**Scaffolding only:** tail majorant for zeros past the table (crude \(a_n=2/t_n\)
extrapolation) — **not** a sharp analytic number theory bound, **not** the full
explicit-formula remainder (prime powers, contours, …).

### Multi-\((T,N)\) evidence (`results/theorem_a_scaffold/remainder/`)

- Full strip \(N_{\mathrm{strip}}=N_{\mathrm{full}}\) on the model sum → residual \(\approx0\)
- M5 bound on stripped block decays as \(T\) grows
- Tail majorant decays in \(T\) but remains a **label: scaffolding_majorant_not_sharp**

Optional light arithmetic peel/weight rows under `arithmetic/` when primes are available.

---

## 3. Open hypotheses for full Theorem A

To upgrade scaffolding → arithmetic Theorem A (still **conditional on RH** for the
classical statement), one still needs:

1. **Infinite zero-sum control** under RH: show the tail \(\sum_{n>N}\) contributes
   \(o(1)\) (or \(O(T^{-2})\)) to \(R_d\) uniformly after admissible weighting — beyond
   the crude majorant here.
2. **Arithmetic remainder** in the explicit formula (\(\psi\) or \(\theta\) vs truncated
   zeros): prime-power and contour terms with explicit constants.
3. **Weight-class theorem:** for \(w\in W_\alpha\), \(R_d(Wq_T^{\mathrm{arith}})\to0\) under RH
   (or an equivalent bulk residual after endpoint removal).
4. **Bridge** from diagnostic \(R_d\) to classical zero-free / \(\psi\)-error language.

Until (1)–(4) exist in writing with proofs, **full Theorem A remains open**.

---

## 4. Relation to open-plateau

| Open-plateau finding | Scaffolding response |
|----------------------|----------------------|
| Taper lowers \(R_d\) | Formal \(W_\alpha\) + \(E_{\mathrm{end}}\) estimator |
| Peel does not kill arithmetic plateau | Remainder path + M5 on finite block + tail majorant labels |
| Dense multi-\(T\) no decay | Arithmetic still open; model CL/EF decay verified |

Pointer: [`RESEARCH_PLATEAU.md`](RESEARCH_PLATEAU.md) · status: [`STATUS.md`](STATUS.md) ·
definitions: [`THEOREMS_AB.md`](THEOREMS_AB.md).

---

## 5. How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/test_weights.py tests/test_remainder.py -v
PYTHONPATH=src python3 experiments/run_theorem_a_scaffold.py \
  --out-dir results/theorem_a_scaffold --workers 32
# smoke (separate out-dir):
PYTHONPATH=src python3 experiments/run_theorem_a_scaffold.py \
  --out-dir results/theorem_a_scaffold_smoke --phases weight,remainder \
  --T-list 15,25 --n-strips 0,5,10 --n-zeros 10 --workers 4 --n-points 1024
```

**RH remains open.** Full Theorem A is **closed conditionally** — see [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md). This scaffold note is historical path documentation.
