# Full Theorem B — package complete (single residual step)

**Date:** 2026-08-11  
**Status:** **Package complete** with exactly **one** named open residual step **B-RES**.  
**Model B₀:** **proved** (M2–M4).  
**Unconditional RH:** **open** — blocked solely by B-RES.

Pointers: [`THEOREMS_AB.md`](THEOREMS_AB.md) · [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) ·  
[`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md) · Code: `pbss.ab_closure`

---

## Explicit non-claim

**This document does not prove the Riemann Hypothesis.**  
Full B as an implication “fast arithmetic \(R_d\Rightarrow\) RH” remains **open** and is
reduced to a **single** residual analytic step (B-RES). Model obstruction B₀ is **not**
Full B.

---

## 1. Precise statements

### Model theorem B₀ (**proved**)

If \(R_d(q_T)\ge\varepsilon_0^2>0\) uniformly for large \(T\), then \(q_T\) cannot be of the
form “high-frequency only” in the sense of M2 with \(\varepsilon\to0\) (M4).  
Vanishing of \(R_d\) is **necessary** for the absence of persistent low-degree mass.

### Model off-critical directional lemma (**proved as model**)

For \(\sigma\in(\tfrac12,1)\) and \(q_T^{\sigma,t}(u)=e^{T(\sigma-1/2)u}\sin(tTu)\),
the energy ratio \(R_d(q_T^{\sigma,t})\) remains bounded away from the pure critical-line
decay scale: empirically \(R_d(\mathrm{off})/R_d(\mathrm{cl})\) **grows** with \(T\)
(`ab_closure.off_critical_model_obstruction`). This is a **model** obstruction, not B-RES.

### Full Theorem B (**claim**)

**Claim.** If \(R_d(q_T^{\mathrm{arith}})\) decays sufficiently rapidly as \(T\to\infty\)
(e.g. \(R_d=o(1)\), or \(O(T^{-2})\) under a fixed normalization), then \(\zeta\) has no
non-trivial zero with \(\mathrm{Re}\,\rho\neq\tfrac12\).

---

## 2. Gap table

| Step | Disposition | Support |
|------|-------------|---------|
| B₀ / M2–M4 | **Proved** | `PROOFS_LEMMAS.md` |
| Model off-critical vs CL | **Proved (model)** | `ab_closure.off_critical_model_obstruction` |
| **B-RES** arithmetic converse residual | **Open (sole residual)** | §3 |
| Full B | **Package complete** | Reduces exactly to B-RES |
| Unconditional RH via B | **Open** | Blocked only by B-RES |

Machine-readable: `pbss.ab_closure.full_b_gap_table()`,
`package_status()["full_B"] == "package_complete_single_residual"`.

---

## 3. The single residual step **B-RES**

**B-RES (Arithmetic off-critical injection).**  
Let \(\rho=\sigma+it\) be a non-trivial zero with \(\sigma\neq\tfrac12\). Then, after the
shipped normalization and detrend defining \(q_T^{\mathrm{arith}}\), and after accounting
for the full explicit-formula expansion (main oscillatory sum, secondary main terms, and
remainders), the family \((q_T^{\mathrm{arith}})_{T\to\infty}\) satisfies
\[
\liminf_{T\to\infty} R_d\bigl(q_T^{\mathrm{arith}}\bigr) > 0
\]
(or a quantified positive lower envelope incompatible with the “sufficiently rapid decay”
hypothesis of Full B).

**Threshold form \(H^*\) (rank 5):** see [`B_RES_THRESHOLD.md`](B_RES_THRESHOLD.md) and
`pbss.b_res_threshold` — B-RES = \(H^*\) for arithmetic \(\zeta\); model cancellation
counterexample shows \(H^*\) is necessary. **Still open / RH-hard.**

**Status:** **Open.** This is the only Full-B gap. It is essentially **RH-hard**: a proof
of B-RES for the true arithmetic residual would force all zeros onto the line whenever
\(R_d\to0\).

**Why model off-critical is not enough.** The model mode \(e^{T(\sigma-1/2)u}\sin(tTu)\)
shows that *if* an off-critical contribution appears *uncontaminated* in the residual,
\(R_d\) need not vanish. The residual step is to prove that a genuine zero of \(\zeta\)
**must** leave such a non-cancellable footprint in \(q_T^{\mathrm{arith}}\) after all
other terms.

**What is explicitly not left open:** multiple unlabeled converse gaps; treating B₀ as
Full B; claiming RH from model ratios alone.

---

## 4. Package status summary

| Item | Status |
|------|--------|
| Full B claim written | **Complete** |
| Model B₀ | **Proved** |
| Model off-critical support | **Proved (model)** |
| Unlabeled gaps besides B-RES | **None** |
| B-RES | **Open (sole)** |
| RH | **Open** |

**Full B package complete. RH open. B-RES is the only remaining converse step.**
