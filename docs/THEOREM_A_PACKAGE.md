# Conditional Theorem A — complete package

**Date:** 2026-07-26  
**Status:** **Package complete** (statement + lemma chain + labeled gaps).  
**Arithmetic Full A:** still **open** (external ANT inputs missing).  
**RH:** **open** — this package does **not** prove RH.

Pointers: [`THEOREMS_AB.md`](THEOREMS_AB.md) · [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md) ·  
[`THEOREM_A_SCAFFOLD.md`](THEOREM_A_SCAFFOLD.md) · [`INFINITE_TAIL_REMAINDER.md`](INFINITE_TAIL_REMAINDER.md) ·  
[`RH_CLOSEOUT_ROADMAP.md`](RH_CLOSEOUT_ROADMAP.md) · Code: `pbss.theorem_a_chain`

---

## Explicit non-claim

**This document does not contain an unconditional proof of the Riemann Hypothesis.**  
It does **not** claim that open-plateau arithmetic multi-\(T\) numerics
(\(R_d\sim 0.15\)–\(0.19\) through \(x\sim5\times10^{10}\)) complete Theorem A.
Those numbers show the **current residual probe has not collapsed**; they neither
prove nor disprove RH.

“Complete Theorem A” **in this package** means: a finished **conditional**
statement under RH, with every step labeled **proved / assumed / open**, and a
checkable model-chain exercise in-repo. It does **not** mean arithmetic A is proved.

---

## 1. Precise statement of Full Theorem A (conditional)

### Definitions (as in THEOREMS_AB)

- Log window: \(x=e^{uT}\), \(u\in[0,1]\).
- Orthonormal basis \(\{\varphi_k\}\) shifted Legendre; \(P_d\) projection onto \(V_d\).
- Energy ratio \(R_d(q)=\|P_d q\|_2^2/\|q\|_2^2\).
- **Arithmetic residual** (shipped):  
  \(q_T^{\mathrm{arith}}=\mathrm{detrend}\bigl((\theta(e^{uT})-e^{uT})/\sqrt{e^{uT}}\bigr)\)  
  with default deg1 detrend (`probes.arithmetic_residual`).
- Optional **admissible weight** \(w\in W_\alpha\) (`pbss.weights`).

### Theorem A (Full, conditional)

**Assume RH:** every non-trivial zero of \(\zeta\) has \(\mathrm{Re}\,\rho=\tfrac12\).

**Claim.** For the arithmetic residual \(q_T^{\mathrm{arith}}\) (or an equivalent
explicit-formula residual on the same log-window, after a fixed admissible weight
\(w\in W_\alpha\) if used),
\[
R_d\bigl(w\,q_T^{\mathrm{arith}}\bigr)\to 0\qquad(T\to\infty)
\]
(for each fixed degree \(d\)).

**Rate (target, not proved for arithmetic):** \(O_d(T^{-2})\) or better under
sufficient zero-density / remainder control — matching model A₀ order.

---

## 2. Lemma chain (what is already proved in-repo)

| Step | Lemma / object | Status | Role in A |
|------|----------------|--------|-----------|
| 1 | M1 pure mode energy | **Proved** | Diagnostic well-posed on \(V_d\) |
| 2 | M2 orthogonal defect | **Proved** | Identifies low-degree mass \(\varepsilon^2\) |
| 3 | M3 pure CL decay | **Proved** | Model A₀: \(R_d(\sin(tTu))=O(T^{-2})\) |
| 4 | M4 defect blocks vanishing | **Proved** | Necessity of \(\varepsilon\to0\) for \(R_d\to0\) |
| 5 | M5 finite CL sum | **Proved** | Truncated EF residual A₀ |
| 6 | M6 weighted model decay | **Proved** | \(R_d(w q)=O(T^{-2})\) for model CL/EF |
| 7 | Finite peel / remainder path | **Shipped + tested** | `peel_via_remainder`, M5 majorant |
| 8 | Infinite zero tail | **Scaffolding only** | `bound_infinite_zero_tail_scaffold` |
| 9 | Arithmetic \(\psi\)/\(\theta\) remainder | **Open** | Not a proved bound in-repo |
| 10 | Arithmetic \(R_d\to0\) under RH | **Open** | Full A — needs 8–9 + EF identification |

### Conditional deduction sketch (if external inputs hold)

Under RH, write (schematically, after smoothing)
\[
q_T^{\mathrm{arith}}
= q_T^{(N)} + r_T^{\mathrm{tail}} + r_T^{\mathrm{arith}},
\]
where \(q_T^{(N)}\) is a finite critical-line mode sum.

1. **M5 + M6:** \(R_d(w q_T^{(N)})=O_d(T^{-2})\).  
2. **Assume (ANT-1):** zero-density / truncation so \(r_T^{\mathrm{tail}}\) contributes
   \(o(1)\) (or \(O(T^{-2})\)) to \(R_d(w\,\cdot)\).  
3. **Assume (ANT-2):** arithmetic EF remainder \(r_T^{\mathrm{arith}}\) similarly controlled.  
4. **Triangle / projection continuity:** then \(R_d(w q_T^{\mathrm{arith}})\to0\).

Steps 2–3 are **not proved here**. Without them, Full A remains open even under RH
as a formal implication from shipped lemmas alone.

---

## 3. Label dictionary

| Label | Meaning |
|-------|---------|
| **Proved** | Continuous \(L^2\) argument in `PROOFS_LEMMAS.md` + discrete tests |
| **Shipped majorant** | Explicit crude bound consistent with a proved lemma (not sharp) |
| **Assumed (ANT)** | Standard-style analytic number theory input **not** proved in-repo |
| **Scaffolding** | Heuristic model (e.g. \(a_n=2/t_n\) tail) — **not** a theorem about \(\zeta\) |
| **Open** | No complete proof path in this repository |
| **Measured** | Numeric campaign only; not a proof |

---

## 4. External analytic inputs still required

| ID | Input | Label | Success criterion for “A closed under RH” |
|----|--------|-------|-------------------------------------------|
| ANT-1 | Control of zeros with \(|\gamma|>G(T)\) in the window residual | Assumed / open | Explicit \(G(T)\) + bound \(\to0\) in \(R_d\) |
| ANT-2 | True \(\psi\) or \(\theta\) explicit-formula remainder (prime powers, contours) | Assumed / open | Constants + \(T\to\infty\) rate |
| ANT-3 | Identification: arithmetic \(q_T\) equals EF sum + remainders up to \(o_{R_d}(1)\) | Open | Norm equivalence on log-window |
| ANT-4 | Optional: weight-class theorem for **arithmetic** \(Wq_T\) (M6 is model-only) | Open | Same decay after \(w\in W_\alpha\) |

See [`INFINITE_TAIL_REMAINDER.md`](INFINITE_TAIL_REMAINDER.md) for (ANT-1)–(ANT-2) detail.

---

## 5. Checkable model chain (in-repo)

API: `pbss.theorem_a_chain.model_chain_report(T, …)`.

For each \(T\), the report returns:

- Empirical \(R_d\) of pure CL, finite EF, and weighted versions (shipped energy + weights).
- **Proved-style majorants** M5 / M6 (must dominate empirical for large \(T\)).
- Scaffolding tail majorant (labeled, not arithmetic).
- Flags: `proved_model_decay_ok`, `full_arithmetic_A_status="open"`, RH non-claim banner.

Multi-\(T\) artifact: `results/theorem_a_model_chain/` (optional campaign).  
Tests: `tests/test_theorem_a_chain.py`.

---

## 6. Relation to numerics (do not misread)

| Campaign | Finding | Does it complete A? |
|----------|---------|---------------------|
| Grand / overnight / extend-\(x\) | Arith \(R_4\sim0.15\)–\(0.19\) | **No** — plateau, not collapse |
| Open-plateau peel | More zeros ≠ lower arith \(R_d\) | **No** — supports remainder gap |
| Arithmetic weights | Tukey lowers \(R_d\) somewhat | **No** — still not A₀; endpoint clue only |
| Model CL / EF / M5–M6 tests | Decay under majorants | **Yes for models only** |

---

## 7. Package status summary

| Item | Status |
|------|--------|
| Full A statement under RH | **Complete** (this document §1) |
| Model lemma chain M1–M6 | **Proved** |
| Conditional deduction sketch | **Complete** (§2) |
| ANT-1…ANT-4 filled with proofs | **Open** |
| Unconditional RH | **Open** (non-goal of this package) |
| RH close-out path | [`RH_CLOSEOUT_ROADMAP.md`](RH_CLOSEOUT_ROADMAP.md) |

**Package complete. Full arithmetic Theorem A incomplete. RH open.**
