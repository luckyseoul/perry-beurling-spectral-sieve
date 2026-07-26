# Roadmap to close out RH (aspirational — not a proof)

**Date:** 2026-07-26  
**Repo:** Perry–Beurling Spectral Sieve (PBSS)  
**Status:** Roadmap only.

## Explicit non-claim

**This roadmap is not a proof of the Riemann Hypothesis.**  
**RH is not closed, not solved, and not claimed.**  
No milestone below is marked “done” for unconditional RH. Completing every
milestone would still require independent verification and community review
before any serious claim about \(\zeta\).

Related: [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) (conditional Full A package).

---

## Current baseline (PBSS, 2026-07-26)

| Asset | Status |
|-------|--------|
| Diagnostic \(R_d\) on log-window | Shipped, tested |
| Model A₀ / finite-mode A₀ / weighted model A₀ | **Proved** (M3, M5, M6) |
| Conditional Full Theorem A statement + labeled gaps | **Package complete** |
| Arithmetic multi-\(T\) through \(5\times10^{10}\) | Soft plateau \(R_d\sim0.15\)–\(0.19\) |
| Beurling ordinary vs defective separation | Holds numerically |
| Unconditional RH | **Open** |

---

## Milestone sequence

### M0 — Freeze claims (done in package)

**Success criteria**

- [x] Conditional Full A written with proved / assumed / open labels  
- [x] Explicit RH non-claim in STATUS and package docs  
- [x] Numerics not misread as A₀ for arithmetic residual  

**Exit:** No status drift claiming “A done” or “RH closed.”

---

### M1 — Arithmetic explicit-formula identification (ANT-3)

**Goal.** Prove (or cite a standard theorem with full constants adapted to PBSS)
that the shipped arithmetic residual differs from a truncated explicit-formula
mode sum by a remainder controlled in the \(R_d\) (or \(L^2\)) metric on the
log-window.

**Success criteria**

- Written theorem: \(w q_T^{\mathrm{arith}} = w q_T^{(N)} + r_{N,T}\) with
  \(\|P_d r_{N,T}\|/\|w q_T^{\mathrm{arith}}\|\to0\) under specified \(N=N(T)\).  
- Constants trackable in code or tables.  
- Unit/integration tests only where the identity is model-level; ANT citations
  labeled **assumed/cited**, not “proved in PBSS.”

**Depends on:** Classical explicit formula literature (e.g. \(\psi\)-form with
smoothing).  
**Does not imply RH.**

---

### M2 — Infinite zero-sum / height truncation (ANT-1)

**Goal.** Under RH, choose \(G=G(T)\) so zeros with \(|\gamma|>G\) contribute
\(o(1)\) (or \(O(T^{-2})\)) to \(R_d(w\,\cdot)\).

**Success criteria**

- Theorem under RH + standard zero-density / large-value estimates (cited).  
- Replace scaffolding `bound_infinite_zero_tail_scaffold` with a bound whose
  hypotheses are standard and fully listed.  
- Numeric checks: model tails decrease in \(T\) (already); arithmetic remains
  diagnostic only until M1 holds.

**Depends on:** RH (hypothesis) + zero-density estimates (external).  
**Does not prove RH.**

---

### M3 — Arithmetic remainder \(R_{\mathrm{arith}}\) (ANT-2)

**Goal.** Bound prime-power / contour / trivial-zero contributions in the same
window after weighting.

**Success criteria**

- Explicit majorant \(\to0\) as \(T\to\infty\) (or \(O(T^{-2})\)).  
- Documented dependence on smoothing parameters matching the residual definition.  

**Depends on:** M1 identification.  
**Does not prove RH.**

---

### M4 — Conditional Full Theorem A closed under RH

**Goal.** Combine M1–M3 + proved M3/M5/M6 into a single theorem:

> **Assume RH (+ listed zero-density inputs). Then**  
> \(R_d(w q_T^{\mathrm{arith}})\to0\).

**Success criteria**

- Proof writeup with no unlabeled gaps.  
- STATUS updates: “Full A **proved conditional on RH + (list)**” — still  
  **not** “RH proved.”  
- Independent proof review (internal or external).

**This still does not prove RH.** It only completes the PBSS conditional A.

---

### M5 — Bridge to classical consequences (optional for “PBSS RH path”)

**Goal.** Show that sufficiently fast decay of \(R_d(w q_T^{\mathrm{arith}})\)
implies a classical zero-free / \(\psi\)-error statement strong enough for RH
(or a known equivalent). This is essentially **Theorem B** direction and is
**RH-hard**.

**Success criteria**

- Precise implication theorem with constants.  
- Clear separation: “decay \(\Rightarrow\) RH” vs “RH \(\Rightarrow\) decay (Full A).”  

**Status today:** **Open** (archive Theorem B).  
**Risk:** May be circular or as hard as RH; treat as research, not a short goal.

---

### M6 — Independent verification & adversarial review

**Goal.** Before any public claim about \(\zeta\):

**Success criteria**

- [ ] Cold read of Full A conditional proof by someone who did not write it  
- [ ] Reproduction of all numeric campaigns from scripts + checkpoints  
- [ ] Explicit search for overclaim language (no “RH solved”)  
- [ ] Optional formalization of model lemmas (Lean/Isabelle) — optional stretch  

---

### M7 — Unconditional RH (external to PBSS alone)

**Goal.** Prove RH by some path (PBSS-mediated or classical).

**Success criteria**

- Peer-reviewed / community-accepted proof.  
- PBSS may supply **motivation and diagnostics** but must not replace peer review.

**Status:** **Open.** This roadmap **does not claim** a path length or ETA.  
**Non-goal of day-to-day PBSS engineering.**

---

## Dependency diagram (ASCII)

```
M0 freeze claims
    │
    ▼
M1 arithmetic EF identification ──► M3 arithmetic remainder
    │                                    │
    ▼                                    │
M2 infinite zeros under RH ──────────────┤
    │                                    │
    └──────────────► M4 conditional Full A (under RH + ANT)
                         │
                         ▼
                    M5 classical bridge (optional, RH-hard)
                         │
                         ▼
                    M6 independent review
                         │
                         ▼
                    M7 unconditional RH  [OPEN — not claimed]
```

---

## Per-milestone “done” vs “RH closed”

| Milestone | Completing it means… | RH status after |
|-----------|----------------------|-----------------|
| M0–M3 | Engineering / conditional ANT packages | Still open |
| M4 | Full A under RH (+ listed inputs) | Still open |
| M5 | Possible equivalence fragment | Still open unless equivalence is full RH |
| M6 | Review quality | Still open |
| M7 | Actual RH proof accepted | Closed only if M7 succeeds |

---

## Immediate next actions (practical)

1. Prioritize **M1** literature alignment with the exact residual
   \((\theta-x)/\sqrt{x}\) on \(x=e^{uT}\).  
2. Keep scaffolding tails clearly labeled until M2 has real citations.  
3. Do **not** treat arithmetic weight decay (mild) as A₀.  
4. Maintain STATUS non-claim on every release.

---

## Final banner

**RH is not closed. This is a roadmap, not a finish line.**  
**Conditional Theorem A package:** see [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md).
