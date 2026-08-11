# Roadmap to close out RH (aspirational — not a proof)

**Date:** 2026-08-11 (Full A conditional closed; Full B → B-RES) / 2026-07-26  
**Repo:** Perry–Beurling Spectral Sieve (PBSS)  
**Status:** Roadmap. Full A **closed conditionally**; Full B packaged to **B-RES** only.

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
| Full Theorem A under RH + cited ANT | **Closed conditionally** ([`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md)) |
| Full Theorem B | **Package complete** — sole gap **B-RES** ([`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md)) |
| Arithmetic multi-\(T\) through \(5\times10^{10}\) | Soft plateau \(R_d\sim0.15\)–\(0.19\) (finite \(T\); not a refutation of conditional A) |
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

### M1 — Arithmetic explicit-formula identification (ANT-3) — **cited / closed as package step**

**Goal.** Prove (or cite a standard theorem with full constants adapted to PBSS)
that the shipped arithmetic residual differs from a truncated explicit-formula
mode sum by a remainder controlled in the \(R_d\) (or \(L^2\)) metric on the
log-window.

**Disposition (2026-08-11):** **Cited** (Davenport / Ingham / Titchmarsh / Ivić) with
hypotheses and adapted conclusion in [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) §3.

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

### M2 — Infinite zero-sum / height truncation (ANT-1) — **cited / closed as package step**

**Goal.** Under RH, choose \(G=G(T)\) so zeros with \(|\gamma|>G\) contribute
\(o(1)\) (or \(O(T^{-2})\)) to \(R_d(w\,\cdot)\).

**Disposition (2026-08-11):** **Cited** truncated-EF under RH; scaffold majorant is
diagnostic only (not sole support).

**Success criteria**

- Theorem under RH + standard zero-density / large-value estimates (cited).  
- Replace scaffolding `bound_infinite_zero_tail_scaffold` with a bound whose
  hypotheses are standard and fully listed.  
- Numeric checks: model tails decrease in \(T\) (already); arithmetic remains
  diagnostic only until M1 holds.

**Depends on:** RH (hypothesis) + zero-density estimates (external).  
**Does not prove RH.**

---

### M3 — Arithmetic remainder \(R_{\mathrm{arith}}\) (ANT-2) — **cited / closed as package step**

**Goal.** Bound prime-power / contour / trivial-zero contributions in the same
window after weighting.

**Disposition (2026-08-11):** **Cited** classical EF remainders + \(\psi-\theta\).

**Success criteria**

- Explicit majorant \(\to0\) as \(T\to\infty\) (or \(O(T^{-2})\)).  
- Documented dependence on smoothing parameters matching the residual definition.  

**Depends on:** M1 identification.  
**Does not prove RH.**

---

### M4 — Conditional Full Theorem A closed under RH — **DONE (package)**

**Goal.** Combine M1–M3 + proved M5/M6/M7 into a single theorem:

> **Assume RH (+ listed cited ANT inputs). Then**  
> \(R_d(w q_T^{\mathrm{arith}})\to0\).

**Disposition (2026-08-11):** **Closed conditionally** — see
[`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md) §4. STATUS:
`full_arithmetic_A = closed_conditional`. Still **not** “RH proved.”

---

### M5 — Bridge to classical consequences / Full B

**Goal.** Show that sufficiently fast decay of \(R_d(w q_T^{\mathrm{arith}})\)
implies RH (Full B).

**Disposition (2026-08-11):** Full B **package complete** with sole residual step
**B-RES** ([`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md)). Model B₀ proved; B-RES open
(RH-hard). Clear separation: “decay \(\Rightarrow\) RH” needs B-RES; “RH \(\Rightarrow\) decay”
is Full A (closed conditional).

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
