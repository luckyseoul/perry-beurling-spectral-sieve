# Roadmap to close out RH (aspirational — not a proof)

**Updated:** 2026-09-05 (original roadmap: 2026-07-26 / 2026-08-11)
**Repo:** Perry–Beurling Spectral Sieve (PBSS)  
**Status:** Exact continuous A is proved under RH; B-ABS is a proved RH
equivalence. The normalized converse and the unconditional bound needed to
deduce RH remain open.

## Explicit non-claim

**This roadmap is not a proof of the Riemann Hypothesis.**  
**RH is not closed, not solved, and not claimed.**  
No milestone below is marked “done” for unconditional RH. Completing every
milestone would still require independent verification and community review
before any serious claim about \(\zeta\).

Current proofs: [`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md),
[`ABSOLUTE_MOMENT_CONVERSE.md`](ABSOLUTE_MOMENT_CONVERSE.md), and
[`NORMALIZED_RECURRENCE_THEOREM.md`](NORMALIZED_RECURRENCE_THEOREM.md).
They use exact continuous moments and continuous orthogonal detrending. The
older milestone wording below does not certify a transfer to fixed sampled
arrays or arbitrary weights.

---

## Current baseline (PBSS, 2026-09-05)

| Asset | Status |
|-------|--------|
| Diagnostic \(R_d\) on log-window | Shipped, tested |
| Model A₀ / finite-mode A₀ / weighted model A₀ | **Proved** (M3, M5, M6) |
| Continuous Theorem A under RH | **Proved:** \(R_d((I-P_1)Q(T\cdot))=O_d(T^{-2})\), fixed \(d\ge2\) |
| Absolute-moment converse B-ABS | **Proved RH equivalence**; its unconditional growth bound remains open |
| Normalized Theorem B | **Open**; see [`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md) |
| Recurrence of small continuous normalized ratios | **Proved unconditionally:** \(\liminf R_d=0\) for the same detrended residual |
| Historical arithmetic multi-\(T\) campaign through \(5\times10^{10}\) | Reported plateau \(R_d\sim0.15\)–\(0.19\); historical finite-sample evidence, not a continuum theorem |
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

### M1 — Arithmetic explicit-formula identification (ANT-3) — **continuous scope proved under RH**

**Goal.** Prove (or cite a standard theorem with full constants adapted to PBSS)
that the shipped arithmetic residual differs from a truncated explicit-formula
mode sum by a remainder controlled in the \(R_d\) (or \(L^2\)) metric on the
log-window.

**Disposition (2026-09-05):** The exact continuous ψ-to-θ calculation is in
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md), using explicit source
equations. The original target of a remainder theorem for a shipped sampled or
arbitrarily weighted residual is broader and remains a separate extension.

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

### M2 — Infinite zero-sum / height truncation (ANT-1) — **continuous moment estimate proved under RH**

**Goal.** Under RH, choose \(G=G(T)\) so zeros with \(|\gamma|>G\) contribute
\(o(1)\) (or \(O(T^{-2})\)) to \(R_d(w\,\cdot)\).

**Disposition (2026-09-05):** Continuous A uses cutoff \(G=e^T\) and bounds the
full retained polynomial moments uniformly through a convergent coefficient
sum. This does not assert the broader arbitrary-weight relative tail target
above. The scaffold majorant remains diagnostic.

**Success criteria**

- Theorem under RH + standard zero-density / large-value estimates (cited).  
- Replace scaffolding `bound_infinite_zero_tail_scaffold` with a bound whose
  hypotheses are standard and fully listed.  
- Numeric checks: model tails decrease in \(T\) (already); arithmetic remains
  diagnostic only until M1 holds.

**Depends on:** RH, the cited explicit formula, and the standard zero-count bound.
**Does not prove RH.**

---

### M3 — Arithmetic remainder \(R_{\mathrm{arith}}\) (ANT-2) — **specified continuous remainder proved under RH**

**Goal.** Bound prime-power / contour / trivial-zero contributions in the same
window after weighting.

**Disposition (2026-09-05):** Continuous A explicitly removes the prime-square
constant and obtains an \(L^1\cap L^2\) remainder in logarithmic coordinates.
Its estimates apply to the stated continuous residual and detrend; other
smoothing and weighting choices require their own argument.

**Success criteria**

- Explicit majorant \(\to0\) as \(T\to\infty\) (or \(O(T^{-2})\)).  
- Documented dependence on smoothing parameters matching the residual definition.  

**Depends on:** M1 identification.  
**Does not prove RH.**

---

### M4 — Continuous Theorem A under RH — **PROVED**

**Result.** For \(Q(y)=e^{-y/2}(\theta(e^y)-e^y)\), continuous orthogonal
degree-one detrending, and every fixed \(d\ge2\):

> **Assume RH (+ listed cited ANT inputs). Then**  
> \(R_d((I-P_1)Q(T\cdot))=O_d(T^{-2})\).

**Disposition (2026-09-05):** See
[`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md) for the full proof,
including a positive limiting denominator. RH is its hypothesis.

---

### M5 — Bridge to classical consequences / Full B

**Goal.** Show that sufficiently fast decay of \(R_d(w q_T^{\mathrm{arith}})\)
implies RH (Full B).

**Disposition (2026-09-05):** B-ABS proves that one fixed absolute continuous
Legendre moment is subexponential exactly when RH holds. A converse from a
specified normalized decay rate remains open. The old positive-floor B-RES
implication is itself RH-equivalent, because the same continuous arithmetic
residual has unconditional \(\liminf R_d=0\). See
[`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md).

A converse alone is insufficient to deduce RH: its matching growth or decay
condition must also be proved unconditionally. This is an additional proof
obligation, even if a normalized equivalence is eventually established.

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
An unconditional argument is required; completing an equivalence does not
establish either side of it.

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
              Unconditional matching growth/decay bound [OPEN]
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
| M4 | Continuous A under RH (+ listed cited inputs) | Still open |
| M5 | A converse or equivalence | Still open; an unconditional matching bound is also required |
| M6 | Review quality | Still open |
| M7 | Actual RH proof accepted | Closed only if M7 succeeds |

---

## Immediate next actions (practical)

1. Pursue an unconditional bound for the absolute criterion in B-ABS, or state
   a different exact hypothesis and prove both its converse and its bound.
2. Keep sampled and weighted transfers separate from the proved continuous scope.
3. Use the phase obstruction and recurrence theorem when formulating normalized targets.
4. Maintain STATUS non-claim on every release.

---

## Final banner

**RH is not closed. This is a roadmap, not a finish line.**  
**Conditional continuous theorem:** see [`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md).
