# Open-ended goal: Why does arithmetic \(R_d\) plateau?

## Intent
Not a checklist of trial counts. **Stay in research mode** until there is a defensible written answer (with multi-T evidence) to:

> Why does the arithmetic residual \(q_T=(\theta-x)/\sqrt{x}\) keep \(R_d\sim O(10^{-1})\) through \(x\sim5\times10^{10}\), while pure/finite CL modes collapse under the same diagnostic?

Finish only when the **research bars** below are met. Do **not** pad wall-clock. Do **not** declare victory after one dense peel/MC batch.

## Explicit non-claim
No unconditional RH proof. No “full Theorem A done.” Open conclusions are allowed if evidence supports them.

## Research bars (all required)

### 1. Intervention catalog (breadth)
Design and run **at least five independent intervention classes** aimed at the plateau, each with a clear hypothesis. Examples (implement what fits; invent others if better):

1. **Zero-fit / peel** — strip or LS-fit more zeros; vary \(N\); residual-minus-modes  
2. **Whitening** — detrend degree, smooth, log-window change, edge taper  
3. **Measure / residual definition** — \(\psi\) vs \(\theta\), different normalizations, weight \(w(u)\)  
4. **Basis / projector** — change \(d\), shifted Legendre vs other polys, high-pass of residual  
5. **Scale** — larger \(x\) if RAM allows, else multi-\(T\) denser near \(T_{\max}\) on existing \(5\times10^{10}\)  
6. **Controls** — Beurling/defective vs ordinary under same intervention (optional 6th)

Each class must produce **durable multi-\(T\) numeric \(R_d\)** (not a single \(T\) anecdote) under `results/open_plateau/` (or clear subdirs).

### 2. Depth (not a 20-minute batch)
For **at least three** of the classes, the experiment must be **expensive enough that a 86-core box spends on the order of hours**, not minutes, e.g.:

- MC-scale randomizations of residual construction / defects with **total trials ≥ 50M**, **or**
- Multi-\(T\) residual rebuilds on the **full \(5\times10^{10}\)** table across **≥20 \(T\)** × **≥5** construction variants with resume, **or**
- Beurling / system battery **≥500** systems × multi-\(T\), **or**
- Equivalent documented heavy compute (show trial counts / \(T\)×variant product in summaries)

Cheap classes can be small; **three deep axes** must show large `elapsed_s` and huge row/trial products in JSON summaries. If a “deep” axis finishes in &lt;30 minutes, **enlarge that axis** (more \(T\), variants, or trials)—not invent empty loops.

### 3. Synthesis (the actual finish line)
Write `docs/RESEARCH_PLATEAU.md` (and STATUS pointer) that:

- States the plateau fact with numbers from shipped campaigns  
- Summarizes **each** intervention class: hypothesis, what was run, what \(R_d(T)\) did  
- Gives a **judgment**: which interventions moved the needle, which failed, best current residual recipe  
- Lists what would be required next for full Theorem A (zero-sum control, etc.)  
- Keeps the RH non-claim  

### 4. Engineering
- Resume stamps under `results/open_plateau/` so sleep/crash doesn’t lose work  
- Multi-core for independent units (never multi-hour single-core)  
- Optional V100/CuPy when free; CPU path must work  
- Full pytest green; no pegged workers at end  

## Stop conditions
**Stop when research bars 1–4 are honestly met.**  
**Do not stop** after only re-running overnight floors (peel 2k / Beurling 100 / MC 200k×8) without new intervention classes and synthesis.

## Non-goals
RH proof; prize posts; token-burning sleep loops; overwriting durable `results/beurling_battery/` with smoke.
