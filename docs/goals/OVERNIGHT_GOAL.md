# Overnight marathon goal

**Activate with something like:**  
`/goal do OVERNIGHT_GOAL.md multi-core GPU resume RH open`

## One sentence
Resume-capable multi-phase PBSS marathon on real compute: dense arithmetic zero-peel on existing 5e10 primes, large Beurling battery, large MC stress, GPU residual multi-T ablations; update STATUS; RH stays open. Finish when the **work floors** below are met with durable artifacts—not when a clock hits a number.

## Do NOT pad wall time
- No minimum hours. No “keep running until 6h.”
- If a phase finishes early with floors met, **move on or stop**—do not invent busywork to burn tokens.
- Wall time is a side effect of large MC / systems / grids, not a success metric.

## Work floors (done only when all are true)
1. **Arithmetic zero-peel** on \(x_{\max}=5\times10^{10}\) (existing prime checkpoint): multi-\((T,N,d,\mathrm{detrend})\) with **≥2000 rows** of numeric \(R_d\); JSON/TXT + plot under marathon results.
2. **Beurling battery:** **≥100** systems, multi-\(T\), scorecard + plot. Write under marathon out-dir (do **not** clobber durable `results/beurling_battery/` with tiny smokes).
3. **MC stress:** **≥200 000 trials per \(T\)**, **≥8** distinct \(T\) values, multi-\(d\); mean/std \(R_d\); multi-core (compute-budget / full workers, not one core).
4. **GPU residual multi-\(T\)** ablations on \(x\le5\times10^{10}\) when V100/CuPy available (CPU fallback OK); numeric \(R_d\) artifacts.
5. **Resume:** phase stamps so a restart does not redo completed phases.
6. **STATUS** marathon section + **explicit RH non-claim**; full `pytest` green; no pegged workers left after exit.

## Non-goals
- Unconditional RH / full A–B proof claims.
- Prize posts.
- Forced sieve to \(10^{12}\) if RAM/disk say no (document stop; enlarge MC/Beurling instead if needed).
- Smoke runs writing over durable scorecards.

## Implementation notes
- Prefer ProcessPool fan-out + CuPy residual/projection when free.
- Evidence under session scratch; durable JSON/TXT/plots under `results/overnight_marathon/` (or clear per-phase subdirs).
- Enlarge grids only to hit **work floors**, not to hit a clock.
