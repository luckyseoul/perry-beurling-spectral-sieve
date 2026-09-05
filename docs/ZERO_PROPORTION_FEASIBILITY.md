# Rank 4 — Zero-proportion (Anthropic-style) feasibility

**Date:** 2026-08-11  
**Code:** `pbss.zero_proportion_feasibility`  
**Upstream:** [`related/anthropic-riemann-zeta/`](related/anthropic-riemann-zeta/)

## Decision: **STOP** (proportion reproduction only)

No incremental inequality is formulated that maps Weil / BGST / Bombieri-style
quadratic-form rank methods onto PBSS \(R_d\) **as a zero-count / proportion engine**.
Reproducing the upstream 41.6%→67.2% proportion bound would be a **separate research
program**, not a continuation of residual diagnostics.

**Proportion-STOP ≠ anti-cancellation assimilate.** The Weil/BGST *anti-cancellation /
coercivity* content toward \(H^*\) is tracked separately in
[`HSTAR_ANTICANCELLATION.md`](HSTAR_ANTICANCELLATION.md) (`pbss.hstar_anticancellation`).
That track must **not** flip these Rank-4 proportion candidates to `ready_to_implement`.
The 2026-09-05 audit retracts its claimed Weil coercivity proxy: see the
[transfer obstruction](WEIL_TRANSFER_OBSTRUCTION.md) and
[exact off-critical phase theorem](OFF_CRITICAL_PHASE_OBSTRUCTION.md).

## Class comparison (summary)

| Axis | Anthropic-style | PBSS |
|------|-----------------|------|
| Object | fraction of zeros on \(\mathrm{Re}s=\tfrac12\) | \(R_d\) of density residual |
| Tools | Weil forms, BGST, Bombieri | Legendre projection, EF peel, Beurling |
| ⇒ RH? | no | no |

## Do not

- Rebuild the killed Spectral Sieve operators for this.  
- Claim proportion methods solve **B-RES**.  
- Launch large zero-count campaigns without a new theorem statement.

## Resume only if

A candidate inequality is written with hypotheses and a PBSS-facing conclusion,
marked `ready_to_implement` in `incremental_inequality_candidates()`.
