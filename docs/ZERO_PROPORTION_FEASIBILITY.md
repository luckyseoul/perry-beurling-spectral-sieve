# Rank 4 — Zero-proportion (Anthropic-style) feasibility

**Date:** 2026-08-11  
**Code:** `pbss.zero_proportion_feasibility`  
**Upstream:** [`related/anthropic-riemann-zeta/`](related/anthropic-riemann-zeta/)

## Decision: **STOP**

No incremental inequality is formulated that maps Weil / BGST / Bombieri-style
quadratic-form rank methods onto PBSS \(R_d\). Reproducing the upstream
41.6%→67.2% proportion bound would be a **separate research program**, not a
continuation of residual diagnostics.

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
