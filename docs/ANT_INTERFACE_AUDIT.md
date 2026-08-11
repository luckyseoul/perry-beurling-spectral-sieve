# Rank 3 — Bounded ANT interface audit (Full A freeze)

**Date:** 2026-08-11  
**Code:** `pbss.ant_audit.ant_interface_audit`  
**Related:** [`THEOREM_A_PACKAGE.md`](THEOREM_A_PACKAGE.md)

## Goal

Verify that every Full-A external input (ANT-1…3) has a **named interface** to PBSS
objects (residual, window, \(R_d\), weight, mode sum), with no unlabeled gaps — then
**stop Full A packaging / re-proof work**.

## Entry

```bash
PYTHONPATH=src python3 -c "from pbss.ant_audit import ant_interface_audit; import json; print(json.dumps(ant_interface_audit(), indent=2, default=str)[:2500])"
```

## Disposition

| ID | Status | Action |
|----|--------|--------|
| ANT-3 EF identification | **matched** | stop re-proof |
| ANT-1 infinite tail under RH | **matched** | stop re-proof; scaffold is diagnostic only |
| ANT-2 arithmetic remainder | **matched** | stop re-proof |
| ANT-4 weight transfer | **not_required** | optional |
| M7 perturbation | **matched** | proved in-repo |

`freeze_full_a_packaging=true` when unlabeled_count=0 and
`full_arithmetic_A=closed_conditional`.

## Non-claims

- Does **not** re-derive Davenport/Ingham/Titchmarsh.  
- Does **not** prove RH.  
- Does **not** reopen Full A as open packaging.
