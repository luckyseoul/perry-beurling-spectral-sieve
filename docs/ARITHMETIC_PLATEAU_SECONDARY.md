# Rank 2 — Arithmetic plateau via secondary EF-style terms

**Date:** 2026-08-11  
**Code:** `pbss.plateau_secondary` · **Tests:** `tests/test_ranks_2_to_5.py`  
**Banner:** **Not a proof of RH. Not a claim that Full A is unconditional.**

## Question

Why does arithmetic \(R_d\) / low-degree remainder energy stay large (soft plateau)
while finite critical-line mode sums decay (M5)?

## Method

Fix residual \(H_{\theta,\sqrt{}}=(\theta-x)/\sqrt{x}\) (deg1 detrend). Enrich only the
model \(m\) using the shipped kill-0.21 columns (`zeros`, `zeros_smooth`, `zeros_Vd`)
via `ef_identify.multi_N_enrich_scan`.

### Predeclared predictions

| ID | Prediction |
|----|------------|
| P1 | `zeros_smooth` mean \(E_d(r)/\|q\|^2\) **strictly below** `zeros` |
| P2 | `zeros_smooth` still has **positive** Ed (not a free \(V_d\) oracle) |
| P3 | `zeros_Vd` drives Ed **near 0** |
| P4 | `zeros`-only Ed is **not** strongly decreasing in \(N\) |

## Checkable entry

```bash
PYTHONPATH=src python3 -c "from pbss.plateau_secondary import plateau_secondary_report; import json; print(json.dumps(plateau_secondary_report(), indent=2, default=str)[:2000])"
```

Default scale: \(T=14\), \(x_{\max}=2\times 10^6\) (large enough for P1).

## Mechanism (if predictions pass)

- Missing **low zeros** are not the story: Ed flat in \(N\) for zeros-only.  
- Remaining mass is largely **\(V_d\)-shaped** (oracle kill by `zeros_Vd`).  
- Classical-ish **exp-decay secondary** columns cut Ed without free Legendre span.  

## Stop conditions

- If P1–P3 fail under moderate \(T/x_{\max}\) changes → do not claim a robust mechanism.  
- Do **not** replace analysis by larger-\(x\) marathons without a transfer theorem.
