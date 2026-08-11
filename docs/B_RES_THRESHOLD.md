# Rank 5 — B-RES as threshold / obstruction (not RH)

**Date:** 2026-08-11  
**Code:** `pbss.b_res_threshold` · **Package:** [`THEOREM_B_PACKAGE.md`](THEOREM_B_PACKAGE.md)

## Explicit non-claims

- **B-RES is not solved.**  
- **RH is not proved.**  
- Model calculations are not arithmetic \(\zeta\).

## Threshold hypothesis \(H^*\)

After all EF main terms, secondary terms, and admissible weights, an off-critical
zero contributes a residual piece \(q_{\mathrm{off}}\) with
\[
\liminf_{T\to\infty}
\frac{\|P_d q_{\mathrm{off}}\|_2}{\|q_T^{\mathrm{arith}}\|_2}
\ge \varepsilon(\sigma,t,d)>0
\]
(no total cancellation into \(V_d^\perp\)).

**B-RES** = “\(H^*\) holds for the true arithmetic residual of \(\zeta\).”

Under \(H^*\), rapid \(R_d\to 0\) forbids off-critical zeros (model Full B).

## Model evidence (checkable)

1. **Pure off-critical mode** \(\mathrm{e}^{T(\sigma-1/2)u}\sin(tTu)\): \(R_d\) stays large
   vs critical-line mode (`off_critical_rd_lower_model`).  
2. **Cancellation counterexample:** subtract \(P_d q\) by force → \(R_d\approx 0\) even
   though \(\sigma\neq\tfrac12\) (`cancelled_off_critical_rd`).  
   ⇒ Without \(H^*\), off-critical *origin* alone does not force nonvanishing \(R_d\).

```bash
PYTHONPATH=src python3 -c "from pbss.b_res_threshold import b_res_threshold_report; import json; print(json.dumps(b_res_threshold_report(), indent=2, default=str)[:2000])"
```

## What to do next (only if pursuing B)

Prove a **weaker unconditional lemma** toward \(H^*\) (e.g. limited cancellation under
explicit secondary control) — do **not** attack full RH head-on.
