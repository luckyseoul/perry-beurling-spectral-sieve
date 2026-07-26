# Kill 0.21: enrich \(m\) only (H_theta_sqrt)

**Date:** 2026-07-26  
**Code:** `pbss.ef_identify.build_m_columns` / `identify_ef(..., m_enrich=...)` / `multi_N_enrich_scan`  
**Campaign:** `results/kill021_enrich_m/`  
**Scope:** residual fixed = H_theta_sqrt \((\theta-x)/\sqrt{x}\) deg1; only enrich the mode model \(m\).

## Metric

\[
\frac{E_d(r)}{\|q\|^2}
=\frac{\|P_d(q-m_{\mathrm{eff}})\|_2^2}{\|q\|_2^2}
\]

Baseline (zeros-only \(m\)): \(\approx 0.21\), flat in \(N\) (prior EF attack).

## Enrichments of \(m\) (joint LS with zero sum)

| `m_enrich` | Extra columns in span of \(m\) |
|------------|--------------------------------|
| `zeros` | \(m_T^{(N)}\) only |
| `zeros_poly1` | \(+\,1,\;(u-\tfrac12)\) |
| `zeros_highleg` | \(+\,\varphi_2,\ldots,\varphi_d\) |
| `zeros_Vd` | \(+\,\varphi_0,\ldots,\varphi_d\) (full \(V_d\)) |
| `zeros_smooth` | \(+\,e^{-uT/2},\;e^{-uT},\;e^{-2uT}\) |
| `zeros_endpoint` | \(+\,u^2(1-u)^2,\;u^2(1-u)^2(u-\tfrac12)\) |

## Multi-\(N\) campaign (primes \(\le 10^{10}\), \(T\in\{14,16,18,20\}\), \(N\in\{5,10,20,40\}\), \(d=4\))

### Mean \(E_d(r)/\|q\|^2\) by enrichment

| enrichment | mean Ed | vs zeros |
|------------|--------:|---------:|
| **zeros_Vd** | **~0** | **killed** |
| **zeros_highleg** | **~0.0002** | **killed (~0.1%)** |
| **zeros_smooth** | **0.052** | **~25% of baseline** |
| zeros_endpoint | 0.142 | ~69% |
| zeros_poly1 | 0.207 | ~100% |
| zeros | 0.207 | 100% |

### Multi-\(N\) (zeros vs smooth vs highleg)

| \(N\) | zeros | zeros_smooth | zeros_highleg |
|------:|------:|-------------:|--------------:|
| 5 | 0.202 | 0.051 | ~0 |
| 10 | 0.206 | 0.051 | ~0 |
| 20 | 0.209 | 0.052 | ~0 |
| 40 | 0.212 | 0.052 | ~0 |

Zeros-only stays flat ~0.21. Smooth enrichment **holds ~0.052 independent of \(N\)** (not zero-summable further). Highleg / Vd zero Ed by spanning the low-degree polynomial space.

## Outcome: **win** (0.21 killed for some enrichments)

1. **Primary kill:** putting \(\varphi_2,\ldots,\varphi_d\) (or full \(V_d\)) into \(m\) drives Ed→0. The 0.21 **is** low-degree Legendre mass of \(q\) after deg1 detrend (degrees ≥2).
2. **Non-tautological partial kill:** classical-ish **exp-decay secondary columns** cut Ed to ~0.05 (~4× improvement) without free \(V_d\).
3. **Still flat in \(N\):** neither zeros-only nor smooth Ed falls as more zeros are added — zeros still do not eat the remaining \(V_d\)-shaped mass.

## What this means for Full A

- Peeling more zeros cannot kill \(R_d\) until \(m\) includes the **secondary / low-degree structure** of \((\theta-x)/\sqrt{x}\).
- Next constructive step: derive EF secondary main terms that reproduce the \(\varphi_{k\ge 2}\) content (or prove that content is \(O(T^{-2})\) under RH after a correct weight).
- `zeros_Vd` is a diagnostic oracle, not an EF theorem; `zeros_smooth` is the interesting quantitative middle ground.

## Reproduce

```bash
PYTHONPATH=src python3 -m pytest tests/test_enrich_m_kill021.py -v
PYTHONPATH=src python3 experiments/run_kill021_enrich_m.py \
  --out-dir results/kill021_enrich_m
```
