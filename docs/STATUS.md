# Status: Perry–Beurling Spectral Sieve

**Date:** 2026-07-26  
**Repo:** `luckyseoul/perry-beurling-spectral-sieve` (private)

## Explicit non-claim

**This repository does not contain an unconditional proof of the Riemann Hypothesis.**

What *is* claimed: precise conditional Theorems A/B, **proved model lemmas M1–M4**
about the diagnostic, and multi-\(T\) numerics supporting the model forms A₀/B₀.

Full writeup: [`THEOREMS_AB.md`](THEOREMS_AB.md) · Proofs: [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md)

---

## Proved in-repo (model diagnostic)

| Lemma | Statement | Code / test |
|-------|-----------|-------------|
| **M1** | \(R_d(\varphi_m)=1_{m\le d}\) | `lemmas.continuous_R_d_pure_mode` · `test_M1_*` |
| **M2** | Orthogonal defect \(R_d=\varepsilon^2\) | `continuous_R_d_orthogonal_defect` · `test_M2_*` |
| **M3** | Critical-line mode \(R_d(\sin(tTu))=O(T^{-2})\) | `test_M3_*` |
| **M4** | Fixed \(\varepsilon>0\) ⇒ \(R_d\not\to0\) | `test_M4_*` |

## Theorems A/B (status split)

| Result | Status |
|--------|--------|
| **A₀** (critical-line pure mode \(R_d\to0\)) | **Proved** (M3) |
| **A** (arithmetic residual under RH) | **Conditional / open** |
| **B₀** (persistent defect blocks \(R_d\to0\)) | **Proved** (M2+M4) |
| **B** (fast decay of prime residual ⇒ RH) | **Open** (as hard as RH) |

---

## Multi-\(T\) campaigns (shipped path)

### Model A0/B0 (`experiments/run_multi_T.py`)

| \(T\) | \(R_4\) critical-line | \(R_4\) persistent defect |
|------:|----------------------:|--------------------------:|
| 3 | 2.28e-02 | 0.250 |
| 80 | 6.13e-05 | 0.250 |

Critical-line **decays**; defect **flat** at \(\varepsilon^2\).

### Arithmetic residual — overnight campaign (`experiments/run_overnight_campaign.py`)

**Scale reached:** \(x_{\max}=10^9\) (50.8M primes, segmented sieve ~5s), 86 workers, clean exit.  
Builder: `arithmetic_residual` — \((\theta(x)-x)/\sqrt{x}\).

**Focus arm** (detrend=`deg1`, smooth=1, d=4):

| \(T\) | \(x_{\max}\) | \(R_4\) arith | notes |
|------:|-------------:|-------------:|:------|
| 8.0 | 3.0e3 | 0.108 | |
| 16.0 | 8.9e6 | 0.193 | peak region |
| 18.7 | 1.3e8 | 0.180 | |
| 20.7 | **1.0e9** | **0.165** | still ≫ A0 mode |

Controls at \(T=20.7\): critical-line \(R_4\approx 1.0\times10^{-3}\), defect \(=0.250\), off-critical \(\approx 2.9\times10^{-3}\).

**Ablations at \(T=20.7\)** (same residual family):

| detrend \ smooth | 1 | 5 | 15 |
|:-----------------|--:|--:|---:|
| none | 0.966 | 0.974 | 0.983 |
| deg0 | 0.200 | 0.249 | 0.344 |
| deg1 | **0.165** | 0.208 | 0.294 |

**Reading:** With linear detrend, arithmetic \(R_d\) **rises then soft-plateaus** near 0.16–0.19 through \(x=10^9\) — **no A0-style \(T^{-2}\) collapse**. Raw residual is almost fully low-degree (\(\sim 0.97\)). Defect floor remains 0.25.  
Plots: `results/overnight_Rd_vs_T.png`, `results/overnight_arith_deg1_linear.png`.  
JSON: `results/overnight_campaign.json`.

**Not a proof of RH or of full Theorem A.**

---

## Progress vs open

**Done:**

1. A0/B0 lemmas M1–M4 proved + tested.  
2. Model multi-\(T\) separation.  
3. Arithmetic multi-\(T\) to \(10^9\) with detrend/smooth ablations + controls + plots.  
4. Segmented sieve for large \(x_{\max}\).

**Still open:**

1. Full Theorem A (arithmetic \(R_d\to 0\) under RH) — **not** seen at \(x\le 10^9\).  
2. Full Theorem B.  
3. Explicit-formula residual / stronger whitening beyond deg1.  
4. Legacy P≈3.92 normalization.

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
PYTHONPATH=src python3 experiments/run_diagnostic.py
```
