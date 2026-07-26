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

### Arithmetic residual (`experiments/run_arithmetic_multi_T.py`)

Builder: `arithmetic_residual` — \((\theta(x)-x)/\sqrt{x}\), detrend=deg1.

| \(T\) | \(x_{\max}\) | \(R_4\) arithmetic | vs defect 0.25 |
|------:|-------------:|-------------------:|:---------------|
| 8 | 3e3 | 0.108 | below |
| 12 | 1.6e5 | 0.189 | below |
| 16 | 8.9e6 | 0.193 | below |

**Reading:** arithmetic \(R_d\) **plateaus** (~0.19), does **not** show A0-style decay at \(T\le 16\). Still separated from defect floor. Full details: `results/arithmetic_multi_T.json`, `docs/THEOREMS_AB.md`.

---

## Progress vs open

**Done:**

1. A0/B0 lemmas M1–M4 proved + tested.  
2. Model multi-\(T\) separation.  
3. **Arithmetic residual multi-\(T\)** measured to \(x_{\max}\approx 9\cdot 10^6\) (honest plateau).  
4. Detrend/smooth options on residual builder.

**Still open:**

1. Full Theorem A (arithmetic \(R_d\to 0\) under RH) — **not** seen yet at accessible \(T\).  
2. Full Theorem B.  
3. Larger \(T\) / better residual whitening / explicit formula.  
4. Legacy P≈3.92 normalization.

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
PYTHONPATH=src python3 experiments/run_diagnostic.py
```
