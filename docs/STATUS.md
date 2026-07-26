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

### Grand campaign (`experiments/run_grand_campaign.py`)

**Scale:** \(x_{\max}=10^{10}\) (455,052,511 primes, segmented sieve + **on-disk checkpoint**).  
**MC:** **20,000** spectral-defect trials **per \(T\)** × 14 windows = **280,000** trials (≫ 2000/T).  
**Ablations:** \(d\in\{2,4,6,8\}\), detrend \(\in\{\mathrm{none},\mathrm{deg0},\mathrm{deg1}\}\), smooth \(\in\{1,5,15\}\).  
**Controls:** critical-line, off-critical \(\sigma=0.75\), persistent defect \(\varepsilon=0.5\).  
**Elapsed:** ~3635 s wall (sieve+arith+controls+MC); clean worker exit.  
**Artifacts:** `results/grand_campaign/` (state JSON, summary, plots); primes under `results/prime_checkpoints/` (gitignored).

**Arithmetic focus** (d=4, deg1, smooth=1):

| \(T\) | \(x_{\max}\) | \(R_4\) |
|------:|-------------:|--------:|
| 10.0 | 2.2e4 | 0.144 |
| 16.0 | 8.9e6 | 0.193 |
| 19.9 | 4.4e8 | 0.174 |
| 23.0 | **1.0e10** | **0.155** |

Controls at \(T=23\): CL \(R_4\approx 4.5\times10^{-4}\), defect \(=0.250\), off-critical \(\approx 6.8\times10^{-4}\).

**MC defect (d=4):** mean \(R_d\approx 0.79\) (std \(\approx 0.16\)) **flat in \(T\)** — controlled defects stay high as expected.

**Reading:** Arithmetic \(R_d\) soft-plateaus ~0.15–0.19 through **ten billion**; still no A0-style collapse. Full Theorem A open. **Not a proof of RH.**

Plots: `results/grand_campaign/grand_Rd_vs_T.png`, `grand_arith_focus_linear.png`.

---

## Progress vs open

**Done:**

1. A0/B0 lemmas M1–M4.  
2. Segmented sieve + prime checkpoint to \(10^{10}\).  
3. Grand campaign: multi-\(T\) arithmetic ablations, controls, MC≥2000/T (20k/T), resume, plots.  
4. Optional CuPy projection backend with NumPy fallback.

**Still open:**

1. Full Theorem A (arithmetic \(R_d\to 0\) under RH) — still not seen at \(x\le 10^{10}\).  
2. Full Theorem B.  
3. Explicit-formula residual / stronger whitening.  
4. Legacy P≈3.92 normalization.

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
PYTHONPATH=src python3 experiments/run_diagnostic.py
```
