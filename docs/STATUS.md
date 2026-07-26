# Status: Perry–Beurling Spectral Sieve

**Date:** 2026-07-26  
**Repo:** `luckyseoul/perry-beurling-spectral-sieve` (private)

## Explicit non-claim

**This repository does not contain an unconditional proof of the Riemann Hypothesis.**

What *is* claimed: precise conditional Theorems A/B, **proved model lemmas M1–M5**
about the diagnostic (including finite-mode A₀), multi-\(T\) numerics supporting
A₀/B₀, truncated explicit-formula residuals, and multi-\((T,N)\) peel scans.

Full writeup: [`THEOREMS_AB.md`](THEOREMS_AB.md) · Proofs: [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md)

---

## Proved in-repo (model diagnostic)

| Lemma | Statement | Code / test |
|-------|-----------|-------------|
| **M1** | \(R_d(\varphi_m)=1_{m\le d}\) | `lemmas.continuous_R_d_pure_mode` · `test_M1_*` |
| **M2** | Orthogonal defect \(R_d=\varepsilon^2\) | `continuous_R_d_orthogonal_defect` · `test_M2_*` |
| **M3** | Critical-line mode \(R_d(\sin(tTu))=O(T^{-2})\) | `test_M3_*` |
| **M4** | Fixed \(\varepsilon>0\) ⇒ \(R_d\not\to0\) | `test_M4_*` |
| **M5** | Finite CL superposition \(R_d=O_d(T^{-2})\) (finite-mode A₀) | `bound_R_d_finite_mode_sum` · `test_M5_*` |

## Theorems A/B (status split)

| Result | Status |
|--------|--------|
| **A₀** (critical-line pure mode \(R_d\to0\)) | **Proved** (M3) |
| **Finite-mode A₀** (finite CL sum \(R_d=O(T^{-2})\)) | **Proved** (M5) |
| **A** (arithmetic residual under RH) | **Conditional / open** |
| **B₀** (persistent defect blocks \(R_d\to0\)) | **Proved** (M2+M4) |
| **B** (fast residual decay ⇒ RH) | **Open** (as hard as RH) |

---

## Explicit-formula residual + peel scan

**Builder:** `pbss.probes.explicit_formula_residual` — truncated sum of first \(N\)
critical-line modes (\(a_n=2/|\rho_n|\), offline ordinates in `pbss.zeros`).  
**Peel:** `peel_residual` strips the first \(N\) modes from a full truncated sum.  
**Experiment:** `experiments/run_explicit_formula_peel.py`  
**Artifacts:** `results/explicit_formula_peel/` (JSON, TXT, PNG).

| \(T\) | \(N\) | \(R_d\) include (d=4) | note |
|------:|------:|----------------------:|------|
| 8 | 1 | \(6.2\times10^{-3}\) | decays with \(T\) |
| 48 | 1 | \(1.7\times10^{-4}\) | |
| 8 | 20 | \(1.3\times10^{-2}\) | more modes, still \(O(T^{-2})\) |
| 48 | 20 | \(3.3\times10^{-4}\) | |

**Reading:** Truncated EF residuals behave like finite-mode A₀ (M5). This does
**not** make the arithmetic residual at \(x\le10^{10}\) collapse, and is **not**
an RH proof.

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
**MC:** **20,000** spectral-defect trials **per \(T\)** × 14 windows = **280,000** trials.  
**Arithmetic focus** (d=4, deg1): soft plateau \(R_4\sim0.15\)–\(0.19\) through \(10^{10}\).

**Reading:** Arithmetic \(R_d\) soft-plateaus; still no A0-style collapse. Full Theorem A open. **Not a proof of RH.**

Plots: `results/grand_campaign/grand_Rd_vs_T.png`, `grand_arith_focus_linear.png`.

---

## Progress vs open

**Done:**

1. A0/B0 lemmas M1–M4; **finite-mode A₀ as M5**.  
2. Explicit-formula truncated residual + peel + multi-\((T,N)\) scan.  
3. Segmented sieve + prime checkpoint to \(10^{10}\).  
4. Grand campaign: multi-\(T\) arithmetic ablations, controls, MC, resume, plots.  
5. Optional CuPy projection backend with NumPy fallback.  
6. Status paper: `docs/paper/`.

**Still open:**

1. Full Theorem A (arithmetic \(R_d\to 0\) under RH) — not seen at \(x\le 10^{10}\).  
2. Full Theorem B.  
3. Zero-peeling the **arithmetic** residual (match primes to EF modes).  
4. Legacy P≈3.92 normalization.

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 4
PYTHONPATH=src python3 experiments/run_explicit_formula_peel.py
PYTHONPATH=src python3 experiments/run_diagnostic.py
```
