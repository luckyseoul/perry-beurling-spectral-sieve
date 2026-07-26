# Status: Perry–Beurling Spectral Sieve

**Date:** 2026-07-26  
**Repo:** `luckyseoul/perry-beurling-spectral-sieve` (private)

## Explicit non-claim

**This repository does not contain an unconditional proof of the Riemann Hypothesis.**

What *is* claimed: model lemmas M1–M5, multi-\(T\) numerics, truncated explicit-formula
residuals, arithmetic zero-peel diagnostics, Beurling battery scorecards, extended-\(x\)
scans, and MC instrument stress. Full Theorems A/B and RH remain open.

Full writeup: [`THEOREMS_AB.md`](THEOREMS_AB.md) · Proofs: [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md)

---

## Marathon legs 1–5 (2026-07-26)

### 1. Arithmetic zero-peel

**Entry:** `experiments/run_arithmetic_zero_peel.py`  
**Artifacts:** `results/arith_zero_peel/` (JSON, TXT, plot)  
**What:** \(q_T=(\theta-x)/\sqrt{x}\) with first \(N\) CL modes stripped (optional LS scale α), multi-\((T,N,d,\mathrm{detrend})\).  
**Reading:** Peeling model modes does **not** collapse arithmetic \(R_d\) to A₀ levels through \(10^{10}\). Diagnostic only.

### 2. Larger-\(x\) stack

**Entry:** `experiments/run_extend_x_scan.py` (+ GPU residual helpers in `pbss.gpu_residual`)  
**Checkpoint:** `results/prime_checkpoints/primes_le_50000000000.*`  
**Artifacts:** `results/extend_x_scan/` (JSON, TXT, plot, csum)  

| Quantity | Value |
|----------|------:|
| Chosen \(x_{\max}\) | \(5\times10^{10}\) |
| \(n_{\mathrm{primes}}\) | 2 119 654 578 |
| Sieve | parallel segment extend from \(10^{10}\), **86 workers**, ~131 s |
| \(10^{11}\) | **Rejected:** est. ~42 GiB primes > ~55% of 60 GiB RAM |
| \(10^{12}\) | Infeasible (disk/RAM) |
| Residual | GPU searchsorted+project when prefix fits V100; mmap for full table |
| Focus \(R_4\) deg1 at \(T\approx24.6\) | ~0.146 (soft plateau continues) |

### 3. Beurling battery

**Entry:** `experiments/run_beurling_battery.py`  
**Artifacts:** `results/beurling_battery/`  
**Systems:** `ordinary_primes` (rh_like), `gapped_gap3`, `thinned_every3` (defective).  
**At \(T=18\):** ordinary \(R_4\sim0.18\); defective \(\sim0.98\)–\(0.99\) — clear separation.

### 4. MC / ablation stress

**Entry:** `experiments/run_mc_stress.py`  
**Artifacts:** `results/mc_stress/`  
**Scale:** **50 000** trials/T × 4 windows = 200 000 trials; degrees \(2,4,6,8\); 86 workers.  
**Result:** mean \(R_d\sim0.79\) flat in \(T\); stamp `MC_STRESS_COMPLETE`.

### 5. Docs / non-claim

This STATUS + theorem pointers. Paper: `docs/paper/`. **No RH proof claim.**

---

## Proved in-repo (model diagnostic)

| Lemma | Statement | Code / test |
|-------|-----------|-------------|
| **M1–M4** | pure mode / defect / CL decay / blocks vanishing | `lemmas` · `test_lemmas` |
| **M5** | Finite CL sum \(R_d=O_d(T^{-2})\) | `bound_R_d_finite_mode_sum` · `test_M5_*` |

| Result | Status |
|--------|--------|
| **A₀ / finite-mode A₀** | **Proved** (M3/M5) |
| **A** arithmetic under RH | **Open** |
| **B₀** | **Proved** (M2+M4) |
| **B** | **Open** |

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_arithmetic_zero_peel.py
PYTHONPATH=src python3 experiments/run_beurling_battery.py
PYTHONPATH=src python3 experiments/run_mc_stress.py --mc-per-t 50000 --workers 86
PYTHONPATH=src python3 experiments/run_extend_x_scan.py --workers 86
PYTHONPATH=src python3 experiments/run_explicit_formula_peel.py
```
