# Status: Perry–Beurling Spectral Sieve

**Date:** 2026-07-29 (tight stress) / 2026-07-26 (marathon)  
**Repo:** `luckyseoul/perry-beurling-spectral-sieve` (private)

## Explicit non-claim

**This repository does not contain an unconditional proof of the Riemann Hypothesis.**

What *is* claimed: model lemmas M1–M5, multi-\(T\) numerics, truncated explicit-formula
residuals, arithmetic zero-peel diagnostics, Beurling battery scorecards, extended-\(x\)
scans, and MC instrument stress. Full Theorems A/B and RH remain open.

Full writeup: [`THEOREMS_AB.md`](THEOREMS_AB.md) · Proofs: [`PROOFS_LEMMAS.md`](PROOFS_LEMMAS.md)


---

## Tight stress session (2026-07-29, Grok Build)

**Entry:** `experiments/run_tight_stress_20260729.py`  
**Artifacts:** `results/tight_stress_20260729/` (`STRESS_REPORT.md`, `formal_diagnostics.json`, MC/Beurling/off-critical JSON)  
**Banner:** **NOT AN UNCONDITIONAL PROOF OF RH.**

### Formalization

Locked working definitions: \(R_d\), \(S_d=T^{2(d+1)}R_d\), \(P(q):=S_d\) (not legacy 3.92). Theorems A₀/M3, finite A₀/M5, B₀/M2–M4 **proved** in-repo; full A (arithmetic under RH) and B (converse) **open**.

### MC ablations (192 000 trials)

| Ablation | mid-T mean \(R_4\) |
|----------|-------------------:|
| baseline | 0.792 |
| heavy_defect | 0.970 |
| light_defect | 0.285 |
| high_freq | 0.792 |
| low_freq | 0.794 |
| high_deg_defect | 0.496 |

- Flat instrument: baseline mean \(R_4\approx0.794\), std across \(T\) \(\approx0.002\).
- **Failure-mode note:** `high_deg_defect` lowers \(R_d\) when defect mass sits **outside** \(V_d\) (degrees \(>d\)) — classifier must fix defect support relative to projection degree.

### Off-critical \(\sigma\) sweep

Model modes at \(\sigma=0.9\) vs \(1/2\): \(R_4\) ratio grows with \(T\) (~3.9 at \(T=8\) → ~15.9 at \(T=32\)). Directional diagnostic only — **not** a zero-free-region theorem.

### Expanded Beurling (35 systems, \(x_{\max}=10^6\))

- Ordinary \(R_4\approx0.189\); defective min/median/max \(\approx0.993/0.995/0.996\).
- **Perfect separation** (no failure systems; no thin-margin cases at this scale).

### Rate check

M3 \(O(T^{-2})\) holds for \(d\ge1\). Stronger archive heuristic \(O(T^{-2(d+1)})\) **not** supported as a proved/empirical sharp bound for \(S_d\).

### Still open

Arithmetic soft plateau, full A/B, legacy \(P\approx3.92\), infinite-zero control. RH open.


---

## Overnight marathon (`OVERNIGHT_GOAL.md`)

**Entry:** `experiments/run_overnight_marathon.py`  
**Artifacts:** `results/overnight_marathon/` (phase stamps + peel / beurling / mc / residual)  
**Work floors (met):** peel **2048** rows on \(x_{\max}=5\times10^{10}\); Beurling **100** systems; MC **200 000**/T × **8** \(T\); residual multi-\(T\) **96** rows.  
**Resume:** `PHASE_{A,B,C,D}_*_COMPLETE` + `MARATHON_COMPLETE`.  
**Compute:** multi-core ProcessPool for MC (86w) and Beurling waves; peel multi-core by \(T\); residual sequential mmap (pool OOM avoided).  
**RH:** open — **not a proof.**

| Phase | Floor | Result |
|-------|-------|--------|
| A peel | ≥2000 rows | 2048 |
| B Beurling | ≥100 systems | 100 (600 rows) |
| C MC | ≥200k/T × ≥8 T | 200k × 8 = 1.6M trials |
| D residual | multi-T large \(x\) | 96 rows on \(5\times10^{10}\) |

Smoke tests must use separate `--out-dir` under `results/overnight_marathon/smoke_*`.

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
**Artifacts:** `results/beurling_battery/` (durable full run; smoke must use `--out-dir` elsewhere)  
**Scale:** \(x_{\max}=10^{8}\), \(T\in\{8,10,12,14,16,18\}\), degree 4, 3 systems, elapsed ~7.3 s.  
**Systems:** `ordinary_primes` (rh_like), `gapped_gap3`, `thinned_every3` (defective).  
**Scorecard at \(T=18\):** ordinary \(R_4\approx0.181\); gapped \(\approx0.983\); thinned \(\approx0.990\) — clear separation.

### 4. MC / ablation stress

**Entry:** `experiments/run_mc_stress.py`  
**Artifacts:** `results/mc_stress/`  
**Scale:** **50 000** trials/T × 4 windows = 200 000 trials; degrees \(2,4,6,8\); 86 workers.  
**Result:** mean \(R_d\sim0.79\) flat in \(T\); stamp `MC_STRESS_COMPLETE`.

### 5. Docs / non-claim

This STATUS + theorem pointers. Paper: `docs/paper/`. **No RH proof claim.**

---

## Open-plateau campaign (`OPEN_GOAL.md`)

**Entry:** `experiments/run_open_plateau.py`  
**Artifacts:** `results/open_plateau/` (per-class JSON + `PHASE_*_COMPLETE` resume stamps)  
**Synthesis:** [`RESEARCH_PLATEAU.md`](RESEARCH_PLATEAU.md) — plateau numbers, per-class judgment, RH non-claim.  
**Classes (≥5):** peel, whiten, measure, basis, scale, mc_rand, beurling.  
**Deep axes (≥3):** MC **50M** trials (~1.7 h, 86w); scale **24 T × 6** variants on \(5\times10^{10}\); Beurling **≥500** systems multi-\(T\) (enlarged battery).  
**Headline judgment:** Peeling zeros does **not** collapse arithmetic \(R_d\); deg1+optional taper is best residual recipe; dense multi-\(T\) shows **no** A₀ decay through \(x\sim5\times10^{10}\). **RH remains open.**

---

## Theorem-A scaffolding (items 1–2)

**Writeup:** [`THEOREM_A_SCAFFOLD.md`](THEOREM_A_SCAFFOLD.md)  
**Code:** `pbss.weights` (admissible \(W_\alpha\), endpoint estimator), `pbss.remainder` (peel-via-remainder, M5 + tail majorant)  
**Entry:** `experiments/run_theorem_a_scaffold.py`  
**Artifacts:** `results/theorem_a_scaffold/` (`PHASE_{WEIGHT,REMAINDER,ARITHMETIC}_COMPLETE`, `SCAFFOLD_COMPLETE`)  
**What it is:** formal weight class for endpoint control + truncated zero-sum remainder path with multi-\((T,N)\) numerics.  
**What it is not:** full Theorem A or RH. Tail majorants are scaffolding labels, not sharp ANT bounds.

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
PYTHONPATH=src python3 experiments/run_open_plateau.py --classes peel,whiten,measure,basis,scale,beurling,mc_rand --workers 86
PYTHONPATH=src python3 experiments/run_theorem_a_scaffold.py --out-dir results/theorem_a_scaffold --workers 32
```
