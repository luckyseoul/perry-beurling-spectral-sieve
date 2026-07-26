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

## Multi-\(T\) campaign (shipped path)

Command:
```bash
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86 --scratch …
```

Representative numbers (`results/multi_T_scan.json`, d=4, n=8192, ε=0.5):

| \(T\) | \(R_4\) critical-line | \(R_4\) persistent defect | gap |
|------:|----------------------:|--------------------------:|----:|
| 3 | 2.28e-02 | 0.250 | 0.227 |
| 20 | 1.00e-03 | 0.250 | 0.249 |
| 80 | 6.13e-05 | 0.250 | 0.250 |

- Critical-line: **decays** \(\sim T^{-2}\) (0.023 → 6e-5).  
- Persistent defect: **flat** at \(\varepsilon^2=0.25\) (Lemma M2).  
- Gap at large \(T\): **≥ 0.249**.

Off-critical envelope modes are implemented (`probe_off_critical_mode`) for further attack; they are **not** yet a complete reduction of arithmetic off-line zeros.

---

## Session progress vs open

**Moved this session (A/B push):**

1. Precise A/B definitions and status split.  
2. Proved M1–M4 with written proofs + unit tests on shipped projection.  
3. Multi-core multi-\(T\) separation: A₀ decay vs B₀ flat defect.  
4. Off-critical probe for next-stage B work.

**Still open:**

1. Full Theorem A for prime residual under RH.  
2. Full Theorem B (obstruction ⇒ RH).  
3. Sharp rate \(O(T^{-2(d+1)})\).  
4. Large-\(T\) arithmetic residual without staircase artifact.  
5. Legacy P≈3.92 normalization.

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
PYTHONPATH=src python3 experiments/run_diagnostic.py
```
