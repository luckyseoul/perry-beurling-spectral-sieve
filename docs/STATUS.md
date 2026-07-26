# Status: Perry–Beurling Spectral Sieve reconstruction

**Date:** 2026-07-26  
**Repo:** `luckyseoul/perry-beurling-spectral-sieve` (private)  
**Author of original research notes:** Nicholas Perry  

## Explicit non-claim

**This repository does not prove or disprove the Riemann Hypothesis.**  
The framework is a **conditional / diagnostic classifier** for RH-like spectral structure on density perturbations of prime (or Beurling generalized prime) systems.

---

## Theorems A / B (conditional, as archived)

### Theorem A (conditional on RH) — diagnostic statement

*Assuming the Riemann Hypothesis*, the degree-\(d\) projection strength of the
normalized density residual \(q_T\) on a logarithmic window of length \(T\)
decays to 0 as \(T\to\infty\):

\[
R_d(q_T)=\frac{\|P_d q_T\|_{L^2}^2}{\|q_T\|_{L^2}^2}\to 0.
\]

The archive notes the rate heuristic \(R_d(q_T)=O(T^{-2(d+1)})\).

**Status:** Conditional on RH; not verified as a full theorem with published
proof in this repo. Used as a design principle for the scaled statistic
\(S_d=T^{2(d+1)}R_d\).

### Theorem B (obstruction / converse direction)

*If* \(R_d(q_T)\) decays sufficiently rapidly as \(T\to\infty\), then there are
no zeros off the critical line (conditional equivalence / obstruction sketch).

**Status:** The converse is essentially as hard as RH itself (README limitation).
**Not** an unconditional proof. Recorded here only as the intended logical
role of the diagnostic.

---

## What this session shipped (needle move)

| Item | Location |
|------|----------|
| Orthonormal shifted Legendre basis on \([0,1]\) | `src/pbss/basis.py` |
| Projection API: \(c_k\), \(E_d\), \(R_d\), \(S_d=P\) | `src/pbss/projection.py` |
| Synthetic + prime residual probes | `src/pbss/probes.py` |
| Unit tests driving real projection path | `tests/test_projection.py` |
| End-to-end experiment (RH-like vs defective) | `experiments/run_diagnostic.py` |
| Numerical results | `results/diagnostic_run.json` |

### Formulas used (reconstruction)

- Basis: \(\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1)\) on \(u\in[0,1]\).
- Coefficients: \(c_k=\langle q,\varphi_k\rangle\) via trapezoid quadrature.
- Energy: \(E_d=\sum_{k=0}^d|c_k|^2\).
- Energy ratio: \(R_d=E_d/\|q\|^2\in[0,1]\).
- Working projection strength: \(P:=S_d=T^{2(d+1)}R_d\).

### Legacy numbers

Archive README quotes **P(q)≈3.92** for zeta and threshold **≈29.5**. Those
came from lost high-precision scripts with incompletely documented normalization.
**This reconstruction does not hard-code or force agreement with 3.92.**  
It reports honest \(R_d\) and \(S_d\) from the shipped code. Matching or
refuting 3.92 is future work once the original normalization is recovered.

### Representative run (d=4, n=4096; see `results/diagnostic_run.json`)

| Probe | \(R_4\) (energy ratio) | Role |
|-------|------------------------:|------|
| High-frequency sinusoid | \(\approx 8.8\times 10^{-4}\) | RH-like synthetic |
| Critical-line mode \(\sin(tTu)\) | \(\approx 1.0\times 10^{-3}\) | RH-consistent form |
| Demeaned prime residual (\(x\le 10^5\)) | \(\approx 0.22\) | real primes; residual still has mid-band mass |
| Defective (HF + degree-1 weight 2.5) | \(\approx 0.93\) | non-RH control |

Classifier: defective \(R_d\) exceeds both RH-like synthetics by \(\approx 0.92\).

---

## Limitations (from README, still in force)

1. Classifier / diagnostic — not a decisive RH proof.  
2. Finite windows cannot exclude zeros at extremely high height.  
3. Not a practical local primality sieve.  
4. Cost grows with degree \(d\) and window size.  
5. Converse (low energy ⇒ RH) is as hard as RH.

---

## Remaining open (next needles)

1. Recover or re-derive the exact historical normalization that produced P≈3.92.  
2. Diamond-system ground-truth battery (Beurling systems with known RH / non-RH).  
3. PSWF basis swap (mentioned in archive notes).  
4. Large-\(T\) prime residual scans and Monte Carlo defect campaigns (MC=4000 lore).  
5. Rigorous statement + proof of Theorems A/B under stated hypotheses (analytic number theory, not code).

---

## How to run

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_diagnostic.py --scratch /path/to/scratch
```
