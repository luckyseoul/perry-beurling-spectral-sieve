# Perry–Beurling Spectral Sieve

Private research archive and **runnable reconstruction** of the **Perry–Beurling Spectral Sieve** (spectral diagnostic / P(q) framework) for testing consistency with the Riemann Hypothesis on Beurling generalized prime systems.

**Author:** Nicholas Perry  
**Status:** Independent research. Reconstruction of core projection diagnostic (2026-07).  
**Not a proof of RH** — see [`docs/STATUS.md`](docs/STATUS.md).

## Overview

A spectral approach combining Beurling’s theory of generalized primes / Beurling zeta functions with a projection-based diagnostic. Analyze density perturbations \(q\) associated with prime systems and test whether their low-degree polynomial energy is consistent with all non-trivial zeros on \(\mathrm{Re}(s)=1/2\).

The framework is a **classifier / diagnostic**, not a full proof of RH.

## Quick start

```bash
pip install -r requirements.txt
PYTHONPATH=src python3 -m pytest tests/ -v
PYTHONPATH=src python3 experiments/run_diagnostic.py
PYTHONPATH=src python3 experiments/run_multi_T.py --workers 86
PYTHONPATH=src python3 experiments/run_arithmetic_multi_T.py --workers 86
```

Arithmetic residual multi-\(T\): real primes \(\theta(x)-x\), compared to A0 critical-line
mode and B0 defect control. At \(T\le 16\) (\(x\le 10^7\)) arithmetic \(R_d\) **plateaus**
near ~0.19 (below defect 0.25) — see `results/arithmetic_multi_T.json`.

## Core math (shipped)

### Projection strength

1. Normalize a density perturbation \(q\) on the unit log-window \(u\in[0,1]\).
2. Project onto orthonormal **shifted Legendre** polynomials
   \(\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1)\).
3. Energy ratio
   \[
   R_d(q)=\frac{\|P_d q\|_{L^2}^2}{\|q\|_{L^2}^2}
   =\frac{\sum_{k=0}^d|\langle q,\varphi_k\rangle|^2}{\|q\|^2}.
   \]
4. Working projection strength (scaled)
   \[
   P(q)\;:=\;S_d(q)=T^{2(d+1)}\,R_d(q),
   \]
   with \(T\) the logarithmic window length (so \(S_d\) is \(O(1)\) under the
   RH decay heuristic \(R_d=O(T^{-2(d+1)})\)).

- **Low** \(R_d\) / controlled \(S_d\): high-frequency content — RH-like signature.  
- **High** \(R_d\): low-degree mass — defective / non-RH-like control.

### Legacy numbers

Earlier notes quoted **P(q)≈3.92** for zeta and threshold **≈29.5**. Those
used lost high-precision scripts. This reconstruction **does not hard-code
those values**; it reports \(R_d\) and \(S_d\) from the shipped path. See
`docs/STATUS.md`.

### Theorems A/B (precise status)

| Result | Status |
|--------|--------|
| **A₀** critical-line mode \(R_d\to0\) at \(O(T^{-2})\) | **Proved** (Lemma M3) |
| **A** arithmetic residual under RH | conditional / open |
| **B₀** persistent defect \(\Rightarrow R_d=\varepsilon^2\not\to0\) | **Proved** (Lemmas M2–M4) |
| **B** fast residual decay \(\Rightarrow\) RH | open |

Details: [`docs/THEOREMS_AB.md`](docs/THEOREMS_AB.md) · Proofs: [`docs/PROOFS_LEMMAS.md`](docs/PROOFS_LEMMAS.md) · Status: [`docs/STATUS.md`](docs/STATUS.md).

**Not an unconditional RH proof.**

## Repository layout

```
src/pbss/           # library: basis, projection, probes
tests/              # pytest driving real projection API
experiments/        # end-to-end RH-like vs defective run
results/            # last diagnostic JSON/TXT
docs/STATUS.md      # theorem status, session progress, open items
```

## Limitations

- Diagnostic, not a decisive RH proof.  
- Finite windows cannot exclude extremely high zeros.  
- Not a practical local primality sieve.  
- Large \(d\) or \(T\) is expensive.  
- Converse (low energy ⇒ RH) is essentially as hard as RH.

## Related private repos

- `perry-spirals`, `wieferich-hunts`, other archives under the same account.

---

*Private research. Not for public distribution without explicit permission.*
