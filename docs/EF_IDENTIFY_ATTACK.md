# EF identification attack: residual = modes + remainder

**Date:** 2026-07-26  
**Code:** `pbss.ef_identify` · **Campaign:** `results/ef_identify_attack/`  
**Entry:** `experiments/run_ef_identify_attack.py`

## Objective

Force Full Theorem A progress by **identifying** the arithmetic residual with the
truncated explicit-formula mode sum already shipped in the repo, quantifying the
remainder, and either closing the gap or naming a **sharp block**.

---

## Decomposition (shipped)

\[
q_H = \alpha\, m_T^{(N)} + r,\qquad
m_T^{(N)}=\sum_{n=1}^{N}a_n\cos(t_n T u-\alpha_n)
\]

with \(m\) from `explicit_formula_residual`, \(\alpha=\arg\min\|q-\alpha m\|\)
(trapezoid \(L^2\)). Optional joint LS with poly bulk \(\beta_0+\beta_1(u-\tfrac12)\).

Metrics: \(\mathrm{corr}(q,m)\), \(\|r\|^2/\|q\|^2\), \(E_d(r)/\|q\|^2\), \(R_d(q)\), \(R_d(r)\),
M5 majorant on \(m\), triangle majorant on \(E_d(r)\).

---

## Hypotheses tested

| ID | Residual object | Norm |
|----|-----------------|------|
| **H_theta_sqrt** | \(\theta(x)-x\) | \(/\sqrt{x}\) |
| **H_psi_sqrt** | \(\psi(x)-x\) (prime powers) | \(/\sqrt{x}\) |
| **H_theta_x** | \(\theta(x)-x\) | \(/x\) |
| **H_psi_x** | \(\psi(x)-x\) | \(/x\) |

Also: detrend `deg1` vs `none`; poly bulk on/off; \(N\in\{5,10,20,40\}\) (and 50);
\(T\in\{12,\ldots,20\}\) on primes to \(10^{10}\).

Model control: if \(q=m_T^{(N)}\), remainder is machine zero (identity holds).

---

## Results (decisive campaign)

**Scale:** 64 arithmetic rows = 5 \(T\) × 4 \(N\) × 4 \(H\), deg1, primes \(\le 10^{10}\).

### Winner among hypotheses

| Hypothesis | mean \|corr\| | mean \(\|r\|^2/\|q\|^2\) | mean \(E_d(r)/\|q\|^2\) |
|------------|--------------:|------------------------:|------------------------:|
| **H_theta_sqrt** | **0.654** | **0.569** | **0.208** |
| H_psi_sqrt | 0.594 | 0.644 | 0.369 |
| H_theta_x | 0.146 | 0.978 | 0.859 |
| H_psi_x | 0.124 | 0.984 | 0.873 |

**Progress:** \(\theta/\sqrt{x}\) is the correct *scale class* for the shipped mode
amplitudes (beats \(\psi\) and \(/x\) by a lot). Modes are not noise: correlation
and L² capture **improve with \(N\)**.

### \(N\)-dependence for H_theta_sqrt (core discovery)

| \(N\) | \|corr\| | \(\|r\|^2/\|q\|^2\) | \(E_d(r)/\|q\|^2\) | \(R_d(r)\) |
|------:|---------:|-------------------:|------------------:|----------:|
| 5 | 0.572 | 0.673 | 0.203 | 0.301 |
| 10 | 0.636 | 0.595 | 0.207 | 0.348 |
| 20 | 0.685 | 0.531 | 0.210 | 0.396 |
| 40 | 0.723 | 0.478 | 0.213 | 0.447 |

- **L² remainder fraction falls** as more zeros enter \(m\) (~0.67 → 0.48).  
- **Low-degree remainder mass \(E_d(r)/\|q\|^2\) is flat** (~0.21, spread \(\sim 0.01\)).  
- Poly bulk joint fit does **not** move \(E_d(r)\) (deg1 already killed linears).  
- No-detrend collapses \(\theta\) match (corr ~0.15); deg1 is required for capture.

Bounds: triangle majorant holds; M5 on modes decays in \(T\) (model-side OK).

---

## Sharp block (named)

### `LOW_DEGREE_MASS_INVARIANT_TO_ZERO_TRUNCATION_N`

**Statement.** Under the best identification (H_theta_sqrt, LS scale, deg1), the
truncated critical-line mode sum captures an increasing share of \(L^2\) mass as
\(N\) grows, but the **Legendre \(V_d\) mass of the remainder stays \(\approx 0.21\)
independent of \(N\in[5,40]\)**. That \(V_d\) mass is therefore **not in the span** of
the first \(N\) shipped EF modes.

**What this is not:** “EF is useless.” Modes *do* align with the residual
(\(\mathrm{corr}\to 0.72\)).  

**What fails for Full A:** the step “subtract zeros ⇒ \(R_d\to 0\)” — low-degree
energy refuses to leave with \(N\).

### What would unblock Full A

1. **Secondary main terms** in \(m\) (smooth EF contributions, not only zeros) so
   \(E_d(q-\alpha m_{\mathrm{full}})/\|q\|^2\to 0\); or  
2. **Change \(q\)** so its \(V_d\) content is proved \(O(T^{-2})\) (weight-in-definition,
   different normalization, smoothed \(\psi\)); or  
3. **Different oscillatory law** for amplitudes/phases matching the exact windowed
   transform of \(\theta\) on \(x=e^{uT}\).

---

## Bounds (shipped)

| Bound | Role |
|-------|------|
| Model identity | \(q=m\Rightarrow r=0\) (proved by construction; tested) |
| Triangle \(E_d(r)\le 2E_d(q)+2\alpha^2 E_d(m)\) | Always true; plumbing check |
| M5 on \(m\) | \(R_d(m)=O(T^{-2})\) majorant; shrinks in \(T\) |

No bound forces \(E_d(r)/\|q\|^2\to 0\) for arithmetic \(q\) under present \(m\).

---

## Status for Full Theorem A

| Piece | Status after attack |
|-------|---------------------|
| Mode object \(m_T^{(N)}\) | Correct for **model**; M5 holds |
| Best residual class for matching | **H_theta_sqrt** (shipped default) |
| L² identification | **Partial — improves with \(N\)** |
| \(V_d\) identification | **Blocked — flat in \(N\)** |
| Full A under RH | Still needs secondary terms or residual redesign |

---

## Reproduce

```bash
PYTHONPATH=src python3 -m pytest tests/test_ef_identify.py -v
PYTHONPATH=src python3 experiments/run_ef_identify_attack.py \
  --out-dir results/ef_identify_attack \
  --primes-path results/prime_checkpoints/primes_le_10000000000.npy
```
