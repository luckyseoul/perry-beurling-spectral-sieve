# Research note: Why does arithmetic \(R_d\) plateau?

**Date:** 2026-07-26  
**Campaign entry:** `experiments/run_open_plateau.py`  
**Artifacts:** `results/open_plateau/` (per-class JSON + `PHASE_*_COMPLETE` resume stamps)  
**Related goal:** [`goals/OPEN_GOAL.md`](goals/OPEN_GOAL.md)

## Explicit non-claim

**This document does not claim a proof of the Riemann Hypothesis.**  
No unconditional RH theorem is established. The work is diagnostic: multi-\(T\) numerics
on the arithmetic residual and controls. Open conclusions are intentional.

---

## Plateau fact (from shipped campaigns)

On the full prime table \(x_{\max}=5\times10^{10}\) (\(n_{\mathrm{primes}}=2\,119\,654\,578\)):

| Setup | Typical \(R_4\) | Notes |
|-------|----------------:|-------|
| Arithmetic \(q_T=(\theta-x)/\sqrt{x}\), detrend deg1 | **\(\sim 0.15\)–\(0.19\)** | Soft plateau through \(T\approx\log(5\times10^{10})\approx 24.6\) |
| Same residual, no detrend (`raw`) | \(\sim 0.96\)–\(0.97\) | Bulk trend dominates |
| Pure / finite CL controls (prior campaigns) | collapse toward small \(R_d\) | Contrast with arithmetic plateau |
| MC random spectral defects (this campaign) | \(\sim 0.795\) flat in \(T\) | Instrument floor well above arithmetic deg1 |

Focus number from extend-\(x\) / overnight residual: \(R_4\) deg1 at \(T\approx 24.6\) is \(\approx 0.146\).
Multi-\(T\) open-plateau **scale** class (24 \(T\) near \(T_{\max}\), deg1): mean \(R_4\approx 0.176\)
with range \(\approx[0.146,0.194]\) — **no clear decay** with denser \(T\) sampling.

**Question:** Is the plateau bulk trend, low zeros, residual definition, basis choice, finite-\(x\)
artifact, or an instrument that is simply insensitive?

---

## Intervention catalog (≥5 classes)

Each class has a hypothesis, multi-\(T\) numeric \(R_d\), and durable JSON under
`results/open_plateau/<class>/`.

### 1. Peel (zero-fit / mode strip) — breadth

- **Hypothesis:** If low zeros dominate low-degree mass, stripping the first \(N\) CL modes
  should reduce \(R_d\) toward A₀ levels as \(N\) grows at large \(T\).
- **What ran:** Cited overnight peel battery — **2048** multi-\((T,N,d,\mathrm{detrend})\) rows on
  \(x_{\max}=5\times10^{10}\); \(T\in[14.6,24.6]\) (16 values); \(N\) up to 50.
- **What \(R_d\) did:**
  - detrend `none`: mean \(R_d\approx 0.983\) (stuck near 1).
  - detrend `deg1`: mean \(\approx 0.303\), range \(\approx[0.034,0.484]\).
  - **Increasing \(N\) does not collapse \(R_d\)**; at deg1, mean \(R_d\) *rises* with \(N\)
    (\(N=0\to\approx0.15\), \(N=50\to\approx0.40\)). Mode stripping removes structure that was
    *helping* the residual look smoother in the Legendre sense, not the plateau mass.
- **Needle:** **Failed** as a path to A₀. Peeling model modes is not enough.

### 2. Whiten (detrend / smooth / taper) — breadth

- **Hypothesis:** Slow bulk (mean/linear trend, light smoothing, edge taper) inflates
  low-degree mass; whitening should lower \(R_d\) if the plateau is bulk not zeros.
- **What ran:** 72 rows, 12 \(T\) near \(T_{\max}\), variants raw / deg0 / deg1 / smooth 5,15 / Hanning taper.
  Elapsed \(\approx 1090\) s on full table.
- **What \(R_d\) did (means over \(T\)):**

  | variant | mean \(R_4\) |
  |---------||-------------:|
  | raw | 0.966 |
  | deg0 | 0.204 |
  | deg1 | 0.168 |
  | deg1_s5 | 0.191 |
  | deg1_s15 | 0.234 |
  | **deg1_taper** | **0.078** |

- **Needle:** Detrend is **necessary** (raw \(\to\) deg1: large drop). Extra smoothing **hurts**.
  Edge taper **moves the needle** most (mean \(\sim0.078\)) but is partly a windowing effect —
  not a proof that arithmetic mass vanished; treat as best *processing* recipe, not asymptotic A₀.

### 3. Measure (residual definition) — breadth

- **Hypothesis:** Plateau depends on normalization of \(\theta-x\); a better \(w(x)\) might reveal A₀-like decay.
- **What ran:** 48 rows, 12 \(T\), norms \(\sqrt{x}\), \(x\), plain, \(\log x\); deg1 detrend; elapsed \(\approx 743\) s.
- **What \(R_d\) did:**

  | norm | mean \(R_4\) |
  |------|-------------:|
  | **sqrt** \((\theta-x)/\sqrt{x}\) | **0.168** |
  | logx | 0.884 |
  | plain | 0.891 |
  | x | 0.917 |

- **Needle:** Classic \(\sqrt{x}\) normalization is **best** among tested; alternatives sit near 0.9.
  Definition choice matters a lot; wrong norm looks like a fake hard plateau.

### 4. Basis / projector — breadth

- **Hypothesis:** Plateau energy is concentrated in lowest Legendre modes; raising \(d\) or
  high-passing \(V_k\) should change \(R_d\) character if the defect is low-degree bulk.
- **What ran:** 120 rows, 12 \(T\), degrees \(0,1,2,4,6,8,12\) + high-pass \(k=0,1,2\) at \(d=8\);
  elapsed \(\approx 1911\) s.
- **What \(R_d\) did:**

  | variant | mean \(R_d\) |
  |---------||-------------:|
  | d0, d1 | 0 (by construction after deg1 detrend) |
  | d2 | 0.079 |
  | d4–d12 | \(\approx 0.168\)–\(0.173\) (flat) |
  | hp2_d8 | 0.100 |

- **Needle:** Mass sits above degree 2 and **plateaus through high \(d\)**. High-pass \(k=2\) only
  modestly helps. Not a pure degree-0/1 bulk problem after deg1 detrend.

### 5. Scale (dense multi-\(T\) on \(5\times10^{10}\)) — **deep**

- **Hypothesis:** Plateau is a finite-\(x\) artifact; denser multi-\(T\) near \(T_{\max}\) across
  construction variants should show decay if asymptotic A₀ applies.
- **What ran:** **24 \(T\) × 6 variants = product 144** residual rebuilds on the full table
  (sequential mmap; pool OOM avoided). Elapsed \(\approx 2203\) s (\(\sim 37\) min).
- **What \(R_d\) did:** Same ordering as whiten; deg1 mean \(\approx 0.176\) across denser \(T\) —
  **no systematic decay** toward 0 as \(T\to T_{\max}\).
- **Depth bar:** product \(\ge 20\times 5\), multi-\(T\) on full \(5\times10^{10}\). Marked `deep: true`.

### 6. MC randomization — **deep**

- **Hypothesis:** Randomized spectral defects keep high \(R_d\) independent of \(T\); large-\(N\) MC
  quantifies the defect floor vs arithmetic plateau levels.
- **What ran:** **50 000 000** trials total (10 \(T\) × 5 000 000), degrees \(2,4,6,8\), \(n_{\mathrm{points}}=8192\),
  **86 workers**. Elapsed \(\approx 6180\) s (\(\sim 1.7\) h). Entry: `experiments/run_mc_stress.py`.
- **What \(R_d\) did:** mean \(R_4\approx 0.7947\) **flat in \(T\)** (range of means \(\sim[0.7945,0.7950]\)).
- **Needle:** Instrument is stable and **insensitive to \(T\)** under random defects. Arithmetic deg1
  plateau (\(\sim0.15\)–\(0.19\)) sits **well below** the MC defect floor (\(\sim0.79\)), so the arithmetic
  signal is not “as random as a defect.” But MC does not collapse either — control only.

### 7. Beurling system battery — **deep**

- **Hypothesis:** If the diagnostic is meaningful, defective Beurling systems keep high \(R_d\)
  while ordinary stays lower; a large battery stress-tests separation vs noise.
- **What ran:** Multi-core ProcessPool (mmap ordinary primes, module-level job).
  **Final expensive run:** **2000** dense systems × 8 \(T\) at \(x_{\max}=10^{8}\),
  **24 workers**, \(n_{\mathrm{points}}=8192\), elapsed \(\approx 2113\) s (\(\sim 35\) min),
  product **16 000** rows. Cheap large-gap pilots at smaller \(x\) finished in minutes even at
  100k–400k systems and were enlarged with denser \(x_{\max}\) for honest elapsed depth.
  Resume stamp `PHASE_BEURLING_COMPLETE`. Durable `results/beurling_battery/` not overwritten.
- **What \(R_d\) did:** rh_like mean \(\approx 0.34\); defective mean \(\approx 0.90\) across
  multi-\(T\) (clear separation). Scorecard in `results/open_plateau/beurling/beurling.json`.
- **Depth bar:** \(n_{\mathrm{systems}}=2000\ge 500\) × multi-\(T\); elapsed \(\ge 30\) min; product 16 000.

---

## Judgment

| Intervention | Moved needle? | Verdict |
|--------------|---------------|---------|
| Peel more zeros | No (often worse) | Reject as path to A₀ at this \(x\) |
| Detrend deg1 | Yes (vs raw) | **Required** preprocessing |
| Extra smooth | No (hurts) | Avoid for low \(R_d\) |
| Edge taper | Yes (down to \(\sim0.08\)) | Best *processing* knobs; interpret carefully |
| Norm \(\neq\sqrt{x}\) | No (worse) | Keep \(\sqrt{x}\) |
| Higher \(d\) / high-pass | Weak | Mass not pure low-\(d\) bulk after deg1 |
| Dense multi-\(T\) scale | No decay | Finite-\(x\) alone does not explain plateau away |
| MC defects | Stable high floor | Instrument OK; arithmetic below defect floor |
| Beurling defective vs ordinary | Separation holds | Diagnostic distinguishes density defects |

### Best current residual recipe (operational, not a theorem)

1. Use primes / \(\theta\) cumulative on largest available table (\(5\times10^{10}\)).
2. Residual \(q_T=(\theta-x)/\sqrt{x}\) with **deg1** detrend (optional light Hanning taper if reporting a
   “windowed” diagnostic).
3. Report multi-\(T\) \(R_d\) for \(d=4\) (and a small \(d\) grid); do **not** claim A₀ collapse from peel alone.
4. Always pair with Beurling ordinary vs defective scorecard and MC floor for instrument health.

### Why the plateau persists (working explanation)

Evidence favors:

1. **Not** removable by fitting the first \(N\le 50\) model zeros (peel fails).
2. **Not** a pure degree-0 bulk after deg1 (basis: \(d\ge 4\) still \(\sim0.17\)).
3. **Partially** window/edge driven (taper helps a lot — residual is sensitive to endpoints).
4. **Not** fixed by denser \(T\) alone near \(T_{\max}\) on this table (scale deep axis).
5. **Consistent with** a genuine arithmetic oscillatory / irregularity content whose low-degree
   Legendre mass stays \(O(10^{-1})\) through \(x\sim5\times10^{10}\), **distinct from** pure CL collapse
   and from fully defective Beurling systems (\(R_d\sim 1\)).

This is a **research judgment from numerics**, not a proof that RH is true or false.

---

## What full Theorem A still needs

1. **Zero-sum control with rigorous truncation error** — peel numerics alone do not bound the
   explicit-formula remainder uniformly enough for an \(A_0\) theorem.
2. **Analytic (not only numeric) control of endpoint / bulk contributions** — taper success
   suggests boundary terms matter; need a clean weight class where endpoint pollution is proven small.
3. **Larger \(x\) or stronger asymptotic diagnostics** if RAM allows (this host rejects \(10^{11}\) primes
   for RAM; \(10^{12}\) infeasible). GPU residual path optional when prefix fits V100.
4. **Bridge from \(R_d\) diagnostic to classical \(\psi\)/\(\theta\) error terms** with explicit constants.
5. Keep RH **open** until a complete written proof exists outside this diagnostic suite.

---

## Engineering notes

- **Resume:** `results/open_plateau/PHASE_{PEEL,WHITEN,MEASURE,BASIS,SCALE,MC_RAND,BEURLING}_COMPLETE`
  + `open_plateau_state.json`.
- **Multi-core:** MC **86** workers; Beurling **86** workers (long-lived pool; \(x_{\max}\) sized so
  full machine fits 60 GiB); residual classes sequential or small waves for full-table mmap safety.
- **Durable Beurling battery** under `results/beurling_battery/` was **not** overwritten by smoke.
- **Non-goal:** wall-clock padding. Depth enforced by trial/row products and enlarge-if-fast policy.

---

## Pointers

- Status board: [`STATUS.md`](STATUS.md)
- **Theorem-A scaffolding (weight class + truncation remainder):** [`THEOREM_A_SCAFFOLD.md`](THEOREM_A_SCAFFOLD.md)
- Overnight floors: [`goals/OVERNIGHT_GOAL.md`](goals/OVERNIGHT_GOAL.md), `results/overnight_marathon/`
- Open goal: [`goals/OPEN_GOAL.md`](goals/OPEN_GOAL.md)
- Lemmas / theorems: `THEOREMS_AB.md`, `PROOFS_LEMMAS.md`
