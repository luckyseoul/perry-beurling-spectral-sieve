# External session check and open-item attempt (2026-09-12)

**Status:** external verification record and negative-results note from an
assistant session on a clean clone at `09fce0f` (up to date with `origin/main`).
**Not a proof of RH.** Not peer review. No new theorem is claimed here.
This file records (1) independent re-checks of shipped claims, (2) the result of
an attempt on the two live open items, and (3) what remains missing for each.

Machine-readable status checked during the session:
`pbss.ab_closure.package_status()` reports `full_B = open_normalized_converse`,
`normalized_unattained_boundary_case = open`,
`absolute_moment_growth_bound = open_unconditional`,
`one_sided_moment_bound = rh_equivalent_open`, `rh = open`.

## 1. Independent re-checks (all with numbers)

| Check | Method | Result |
|---|---|---|
| Test suite | `PYTHONPATH=src python3 -m pytest tests/ -q -n 40` | **373 passed** in 10.5 s |
| Shipped moment evaluator | own mpmath 45-digit interval implementation (own formulas, complete prime prefixes from `results/prime_checkpoints/primes_le_10000000000.npy`), compared with `pbss moments --max-prime-limit 50000000000` | T = 12, 14, 16, 18, 20 (up to 25 614 562 primes): max relative difference **7.5e-12** on coefficients, **3.5e-15** on the norm, **2.5e-14** on `R_d`. Extends the in-repo corroboration (`experiments/check_arithmetic_moments.py`, x <= 1e4) to T <= 20 |
| `pbss moments` CLI | T = 8, 10, 12, 14, 16, 18, 20 with explicit `--max-prime-limit` | runs; T = 10 reproduces the shipped artifact decimals |
| Phase-trough constants (A=0.5, t=3, fresh parameter point) | own exact closed-form scan of trough windows in T in [1000, 1300]; compared with `A d^2(d+1)^2(d+2)^2/(3(A^2+t^2)^2)` and with `pbss.off_critical_exact.asymptotic_constants` | d = 1..4 observed minima exceed the prediction uniformly by **1.00025** (O(1/T) correction); repo constants equal the closed form to all printed digits |
| Unattained-spectrum example | own closed-form evaluation of the explicit series (`e^{-n}`, `a_n = 1/4 - 1/n`, `gamma_n = e^n`) | T = 100/400/1600/6400: measured decay slope of `log R_d` vs `sqrt(T)` is **-1.646** against the proved **-1.6569**; two-sided bound (eq. 12) holds with `V/(e^{2AT}S2) in [1.05, 1.66]`; `log(V/e^{2AT})/(-4 sqrt(T)) -> 0.99` |
| Cross-node reproducibility | identical sweep jobs re-run on soulkiller, jellyfin, orin, nuka | 24 replicated jobs **bit-identical** (difference 0.0) |

These checks support the shipped arithmetic and model claims above; they say
nothing about the unresolved limit statements.

## 2. Attempt on open item 1: normalized converse, unattained-boundary case

**Where the proof stops.** The attained-boundary theorem uses two matched
ingredients: a signed-moment lower bound of size `A_{k,rho} e^{aT}/T` for one
zero `rho` with `a = Re rho - 1/2`, and a total-energy upper bound
`V(T) <= K e^{2aT}` valid for all large T. Under an unattained supremum
`Theta > 1/2`, every fixed zero has `a < Theta - 1/2` while only
`V(T) = o(e^{2(Theta-1/2)T})` holds. The two rates no longer match. The
**exact missing lemma** is a matching-rate statement at the moment-peak times:

> (L*) there is a zero `rho` and a sequence `T_n -> inf` on which the signed
> moment attains its scale `e^{aT_n}/T_n` and simultaneously
> `V(T_n) <= K e^{2aT_n}`; equivalently, no competitor zero `rho_i` with
> `Re rho_i > Re rho` has `(Re rho_i - Re rho) T_n - log|rho_i|` unbounded
> along that sequence.

Such a statement is ordinary-prime arithmetic input: it constrains the actual
configuration of zeros near the boundary. It is not implied by discreteness, the
local counting upper bound, the residue coefficients `-1/rho`, or a quantified
vanishing prefactor, as `docs/UNATTAINED_SPECTRUM_OBSTRUCTION.md` proves.
Attainment of `Theta` would give (L*) trivially; whether `Theta` is attained is
itself open, and no off-critical zero of `zeta` is known, so no finite
computation can supply the missing input.

**What this session added (exact-arithmetic sweep, 9 families, 26 cases,
sharded across the four-node mesh; see section 4).** Family:
`Q(y) = sum_n gamma_n^{-q} e^{a_n y} cos(gamma_n y)`, `gamma_n = exp(n^beta)`,
`a_n = A - delta_n`, `A = 1/4` (`q = 1` gives arithmetic-size residues
`1/|rho|`). Measured decay slopes `d log R_d / d sqrt(T)` (natural log):

| family | beta | amplitude q | deficit delta_n | slope |
|---|---|---|---|---|
| reference (= repo example) | 1 | 1 | 1/log gamma | **-1.646** (proved -1.657) |
| denser ordinates | 3/4 | 1 | 1/log gamma | -1.739 |
| denser ordinates | 1/2 | 1 | 1/log gamma | -1.704 |
| smaller residues | 1 | 1/2 | 1/log gamma | -2.010 |
| larger residues | 1 | 3/2 | 1/log gamma | -1.408 |
| faster deficit decay | 1 | 1 | 0.5 n^{-1} | -1.159 |
| slower deficit decay | 1 | 1 | n^{-1/2} | -3.485 |
| faster deficit decay | 1 | 1 | n^{-2} | -0.353 |

Observations: the obstruction is robust in `R_d -> 0` for every family; density
changes the slope only mildly; the deficit profile moves it most. In every case
the moment-dominant and energy-dominant index scales coincide, i.e. the models
are "balanced": the failure is the missing `1/gamma` dilution of the moments
against the undiluted energy, exactly the mechanism the repo identifies.

**A structural limit of this model class.** An absolutely convergent series with
residue-size amplitudes needs `sum_n 1/gamma_n < inf`, which forbids
Riemann-von-Mangoldt density `N(T) ~ T log T` (there `1/gamma_n ~ 1/(n log n)`
and the sum diverges). `zeta`'s zero sum is only conditionally convergent. A
saddle-point estimate with `zeta`-density and `delta log gamma = 1` gives
`R_d ~ exp(-2 sqrt(2) sqrt(T))` times powers of `T`, i.e. still vanishing, but
that estimate is heuristic and is not a construction. So the model evidence
blocks every argument that uses only the listed spectral hypotheses; it cannot
decide the ordinary-prime case.

**A sharper finding (this session): the dichotomy is attainment, not deficit size.**
A second exact-arithmetic campaign (59 cases; deficit profiles `delta_n = c n^-p`,
`e^-n`, `e^-n^2`, and attained families with the top mode exactly on the
boundary; T up to 409600; sharded across the mesh) gives:

| configuration | `T R_d` at T = 100 | 6400 | 409600 |
|---|---|---|---|
| attained (`a_8 = A`, others below) | 1.26e-7 | 2.91e-7 | 4.53e-7 |
| unattained, `delta = 0.25 n^-1` | 5.84e-8 | 1.72e-28 | — |
| unattained, `delta = n^-1` | 2.56e-9 | 1.33e-57 | — |
| unattained, `delta = 0.25 n^-3` (slowest onset) | 8.55e-8 | 1.59e-7 (plateau) | 8.9e-11 (at T=25600) |

The attained family sits at the predicted single-mode scale
`2A/gamma_top^2 = 0.5 e^-16 = 5.6e-8` (same order). Every unattained family
decays; the decay onset is delayed by faster-vanishing deficits but always
arrives. Uniformly positive deficits do not create a floor
either: `delta_n = 0.05(1 + 1/n)` keeps the supremum unattained and
`T R_d` falls from `6.4e-8` (T=100) to `7.6e-14` (T=6400). A uniform deficit
only rescales the effective time: `delta_n = 0.01(1 + 1/n)` is still near-flat
at T=6400, matching the reference family at `T ~ 64`.

Consequently, **in this model class the positive floor holds if and only if the
supremum of the pole real parts is attained.** The decisive input for the
ordinary-prime case is therefore not an estimate but the arithmetic fact
"the supremum of zero real parts is attained" (or a proof that the approach to
it forces the alignment needed by (L*)). No such theorem is known; the
literature confirms that not even a fixed zero-free strip `Re s > 1/2 + c`
is known (the "weak Riemann hypothesis").

**Verdict for item 1:** open, and the model evidence says the unattained case
should not be expected to be closable by an estimate: if the ordinary zeros
approach an unattained boundary the way every tested model does, the normalized
converse fails there. Record it as "open; model-supported counterexample class"
rather than "open pending an estimate".

## 2b. Generalized obstruction theorem (proved this session)

The repo proves the obstruction for one deficit profile (`delta_n = 1/n`).
Assembling the note's own lemmas with an elementary dominated-sum estimate
extends it to a two-parameter family. **Statement proved here:**

> Let `A = 1/4`, `c > 0`, `p > 0`, `gamma_n = e^n`, `a_n = A - c n^-p`, and
> `Q_{c,p}(y) = sum_{n>=8} e^-n e^{a_n y} cos(gamma_n y)`. Then `Q_{c,p}` has
> genuine simple Laplace poles at `a_n +- i gamma_n` with arithmetic-size
> residues, the real-part supremum `A` is unattained, and for every fixed
> `d >= 2` the continuous detrended ratio satisfies
> `R_d(T) <= C_{c,p,d} T^theta exp( -K(c,p) T^{1/(p+1)} ) = O_N(T^-N)` for
> every `N`, where
> `K(c,p) = 2(1 + 1/p) (cp/2)^{1/(p+1)} (2 - 2^{1/(p+1)})`.

Ingredients: (i) the note's two-sided Gram bounds (eqs. (11)-(12)), which hold
for `gamma_n = e^n` and arbitrary coefficients; (ii) the note's integration-
by-parts moment bound (eqs. (13)-(14)) for every fixed `g in C^1`; (iii) a
dominated-sum lemma: for `f(n) = 2n + b t n^-p`, `sum_n e^{-f(n)}` is at most
`poly(t) exp(-min_n f)` (split at `n*/2` and `2n*`, expand quadratically near
`n* = (bpt/2)^{1/(p+1)}`). The assembly gives
`T R_d <= C T (sigma(T))^2 / Sigma(T)` with
`sigma(t) = sum_n e^{-2n - c t n^-p}` and `Sigma(t) = sum_n e^{-2n - 2c t n^-p}`,
whose exponents differ by exactly the factor `2 - 2^{1/(p+1)}`.

**Check of the constant.** `K(1,1) = 4(sqrt 2 - 1) = 1.6568542494923801`
exactly reproduces the repo's proved constant for the original example.
Measured slopes of `ln(T R_d)` against `T^{1/(p+1)}` (T up to 6400):

| c | p | measured slope | predicted `-K(c,p)` |
|---|---|---|---|
| 4 | 1 | -3.3402 | -3.3137 |
| 1 | 1 | -1.5907 | -1.6569 |
| 0.25 | 1 | -0.6922 | -0.8284 |
| 4 | 2 | -3.2268 | -3.5244 |
| 1 | 2 | -1.8504 | -2.2202 |
| 4 | 3 | -2.3554 | -3.3839 |

Every measured slope is below the predicted magnitude and converges upward in
T, as expected for a one-sided bound whose prefactor `T^theta` still matters in
the tested range. The `p = 1` family (the repo's case) is the converged one.

**Dichotomy, now sharp.** Combining this theorem with the single-mode
dominance argument (a mode exactly on the boundary keeps `T R_d` bounded below
by `2a_top/gamma_top^2` whenever the top mode has the smallest ordinate among
attained modes) gives, for this class: the positive floor holds if and only if
the real-part supremum is attained. Numerically the attained family sits at
`T R_d ~ 1.3e-7 .. 4.5e-7` over `T = 100 .. 409600`, while every unattained
profile decays (fastest `e^{-1.66 sqrt(T)}` for `1/n`, slowest a delayed onset
for `n^-3`).

**Why this does not transfer to zeta by itself.** Both directions of the
transfer stall on the same missing input: the energy lower bound needs
diagonal dominance in the zero Gram matrix, whose off-diagonal row sums
`sum_j min(1, 2/|gamma_i - gamma_j|)` diverge for zeta's dense ordinates
(they converge for the model's exponential gaps); and the energy upper bound in
the unattained case carries no rate. So the ordinary-prime unattained case is
blocked in both directions by quantitative control of the oscillatory zero sums
near the boundary -- pair-correlation-type arithmetic input -- not by any
bookkeeping gap.

## 3. Attempt on open item 2: unconditional one-sided moment bound

The stated target (for one fixed `k`), `c_k(T) <= C_eps e^{eps T}` eventually
for every `eps > 0`, is **exactly RH** (B-ABS plus the signed oscillation
theorem; the pole-order bookkeeping was re-derived and is correct). Two
consequences recorded for later use:

- **Sharpened equivalent form.** Because continuous Theorem A gives
  `R_d = O(T^{-2})` under RH and the `beta = 1/2` case of the unit-interval Gram
  bound in `docs/ATTAINED_BOUNDARY_CONVERSE.md` gives `V(T) = O(T)` under RH,
  one has `|c_k(T)| <= ||q_T|| sqrt(R_d) = O(1/T)` under RH. Hence, for each
  fixed `k >= 2`,
  `RH  <=>  (exists delta > 0 : |c_k(T)| = O(T^{-delta}))`,
  and likewise with either one-sided sign `+/- c_k(T) <= C T^{-delta}`.
  Any inverse-power decay, not just subexponential growth, is RH-equivalent.
- **No intermediate rate is reachable.** `c_k(T) = O(e^{cT})` for a fixed
  `c < 1/2` would imply (by the signed oscillation theorem) that no zero has
  real part `> 1/2 + c`, i.e. a fixed zero-free strip. No zero-free strip for
  `zeta` is known unconditionally; the classical zero-free region is a curve
  `1 - c/log t`. So even the weakest fixed-rate improvement below `1/2` is a
  famous open problem, and the stated target is RH itself.

**Verdict for item 2:** not solvable without proving RH; no unconditional
partial rate is within reach of current zero-free-region technology. The
obstacle for any partial result is exactly the unknown zero-free strip
("weak Riemann hypothesis"): it is not known unconditionally that
`zeta(s) != 0` for `Re s > sigma_0` for any `sigma_0 < 1`.

## 4. Compute mesh used

Inventory (`~/.grok/skills/use-available-compute/scripts/mesh-inventory.sh`):

| node | CPU | RAM | accelerator | used for |
|---|---|---|---|---|
| soulkiller (local) | 2x Xeon E5-2696 v4, 88 threads | 60 GiB | Tesla V100-SXM2 16 GB (CUDA) | sweep shard (26 jobs), pytest, 44-worker mpmath moment check, phase-trough scan |
| nuka | Ryzen 7 5700X3D, 16 threads (3D V-Cache) | 15.6 GiB | Radeon RX 9070 XT (ROCm 7.15) | heavy beta=1/2 shard |
| jellyfin | Ryzen 7 5700G, 16 threads | 14.5 GiB | Intel Arc A380 (no NVIDIA driver) | sweep shard |
| orin | Jetson Orin Nano Super, 6x A78AE | 7.4 GiB | Ampere sm_87, CUDA 13.2, DLA (offline) | sweep shard |

`lucky` is DNS only and was not used. The workloads here are high-precision
serial arithmetic; GPUs were therefore not used (no dense batched kernel), and
each node ran shards sized to its own core count. Results were replicated on two
nodes per job; all 23 replications agreed to the last printed digit.

## 5. Reproduce

Scripts and raw outputs: `~/scratch/pbss_attempt/` on soulkiller
(`model_family.py`, `run_jobs.py`, `my_check.py`, `phase_trough_check.py`,
`out_full/`). Key commands:

```bash
python3 my_check.py --T 12 14 16 18 20 --workers 44 --digits 45
python3 phase_trough_check.py
python3 run_jobs.py jobs_all.json out_full 12     # 26-job family sweep
cd /home/nick/perry-beurling-spectral-sieve && \
  PYTHONPATH=src python3 -m pbss moments --T 20 --degree 4 --max-prime-limit 50000000000
```

On these scripts, note: `pbss moments` needs `--max-prime-limit` above
`exp(T)` (default 1e6 caps T at 13.8).

## 6. Non-claims

No unconditional RH claim, no new theorem, no claim that the sweeps prove
anything about `zeta`. The sweep families have no Euler product, no functional
equation, and at most polylogarithmic pole counts. The heuristic density
estimate is labeled as such. Item 1 remains open; item 2 is RH-equivalent and
open.
