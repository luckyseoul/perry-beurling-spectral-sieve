# External session record: checks, sweep, and obstruction-family theorem

**NOT AN UNCONDITIONAL PROOF OF RH.** This directory is an external session
record: independent re-checks of shipped claims, a mesh-run exact-arithmetic
sweep over obstruction families, and a generalized obstruction theorem. It is
kept on a branch, not on `main`, pending review.

Write-up: [`docs/EXTERNAL_CHECK_20260912.md`](../../docs/EXTERNAL_CHECK_20260912.md).

## Contents

| File | Role |
|---|---|
| `model_family.py` | Exact closed-form `R_d(T)`, moments and energy for explicit exponential-spectrum families `Q(y) = sum A_n e^{a_n y} cos(gamma_n y)`; `--selftest` validates every identity against brute-force quadrature |
| `run_jobs.py` | Parallel (ProcessPoolExecutor) runner over a job list; writes one JSON per case |
| `my_check.py` | Independent 45-digit mpmath re-implementation of `theta` moments, compared against the shipped `pbss moments` evaluator for T = 12..20 |
| `phase_trough_check.py` | Independent scan of the pure-mode phase-trough constant `liminf T^3 R_d` |
| `jobs_*.json` | Job definitions for each campaign |
| `results/family_sweep/` | 26 jobs: ordinate density, amplitude, and deficit-profile families |
| `results/deficit_grid/` | 44 jobs: `delta_n = c n^-p` over c in {0.25,1,4}, p in {1,2,3}, plus attained families |
| `results/dichotomy/` | 18 jobs: super-fast vanishing deficits and uniformly-positive deficits |
| `results/large_T/` | 15 jobs: the same families at T = 25600..409600 |
| `results/my_check.json`, `results/phase_trough_check.json` | Independent-check outputs |

## Headline results

1. **Generalized obstruction theorem.** For `A = 1/4`, `c > 0`, `p > 0`,
   `gamma_n = e^n`, `a_n = A - c n^-p`, the explicit series has genuine simple
   Laplace poles with arithmetic-size residues, an unattained real-part
   supremum, and `R_d(T) <= C T^theta exp(-K(c,p) T^{1/(p+1)}) = O_N(T^-N)`,
   where `K(c,p) = 2(1+1/p)(cp/2)^{1/(p+1)}(2 - 2^{1/(p+1)})`.
   `K(1,1) = 4(sqrt2 - 1)` reproduces the repo's proved constant exactly.
2. **Attainment dichotomy.** In this class the positive floor
   `limsup T R_d > 0` holds iff the real-part supremum is attained:
   attained families sit at `T R_d ~ 1.3e-7..4.5e-7` for T up to 409600; every
   unattained family decays.
3. **Transfer block for zeta.** Both directions stall on the energy `V(T)`:
   the lower bound needs zero-Gram diagonal dominance (row sums
   `sum_j min(1, 2/|gamma_i - gamma_j|)` diverge for dense ordinates); the
   unattained upper bound `o(e^{2(Theta-1/2)T})` carries no rate.

## Reproduce

```bash
python3 model_family.py --selftest                      # identity checks vs quadrature
python3 run_jobs.py jobs_grid.json out_grid 12          # any job list
python3 my_check.py --T 12 14 16 18 20 --workers 44 --digits 45
python3 phase_trough_check.py
```

`my_check.py` reads the prime checkpoint at
`results/prime_checkpoints/primes_le_10000000000.npy`. Higher T values need
`--max-prime-limit` above `exp(T)` when using `pbss moments` (its default 1e6
caps T at 13.8).

## Scope limits

The families have no Euler product, no functional equation, no prime
realization; at most polylogarithmic pole counts. The theorem above is a
statement about the explicit family. The ordinary-prime unattained case and
the unconditional moment bound remain open.
