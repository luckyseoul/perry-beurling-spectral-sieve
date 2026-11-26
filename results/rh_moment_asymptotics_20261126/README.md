# RH moment asymptotics and ramp statistic (2026-11-26)

`rh_moment_asymptotics.json` is from `experiments/check_rh_moment_asymptotics.py` (defaults: 200 zeros,
X = 10^3..10^9, ramp scan to 10^8); `run.log` is its console output. The explainer is
[`docs/RH_PROOF_SIMPLIFICATIONS.md`](../../docs/RH_PROOF_SIMPLIFICATIONS.md).

| Part | Content |
|------|---------|
| `part_a_ramp` | published G values, explicit formula vs primes, r(x) scan, windows [Y, Y^1.5] |
| `part_b_mellin` | mu_0..mu_8 and their ratio to the s = 0 pole term |
| `part_c_moments` | exact T c_k (psi, theta, T Delta_k) vs S_k + O_k, tail bounds and RH envelopes |
| `part_d_r4_plateau` | predicted vs measured R4 (normalizers from the 2026-11-20 artifact), crossover windows, T^2 R4 limit |
| `predictions` | continuous theta moments at 5e10 and 1e12 with rigorous zero-tail bars |

S_k, O_k and all envelopes assume RH. Nothing here proves RH.
