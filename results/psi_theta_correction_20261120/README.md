# ψ − θ correction run (2026-11-20)

Produced by

```bash
PYTHONPATH=src python3 experiments/check_psi_theta_correction.py --r4-x 1e3,1e4,1e5,1e6,1e7,1e8
```

`psi_theta_correction.json` contains:

| Key | Content |
|---|---|
| `part_a_exact_corrections` | Exact \(\Delta_k=c_k[F]-c_k[Q]\), \(k=1..4\), at each endpoint; \(K_m(T)\) and the three-term split of \(\Delta_2\) |
| `part_b_limits` | \(B\), \(C_3\), closed-form \(K_0(\infty)\) and \(\sqrt5K_0(\infty)\); Cauchy-formula \(K_0,K_1,K_2\) at infinity; sign changes of \(\Delta_2\) |
| `part_c_self_check` | θ path versus `pbss.arithmetic_moments` at \(X=10^5\) |
| `part_c_continuous_R4_theta_vs_psi` | Continuous \(c_0..c_4\), norms, affine-detrended \(R_4\) and \(T^2R_4\) for both residuals |

Floating-point evaluations of exact identities: `certified_interval_result` is
`false`, and no RH claim is made. Explanation and tables:
[`docs/PSI_THETA_CORRECTION.md`](../../docs/PSI_THETA_CORRECTION.md).
