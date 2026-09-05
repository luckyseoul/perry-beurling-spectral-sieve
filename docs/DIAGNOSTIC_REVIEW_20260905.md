# Diagnostic review — 2026-09-05

Reviewed and corrected with four agents in `/home/nick/perry-beurling-spectral-sieve`,
on `main` at base `cb7cce320148e72eba22b9b71ee273fd7b29cca0`.
Validation was performed before the diagnostic checkpoint commit, including the
pre-existing H* working-tree tests. Existing H* work and historical campaign artifacts
were preserved. This review improves the finite diagnostic; it does not close RH or B-RES.

## Corrected findings

| Finding | Resulting behavior |
|---|---|
| Raw quadrature moments were treated as orthonormal coefficients on every sampled grid. A constant on five uniform points at degree 4 gave `R_d=1.3767242431640625`. | Weighted QR now gives `R_d=1`; fits satisfy weighted orthogonality and Pythagoras. Coefficients reconstruct the polynomial; energy uses `c^T G c`. |
| Shapes were flattened, invalid grids/weights could enter the norm, and `project` accepted invalid `T`. | Real finite vectors, increasing grids, sufficient numerical rank, and positive finite scaling are required. Ratio-only calculations withstand extreme amplitude/weight scaling; unrepresentable absolute outputs fail clearly. |
| `diagnose --input` used 2048 grid nodes regardless of input length. Conflicting flags could silently ignore data. | Grid length comes from the input, `--u` supports nonuniform grids, and incompatible flags fail. Reports include the method, parameters, grid/sample hashes, and finite JSON values. |
| The public scorecard called the Beurling builder incorrectly, unpacked the residual incorrectly, and still returned success after errors. | It uses supported system specs and a shared prime table, retains per-system errors, and returns a nonzero failure status. A real three-system regression verifies the success path. |
| CLI and experiment used different separation thresholds and omitted resolution/support checks. | Both require the specified contamination to exceed both oscillatory controls by 0.5; unsupported degree or aliasing is inconclusive. This is an instrument check for these controls only. |
| A refined grid could cross the separation threshold without changing the overall pass result; three-sample smoothing changed physical width on refinement. | The experiment assesses both grids, treats disagreement as inconclusive, and disables prime smoothing on both. |
| M3/M5 used a false degree-independent coefficient constant; M6 used a false bound on a weighted inner product. | Integration by parts uses `C_k(W,V)=sqrt(2k+1)[2(k+1)W+V]`. Fixed absolute continuity/variation and positive norm hypotheses are explicit; default M5 floors are checked at finite T. |
| Weight admissibility ignored endpoint and upper-bound conditions; endpoint attribution omitted projected cross terms. | Necessary sampled conditions are checked and distinguished from unverified continuous regularity. Endpoint energies include their cross term; undefined zero ratios are rejected. |

The mathematical correction and a concrete counterexample to the old M6
inequality are recorded in [PROOFS_LEMMAS.md](PROOFS_LEMMAS.md).
The continuous [Legendre orthogonality](https://dlmf.nist.gov/18.3) formula
must be distinguished from a sampled [least-squares fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html).

## Verified results

`python3 -m pytest -n 86 --dist loadscope -q`: **264 passed in 21.44 seconds**.
The full suite includes the pre-existing H* tests. The exact staged diagnostic-only
snapshot separately passed **253 tests in 14.35 seconds**, excluding those
pre-existing uncommitted additions. Python multiprocessing needs
local socket access; the initial sandbox-only run blocked two such tests.

The corrected default experiment uses four CPU workers for four independent
probes, `d=4`, `x_max=100000`, and grids of 4096 and 8191 points:

| Probe | Coarse R_d | Refined R_d | Absolute change |
|---|---:|---:|---:|
| high_frequency | 0.000877916 | 0.000878513 | 5.97e-07 |
| critical_line | 0.000999393 | 0.000999989 | 5.96e-07 |
| prime_residual | 0.177455915 | 0.177810159 | 0.000354 |
| defective | 0.925356170 | 0.925355991 | 1.79e-07 |

Both synthetic separation checks pass. All four refinement changes satisfy the
explicit `0.001 + 0.01 * abs(R_d_fine)` tolerance. The initial run with 3-sample
prime smoothing failed this test; holding smoothing disabled resolves that
specific discrepancy without changing the tolerance.

NumPy/CuPy parity passed on the Tesla V100 for 8192/16384/65536 points and
`d=4/8/12`, including nonuniform grids and custom weights. Maximum absolute
ratio disagreement: **1.8041124150158794e-16**; polynomial coefficients agree
within the replay assertions. All GPU cases report `cupy`, not fallback.
This is a correctness check, not a throughput benchmark.

An independent run on Nuka's AMD Radeon RX 9070 XT (CuPy 13.5.1/HIP) also
passed all three cases with actual GPU backends. Maximum absolute ratio
disagreement was **4.3021142204224816e-16**, and coefficient assertions passed.
All 25 imported source-file hashes match this checkout. Nuka used one GPU
context and one BLAS thread; Orin and Jellyfin were inventoried and available.

The installed CuPy initially mixed CUDA 13 headers with NVRTC 12.8. Setting
`CUDA_PATH=/usr/local/cuda-12.8` in the replay process fixed it; no installation
or global configuration change was needed. See CuPy's
[environment-variable documentation](https://docs.cupy.dev/en/stable/reference/environment.html).

## Replay and provenance

Run from the repository root:

```bash
python3 -m pytest -n 86 --dist loadscope -q
python3 experiments/run_diagnostic.py --output-dir results/diagnostic_review_20260905 --workers 4
python3 results/diagnostic_review_20260905/check_gpu_parity.py
```

The GPU replay is a machine-specific snapshot: it selects the local CUDA 12.8
installation before importing CuPy and requires a real GPU backend. The CPU
experiment is independent of that configuration.

- [Diagnostic JSON](../results/diagnostic_review_20260905/diagnostic_run.json)
- [Readable report](../results/diagnostic_review_20260905/diagnostic_run.txt)
- [GPU evidence](../results/diagnostic_review_20260905/gpu_parity.json)
- [Independent AMD evidence](../results/diagnostic_review_20260905/nuka_gpu_parity.json)
- [Standalone Nuka replay snapshot](../results/diagnostic_review_20260905/nuka_replay.tar.gz)
- [Validation and SHA-256 manifest](../results/diagnostic_review_20260905/validation.json)

The Nuka archive contains the exact source, script, and source hash manifest.
Extract it into a fresh directory on Nuka, then run
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/nick/.venvs/rocm72/bin/python parity.py`
there. The script verifies its source hashes before testing and writes
`nuka_gpu_parity.json`. The reviewed run used
`/tmp/pbss-nuka-parity-20260905-01/` over
`ssh -F /home/nick/.ssh/mesh.config nuka`.

## Remaining limits

A sampled fit and one grid refinement do not certify continuum accuracy or
exclude unseen oscillations, out-of-subspace defects, or cancellation. `S_d`
remains a scaling convention; generic model decay is `O(T^-2)`, and a universal
RH threshold is unavailable. Weighted norm floors are hypotheses, and sampled
admissibility cannot prove regularity between nodes. Low-level arithmetic
residual APIs still trust the provenance and coverage of caller-supplied prime
tables; the public scorecard generates its own complete table.

The scoped changed-file whitespace check passes. The repository-wide check also
sees two pre-existing Markdown trailing spaces in `docs/goals/README.md`; that
unrelated dirty file was preserved.
