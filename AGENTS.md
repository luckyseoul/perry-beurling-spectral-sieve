# PBSS project rules

## Compute (always on)

Follow `use-available-compute` for every campaign, experiment, pytest fan-out, or
subagent that burns CPU/GPU. Home rule: `~/.grok/rules/use-local-compute.md`.

```bash
~/.grok/skills/use-available-compute/scripts/compute-budget.sh
```

GPU if dense fit; else `full_workers` (nproc−2). Serial only with an explicit reason.

## Product surface

CLI: `PYTHONPATH=src python3 -m pbss …` — see `docs/TOOL.md`.

## Non-claim

Not an unconditional proof of RH. Conditional A is proved for the exact continuous
degree-one-detrended theta residual; see `docs/CONTINUOUS_THEOREM_A.md`.
`docs/ABSOLUTE_MOMENT_CONVERSE.md` proves an RH equivalence whose unconditional
growth bound remains open. The normalized converse remains open, and
`docs/NORMALIZED_RECURRENCE_THEOREM.md` proves unconditional liminf zero for the
same continuous residual. Do not identify these theorems with fixed sampled
arrays, arbitrary weights, or an unconditional proof of their RH hypotheses.
