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

Not an unconditional proof of RH. Full A is closed conditionally; Full B reduces to B-RES.
