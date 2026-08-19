# PBSS as a general mathematics tool

**Not a proof of RH.** The same projection diagnostic is useful for any 1D residual
on a unit interval: measure how much energy sits in a low-degree polynomial subspace.

## Install

```bash
# Prefer if venv available:
# python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pbss version
```

Without a venv (editable path):

```bash
PYTHONPATH=src python3 -m pbss version
```

## Commands

| Command | Purpose |
|---------|---------|
| `pbss project -i q.npy -d 4 --T 20` | Project arbitrary \(q(u)\) → \(R_d\), \(S_d\), coeffs |
| `pbss diagnose --demo` | Built-in HF / critical-line / defective scorecard |
| `pbss diagnose -i q.csv` | Diagnose your residual |
| `pbss sensitivity` | Flat vs Gamma-weight Fisher \(d'\) |
| `pbss sensitivity --confirm-53` | Confirm ≥53% relative gain claim |
| `pbss scorecard --x-max 1e6` | Ordinary vs gapped/thinned Beurling systems |

## Library API (general use)

```python
import numpy as np
from pbss import project, sample_grid, confirm_sensitivity_claim

u = sample_grid(2048)
q = np.sin(80 * np.pi * u)          # any residual on [0,1]
r = project(q, u, degree=4, T=20.0)
print(r.energy_ratio, r.scaled_strength)

rep = confirm_sensitivity_claim()   # Gamma weight discriminability
print(rep["verdict"], rep["noisy"]["relative_gain_percent"])
```

## Gamma-weight sensitivity (≥53%)

Project record: Gamma bump \(w(s)=s^{k-1}e^{-\sigma s}\) improved offline/online
discriminability by **~53%** on a noisy Beurling-like ensemble (idealized run often
larger). Shipped confirmation:

```bash
pbss sensitivity --confirm-53 --json-out /tmp/sens.json
```

See [`MEASURE_SENSITIVITY.md`](MEASURE_SENSITIVITY.md).

## Typical non-RH uses

- Detect **low-mode contamination** in a high-frequency residual  
- Compare **weighted** vs flat energy for matched-filter style detection  
- Scorecard **ordinary vs defective** counting systems (Beurling-style)  
- Teaching / demos of orthonormal projection energy ratios  

## Non-claims

- Does not prove RH or close B-RES.  
- Sensitivity gain is **ensemble-dependent**; 53% is a **lower-bound claim** that the
  shipped noisy ensemble meets (often exceeds).
