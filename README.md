# Perry–Beurling Spectral Sieve

Private research archive and reconstruction of the **Perry–Beurling Spectral Sieve** (also referred to as the spectral diagnostic / P(q) framework) for testing consistency with the Riemann Hypothesis on Beurling generalized prime systems.

**Author:** Nicholas Perry  
**Status:** Independent research notes reconstructed primarily from November–December 2025 technical sessions (and related follow-ups). Full original scripts, high-precision code, and complete derivations may reside in local archives or earlier session logs.

## Overview

A spectral approach combining Beurling’s theory of generalized primes / Beurling zeta functions with a projection-based “sieve-like” diagnostic. The goal is to analyze density perturbations or spectral measures associated with prime systems and test whether their behavior is consistent with all non-trivial zeros lying on the critical line Re(s) = 1/2.

The framework is numerically stable, basis-invariant, and designed as an exploratory diagnostic / classifier rather than a full proof of RH.

## Core Components

### 1. Projection Strength Metric P(q)

- Construct a normalized density perturbation function *q*.
- Project *q* onto a finite-dimensional space of orthonormalized shifted Legendre polynomials.
- **P(q)** measures the projection strength (energy in the low-degree polynomial subspace).
- Low P(q) indicates high-frequency oscillatory behavior orthogonal to low-degree polynomials — characteristic of RH-like systems.
- Empirical threshold ≈ 29.5 for classification.
- Computed value for the Riemann zeta function: **P(q) ≈ 3.92** (lowest among systems tested).

### 2. Theorems (Conditional Framework)

- **Theorem A**: Assuming the Riemann Hypothesis, the degree-*d* projection strength P_d(q_T) decays to 0 as the logarithmic window size T → ∞.
- **Theorem B**: If P_d(q_T) decays sufficiently rapidly, then there are no zeros off the critical line (provides a conditional equivalence / obstruction argument, not an unconditional proof).

### 3. L² Energy Ratio Diagnostic

- Alternative / complementary formulation using a finite-window L² energy ratio of the polynomial projection component.
- Under RH the projection energy decays as O(T^{-2(d+1)}).
- Ground-truth validation performed on Diamond’s constructed Beurling systems.
- Confirms the diagnostic correctly separates RH-consistent from non-consistent generalized prime systems in tested cases.

### 4. Numerical Experiments & Validation

- High-precision Monte Carlo simulations (including runs with MC = 4000 samples) examining baseline proportionality and behavior under controlled spectral defects / perturbations.
- Tests against real primes (up to 10^5): observed reduced defect spread relative to theoretical baseline, consistent with spectral rigidity expected under RH-like statistics.
- Results supported claims of spectral equilibrium for the sieve diagnostic in the context of RH testing.
- Framework shown to be robust for distinguishing systems while remaining computationally intensive for very large windows.

## Limitations (Explicitly Noted in Development)

- Functions as a **classifier / diagnostic**, not a decisive proof of the Riemann Hypothesis.
- Finite windows cannot rigorously exclude the possibility of zeros at extremely high heights.
- Global / continuous character makes it less suitable as a practical local primality sieve.
- High computational cost for large degree or large T.
- Converse direction (low energy ⇒ RH) is essentially as hard as RH itself.

## Related Private Repositories

- `perry-spirals` — related spiral geometry / Perry’s Law work
- `wieferich-hunts` — parallel prime-search experiments
- Other space-tech / math archives under the same account

## Notes on Reconstruction

This repository was created to centralize the work after it was confirmed that no dedicated GitHub archive existed. Content is synthesized from detailed technical discussions (primarily late 2025). If original Python scripts, notebooks, exact Legendre projection implementations, or full theorem write-ups are recovered from local storage or session exports, they should be added here.

Further extensions discussed in later notes include PSWF (Prolate Spheroidal Wave Function) basis swaps for the zeta sieve.

---

*Private research. Not for public distribution without explicit permission.*
