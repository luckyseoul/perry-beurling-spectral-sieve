# Unattained-boundary milestone and remaining arithmetic target

**Date:** 2026-09-05. **Base:** `9819f4c72fef1746dcbdf4ef3544572fe9bb2325`.
**RH and the ordinary-prime unattained-boundary converse:** open.

The [preceding proof milestone](SIGNED_BOUNDARY_MILESTONE_20260905.md) proved
the normalized converse when a rightmost off-line zero exists. The next attack
was to extend its energy-to-pole comparison to an unattained supremum using
the established vanishing energy prefactor and discrete meromorphic poles.
The [new explicit example](UNATTAINED_SPECTRUM_OBSTRUCTION.md) proves that
these inputs are insufficient.

For a_n=1/4-1/n and gamma_n=exp(n), n>=8, a uniformly convergent real
exponential series has genuine simple Laplace poles at a_n ± i gamma_n. Its
energy satisfies explicit two-sided bounds at scale

    V(T) = exp(T/2 - 4 sqrt(T)) times a factor between c and C(T+2).

Yet every fixed-degree continuous affine-detrended projection ratio satisfies

    R_d(T) = O_d(T exp(-4(sqrt(2)-1)sqrt(T))) = O_{d,N}(T^(-N))

for every fixed N>0 and all sufficiently large T. The result is uniform in
the endpoint phases. A corollary retains these conclusions with the exact
arithmetic residue coefficients -1/rho_n at the poles rho_n-1/2.

This is an analytic spectrum example. It asserts no ordinary-prime or positive
generalized-prime realization, Euler product, or functional equation. Its
local pole count satisfies the upper bound used in the preceding proof, but
it does not assert the full ordinary-zeta counting asymptotic. Those scope
limits prevent confusing a counterexample to the spectral inference with a
counterexample for ordinary zeta.

## Review and evidence

The root and three agents checked the proof. Independent reviews covered
local uniform convergence, termwise Laplace integration, meromorphic
continuation and nonzero residues, the phase-uniform Gram constants, the two
sum estimates for every real T>=64, preservation under affine detrending, and
the arithmetic-residue variant. The note contains every required estimate.
These are internal mathematical reviews, not external peer review or formal
verification. No novelty claim is made.

This milestone changes proof and status documentation only. Its source hashes
and review dispositions are recorded in
[`results/unattained_milestone_20260905/`](../results/unattained_milestone_20260905/).
The last implementation milestone passed 373 tests; all 75 tested Python
files remain unchanged. No repeated numerical campaign was needed to prove
the explicit infinite-series limit statements.

## Exact remaining implication

For the ordinary-prime residual Q(y)=exp(-y/2)(theta(exp(y))-exp(y)), fix

    c_2(T) = (1/T) integral_0^T Q(y) sqrt(5) L_2(2y/T-1) dy.

It would suffice to prove, without RH,

    for every epsilon>0, c_2(T) <= C_epsilon exp(epsilon T)
    for every sufficiently large T.

The corresponding lower bound would also suffice. The signed oscillation
theorem proves this equivalence, not either unconditional bound. Finite
evaluations from `pbss moments` do not establish it.

For the normalized route, a new implication must use arithmetic information
beyond the spectral assumptions realized by the example, or a stronger
uniform estimate that excludes it. Neither that implication nor the one-sided
growth bound has been established in this milestone.
