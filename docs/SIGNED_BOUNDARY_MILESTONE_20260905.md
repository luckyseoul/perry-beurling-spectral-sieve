# Successor proof milestone: signed moments and the rightmost boundary

**Date:** 2026-09-05. **Base:** `8f36de49256c347c9e4af86ed54bef5be0d847a3`.
**RH:** open. The unconditional one-sided moment bound and the ordinary-zeta
normalized converse with an unattained boundary remain open.

## Verified progress

| Target | Result | Evidence |
|---|---|---|
| Weaken the sufficient absolute growth estimate | Either one-sided subexponential bound for one fixed signed moment implies RH | [Signed moment proof](SIGNED_MOMENT_OSCILLATION.md) |
| Prevent cancellation from erasing a genuine arithmetic pole | Every off-line zero forces positive and negative excursions at scale exp(aT)/T, with explicit constant | Same proof; Landau lemma included |
| Recover normalized visibility with a controlled denominator | A rightmost off-critical zero forces limsup T R_d>0, without a gap or finite block | [Attained-boundary proof](ATTAINED_BOUNDARY_CONVERSE.md) |
| Determine whether positive Beurling primes suffice for the full normalized converse | No: a discrete system has genuine off-line zeros and R_d=O(exp(-T/10)/T) for every fixed d>=2 | [Constructive counterexample](BEURLING_NORMALIZED_COUNTEREXAMPLE.md) |
| Establish the unconditional RH-equivalent bound | **Not established** | Explicit remaining target below |

The root and two agents derived and checked the load-bearing implications
independently. The signed theorem has a self-contained proof of the integral
Landau lemma. The arithmetic denominator estimate uses the unconditional
explicit formula and the standard local zero count, with its complete
absolutely summable two-index Gram estimate written out. The Beurling proof
includes the prime-floor discretization error and a single-valued primitive
construction of the zeta continuation, so it proves genuine zeros rather than
formal logarithmic-derivative poles.

Primary sources were checked directly: Akbary–Ng–Shahabi p. 2 before equation
(1.5), Kedlaya Lemma 9.4 and Theorem 2.4, Johnston–Yang Theorem 1.4(1.8), and
Hilberdink–Lapidus §1.1. Precise links are in the individual proofs. These are
internal derivation reviews, not external peer review or formal verification.
No historical novelty claim is made.

## Validation scope

The changes to the theorem registry preserve `rh="open"` and
`full_B="open_normalized_converse"`, while separately recording signed
oscillation, the conditional attained-boundary converse, and the remaining
unattained case. The two affected registry tests passed in 1.28 seconds with
two independent workers. This verifies registry behavior, not the analytic
theorems.

No numerical campaign was needed for these existence and limit proofs. In
particular, finite ratios would not prove the one-sided bound or the Beurling
asymptotic. Source hashes, the registry report, and the focused replay command
are recorded in
[`results/signed_boundary_milestone_20260905/`](../results/signed_boundary_milestone_20260905/).

## Remaining target

Fix k=2 and use the exact continuous theta moment c_2(T). It now suffices to
prove just

    for every epsilon>0, c_2(T) <= C_epsilon exp(epsilon T)
    for every sufficiently large T,

without RH. The corresponding lower bound would also suffice. A published
unconditional PNT estimate still allows exponential rate 1/2 for this moment;
improving constants within that error-term shape does not reach the target.

For the normalized approach, the unresolved case is
Theta=sup Re(rho)>1/2 with no zero at Re(rho)=Theta. The energy estimate improves
to V(T)=o(exp(2(Theta-1/2)T)), but its unspecified vanishing factor cannot be
compared with a pole at a strictly smaller exponent. New uniform control of
zeros approaching that supremum, or a different arithmetic argument, is needed.
Supremum attainment is not a consequence of meromorphic continuation on an
unbounded strip.

The existing absolute-arithmetic-moment evaluator work is a separate numerical
milestone. It can expose the signed finite-window quantities directly; it
cannot establish the all-window bound by sampling.
