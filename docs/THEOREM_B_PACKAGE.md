# Theorem B: arithmetic converse and remaining growth estimate

**Date:** 2026-09-05. **RH:** open.
The earlier label “package complete, blocked solely by B-RES” conflated a
conditional implication with an unconditional proof. This document supersedes
that label and the claimed positive-floor support from intact off-critical modes.

## Current results

| Statement | Status | Proof |
|---|---|---|
| Persistent low-degree mass gives a positive ratio (model B₀) | Proved projection statement | [M2–M4](PROOFS_LEMMAS.md) |
| Real off-critical mode has a persistent positive ratio | False for the stated continuous model | [Exact phase obstruction](OFF_CRITICAL_PHASE_OBSTRUCTION.md) |
| RH implies \(R_d=O(T^{-2})\) for the exact continuously detrended arithmetic residual | Proved conditional theorem | [Continuous A](CONTINUOUS_THEOREM_A.md) |
| Arithmetic \(\liminf R_d=0\) for that residual | Proved unconditionally | [Recurrence theorem](NORMALIZED_RECURRENCE_THEOREM.md) |
| An off-critical zero forces \(\liminf R_d>0\) | RH-equivalent; open | Recurrence corollary |
| A fixed absolute moment is subexponential iff RH | Proved equivalence (B-ABS) | [Absolute-moment converse](ABSOLUTE_MOMENT_CONVERSE.md) |
| Either one-sided subexponential absolute bound iff RH | Proved equivalence; bound open | [Signed oscillation](SIGNED_MOMENT_OSCILLATION.md) |
| A rightmost off-line zero forces limsup T R_d>0 | Proved if that maximum is attained; no gap needed | [Attained boundary](ATTAINED_BOUNDARY_CONVERSE.md) |
| A specified normalized decay rate implies RH | Unattained-supremum case open | Scope below |
| The same normalized converse for general positive Beurling primes | False | [Discrete counterexample](BEURLING_NORMALIZED_COUNTEREXAMPLE.md) |
| Unconditional bound needed to deduce RH | Open | B-ABS growth condition |

Machine-readable: `pbss.ab_closure.package_status()` reports
`full_B="open_normalized_converse"`,
`absolute_moment_converse="proved_rh_equivalence"`, and
`normalized_recurrence="proved_unconditional_continuous"`.

## Definitions and B-ABS

Let \(Q(y)=e^{-y/2}(\theta(e^y)-e^y)\). Fix one integer \(k\ge0\), and set
\[
c_k(T)=\int_0^1Q(Tu)\varphi_k(u)\,du.
\]
Then RH is equivalent to the assertion that for every \(\epsilon>0\),
\[
|c_k(T)|=O_\epsilon(e^{\epsilon T})\qquad(T\to\infty).
\]
For \(k\ge2\), continuous orthogonal degree-one detrending leaves this moment
unchanged. The statement concerns all sufficiently large windows, not a
subsequence or a finite collection of measurements.

The proof keeps the arithmetic denominator out of the criterion. The Laplace
transform of \(Q\) has residue \(-m_\rho/\rho\) at every off-critical zero
\(s=\rho-1/2\) in the right half-plane. The transform of
\(T^{k+1}c_k(T)\) retains a nonzero pole of order \(k+1\). Subexponential
growth would make that transform holomorphic there, excluding such zeros.
No isolated rightmost zero, finite truncation, simplicity, or independence
hypothesis is assumed.

**What is missing:** an unconditional proof of the displayed absolute-growth
bound. The equivalence itself does not supply it.
The [signed oscillation theorem](SIGNED_MOMENT_OSCILLATION.md) strengthens this:
either an eventual upper or an eventual lower subexponential bound for the
signed moment suffices. A zero rho=1/2+a+i gamma forces excursions of both signs
at scale exp(aT)/T with an explicit positive constant. The one-sided bounds
remain unproved without RH.

## The normalized positive-floor B-RES formulation

With \(q_T=(I-P_1)Q(T\cdot)\), fixed \(d\ge2\), and the same continuous
inner product, the recurrence theorem gives
\[
\liminf_{T\to\infty}R_d(q_T)=0
\]
unconditionally. Its proof shows that an eventual positive raw projection
floor would force polynomial cumulative energy, hence RH by the arithmetic
Laplace transform. Continuous A then contradicts the assumed detrended floor.

Thus the particular B-RES implication
\[
\text{off-critical zero exists}\ \Longrightarrow\ \liminf R_d(q_T)>0
\]
is itself equivalent to RH. Under RH it is vacuously true; combined with
unconditional recurrence it implies RH. This does not prove the implication.

The old component-only \(H^*\) was insufficient even as a projection argument:
other projected components may cancel \(P_dq_{\rm off}\). A sufficient margin
must also control the projected remainder. See [B_RES_THRESHOLD.md](B_RES_THRESHOLD.md).

## Why the model route did not establish B-RES

For the intact mode \(e^{ATu}\sin(tTu)\), \(A,t>0\), the exact continuous
ratio tends to zero. Its generic upper scale is \(T^{-1}\); phase troughs
have scale \(T^{-3}\) for \(d\ge1\). Hence neither nonvanishing nor a
phase-independent \(T^{-1}\) lower bound follows from off-critical origin.
Finite sampled comparisons against a critical-line mode establish only those
measured comparisons. Forced projection cancellation is a useful control,
but is not needed to make this intact mode's ratio vanish.

The localized Weil proxy also lacks the positive coercivity estimate needed
to turn an indefinite quadratic form into a lower bound for projection energy;
[WEIL_TRANSFER_OBSTRUCTION.md](WEIL_TRANSFER_OBSTRUCTION.md) proves that an
injective restriction alone cannot provide it.

## Scope of an unresolved normalized converse

The [attained-boundary theorem](ATTAINED_BOUNDARY_CONVERSE.md) closes a precise
arithmetic case. If some off-critical zero has maximal real part beta among
all nontrivial zeros, the explicit formula and a summable Gram bound give
V(T)=integral_0^T Q(y)^2 dy=O(exp(2(beta-1/2)T)). Signed moment oscillation
then forces limsup T R_d>0. This excludes O(T^-2) without requiring a spectral
gap or a finite zero block. The case in which the supremum is unattained,
approached only at unbounded heights, remains open.

A positive discrete Beurling system can have actual simple zeta zeros at
11/20 ± i gamma while R_d=O(exp(-T/10)/T) for every fixed d>=2. Its
[construction](BEURLING_NORMALIZED_COUNTEREXAMPLE.md) proves that positivity,
discreteness, and an Euler product alone do not supply the ordinary-zeta
energy control. No Riemann-type functional equation is asserted for it.

The statements “\(R_d=o(1)\)” and “\(R_d=O(T^{-2})\)” are different hypotheses;
the former does not imply the latter. A valid normalized converse must state
one precise rate for one precise residual, degree, measure, and preprocessing.
It must also control the absolute size of the residual or otherwise retain the
arithmetic information lost on normalization:
\[
|c_k(T)|\le\|q_T\|_2\sqrt{R_d(q_T)}.
\]
A decaying ratio alone does not bound an exponentially growing absolute norm.
The recurrence theorem establishes neither a full limit nor a decay rate.

Even conditional A together with a valid converse for its matching rate would
give an equivalence; an unconditional proof of the decay property would still
be required to deduce RH. Fixed sample arrays, injected model signals, and
finite GPU or high-precision checks do not supply that proof.
