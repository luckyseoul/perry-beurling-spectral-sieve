# B-RES: corrected threshold and arithmetic obstruction

**Date:** 2026-09-05. **RH and normalized B-RES:** open.
**Code:** `pbss.b_res_threshold` · [Theorem B package](THEOREM_B_PACKAGE.md).

## The exact arithmetic result

For the exact continuous residual
\[
Q(y)=e^{-y/2}(\theta(e^y)-e^y),\qquad
q_T=(I-P_1)Q(T\cdot),
\]
with orthogonal projection in \(L^2([0,1],du)\), and every fixed \(d\ge2\),
\[
\boxed{\liminf_{T\to\infty}R_d(q_T)=0\quad\text{unconditionally}.}
\]
The proof is in [NORMALIZED_RECURRENCE_THEOREM.md](NORMALIZED_RECURRENCE_THEOREM.md).
It establishes recurrent small values, not a full limit or a rate.
Consequently the particular assertion
\[
\text{off-critical zeta zero exists}\ \Longrightarrow\ \liminf R_d(q_T)>0
\]
is **equivalent to RH**, including its vacuous direction under RH. It has not
been proved here. A quantified vanishing lower envelope is a different converse
problem and must specify its rate and hypotheses separately.

## The omitted cancellation term in the old H* statement

Write \(q_T^{\rm arith}=q_{\rm off}+q_{\rm rest}\). The former condition
\[
\liminf\frac{\|P_dq_{\rm off}\|}{\|q_T^{\rm arith}\|}>0
\]
does **not** imply a lower bound for \(\|P_dq_T^{\rm arith}\|\). For example,
\(q_{\rm rest}=-P_dq_{\rm off}+h\), where \(h\in V_d^\perp\), cancels the
entire projected component. This is an abstract projection example, not an
assertion about the actual arithmetic remainder.

A valid sufficient condition is the stronger margin
\[
\liminf_{T\to\infty}
\frac{\|P_dq_{\rm off}\|-\|P_dq_{\rm rest}\|}
     {\|q_T^{\rm arith}\|}\ge\varepsilon>0.
\]
The reverse triangle inequality then gives \(\liminf R_d\ge\varepsilon^2\).
This margin is sufficient; it is not the weakest condition and is not necessary.
Its arithmetic establishment is missing.

## Corrected model evidence

For the intact real mode \(e^{ATu}\sin(tTu)\), fixed \(A,t>0\) and fixed
\(d\), the continuous ratio satisfies \(R_d=O(T^{-1})\) uniformly in phase.
Moreover,
\[
\limsup TR_d=\frac{4A(d+1)^2}{A^2+t^2},\qquad
\liminf TR_d=0.
\]
For \(d\ge1\) the minimum scale is \(T^{-3}\), with an explicit positive
constant. See [OFF_CRITICAL_PHASE_OBSTRUCTION.md](OFF_CRITICAL_PHASE_OBSTRUCTION.md).
The former phase-independent expression is an upper envelope, not a lower bound.
Even the intact mode has no persistent positive self-normalized floor.

`off_critical_rd_lower_model` retains its historical function name and returns
finite sampled values; its name does not certify an asymptotic lower bound.
Subtracting the sampled projection still gives a valid finite cancellation
control. Neither that control nor an injected bump establishes an arithmetic
explicit-formula contribution.

The former Weil-visible coercivity claim also lacks a valid transfer from an
indefinite form to projection energy. The precise counterexample and needed
positive-Gram estimates are in [WEIL_TRANSFER_OBSTRUCTION.md](WEIL_TRANSFER_OBSTRUCTION.md).
The corrected campaign interpretation is in [HSTAR_ANTICANCELLATION.md](HSTAR_ANTICANCELLATION.md).

## A converse that is proved

[ABSOLUTE_MOMENT_CONVERSE.md](ABSOLUTE_MOMENT_CONVERSE.md) proves that, for any
one fixed Legendre degree, subexponential growth of the **absolute** continuous
moment for all sufficiently large windows is equivalent to RH. Its pole argument
uses the complete arithmetic residual. The equivalent growth estimate remains
unproved unconditionally. This closes a precise converse, but does not close RH
or justify a normalized-ratio classifier.
