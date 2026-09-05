# Continuous Theorem A under RH

**Date:** 2026-09-05
**Status:** Conditional theorem for the exact continuous arithmetic residual.
**Unconditional RH:** Open. RH is a hypothesis throughout this deduction.

This supplies a precise arithmetic implication using a cited mean-square theorem
and explicit estimates. It does not identify a fixed sampled array with a
continuous function, and it does not derive RH from a normalized projection ratio.

## Statement and definitions

Let
\[
\theta(x)=\sum_{p\le x}\log p,\qquad
\psi(x)=\sum_{p^j\le x}\log p,
\]
and, for \(y\ge0\), define
\[
F(y)=e^{-y/2}\bigl(\psi(e^y)-e^y\bigr),\qquad
Q(y)=e^{-y/2}\bigl(\theta(e^y)-e^y\bigr).
\]
Values at jumps do not affect the integrals below. In \(L^2([0,1],du)\), let
\(\varphi_k(u)=\sqrt{2k+1}L_k(2u-1)\), let \(P_d\) be the orthogonal projector
onto \(V_d=\operatorname{span}(\varphi_0,\ldots,\varphi_d)\), and set
\[
q_T(u)=(I-P_1)\bigl[Q(Tu)\bigr].
\]
Thus the degree-one detrend is the **continuous orthogonal projection** with
respect to \(du\).

**Theorem.** Assume RH. For every fixed integer \(d\ge2\),
\[
\boxed{R_d(q_T)=\frac{\|P_dq_T\|_2^2}{\|q_T\|_2^2}
       =O_d(T^{-2})\qquad(T\to\infty).}
\]
More precisely, with distinct positive ordinates \(\gamma\) and their zero
multiplicities \(m_\gamma\),
\[
\|q_T\|_2^2\longrightarrow
V:=2\sum_{\gamma>0}\frac{m_\gamma^2}{1/4+\gamma^2}\in(0,\infty),
\qquad
\langle q_T,\varphi_k\rangle=O_k(T^{-1})\quad(k\ge2).
\]
For \(d=0,1\), the numerator vanishes identically by construction; those degrees
provide no converse information after this detrend.

## Cited analytic input

The source is Amir Akbary, Nathan Ng, and Majid Shahabi,
*Limiting distributions of the classical error terms of prime number theory*,
Quarterly Journal of Mathematics **65** (2014), 743–780,
[DOI 10.1093/qmath/hat059](https://doi.org/10.1093/qmath/hat059).
The equation and result numbers used here are those of
[arXiv:1306.1657v1](https://arxiv.org/abs/1306.1657v1), also available as the
[author-hosted PDF](https://www.cs.uleth.ca/~akbary/Akbary-Ng-Shahabi.pdf).

We use two explicit statements from that source:

1. Equation **(1.5)**, p. 2, gives under RH, for \(y\ge\log2\) and \(G>1\),
   \[
   F(y)=\operatorname{Re}\sum_{0<\gamma\le G}
       \frac{-2m_\gamma e^{i\gamma y}}{1/2+i\gamma}
       +O\!\left(\frac{e^{y/2}(y+\log G)^2}{G}+ye^{-y/2}\right).
   \tag{EF}
   \]
2. Corollary **1.15(i)**, p. 11, specialized to \(L(s)=\zeta(s)\), gives
   \[
   \lim_{Y\to\infty}\frac1Y\int_0^Y F(y)^2\,dy=V.
   \tag{MS}
   \]
   This specialization assumes RH; it does **not** assume simplicity or linear
   independence of the zero ordinates. Multiplicities therefore remain in \(V\).

There is no central-zero contribution since \(\zeta(1/2)\ne0\). The standard
zero-count estimate \(N(H)=O(H\log(H+2))\), counting multiplicity, will also be used
to bound convergent coefficient sums. These are external classical inputs;
the PBSS deductions follow below.

## 1. Uniform control of the polynomial moments of the ψ residual

Fix a polynomial \(g\) on \([0,1]\). For \(T\ge\log2\), use (EF) with the
single cutoff \(G=e^T\) throughout \(\log2\le y\le T\).
Integration by parts gives
\[
\left|\frac1T\int_{\log2}^{T}e^{i\gamma y}g(y/T)\,dy\right|
\le \frac{2\|g\|_\infty+\|g'\|_1}{\gamma T}.
\tag{1}
\]
Consequently the complete retained zero block contributes at most
\[
\frac{2(2\|g\|_\infty+\|g'\|_1)}{T}
\sum_{\gamma>0}\frac{m_\gamma}{\gamma\sqrt{1/4+\gamma^2}}
=O_g(T^{-1}).
\tag{2}
\]
The series converges by the zero-count estimate and partial summation, since
its summands are \(O(m_\gamma/\gamma^2)\). In particular, the constant in (2)
is independent of the growing cutoff. No fixed-block M5 constant is silently
applied to an arbitrary increasing number of modes.

The integrated error in (EF) is bounded by a constant times
\[
\frac{\|g\|_\infty}{T}
\left[e^{-T}\int_{\log2}^{T}e^{y/2}(y+T)^2\,dy
      +\int_{\log2}^{T}ye^{-y/2}\,dy\right]
=O_g(Te^{-T/2}+T^{-1})=O_g(T^{-1}).
\tag{3}
\]
On the fixed initial interval \([0,\log2]\), \(F\) is integrable, so its
contribution after division by \(T\) is also \(O_g(T^{-1})\). Thus
\[
\boxed{\int_0^1 F(Tu)g(u)\,du=O_g(T^{-1}).}
\tag{4}
\]
This argument bounds the full arithmetic moment directly. It does not require
an asserted relative \(L^2\) tail estimate with an unspecified denominator.

## 2. Remove the prime-square constant explicitly

The identity
\[
\psi(x)-\theta(x)=\sum_{j\ge2}\theta(x^{1/j})
\]
contains the main term \(\theta(\sqrt{x})\); it cannot be discarded merely by
dividing an \(O(\sqrt{x}\log x)\) bound by \(\sqrt{x}\).

For completeness, (EF) with \(G=x\), together with partial summation of the
zero count, yields the usual RH estimate
\[
\psi(x)-x=O(\sqrt{x}\log^2x).
\]
It follows that \(\psi(x)=O(x)\), hence \(\theta(x)=O(x)\). The prime-power
identity first gives \(\psi(x)-\theta(x)=O(\sqrt{x}\log x)\), and therefore
\(\theta(x)-x=O(\sqrt{x}\log^2x)\). Isolating the square term now gives
\[
\begin{aligned}
\psi(x)-\theta(x)
 &=\theta(\sqrt{x})+O(x^{1/3}\log x)\\
 &=\sqrt{x}+O(x^{1/4}\log^2x+x^{1/3}\log x).
\end{aligned}
\tag{5}
\]
In logarithmic coordinates, (5) states
\[
Q(y)=F(y)-1+h(y),\qquad
h(y)=O(y^2e^{-y/4}+ye^{-y/6})\quad(y\to\infty).
\tag{6}
\]
Define \(h=Q-F+1\) on the initial bounded interval as well. Then
\(h\in L^1([0,\infty))\cap L^2([0,\infty))\). For fixed \(g\),
\[
\left|\int_0^1 h(Tu)g(u)\,du\right|
\le\frac{\|g\|_\infty\|h\|_1}{T},\qquad
\|h(T\cdot)\|_2^2\le\frac{\|h\|_2^2}{T}.
\tag{7}
\]
The degree-one detrend removes the constant \(-1\) in (6) exactly.

## 3. Projection numerator and nonvanishing denominator

Write \(f_T(u)=F(Tu)\) and \(h_T(u)=h(Tu)\). Then
\[
q_T=(I-P_1)(f_T+h_T).
\]
Equations (4) and (7), with \(g=\varphi_k\), show that every fixed coefficient
of \(f_T+h_T\) is \(O_k(T^{-1})\). Orthogonality therefore gives
\[
\|P_dq_T\|_2^2
=\sum_{k=2}^{d}|\langle f_T+h_T,\varphi_k\rangle|^2
=O_d(T^{-2}),\qquad
\|P_1(f_T+h_T)\|_2^2=O(T^{-2}).
\tag{8}
\]
By (MS), \(\|f_T\|_2^2\to V\). By (7), \(\|h_T\|_2\to0\), so
Cauchy–Schwarz implies
\[
\|f_T+h_T\|_2^2\longrightarrow V.
\]
Subtracting the orthogonal component in (8) yields
\[
\|q_T\|_2^2
=\|f_T+h_T\|_2^2-\|P_1(f_T+h_T)\|_2^2
\longrightarrow V>0.
\tag{9}
\]
Positivity follows because the zeta function has nontrivial zeros and each term
in the mean-square sum is nonnegative. Equations (8)–(9) prove the theorem. □

## Scope and relation to the diagnostic

- The theorem concerns the exact continuous θ residual, without an additional
  smoothing operator, and with continuous degree-one detrending. The optional
  weighted version is not asserted in this document.
- The pointwise convention \(x\mapsto\max(x,2)\) changes only a fixed bounded
  interval in \(y\). Such a change has \(O(T^{-1})\) moments and
  \(O(T^{-1})\) squared norm, and does not alter the proof.
- Fixed sample counts cannot justify a \(T\to\infty\) statement. To transfer this
  result to sampled calculations, the grid and detrend must converge to the
  stated integrals with errors controlled on the required scale. The code's
  unweighted sample least-squares detrend is not this continuous operator on an
  arbitrary grid.
- The negative sign and complex phase in (EF), the prime-square constant, the
  cutoff scale, and the positive denominator are explicit here. These replace
  the former unsupported adaptations of generic ANT references.
- Together with any valid converse for the *same* normalized decay property,
  this would give an RH equivalence. It would not supply an unconditional proof
  of that decay property or of RH.

The separate absolute-moment converse is discussed in
[`ABSOLUTE_MOMENT_CONVERSE.md`](ABSOLUTE_MOMENT_CONVERSE.md).
