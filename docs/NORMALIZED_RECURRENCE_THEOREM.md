# Unconditional recurrence of small normalized projection ratios

**Date:** 2026-09-05
**Status:** Theorem for the exact continuous arithmetic residual.
**RH:** Open. This theorem proves neither RH nor a full limit or decay rate.

## Statement

Let
\[
Q(y)=e^{-y/2}\bigl(\theta(e^y)-e^y\bigr),\qquad y\ge0,
\]
and let \(P_d\) be the orthogonal projector in \(L^2([0,1],du)\) onto shifted
Legendre polynomials of degrees at most \(d\). For each \(T>0\), set
\[
f_T(u)=Q(Tu),\qquad q_T=(I-P_1)f_T.
\]
The detrend is the exact continuous orthogonal projection with respect to
\(du\). For every fixed integer \(d\ge2\),
\[
\boxed{\liminf_{T\to\infty}R_d(q_T)=0.}
\tag{NR}
\]
No assumption about the location of zeta zeros is needed for (NR).

Both \(\|f_T\|_2\) and \(\|q_T\|_2\) are positive for every \(T>0\).
Indeed, \(Q(y)=-e^{y/2}\) on \(0\le y<\log2\). Its restriction to any nonempty
initial interval is nonzero and is not affine. Consequently \(f_T\notin V_1\),
so none of the displayed ratios has a zero denominator.

## 1. A positive projection floor bounds cumulative energy growth

First consider any real locally square-integrable function \(Q\) on
\([0,\infty)\), and write
\[
U(T)=\int_0^T|Q(y)|\,dy,\qquad
V(T)=\int_0^TQ(y)^2\,dy.
\]
Assume \(V(T)>0\) for all sufficiently large \(T\). Denote the raw ratio by
\[
r_d(T)=R_d\bigl(Q(T\cdot)\bigr).
\]
Since \(\|\varphi_k\|_\infty=\sqrt{2k+1}\),
\[
\begin{aligned}
\left|\langle Q(T\cdot),\varphi_k\rangle\right|
&\le\frac{\sqrt{2k+1}}{T}U(T),\\
r_d(T)&\le\frac{(d+1)^2U(T)^2}{T V(T)}.
\end{aligned}
\tag{1}
\]

Suppose \(\liminf_{T\to\infty}r_d(T)>0\). Choose \(0<\delta\le1\) and
\(T_0>0\) such that \(r_d(T)\ge\delta\) for every \(T\ge T_0\), with
\(V(T_0)>0\). Let
\[
c=\frac{\sqrt\delta}{d+1},\qquad
\alpha=1-\frac{c^2}{4}\in(0,1).
\]
Equation (1) implies \(U(T)/\sqrt{TV(T)}\ge c\). Split the integral for \(U\)
at \(\alpha T\) and use Cauchy–Schwarz on both pieces:
\[
\begin{aligned}
\frac{U(T)}{\sqrt{TV(T)}}
&\le
\sqrt{\frac{\alpha V(\alpha T)}{V(T)}}
+\sqrt{\frac{(1-\alpha)(V(T)-V(\alpha T))}{V(T)}}\\
&\le\sqrt{\frac{\alpha V(\alpha T)}{V(T)}}+\sqrt{1-\alpha}.
\end{aligned}
\tag{2}
\]
As \(\sqrt{1-\alpha}=c/2\), equations (1)–(2) give
\[
V(T)\le K V(\alpha T),\qquad K=\frac{4\alpha}{c^2}>1,
\qquad T\ge T_0.
\tag{3}
\]
This is a fixed-ratio energy bound. For \(T>T_0\), take
\(n=\lceil\log(T/T_0)/\log(1/\alpha)\rceil\). Each application of (3)
has argument at least \(T_0\); the final argument satisfies
\(\alpha^nT\le T_0\). Monotonicity of \(V\) then yields
\[
V(T)\le K^nV(T_0)\le K(T/T_0)^pV(T_0)=O(T^p),\qquad
p=\frac{\log K}{\log(1/\alpha)}<\infty.
\tag{4}
\]
Thus an eventual positive raw projection floor forces at most polynomial
cumulative \(L^2\) energy. This step is purely a statement about functions and
finite-dimensional polynomial projection.

## 2. Polynomial cumulative energy forces RH for the arithmetic residual

Now specialize again to the exact arithmetic \(Q\). If (4) holds, then for
each \(\sigma>0\), integration by parts gives
\[
\int_0^\infty e^{-\sigma y}Q(y)^2\,dy
=\sigma\int_0^\infty e^{-\sigma y}V(y)\,dy<\infty.
\]
Cauchy–Schwarz therefore implies
\[
\int_0^\infty |Q(y)|e^{-\sigma y}\,dy<\infty.
\tag{5}
\]
The same bounds on every closed sub-half-plane imply that the Laplace transform
\[
F(s)=\int_0^\infty Q(y)e^{-sy}\,dy
\]
is holomorphic on \(\operatorname{Re}s>0\).

Initially for \(\operatorname{Re}s>1/2\), integrating the prime-counting step
function and using the Euler-product logarithmic derivative gives
\[
F(s)=
\frac{-\zeta'/\zeta(s+1/2)-A(s)}{s+1/2}
-\frac1{s-1/2},
\tag{6}
\]
where
\[
A(s)=\sum_p\sum_{j\ge2}\frac{\log p}{p^{j(s+1/2)}}.
\]
The series for \(A\) converges absolutely and locally uniformly on
\(\operatorname{Re}s>0\): its terms are bounded by a constant times
\((\log p)p^{-1-2\varepsilon}\) on any region
\(\operatorname{Re}s\ge\varepsilon>0\). Hence it is holomorphic there.

The apparent pole at \(s=1/2\) in (6) cancels, because the pole of \(\zeta\)
at 1 gives residue 1 to the first term. If a nontrivial zero
\(\rho\) with multiplicity \(m_\rho\) had \(\operatorname{Re}\rho>1/2\),
then (6) would have a pole at \(s_0=\rho-1/2\), with nonzero residue
\[
-\frac{m_\rho}{\rho}.
\]
Neither the prime-power term nor the main-term subtraction can cancel this
pole. By the meromorphic identity theorem, it contradicts the holomorphy
established in (5). There are therefore no such zeros. The functional equation
then excludes nontrivial zeros to the left of the critical line as well.

We have proved
\[
\boxed{\liminf_{T\to\infty}r_d(T)>0\quad\Longrightarrow\quad\mathrm{RH}.}
\tag{7}
\]
The same Laplace identity is used in the separate
[`ABSOLUTE_MOMENT_CONVERSE.md`](ABSOLUTE_MOMENT_CONVERSE.md).

## 3. A positive floor survives removal of lower polynomial degrees

Return to \(q_T=(I-P_1)f_T\). Since \(V_1\subseteq V_d\), orthogonality gives
\[
R_d(q_T)=\frac{r_d(T)-r_1(T)}{1-r_1(T)}\le r_d(T).
\tag{8}
\]
The denominator is positive by the non-affinity observation in the statement.
The inequality follows from \(0\le r_1(T)\le r_d(T)\le1\).

If \(\liminf R_d(q_T)>0\), equation (8) gives a positive raw floor, so (7)
forces RH. But RH implies
\[
R_d(q_T)=O_d(T^{-2})\longrightarrow0
\]
by [`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md), contradicting the
assumed positive lower limit. Since the ratios are nonnegative, this proves
\(\liminf R_d(q_T)=0\) unconditionally. □

This is a proof by contradiction, not an assumption that RH is already known:
the hypothetical positive floor first implies RH, and the conditional theorem
then contradicts that same floor.

## Corollary: the positive-floor version of B-RES is equivalent to RH

Fix the same degree, continuous residual, and detrend. Let
\(\mathrm{B\text{-}RES}_{\mathrm{floor}}\) denote the precise assertion
\[
\text{if an off-critical zeta zero exists, then }
\liminf_{T\to\infty}R_d(q_T)>0.
\]
Then
\[
\boxed{\mathrm{B\text{-}RES}_{\mathrm{floor}}\quad\Longleftrightarrow\quad\mathrm{RH}.}
\]
If RH holds, the implication defining \(\mathrm{B\text{-}RES}_{\mathrm{floor}}\)
is vacuously true because there is no off-critical zero. Conversely,
\(\mathrm{B\text{-}RES}_{\mathrm{floor}}\), together with (NR), excludes the
existence of an off-critical zero and therefore implies RH.

Consequently, a persistent positive floor is not a separately available
arithmetic signature that model calculations can certify. Proving that
particular B-RES implication would itself prove RH. This corollary neither
proves nor disproves the implication without resolving RH.

## Exact scope

- (NR) says that arbitrarily large windows have arbitrarily small normalized
  ratios. It does **not** prove \(R_d(q_T)\to0\), an \(O(T^{-2})\) rate, or any
  other decay rate unconditionally.
- The result concerns the exact continuous θ residual and continuous
  degree-one detrend. Fixed sample arrays and their discrete detrends require
  separate resolution estimates; the theorem makes no assertion about their
  large-\(T\) behavior.
- No conclusion about every off-critical mode in isolation is needed. The
  pole argument uses the complete arithmetic residual and its prime-power
  identity.
- Only the positive-floor version of B-RES is covered by the equivalence.
  A converse based on a vanishing, phase-aware lower envelope is a different
  statement and remains a separate problem.
