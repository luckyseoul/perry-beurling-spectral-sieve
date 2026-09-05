# Absolute Legendre moments: a rigorous arithmetic converse

**Date:** 2026-09-05. **Status:** proved equivalence for exact continuous moments.
**RH:** not proved. The equivalent absolute-growth estimate remains unproved
without assuming RH. This theorem replaces an unsupported anti-cancellation
argument; it does not establish a converse from a normalized energy ratio.

## Statement

For y >= 0 define the exact arithmetic residual

\[
 Q(y)=e^{-y/2}\bigl(\theta(e^y)-e^y\bigr),\qquad
 \theta(x)=\sum_{p\le x}\log p.
\]

Let k >= 0 be one fixed integer and let
\(\varphi_k(u)=\sqrt{2k+1}L_k(2u-1)\). Define the **absolute** moment

\[
 c_k(T)=\int_0^1 Q(Tu)\varphi_k(u)\,du
       =\frac1T\int_0^T Q(y)\varphi_k(y/T)\,dy.
 \tag{1}
\]

**Theorem B-ABS.** For every fixed k, the following are equivalent:

1. RH holds.
2. For every epsilon > 0 there are constants C_epsilon and T_epsilon such that
   \(|c_k(T)|\le C_\epsilon e^{\epsilon T}\) for every T >= T_epsilon.

The second statement means subexponential growth; it is an upper bound for
all sufficiently large T, not for a selected subsequence. No independence,
simplicity, isolated rightmost zero, or non-cancellation hypothesis is needed.

For k >= 2, (1) is unchanged if Q(Tu) is first detrended by the **continuous
orthogonal** degree-one projector P_1. This follows from orthogonality to V_1.
The same assertion does not hold automatically for sampled regression or a
smoothed residual.

## 1. Arithmetic Laplace transform and its poles

Put z = s + 1/2. For Re s > 1/2, absolute convergence and summation under the
integral give

\[
\begin{aligned}
 F(s)&:=\int_0^\infty Q(y)e^{-sy}\,dy\\
 &=\frac1z\sum_p\frac{\log p}{p^z}-\frac1{z-1}\\
 &=\frac{-\zeta'(z)/\zeta(z)-H_{\rm pp}(z)}{z}-\frac1{z-1},
 \qquad
 H_{\rm pp}(z)=\sum_p\sum_{m\ge2}(\log p)p^{-mz}.
\end{aligned}
\tag{2}
\]

Indeed each prime contributes
\(\log p\int_{\log p}^\infty e^{-zy}dy=(\log p)p^{-z}/z\).
The logarithmic derivative identity in (2) follows by differentiating the
absolutely convergent Euler product in Re z > 1.

The prime-power term H_pp is holomorphic in Re z > 1/2. On a compact subset
with Re z >= sigma > 1/2, its absolute sum is bounded by

\[
 \frac1{1-2^{-\sigma}}\sum_{n\ge2}(\log n)n^{-2\sigma}<\infty.
\]

Consequently the right side of (2) extends meromorphically to Re s > 0.
The apparent pole at z = 1 cancels: the singular part is
\(1/[z(z-1)]-1/(z-1)=-1/z\).

If rho is a nontrivial zero with Re rho > 1/2 and multiplicity m_rho, then
at \(s_0=\rho-1/2\) the continued F has a simple pole with residue

\[
 r_\rho=-\frac{m_\rho}{\rho}\ne0.
 \tag{3}
\]

Here s_0 != 0 and rho != 0. No other term in (2) has a pole that can cancel (3).

## 2. The moment transform preserves every off-line pole

Write the polynomial \(\varphi_k(u)=\sum_{j=0}^k a_j u^j\) and let
\(A_j(T)=\int_0^T y^j Q(y)dy\). Then

\[
 H_k(T):=T^{k+1}c_k(T)
       =\sum_{j=0}^k a_j T^{k-j}A_j(T).
 \tag{4}
\]

For Re s > 1/2, Fubini and differentiation of absolutely convergent Laplace
integrals yield, with D = d/ds,

\[
 \mathcal L H_k(s)
 =\sum_{j=0}^k a_j(-D)^{k-j}
       \left[\frac{(-D)^jF(s)}s\right].
 \tag{5}
\]

On expanding the product derivatives, the coefficient of the highest
F derivative is

\[
 \frac{(-1)^k\sum_{j=0}^k a_j}{s}F^{(k)}(s)
 =\frac{(-1)^k\varphi_k(1)}sF^{(k)}(s).
 \tag{6}
\]

Every remaining term contains F^(r) with r < k, multiplied by a function
holomorphic near s_0. Since \(\varphi_k(1)=\sqrt{2k+1}\ne0\), (3) and (6)
show that (5) has a pole of order k+1 at s_0 with leading Laurent coefficient

\[
 \boxed{\frac{k!\varphi_k(1)r_\rho}{s_0}
 =-\frac{m_\rho k!\varphi_k(1)}{\rho(\rho-1/2)}\ne0.}
 \tag{7}
\]

The lower derivatives have pole order at most k and therefore cannot cancel it.
This conclusion concerns the actual arithmetic transform, including all zeros
and prime-power terms. It does not isolate a zero by deleting other terms.

## 3. Subexponential absolute moments force RH

Suppose condition 2 holds. Multiplication by T^(k+1) preserves
subexponential growth, so the Laplace transform of H_k is holomorphic throughout
Re s > 0. Differentiation under the integral is justified uniformly on compact
subsets by choosing epsilon smaller than their minimum real part. Local
integrability near T = 0 follows directly from (1).

This holomorphic transform equals (5) on Re s > 1/2. The meromorphic identity
theorem extends that identity throughout the connected half-plane Re s > 0.
The nonzero pole (7) is impossible. Thus zeta has no zero with Re rho > 1/2.
The functional equation reflects a nontrivial zero on the left of the line to
one on the right, so all nontrivial zeros lie on Re rho = 1/2. This proves RH.

## 4. Converse implication

Under RH, the classical estimate
\(\theta(x)-x=O(\sqrt{x}\log^2x)\) gives
\(Q(y)=O(1+y^2)\). Since phi_k is bounded on [0,1], (1) gives
\(|c_k(T)|=O_k(1+T^2)\), which is subexponential. This proves the equivalence. □

The RH error estimate also follows from the explicit formula used in
[CONTINUOUS_THEOREM_A.md](CONTINUOUS_THEOREM_A.md), which proves the stronger
\(c_k(T)=O_k(T^{-1})\) for k >= 1 after accounting for the theta prime-square
constant. Its primary source is Akbary–Ng–Shahabi,
[equation (1.5)](https://www.cs.uleth.ca/~akbary/Akbary-Ng-Shahabi.pdf).
The Euler product, meromorphic continuation, and functional equation used above
are classical properties; see
[DLMF, Riemann zeta function](https://dlmf.nist.gov/25.2) and
[reflection formula](https://dlmf.nist.gov/25.4).

## What this closes, and what remains

B-ABS is an arithmetic converse with a complete pole argument. Its equivalent
absolute-growth condition has not been established unconditionally. It does
not prove the old normalized B-RES claim.

For a residual q_T, a normalized measurement gives only

\[
 |\langle q_T,\varphi_k\rangle|
 \le\|q_T\|_2\sqrt{R_d(q_T)}\qquad(k\le d).
\]

Polynomial decay of R_d can coexist with exponential growth of the absolute
moment if the denominator grows exponentially. Discarding that denominator
loses the information used in (2)–(7). Conversely, assuming subexponential
absolute moments is substantive; the theorem identifies it as an RH-equivalent
condition rather than supplying an independent proof of it.

This is a derived criterion using classical complex analysis. No claim of
historical novelty is made. Finite samples, successful model tests, and finitely
many values of c_k cannot establish the all-T growth condition.
