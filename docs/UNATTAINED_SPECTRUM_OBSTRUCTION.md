# An unattained exponential spectrum can pass every polynomial rate test

**Date:** 2026-09-05
**Status:** Proved analytic spectrum example. This is a real continuous
exponential series with a meromorphic Laplace transform. It is not an
ordinary-prime example, a positive generalized-prime construction, or a zeta
function with an asserted functional equation.

The [attained-boundary converse](ATTAINED_BOUNDARY_CONVERSE.md) leaves open
the case in which off-critical real parts approach a supremum at unbounded
heights. Its energy estimate then has a vanishing prefactor. The example below
shows that this prefactor, discrete poles, and the local zero-count upper
bound used in that proof do not by themselves close the remaining case.
In fact the normalized ratios in this example decay faster than every inverse
power of the window length.

## Statement

For integers \(n\ge8\), set
\[
A=\frac14,\qquad a_n=A-\frac1n,\qquad
\gamma_n=e^n,\qquad \lambda_n=a_n+i\gamma_n,
\]
and define, for \(y\ge0\),
\[
Q(y)=\sum_{n\ge8}e^{-n}e^{a_n y}\cos(\gamma_n y).
\tag{1}
\]
Thus \(1/8\le a_n<A\), \(a_n\uparrow A\), and the amplitudes are exactly
\(1/\gamma_n\), with \(\sum_n1/\gamma_n<\infty\).

For \(T>0\), let
\[
V(T)=\int_0^TQ(y)^2\,dy,\qquad
c_k(T)=\int_0^1Q(Tu)\varphi_k(u)\,du,
\]
\[
\varphi_k(u)=\sqrt{2k+1}L_k(2u-1),\qquad
q_T=(I-P_1)Q(T\,\cdot),\qquad
R_d(q_T)=\frac{\|P_dq_T\|_2^2}{\|q_T\|_2^2}.
\tag{2}
\]
Here \(L_k\) is the Legendre polynomial, and \(P_d\) is the continuous
orthogonal projector onto degree at most \(d\) in \(L^2([0,1],du)\).

**Theorem.** The function \(Q\) is real, continuous, and of exponential
order at most \(A\). Its Laplace transform continues meromorphically to the
whole complex plane, with exactly the simple poles
\(\lambda_n,\bar\lambda_n\), each of residue \(e^{-n}/2\). In particular
the supremum of the pole real parts is positive and unattained. Moreover,
\[
V(T)=o(e^{2AT}),\qquad
\frac{V(T)}{e^{2aT}}\longrightarrow\infty
\quad\text{for every fixed }a<A.
\tag{3}
\]
Nevertheless, for every fixed integer \(d\ge2\), the detrended denominator
is positive for all sufficiently large \(T\), and
\[
\boxed{
R_d(q_T)=O_d\!\left(T e^{-4(\sqrt2-1)\sqrt T}\right)
       =O_{d,N}(T^{-N})\quad\text{for every fixed }N>0.
}
\tag{4}
\]
The estimates hold for all sufficiently large real \(T\), with no phase
subsequence and no growing projection degree.

## 1. Regularity and meromorphic continuation

On every bounded interval in \(y\ge0\), the series (1) converges absolutely
and uniformly, since its terms have absolute values at most \(e^{-n}e^{Ay}\).
Thus \(Q\) is continuous, and
\[
|Q(y)|\le e^{Ay}\sum_{n\ge8}e^{-n}.
\tag{5}
\]
No derivative of the whole series is required or asserted. The subsequent
integrations by parts apply to individual exponential terms, and the stated
summable bounds justify summing their integrals.

For \(\Re s>A\), dominated convergence gives
\[
F(s):=\int_0^\infty Q(y)e^{-sy}\,dy
 =\frac12\sum_{n\ge8}e^{-n}
 \left(\frac1{s-\lambda_n}+\frac1{s-\bar\lambda_n}\right).
\tag{6}
\]
On a fixed compact set of \(s\)-values, for all sufficiently large \(n\)
the denominators in (6) have modulus at least \(e^n/2\). The corresponding
terms are \(O(e^{-2n})\) uniformly on that set. Therefore the series is
locally normally convergent away from the displayed poles. The poles have
no finite accumulation point, their residues are nonzero, and no two of them
coincide. This proves the claimed meromorphic continuation and exact pole set.

The positive ordinates have gaps greater than one, so each unit ordinate
interval contains at most one pole. In particular they obey the local
\(O(\log(t+2))\) counting upper bound used in the arithmetic Gram estimate.
This does not assert the full Riemann–von Mangoldt counting formula.

## 2. A uniform lower and upper Gram bound

Use signed indices \((n,\varepsilon)\) with \(n\ge8\) and
\(\varepsilon\in\{-1,1\}\), and put
\[
f_{n,\varepsilon}(v)=e^{a_n v+i\varepsilon\gamma_n v}.
\]
For \(T\ge1\), the Gram matrix on \([-T,0]\) has diagonal entries
\[
G_{n\varepsilon,n\varepsilon}
 =\int_{-T}^0e^{2a_n v}\,dv,
\qquad
\frac12<e^{-1/2}\le G_{n\varepsilon,n\varepsilon}
 \le\frac1{2a_n}\le4.
\tag{7}
\]
For distinct signed indices, integration of the exponential gives
\[
|G_{n\varepsilon,m\eta}|
 \le\frac2{|\varepsilon e^n-\eta e^m|}.
\tag{8}
\]
The numerator in the antiderivative has modulus at most two, because
\(a_n+a_m>0\).

Here is an explicit row estimate. By sign symmetry take \(\varepsilon=1\).
Using \(1/(1-e^{-1})<2\) and
\(\sum_{m>n}e^{-m}=e^{-n}/(e-1)<e^{-n}\), split the other indices into
four sets:
\[
\begin{aligned}
\sum_{8\le m<n}\frac1{e^n-e^m}
 &\le2(n-8)e^{-n},\\
\sum_{m>n}\frac1{e^m-e^n}
 &\le2e^{-n},\\
\sum_{8\le m\le n}\frac1{e^n+e^m}
 &\le(n-7)e^{-n},\\
\sum_{m>n}\frac1{e^n+e^m}
 &\le e^{-n}.
\end{aligned}
\tag{9}
\]
The first two sets have the same sign as the row and exclude its diagonal;
the last two have the opposite sign and include \(m=n\). Equations (8)–(9)
give
\[
\sum_{(m,\eta)\ne(n,\varepsilon)}|G_{n\varepsilon,m\eta}|
 \le6n e^{-n}\le48e^{-8}<\frac14.
\tag{10}
\]
The middle inequality uses that \(xe^{-x}\) decreases for \(x\ge1\).

For any finite coefficient vector \(d\), the elementary inequality
\(2|d_i d_j|\le |d_i|^2+|d_j|^2\), followed by (7) and (10), therefore
gives
\[
\frac14\sum_i|d_i|^2
 \le\left\|\sum_i d_i f_i\right\|_{L^2(-T,0)}^2
 \le\frac{17}4\sum_i|d_i|^2.
\tag{11}
\]
These constants are independent of the phases, \(T\ge1\), and the number
of retained terms.

Apply (11) to
\[
Q(T+v)=\sum_{n,\varepsilon}d_{n,\varepsilon}(T)
                                      f_{n,\varepsilon}(v),\qquad
d_{n,\varepsilon}(T)=\frac12e^{-n}e^{a_nT}
                                      e^{i\varepsilon\gamma_nT}.
\]
The series converges uniformly on \([-T,0]\), so the finite-sum inequalities
pass to the limit. Define
\[
S_2(T)=\sum_{n\ge8}e^{-2n-2T/n}.
\]
Since \(\sum_{n,\varepsilon}|d_{n,\varepsilon}(T)|^2
=\tfrac12e^{2AT}S_2(T)\), we obtain
\[
\boxed{
\frac18e^{2AT}S_2(T)\le V(T)\le\frac{17}8e^{2AT}S_2(T)
\qquad(T\ge1).
}
\tag{12}
\]

## 3. Every fixed polynomial moment gains another inverse ordinate

Let \(g\in C^1([0,1])\) be fixed. For \(T\ge1\), integration by parts gives
\[
\int_0^T e^{\lambda_n y}g(y/T)\,dy
 =\frac{e^{\lambda_n T}g(1)-g(0)}{\lambda_n}
  -\frac1{T\lambda_n}\int_0^T e^{\lambda_n y}g'(y/T)\,dy.
\]
Since \(a_n\ge1/8\) and \(|\lambda_n|\ge\gamma_n\), its absolute value
is at most
\[
\left(2\|g\|_\infty+8\|g'\|_\infty\right)
                              \frac{e^{a_nT}}{\gamma_n}.
\tag{13}
\]
Multiplying by \(e^{-n}=1/\gamma_n\), summing, taking real parts, and
dividing by \(T\) gives
\[
\left|\int_0^1Q(Tu)g(u)\,du\right|
 \le\frac{C_g e^{AT}}T S_1(T),\qquad
S_1(T)=\sum_{n\ge8}e^{-2n-T/n}.
\tag{14}
\]
In particular (14) applies to every fixed Legendre polynomial.

## 4. Explicit bounds for the two sums

All estimates in this section hold for \(T\ge64\); this ensures that the
index chosen for the lower bound lies in \(n\ge8\).
The arithmetic–geometric mean inequality yields
\[
2n+T/n\ge2\sqrt{2T},\qquad
2n+2T/n\ge4\sqrt T.
\]
There are at most \(T\) integers \(8\le n\le T\), and
\[
\sum_{n>T}e^{-2n}<\frac{e^{-2T}}{1-e^{-2}}<2e^{-2T}.
\]
Using \(2T\ge2\sqrt{2T}\) and \(2T\ge4\sqrt T\) in this range gives
\[
S_1(T)\le(T+2)e^{-2\sqrt{2T}},\qquad
S_2(T)\le(T+2)e^{-4\sqrt T}.
\tag{15}
\]
For the lower bound choose \(n_T=\lceil\sqrt T\rceil\ge8\). Since
\(n_T\le\sqrt T+1\) and \(T/n_T\le\sqrt T\),
\[
\boxed{S_2(T)\ge e^{-2n_T-2T/n_T}
                    \ge e^{-2-4\sqrt T}.}
\tag{16}
\]
Consequently
\[
\frac{S_1(T)^2}{T S_2(T)}
 \le\frac{e^2(T+2)^2}{T}
                e^{-4(\sqrt2-1)\sqrt T}
 =O\!\left(T e^{-4(\sqrt2-1)\sqrt T}\right).
\tag{17}
\]

## 5. Detrending and the normalized bound

Let \(C_k\) denote the constant in (14) for \(g=\varphi_k\), and set
\(C_{r:s}=\sum_{k=r}^sC_k^2\). Orthogonality gives
\[
\|P_1Q(T\,\cdot)\|_2^2
 \le C_{0:1}\frac{e^{2AT}S_1(T)^2}{T^2},\qquad
\|Q(T\,\cdot)\|_2^2=\frac{V(T)}T.
\]
By (12) and (17), the ratio of the first quantity to the second tends to
zero. Therefore for all sufficiently large \(T\),
\[
\|q_T\|_2^2\ge\frac{V(T)}{2T}
                 \ge\frac{e^{2AT}S_2(T)}{16T}>0.
\tag{18}
\]
For fixed \(d\ge2\), affine detrending leaves \(c_k\) unchanged for
\(2\le k\le d\), so
\[
\|P_dq_T\|_2^2=\sum_{k=2}^d|c_k(T)|^2
 \le C_{2:d}\frac{e^{2AT}S_1(T)^2}{T^2}.
\]
Dividing by (18) and using (17) proves (4). This argument controls the
complete series uniformly in its endpoint phases; it does not select a
dominant term and discard a potentially cancelling remainder.

## 6. What the vanishing energy prefactor does not imply

The upper bounds (12) and (15) prove \(V(T)=o(e^{2AT})\). More explicitly,
for \(T\ge64\),
\[
\frac{e^{-2}}8e^{2AT-4\sqrt T}
 \le V(T)\le\frac{17}8(T+2)e^{2AT-4\sqrt T}.
\tag{19}
\]
For any fixed \(a<A\), choose an index \(n\) with \(a_n>a\). Keeping
just its positive summand in the lower bound (12) gives
\[
\frac{V(T)}{e^{2aT}}
 \ge\frac18e^{-2n}e^{2(a_n-a)T}\longrightarrow\infty.
\tag{20}
\]
In particular no pole provides an exponent \(a_n\) for an estimate
\(V(T)=O(e^{2a_nT})\). Even the quantified prefactor in (19) does not
produce the energy-to-pole comparison needed by the attained-boundary proof.

Here the pole gap and height have the explicit relation
\[
(A-a_n)\log\gamma_n=1.
\tag{21}
\]
The indices near \(\sqrt T\) supplying the energy in (19) have ordinates
on an exponential scale in \(\sqrt T\). A fixed polynomial moment incurs
an additional factor \(1/\gamma_n\); estimates (14)–(17) quantify the
resulting loss. Thus further information relating near-boundary real parts,
heights, and the complete arithmetic moments is needed to exclude this
mechanism for ordinary zeta.

This example supplies no monotonicity statement for
\(e^y+e^{y/2}Q(y)\), no prime counting interpretation, and no Euler product.
It refutes an inference from the stated spectral and energy estimates alone.
The unattained-boundary case for the exact ordinary-prime diagnostic remains
open.

## 7. The same obstruction with the arithmetic residue coefficients

**Corollary.** Let
\[
\rho_n=\frac12+\lambda_n=\frac34-\frac1n+i e^n,
\qquad
Q_*(y)=-2\Re\sum_{n\ge8}\frac{e^{\lambda_n y}}{\rho_n}.
\tag{22}
\]
The function \(Q_*\) has all the conclusions (3)–(4) with its own energy,
moments, and continuous affine detrend. Its meromorphic Laplace transform
has exactly the simple poles \(\lambda_n,\bar\lambda_n\), with respective
residues \(-1/\rho_n,-1/\bar\rho_n\). These are the coefficient formulas
for simple zeros in the [arithmetic residual transform](ABSOLUTE_MOMENT_CONVERSE.md).
Here they are prescribed coefficients of an analytic spectrum, rather than
residues derived from any zeta function.

**Proof.** Since \(0<\Re\rho_n<3/4\) and \(e^n>1\),
\[
e^n\le|\rho_n|\le\sqrt2\,e^n,\qquad
\frac12e^{-2n}\le\frac1{|\rho_n|^2}\le e^{-2n}.
\tag{23}
\]
The series (22) is consequently absolutely and locally uniformly convergent
on \(y\ge0\), with \(|Q_*(y)|\le2e^{Ay}\sum_{n\ge8}e^{-n}\).
As in section 1, its Laplace transform for \(\Re s>A\) is
\[
F_*(s)=-\sum_{n\ge8}
 \left(\frac1{\rho_n(s-\lambda_n)}
       +\frac1{\bar\rho_n(s-\bar\lambda_n)}\right).
\tag{24}
\]
Each tail summand is \(O(e^{-2n})\) uniformly on a fixed compact set, so
(24) continues meromorphically to the whole plane, with exactly the stated
nonzero residues.

The Gram estimate (11) applies to arbitrary complex coefficients and phases.
For (22), the two signed coefficients at the endpoint have squared moduli
\(e^{2a_nT}/|\rho_n|^2\). Therefore (11) and (23) give the explicit bounds
\[
\frac14e^{2AT}S_2(T)
 \le V_*(T):=\int_0^TQ_*(y)^2\,dy
 \le\frac{17}2e^{2AT}S_2(T),\qquad T\ge1.
\tag{25}
\]
Equation (13), with \(|\rho_n|^{-1}\le e^{-n}\), similarly gives
\[
\left|\int_0^1Q_*(Tu)g(u)\,du\right|
 \le\frac{2C_g e^{AT}}T S_1(T).
\tag{26}
\]
Thus the sum bounds in section 4 and the affine-projection argument in
section 5 apply with changed constants and exactly the same rates. The
energy conclusions (3) follow from (25) in the same way as (19)–(20). □

This corollary imposes the arithmetic residue coefficients alone. It asserts
no Euler product, functional equation, positive-prime realization, or
ordinary-prime interpretation for \(Q_*\).
