# A discrete Beurling counterexample to the normalized rate converse

**Date:** 2026-09-05
**Status:** Proved construction of a discrete generalized-prime system. Its
Euler-product zeta function has a conjugate pair of zeros to the right of
\(\Re s=1/2\), while every fixed-degree continuously detrended projection
ratio is \(O(e^{-T/10}/T)\). This is a counterexample in the Beurling class;
it does not settle the converse for the ordinary primes or prove anything
about off-critical zeros of the Riemann zeta function. No functional equation
is asserted for the constructed zeta function.

The obstruction is visible directly: a rapidly oscillating term supplies
exponentially more total energy than any fixed polynomial projection can
detect. The absolute moments still detect the off-critical zeros.

## Statement and definitions

Fix any \(\gamma>0\). There is a strictly increasing sequence of real
generalized primes \(1<p_1<p_2<\cdots\to\infty\) with
\[
\theta_g(x):=\sum_{p_j\le x}\log p_j\sim x
\]
such that the absolutely convergent Euler product
\[
\zeta_g(s)=\prod_{j\ge1}(1-p_j^{-s})^{-1},\qquad \Re s>1,
\tag{1}
\]
continues meromorphically to \(\Re s>1/2\), with a simple pole at \(1\)
and simple zeros at
\[
\rho=\frac{11}{20}+i\gamma,\qquad \bar\rho=\frac{11}{20}-i\gamma.
\tag{2}
\]
These are its only pole and zeros in that open half-plane.

Use the same continuous diagnostic as in
[continuous Theorem A](CONTINUOUS_THEOREM_A.md):
\[
Q_g(y)=e^{-y/2}\bigl(\theta_g(e^y)-e^y\bigr),\qquad
\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1),
\]
\[
q_T=(I-P_1)[Q_g(T\,\cdot)],\qquad
R_d(q_T)=\frac{\|P_dq_T\|_{L^2(0,1)}^2}{\|q_T\|_{L^2(0,1)}^2},
\tag{3}
\]
where \(L_k\) is the Legendre polynomial and \(P_d\) projects onto the
polynomials of degree at most \(d\) with respect to \(du\). For every fixed
integer \(d\ge2\),
\[
\boxed{R_d(q_T)=O_{d,\gamma}\!\left(\frac{e^{-T/10}}T\right)
       =O_{d,\gamma}(T^{-2}).}
\tag{4}
\]
The norm in (3) is nonzero for every \(T>0\). All limits hold with
\(\gamma\) and the polynomial degree fixed; no assertion about increasing
degree or a sampled projection is needed.

Here generalized integers are finite formal products of generalized primes,
counted with multiplicity when distinct exponent vectors give the same
numerical product. This is the standard positive real Beurling-prime and
Euler-product convention: see Hilberdink and Lapidus, *Beurling Zeta Functions,
Generalised Primes, and Fractal Membranes*, §1.1, pp. 1–2, especially (1.3)
and (1.3′), [arXiv:math/0410270](https://arxiv.org/abs/math/0410270).
All estimates establishing this particular construction are given below.

## 1. A positive continuous prime density

Set
\[
a=\frac1{10},\qquad b=\frac15=2a,\qquad c=\frac1{20},\qquad
\alpha=\frac12+a=\frac35,\qquad X_0=2^{20},
\]
and retain \(\rho=1/2+c+i\gamma\). Powers of positive real numbers use
their real logarithms. For \(x\ge X_0\), define the real function
\[
H(x)=x+x^\alpha\sin(x^b)-\frac{x^\rho}{\rho}
                                      -\frac{x^{\bar\rho}}{\bar\rho}.
\tag{5}
\]
Its derivative satisfies
\[
\begin{aligned}
|H'(x)-1|
&\le \frac35x^{-2/5}+\frac15x^{-1/5}+2x^{-9/20}\\
&\le \frac35\,2^{-8}+\frac15\,2^{-4}+2\,2^{-9}
 =\frac3{160}.
\end{aligned}
\tag{6}
\]
The bound is uniform in \(\gamma>0\), since differentiating
\(x^\rho/\rho\) cancels its denominator. Therefore
\[
\frac{157}{160}\le H'(x)\le\frac{163}{160}\quad(x\ge X_0).
\tag{7}
\]
Define \(\theta_c(x)=0\) for \(1\le x<X_0\), and for \(x\ge X_0\) set
\[
\theta_c(x)=H(x)-H(X_0),\qquad
\pi_c(x)=\int_{X_0}^x\frac{H'(t)}{\log t}\,dt.
\tag{8}
\]
Also set \(\pi_c(x)=0\) below \(X_0\). The function \(\pi_c\) is continuous,
strictly increasing on \([X_0,\infty)\), and unbounded. Moreover,
\(\theta_c(x)=x+O_\gamma(x^{3/5})\).

## 2. Discretization with a logarithmic error

For every positive integer \(j\), define uniquely
\[
p_j=\pi_c^{-1}(j),\qquad
\pi_g(x)=\#\{j:p_j\le x\}=\lfloor\pi_c(x)\rfloor.
\tag{9}
\]
This gives real primes \(p_j>X_0>1\), with no finite accumulation point.
Let \(E(x)=\pi_g(x)-\pi_c(x)\). Then \(-1<E(x)\le0\) and \(E(X_0)=0\).
Stieltjes integration by parts gives, for \(x\ge X_0\),
\[
\begin{aligned}
\Delta(x):=\theta_g(x)-\theta_c(x)
&=\int_{X_0}^x\log t\,dE(t)\\
&=\log x\,E(x)-\int_{X_0}^x\frac{E(t)}t\,dt.
\end{aligned}
\tag{10}
\]
Consequently
\[
-\log x\le\Delta(x)\le\log(x/X_0),\qquad
\boxed{|\Delta(x)|\le\log x.}
\tag{11}
\]
In particular \(\theta_g(x)\sim x\).

There are finitely many generalized integers below any fixed bound: each
uses only the finitely many primes below that bound, and its total exponent
is bounded because \(p_1>1\). With the multiplicity convention stated above,
the Euler product (1) equals the Dirichlet series over those generalized
integers. Since \(\theta_g(x)=O(x)\),
\(\sum_j(\log p_j)p_j^{-\sigma}<\infty\) for every \(\sigma>1\).
Thus (1) and its logarithm converge absolutely and locally uniformly on
\(\Re s>1\), where \(\zeta_g\) is holomorphic and nonvanishing.

## 3. Meromorphic continuation of the prime transform

Initially for \(\Re s>1\), write
\[
D_g(s)=\sum_j(\log p_j)p_j^{-s},\qquad
D_c(s)=\int_{X_0}^{\infty}x^{-s}H'(x)\,dx.
\]
By (10)–(11) and integration by parts,
\[
D_g(s)-D_c(s)=s\int_{X_0}^{\infty}\Delta(x)x^{-s-1}\,dx.
\tag{12}
\]
The right side is holomorphic on \(\Re s>0\), by locally uniform absolute
convergence, including after differentiation in \(s\).

The three power terms in (5) give
\[
D_c(s)=\frac{X_0^{1-s}}{s-1}
       -\frac{X_0^{\rho-s}}{s-\rho}
       -\frac{X_0^{\bar\rho-s}}{s-\bar\rho}
       +C_\chi(s),
\tag{13}
\]
where the remaining function is entire. Here is a direct verification of
that last assertion. Put \(v_0=X_0^b=16\) and, initially for \(\Re w<0\),
\[
S(w)=\int_{v_0}^{\infty}v^{w-1}\sin v\,dv.
\]
Two integrations by parts give
\[
S(w)=v_0^{w-1}\cos v_0-(w-1)v_0^{w-2}\sin v_0
                         -(w-1)(w-2)S(w-2).
\tag{14}
\]
Repeated use of (14) continues \(S\) to every half-plane \(\Re w<2N\),
using finitely many entire boundary terms and the absolutely convergent
integral for \(S(w-2N)\). These continuations agree on overlaps, so \(S\)
is entire. An integration by parts in the oscillatory term of (5), followed
by \(v=x^b\), gives
\[
C_\chi(s)=-X_0^{\alpha-s}\sin(X_0^b)
                +\frac{s}{b}S\!\left(\frac{\alpha-s}{b}\right).
\tag{15}
\]
This proves the asserted entirety without interpreting a divergent integral
as an ordinary convergent Mellin transform.

Equations (12)–(15) continue \(D_g\) meromorphically to \(\Re s>0\).
Its only poles there are simple poles at \(1,\rho,\bar\rho\), with residues
\(+1,-1,-1\), respectively.

## 4. Actual zeros of the Euler-product continuation

The higher prime powers contribute
\[
H_{\rm pp}(s)=\sum_j\sum_{m\ge2}(\log p_j)p_j^{-ms}
 =\sum_j\frac{(\log p_j)p_j^{-2s}}{1-p_j^{-s}}.
\tag{16}
\]
This is holomorphic on \(\Re s>1/2\). Indeed, on any compact subset choose
\(\sigma>1/2\) below all its real parts; the absolute values of the summands
are bounded by
\[
\frac{(\log p_j)p_j^{-2\sigma}}{1-p_1^{-\sigma}},
\]
a summable sequence. Therefore
\[
L_g(s)=D_g(s)+H_{\rm pp}(s)
\]
is meromorphic on that half-plane and agrees with
\(-\zeta_g'(s)/\zeta_g(s)\) on \(\Re s>1\).

To pass from this logarithmic derivative to a single-valued zeta function,
define
\[
A(s)=L_g(s)-\frac1{s-1}+\frac1{s-\rho}+\frac1{s-\bar\rho}.
\tag{17}
\]
The function \(A\) is holomorphic throughout the simply connected domain
\(\Re s>1/2\). It has the path-independent primitive
\[
G(s)=\int_2^s A(w)\,dw.
\]
Then
\[
Z(s)=C\,\frac{(s-\rho)(s-\bar\rho)}{s-1}\,e^{-G(s)},\qquad
C=\frac{\zeta_g(2)}{(2-\rho)(2-\bar\rho)}>0,
\tag{18}
\]
is meromorphic on \(\Re s>1/2\), satisfies \(-Z'/Z=L_g\), and has
\(Z(2)=\zeta_g(2)\). On the connected half-plane \(\Re s>1\), the quotient
\(Z/\zeta_g\) has zero derivative and equals one at \(2\). Thus \(Z\)
is a meromorphic continuation of the original Euler product.

Since the exponential in (18) never vanishes, it gives exactly the simple
zeros and pole stated in (2). This argument establishes actual zeros of a
single-valued continuation, rather than only formal poles of a logarithmic
derivative.

## 5. The normalized ratios decay exponentially

Write \(y_0=\log X_0\), \(\lambda=c+i\gamma\), and
\[
\chi(y)=e^{ay}\sin(e^{by}),\qquad
z(y)=-2\Re\frac{e^{\lambda y}}{\rho}.
\]
Equations (5), (8), and (11) give the exact decomposition on \([0,\infty)\)
\[
Q_g(y)=\chi(y)+z(y)+h(y),
\tag{19}
\]
where on \(y\ge y_0\),
\[
h(y)=e^{-y/2}\bigl(\Delta(e^y)-H(X_0)\bigr)
     =O_\gamma((1+y)e^{-y/2}).
\]
On the fixed initial interval \([0,y_0]\), the function \(h\) is bounded.
In particular \(h\in L^1(0,\infty)\cap L^2(0,\infty)\).

For any fixed polynomial \(g\) on \([0,1]\), integrate using
\(d\cos(e^{by})=-b e^{by}\sin(e^{by})\,dy\):
\[
\begin{aligned}
\int_0^T\chi(y)g(y/T)\,dy
={}&-\left[\frac{e^{(a-b)y}}b g(y/T)\cos(e^{by})\right]_0^T\\
&+\frac1b\int_0^T e^{(a-b)y}
       \left((a-b)g(y/T)+\frac1Tg'(y/T)\right)\cos(e^{by})\,dy.
\end{aligned}
\tag{20}
\]
Because \(b>a\), the right side is \(O_g(1)\) uniformly for \(T\ge1\).
The corresponding integrals for \(z\) and \(h\) are respectively
\(O_{g,\gamma}(e^{cT})\) and \(O_{g,\gamma}(1)\). Hence, for every fixed
\(k\ge0\),
\[
c_k(T):=\langle Q_g(T\,\cdot),\varphi_k\rangle
       =O_{k,\gamma}(e^{cT}/T).
\tag{21}
\]

The choice \(b=2a\) makes the oscillatory norm particularly transparent:
\[
\begin{aligned}
\|\chi(T\,\cdot)\|_2^2
&=\frac{e^{2aT}-1}{4aT}
  -\frac{\sin(2e^{bT})-\sin 2}{4bT}\\
&\sim\frac{e^{2aT}}{4aT}.
\end{aligned}
\tag{22}
\]
Meanwhile
\[
\|z(T\,\cdot)\|_2=O_\gamma(e^{cT}/\sqrt T),\qquad
\|h(T\,\cdot)\|_2=O_\gamma(T^{-1/2}).
\]
Since \(c<a\), the triangle inequality and its reverse show
\[
\|Q_g(T\,\cdot)\|_2^2\sim\frac{e^{2aT}}{4aT}.
\tag{23}
\]
By orthogonality and (21), removing the first two coefficients changes this
squared norm by at most \(O_\gamma(e^{2cT}/T^2)\). Therefore
\[
\boxed{\|q_T\|_2^2\sim\frac{e^{2aT}}{4aT}},\qquad
\|P_dq_T\|_2^2=\sum_{k=2}^d|c_k(T)|^2
               =O_{d,\gamma}(e^{2cT}/T^2).
\tag{24}
\]
Division proves (4), because \(2(a-c)=1/10\). For any \(T>0\),
\(Q_g(Tu)=-e^{Tu/2}\) on a nonempty initial interval before the first prime.
It cannot agree almost everywhere with an affine function there. Thus its
continuous degree-one detrend has strictly positive norm, as required for (3).

## 6. The absolute moments retain the zero contribution

One more integration by parts in the \(z\) contribution, together with
(20) and \(h\in L^1\), gives for every fixed \(k\ge0\)
\[
c_k(T)=-\frac{2\varphi_k(1)}T
          \Re\frac{e^{\lambda T}}{\rho\lambda}
       +O_{k,\gamma}(e^{cT}/T^2)+O_{k,\gamma}(T^{-1}).
\tag{25}
\]
Choose the positive terms of the sequence
\[
T_n=\frac{2\pi n-\arg(1/(\rho\lambda))}{\gamma}.
\]
Then
\[
|c_k(T_n)|\sim
  \frac{2\varphi_k(1)}{|\rho\lambda|}\frac{e^{cT_n}}{T_n}.
\tag{26}
\]
Every fixed absolute moment has an exponentially growing subsequence. For
\(k\ge2\), it is unchanged by the continuous detrend. Thus the fast decay
of the normalized ratios comes from their denominator and does not imply
subexponential growth of the absolute moments.

## Consequence and scope

This example meets positivity, discreteness, a prime asymptotic
\(\theta_g(x)\sim x\), an absolutely convergent positive Euler-product
Dirichlet series on \(\Re s>1\), and meromorphic continuation with genuine
off-critical zeros. Those properties alone do not make a fixed-degree
\(O(T^{-2})\) normalized ratio exclude such zeros.

The construction has no asserted Riemann-type functional equation or the
remaining structure of the ordinary primes. A converse using that additional
structure requires its own argument. The result here is an analytic
counterexample for generalized primes; it is not a finite numerical
certificate or a counterexample to the ordinary-prime converse.
