# Signed absolute moments: oscillation and a one-sided RH criterion

**Date:** 2026-09-05. **Status:** proved for exact continuous arithmetic moments.
**RH:** not proved. Neither of the equivalent one-sided bounds below has been
established unconditionally. This strengthens the absolute-moment criterion;
it does not supply a normalized-ratio converse.

Use the exact residual and moments from
[B-ABS](ABSOLUTE_MOMENT_CONVERSE.md):

\[
 Q(y)=e^{-y/2}(\theta(e^y)-e^y),\qquad
 c_k(T)=\frac1T\int_0^T Q(y)\varphi_k(y/T)\,dy,
 \qquad \varphi_k(u)=\sqrt{2k+1}L_k(2u-1).
 \tag{1}
\]

Here k >= 0 is fixed. When k >= 2, the moment is unchanged by continuous
orthogonal affine detrending.

## 1. Quantitative signed oscillation

**Theorem.** Suppose rho = 1/2 + a + i gamma is a nontrivial zeta zero, where
a > 0, of multiplicity m. Set

\[
 A_{k,\rho}=\frac{m\sqrt{2k+1}}{|\rho|\,|\rho-1/2|}>0.
 \tag{2}
\]

Then, with the limits understood in the extended real numbers,

\[
 \boxed{\limsup_{T\to\infty}T e^{-aT}c_k(T)\ge A_{k,\rho},\qquad
 \liminf_{T\to\infty}T e^{-aT}c_k(T)\le-A_{k,\rho}.}
 \tag{3}
\]

No simplicity, independence, isolated rightmost zero, or finite zero-sum
assumption is used. The constant is a lower bound, not an asserted exact
limsup or liminf. In particular the moments satisfy
\(c_k(T)=\Omega_\pm(e^{aT}/T)\) if such a zero exists.

## 2. Positivity lemma for Laplace transforms

We use the following integral version of Landau's theorem. Its short proof is
included to specify exactly the positivity required. The analogous Dirichlet
series theorem and Taylor-series argument are recorded in
[Kedlaya, Theorem 2.4](https://kskedlaya.org/ant/part-1-3.html).

**Lemma.** Let g be nonnegative and locally integrable on [0,infinity), and
suppose its Laplace integral converges for some real s. If its abscissa of
convergence sigma_c is finite, its Laplace transform cannot have a holomorphic
continuation to a neighborhood of the real point sigma_c.

**Proof.** If such a continuation existed on a disk of radius r centered at
sigma_c, choose a real b > sigma_c sufficiently close to sigma_c. The transform
would be holomorphic on a disk centered at b of radius greater than
b - sigma_c. For every integer n >= 0,

\[
 (-1)^nG^{(n)}(b)=\int_0^\infty t^n g(t)e^{-bt}\,dt\ge0.
\]

Choose h > b - sigma_c within that Taylor disk. Its convergent Taylor series
and Tonelli's theorem imply

\[
 G(b-h)
 =\sum_{n=0}^\infty\frac{h^n}{n!}
     \int_0^\infty t^n g(t)e^{-bt}\,dt
 =\int_0^\infty g(t)e^{-(b-h)t}\,dt<\infty.
\]

This contradicts b - h < sigma_c. Derivative integrability follows by using
an intermediate real convergence point between sigma_c and b. This proves
the lemma. If sigma_c = -infinity there is no finite boundary to consider. □

An eventually nonnegative function can be replaced on a compact initial
interval by zero. This changes its Laplace transform by an entire function
and leaves its convergence abscissa unchanged.

## 3. The two analytic inputs

Put \(H_k(T)=T^{k+1}c_k(T)\) and write M_k for its Laplace transform.
The [B-ABS pole calculation](ABSOLUTE_MOMENT_CONVERSE.md) proves that M_k
continues meromorphically to Re s > 0. At \(s_0=\rho-1/2=a+i\gamma\),
its leading Laurent term is

\[
 M_k(s)=\frac{L_{k,\rho}}{(s-s_0)^{k+1}}
           +O((s-s_0)^{-k}),\qquad
 L_{k,\rho}=-\frac{m k!\sqrt{2k+1}}{\rho(\rho-1/2)}.
 \tag{4}
\]

M_k is holomorphic near **every positive real point**. To check this additional
fact, the underlying transform is

\[
 F(s)=\frac{-\zeta'(s+1/2)/\zeta(s+1/2)-H_{\rm pp}(s+1/2)}{s+1/2}
       -\frac1{s-1/2}.
 \tag{5}
\]

The main pole at s = 1/2 cancels, and H_pp is holomorphic in this half-plane.
For a real z > 1, the Euler product is positive and nonzero. For 0 < z < 1,
the paired alternating series

\[
 \eta(z)=\sum_{n\ge1}\bigl((2n-1)^{-z}-(2n)^{-z}\bigr)>0
\]

and \(\zeta(z)=\eta(z)/(1-2^{1-z})\) show that zeta has no real zero there.
The formula expressing M_k through F and its derivatives involves only
additional factors 1/s, so it preserves holomorphy on the positive real axis.
In particular gamma != 0 in (4).

These transforms start in a right half-plane of absolute convergence:
the elementary estimate theta(x) <= x log x gives
\(Q(y)=O((1+y)e^{y/2})\), and H_k has finite exponential order.

## 4. Proof of the oscillation theorem

Suppose the first inequality in (3) fails. Choose a real constant C such that
0 < C < A_(k,rho) and

\[
 c_k(T)\le C e^{aT}/T
\]

for all sufficiently large T. The function

\[
 g(T)=C T^k e^{aT}-H_k(T)
\]

is then eventually nonnegative. Replace it by zero on a compact initial
interval so that it is nonnegative everywhere. Its transform has the
meromorphic continuation

\[
 G(s)=\frac{C k!}{(s-a)^{k+1}}-M_k(s)+E(s),
 \tag{6}
\]

where E is entire. Its finite upper bound on exponential order ensures that
the convergence abscissa sigma_c is less than +infinity. If sigma_c > a,
(6) is holomorphic near the real point sigma_c, contradicting the positivity
lemma. Thus sigma_c <= a, including the possible value -infinity.

The Laplace integral for G consequently converges absolutely in Re s > a.
Its equality with (6) first holds far to the right and extends by the
meromorphic identity theorem. Positivity yields, for every real sigma > a,

\[
 |G(\sigma+i\gamma)|\le G(\sigma).
 \tag{7}
\]

Multiply by (sigma - a)^(k+1) and let sigma decrease to a. Since M_k is
holomorphic near the real point a, the right side tends to C k!. At
a + i gamma, the pole term in (6) contributed by its first summand is
absent, while (4) gives the limit |L_(k,rho)| on the left. Therefore

\[
 |L_{k,\rho}|\le C k!,\qquad A_{k,\rho}\le C,
\]

contradicting the choice of C. For the second inequality in (3), a failed
liminf bound makes \(C T^k e^{aT}+H_k(T)\) eventually nonnegative with
0 < C < A_(k,rho). The same proof, with +M_k in (6), gives the identical
contradiction. □

This argument does not assume that no zeros lie to the right of rho. If such
poles obstruct the claimed continuation from a convergent Laplace integral,
the supposed eventual one-sided bound has already contradicted the identity
theorem.

## 5. One-sided equivalent criteria

For any one fixed k >= 0, each of the following conditions is equivalent to RH:

1. For every epsilon > 0, eventually
   \(c_k(T)\le C_\epsilon e^{\epsilon T}\).
2. For every epsilon > 0, eventually
   \(c_k(T)\ge-C_\epsilon e^{\epsilon T}\).

The constants can be taken positive. Either condition contradicts the
corresponding half of (3) if a zero with a > 0 exists: take epsilon < a and
observe \(T e^{-(a-\epsilon)T}\to0\). The functional equation then gives RH.
Conversely RH gives the two-sided subexponential bound of B-ABS and therefore
both one-sided conditions.

For k = 2 the concrete open target can therefore be reduced to just

\[
 \forall\epsilon>0\quad
 \frac{\sqrt5}{T}\int_0^T
 e^{-y/2}(\theta(e^y)-e^y)
 \left(6\frac{y^2}{T^2}-6\frac yT+1\right)dy
 \le C_\epsilon e^{\epsilon T}
 \quad(T\ge T_\epsilon),
 \tag{8}
\]

or the reversed bound with a minus sign on its right side. This is a weaker
sufficient estimate than a two-sided bound, but is still RH-equivalent and
unproved without RH. Positivity is used on the hypothetical dominating
function g, not asserted for the arithmetic integrand or Legendre kernel.

## 6. What an unconditional PNT estimate actually supplies

A published explicit estimate is

\[
 |\theta(x)-x|\le0.027x(\log x)^{1.801}
 \exp\!\left[-0.1853\frac{(\log x)^{3/5}}{(\log\log x)^{1/5}}\right]
 \qquad(x\ge23).
 \tag{9}
\]

See [Johnston--Yang, Theorem 1.4, equation (1.8)](https://arxiv.org/pdf/2204.01980).
Set h(T) = T^(3/5)/(log T)^(1/5). Boundedness of phi_k on [0,1] and (9) give

\[
 |c_k(T)|=O_k\!\left(T^{0.801}
             \exp(T/2-0.1853h(T))\right).
 \tag{10}
\]

Indeed the integrand's majorant is
\(y^{1.801}\exp(y/2-0.1853h(y))\); its logarithmic derivative tends to 1/2.
For large y that derivative is at least 1/4, so its integral up to T is at
most a constant times its endpoint value. The fixed initial interval
contributes O_k(1/T), and the 1/T in (1) gives (10).

The exponential rate allowed by this envelope is still 1/2 because h(T) = o(T).
It does not establish (8) for epsilon < 1/2. Sharpening constants within this
PNT error shape does not remove that gap. The signed theorem identifies a
precise weaker target; it supplies no unconditional estimate at that target.

The oscillation argument is an application of classical transform positivity.
No historical novelty is claimed. Numerical moment evaluations can check
finite ranges but cannot establish the eventual all-T conditions above.
