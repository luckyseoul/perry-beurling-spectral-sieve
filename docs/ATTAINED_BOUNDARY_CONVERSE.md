# A normalized converse when the rightmost zero boundary is attained

**Date:** 2026-09-05. **Status:** conditional arithmetic converse proved.
**RH:** open. The hypothesis that an off-critical rightmost real part is
attained is not proved and is not imposed on zeta in the unconditional claims.

Let
\[
 Q(y)=e^{-y/2}(\theta(e^y)-e^y),\qquad
 q_T=(I-P_1)Q(T\cdot),\qquad
 V(T)=\int_0^T Q(y)^2\,dy.
\]
All projectors use continuous Lebesgue measure on [0,1]. Fix d>=2.
Here P_d is the orthogonal projector onto polynomials of degree at most d, and
\(R_d(q_T)=\|P_dq_T\|_2^2/\|q_T\|_2^2\).
The denominator is positive for T>0: before the first prime, Q(y)=-exp(y/2)
on a nonempty interval, so Q(Tu) cannot agree almost everywhere with an affine
polynomial.

**Theorem.** Suppose there is a nontrivial zeta zero rho with
\[
 \beta=\operatorname{Re}\rho>1/2,
 \qquad \operatorname{Re}\rho'\le\beta
 \quad\text{for every nontrivial zero rho'.}
 \tag{1}
\]
Then
\[
 \boxed{\limsup_{T\to\infty}T R_d(q_T)>0.}
 \tag{2}
\]
In particular R_d(q_T)=O(T^-2), and even R_d(q_T)=o(T^-1), is impossible
under (1). This does not assume a finite dominant zero block, a gap below the
rightmost real part, simplicity, or linear independence. Infinitely many zeros
may have real part beta, and other zeros may approach it from the left.

## 1. A pole becomes a normalized obstruction if the energy is controlled

First suppose only that rho=1/2+a+i gamma is a zero with a>0 and multiplicity m,
and that for some K>0,
\[
 V(T)\le K e^{2aT}\qquad(T\ge T_0).
 \tag{3}
\]
The [signed-moment theorem](SIGNED_MOMENT_OSCILLATION.md), applied to any fixed
k with 2<=k<=d, gives
\[
 \limsup T e^{-aT}c_k(T)\ge
 A_{k,\rho}:=\frac{m\sqrt{2k+1}}{|\rho|\,|\rho-1/2|}>0.
\]
The same c_k is a coefficient of q_T because affine detrending leaves degrees
k>=2 unchanged. Since ||q_T||_2^2<=V(T)/T,
\[
 T R_d(q_T)\ge\frac{T^2c_k(T)^2}{V(T)}.
\]
On sequences approaching the displayed positive limsup, (3) therefore gives
\[
 \boxed{\limsup T R_d(q_T)\ge A_{k,\rho}^2/K>0.}
 \tag{4}
\]
This bridge uses the complete arithmetic pole and an upper energy bound. It
does not require a lower bound on an individually selected zero's projection
relative to the remaining zeros.

## 2. Uniform energy on unit logarithmic intervals

We prove the energy bound needed in (3) from the zero-location hypothesis
Re rho'<=beta, for any fixed beta>1/2. Put
\[
 E_\beta(y)=e^{-\beta y}(\psi(e^y)-e^y).
\]
Then
\[
 \sup_{Y\ge1}\int_Y^{Y+1}|E_\beta(y)|^2\,dy<\infty.
 \tag{5}
\]

The classical truncated explicit formula, stated unconditionally immediately
before equation (1.5) on p. 2 of
[Akbary–Ng–Shahabi](https://www.cs.uleth.ca/~akbary/Akbary-Ng-Shahabi.pdf), is
\[
 \psi(x)=x-\sum_{|\operatorname{Im}\rho|\le G}\frac{x^\rho}{\rho}
 +O\!\left(\frac{x\log^2(xG)}G+\log x\right),
 \qquad x\ge2,\quad G>1.
 \tag{6}
\]
Zeros in the sum are counted with multiplicity. No RH specialization is made
here. Choose G=exp(2(Y+1)) throughout y in [Y,Y+1]. After multiplication by
exp(-beta y), the error in (6) is
\[
 O\bigl((Y+1)^2e^{-(1+\beta)Y}+(Y+1)e^{-\beta Y}\bigr),
\]
which is uniformly bounded and tends to zero. It suffices to bound the
retained zero sum independently of G and Y.

For a positive-ordinate zero rho=beta_rho+i gamma, write
\[
 b_\rho(Y)=\frac{e^{(\rho-\beta)Y}}\rho,
 \qquad |b_\rho(Y)|\le\frac1{|\rho|}.
\]
With y=Y+v, the positive-ordinate sum is
\(\sum b_\rho(Y)e^{(\rho-\beta)v}\), 0<=v<=1. In its squared norm, a pair
rho,rho' has the Gram factor
\[
 \left|\int_0^1e^{\kappa v+i\Delta v}\,dv\right|
 \le\min\left(1,\frac2{|\Delta|}\right),
 \quad
 \kappa=\beta_\rho+\beta_{\rho'}-2\beta\le0,
 \quad \Delta=\gamma-\gamma'.
 \tag{7}
\]
At Delta=0 the first bound applies. The second follows directly from the
antiderivative; its numerator has modulus at most two. Consequently arbitrarily
close zero ordinates cause no divergence.

Group positive ordinates into unit bins and put
\[
 B_n=\sum_{n\le\gamma<n+1}\frac1{|\rho|},\qquad n=0,1,2,\ldots.
\]
The standard local zero count N(t+1)-N(t)=O(log(t+2)), with multiplicities,
gives B_n=O(log(n+2)/(n+1)); finitely many small bins are harmless. For bins n,m,
(7) is at most C/(1+|n-m|). The complete absolute Gram majorant is therefore
\[
 C\sum_{n,m\ge0}
 \frac{\log(n+2)\log(m+2)}{(n+1)(m+1)(1+|n-m|)}<\infty.
 \tag{8}
\]
For clarity, split the sum into n/2<=m<=2n and its complement. In the first
region the nth row is O(log^3(n+2)/(n+1)^2), using the harmonic sum in |n-m|.
In the region m>2n the row is
\[
 O\!\left(\frac{\log(n+2)}{n+1}
           \sum_{m>2n}\frac{\log(m+2)}{(m+1)^2}\right)
 =O\!\left(\frac{\log^2(n+2)}{(n+1)^2}\right).
\]
The other region is its transpose. All three sums converge. Negative-ordinate
terms are conjugates, so multiplying the positive-block norm bound by four
bounds the complete real sum. This proves (5).

The local zero-count input is stated in
[Kedlaya, Lemma 9.4](https://kskedlaya.org/ant/chap-von-mangoldt.html);
Remark 9.7 there also records the Riemann–von Mangoldt formula with O(log t)
error. No positive minimum spacing between zeros is assumed.

## 3. Transfer to theta and sum the energy

The elementary Chebyshev bound theta(x)=O(x) and the exact prime-power identity
give, for x>=2,
\[
 \psi(x)-\theta(x)
 =\theta(\sqrt x)+\sum_{j\ge3}\theta(x^{1/j})
 =O(\sqrt x+x^{1/3}\log x)=O(\sqrt x).
\]
Thus for beta>1/2,
\[
 e^{-\beta y}(\psi(e^y)-\theta(e^y))
 =O(e^{-(\beta-1/2)y}),
\]
and (5) holds with psi replaced by theta as well.

Set a=beta-1/2>0. On each unit interval [j,j+1],
\[
 \int_j^{j+1}Q(y)^2\,dy
 \le e^{2a(j+1)}
      \int_j^{j+1}|e^{-\beta y}(\theta(e^y)-e^y)|^2\,dy
 \le C_\beta e^{2a(j+1)}.
\]
The bounded initial interval is harmless. Summing the geometric progression up
to ceil(T) gives
\[
 \boxed{V(T)=O_\beta(e^{2(\beta-1/2)T}).}
 \tag{9}
\]
The implied constant may depend on beta and diverge as beta decreases to 1/2.
At beta=1/2 the same unit-interval argument only gives V(T)=O(T), so a positive
exponent must not be silently replaced by zero in the geometric sum.

Under (1), the zero rho lies exactly at the exponent in (9), and (4) proves
the theorem. □

## 4. The remaining normalized case is now precise

Let
\[
 \Theta=\sup\{\operatorname{Re}\rho:\zeta(\rho)=0,\ 0<\operatorname{Re}\rho<1\}.
\]
If the exact continuous diagnostic satisfies R_d=O(T^-2), then either:

1. RH holds; or
2. Theta>1/2 and **no zero has real part Theta**.

In the second case, real parts approaching Theta must occur at unbounded
heights, because a bounded region contains only finitely many zeros. This case
has not been excluded. Supremum attainment cannot be inferred from meromorphic
continuation on an unbounded strip.

The bound (9) remains valid with beta=Theta>1/2 even when Theta is not attained,
but each particular zero then has a strictly smaller exponential rate.
The signed-moment lower bound and total-energy upper bound no longer match;
(4) would acquire an exponentially vanishing factor. The theorem does not
justify discarding that factor or choosing a zero depending on T without new
uniform estimates.

In fact an unattained Theta>1/2 gives the slightly stronger energy conclusion
\[
 \boxed{V(T)=o(e^{2(\Theta-1/2)T}).}
 \tag{10}
\]
To see this, use beta=Theta in section 2. Every fixed coefficient
\(b_\rho(Y)=e^{(\rho-\Theta)Y}/\rho\) tends to zero. Each term in the absolute
Gram bound therefore tends to zero, and the summable majorant (8) allows
dominated convergence even with the changing cutoff G(Y). The explicit-formula
error and prime-power correction also tend to zero. Thus
\[
 \varepsilon_j:=\int_j^{j+1}
 |e^{-\Theta y}(\theta(e^y)-e^y)|^2\,dy\longrightarrow0.
\]
For a=Theta-1/2>0 and N=floor(T), section 3 now yields
\[
 e^{-2aT}V(T)
 \le e^{-2aT}V(1)+e^{2a}\sum_{j=1}^{N}e^{-2a(N-j)}\varepsilon_j
 \longrightarrow0.
\]
The last limit follows by splitting off finitely many j and bounding the
remaining epsilon_j by an arbitrarily small constant against a summable
geometric progression. This proves (10). It does not supply a rate for the
little-o term, so it still cannot compare that term with the smaller
exponential rate of any fixed zero or settle the unattained-boundary case.

This closes the attained-boundary case of the normalized converse. It does
not prove unconditional decay, RH, or the full normalized converse. The
independent [Beurling counterexample](BEURLING_NORMALIZED_COUNTEREXAMPLE.md)
shows why an Euler product, positivity, and poles alone do not supply the
missing norm control; zeta's explicit formula and zero-count bounds are
essential inputs to the arithmetic estimate above.
