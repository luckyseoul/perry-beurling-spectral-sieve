# Absolute arithmetic moments from finite prime intervals

**Scope:** numerical evaluation of exact continuous integral identities over a
verified finite ordinary-prime prefix. The returned floating-point values and
error indicators are **not interval certificates**. No finite evaluation proves
the all-window subexponential growth estimate in
[B-ABS](ABSOLUTE_MOMENT_CONVERSE.md), and this evaluator does not prove RH.

## Object and finite-prime formula

For \(T>0\), use the actual arithmetic residual
\[
 Q(y)=e^{-y/2}\bigl(\theta(e^y)-e^y\bigr),\quad
 \theta(x)=\sum_{p\le x}\log p,\qquad 0\le y\le T.
\]
There is no interpolation between sampled residuals, smoothing, or clamping
\(e^y\) to 2. On \(0\le y<\log2\), in particular, \(Q(y)=-e^{y/2}\).
The absolute Legendre moments and norm are
\[
 c_k(T)=\frac1T\int_0^T Q(y)\varphi_k(y/T)\,dy,\qquad
 N(T)=\frac1T\int_0^T Q(y)^2\,dy.
\]
Interchanging a **finite** sum with the integral gives
\[
 c_k(T)=\frac1T\left[
 \sum_{p\le e^T}\log p\int_{\log p}^T e^{-y/2}\varphi_k(y/T)\,dy
 -\int_0^T e^{y/2}\varphi_k(y/T)\,dy\right].
 \tag{1}
\]
If \(\varphi_k(y/T)=\sum_{j=0}^k a_{k,j}y^j/T^j\), where
\[
 a_{k,j}=\sqrt{2k+1}(-1)^{k+j}\binom{k}{j}\binom{k+j}{j},
\]
then every integral in (1) has a terminating elementary antiderivative.
For \(s\ne0\),
\[
 F_j(s,y)=e^{sy}\sum_{r=0}^j
 \frac{(-1)^r j!}{(j-r)!}\frac{y^{j-r}}{s^{r+1}},
 \qquad \frac{\partial F_j}{\partial y}=e^{sy}y^j.
 \tag{2}
\]
The independent optional `mpmath` checker evaluates (1) and (2) at arbitrary
precision. This global formula can involve substantial cancellation between
the prime sum and the exponential main term. The library therefore evaluates
the same integrals locally, between consecutive prime jumps.

## Stable local moment identity

Partition \([0,T]\) at the prime logarithms. On an interval \([a,a+h]\),
write \(x=e^a\), \(C=\theta(e^a)\), and \(\delta=C/x-1\). Except at
irrelevant measure-zero endpoints, \(C\) is constant there, and
\[
 Q(a+s)=\sqrt{x}\left[\delta e^{-s/2}-2\sinh(s/2)\right].
 \tag{3}
\]
This form computes the local residual before integration, rather than
subtracting its two large integrated main terms. At a prime boundary, \(x\)
is the integer prime itself, avoiding an unnecessary exp(log(p)) round trip.

For \(v=s/h\), expand the polynomial exactly as
\[
 \varphi_k((a+hv)/T)=\sum_{j=0}^k b_{k,j}v^j.
\]
The implementation generates these coefficients by the three-term Legendre
recurrence in the affine argument \((2a/T-1)+(2h/T)v\). Define
\[
 f_j(h)=\int_0^1 v^j e^{-hv/2}\,dv,
 \qquad
 g_j(h)=2\int_0^1 v^j\sinh(hv/2)\,dv.
\]
Then the interval contribution is exactly
\[
 \frac{\sqrt{x}h}{T}\sum_{j=0}^k b_{k,j}
 \left[\delta f_j(h)-g_j(h)\right].
 \tag{4}
\]
The entire-function expansions
\[
 f_j(h)=\sum_{n\ge0}\frac{(-h/2)^n}{n!(j+n+1)},\qquad
 g_j(h)=2\sum_{\substack{n\ge1\\n\text{ odd}}}
       \frac{(h/2)^n}{n!(j+n+1)}
\]
are stable for these short prime intervals. By Bertrand's postulate and the
initial interval ending at 2, their widths satisfy \(h\le\log2\). The
implementation verifies this bound after verifying the finite prime prefix.

The series are truncated after \(N=40\). In exact arithmetic, with
\(z=h/2\), the absolute tail in \(f_j\) is at most
\[
 \frac{e^z z^{N+1}}{(N+1)!(j+N+2)},
\]
and twice this bound suffices for \(g_j\). This follows directly by
integrating the Taylor remainder for the exponential. The report propagates
these estimates through the absolute polynomial coefficients in (4).
The estimates themselves are evaluated in working precision; they do not
include transcendental-function or arithmetic roundoff and are not certified
interval endpoints.

## True integrated norm and continuous affine detrending

On one prime interval, direct integration gives
\[
 \int_a^{a+h} Q(y)^2\,dy
 =C^2(e^{-a}-e^{-a-h})-2Ch+(e^{a+h}-e^a).
 \tag{5}
\]
For stable evaluation, let
\[
 a_0=1-e^{-h},\qquad b_0=h-1+e^{-h},\qquad
 c_0=2(\sinh h-h).
\]
Equivalently, (5) is
\[
 x\left[a_0\delta^2-2b_0\delta+c_0\right]
 =x\left[a_0\left(\delta-\frac{b_0}{a_0}\right)^2
             +c_0-\frac{b_0^2}{a_0}\right].
 \tag{6}
\]
The final difference is nonnegative by Cauchy-Schwarz for the measure
\(e^{-s}ds\): \(a_0\), \(b_0\), and \(c_0\) are respectively the
integrals of \(1\), \(e^s-1\), and \((e^s-1)^2\) in that measure.
Its small-h leading term is \(h^3/12\). Evaluating \(b_0,c_0\) by their
short power series avoids cancellation of their leading h terms. Summing (6)
and dividing by T yields the actual continuous norm \(N(T)\).

The continuous orthogonal affine projector is
\[
 P_1Q(Tu)=c_0(T)\varphi_0(u)+c_1(T)\varphi_1(u).
\]
It removes precisely the first two coefficients. All \(c_k\) for \(k\ge2\)
are unchanged, and
\[
 N_{\rm aff}(T)=N(T)-c_0(T)^2-c_1(T)^2.
 \tag{7}
\]
Raw and detrended projection energies use the corresponding sums of squared
**continuous** orthonormal coefficients. These energies differ from the
coefficient Euclidean norm in the separately implemented sampled weighted QR
projector.

If subtraction in (7) falls below the working roundoff indicator, the library
returns the affine norm and its ratio as `None`, with an explicit unresolved
precision status. It retains the absolute moments. Higher precision is then
needed to resolve the denominator, rather than replacing a small negative
numerical result by an invented norm.

## API and precision provenance

```python
from pbss.arithmetic_moments import arithmetic_moments

result = arithmetic_moments(T=6.0, degree=4)
report = result.to_dict()
```

`max_prime_limit=1_000_000` is the default explicit finite-prefix budget.
Callers may supply an integer prime table, but its required prefix is compared
with the complete existing sieve. Missing, duplicate, unsorted or composite
entries within the window are rejected. A larger table may supply the needed
prefix. These are ordinary primes; a thinned Beurling table is a different
object and is not accepted as arithmetic theta.

The library keeps NumPy `longdouble` absolute coefficients, energies and norms.
Its report records the actual mantissa size on the executing platform; the
type does not promise the same precision on every machine. Decimal strings
in the JSON-compatible report retain these absolute values. Signs,
\(\log|c_k|\), and \(\log|c_k|/T\) remain available without discarding the
amplitude by division by \(\|Q\|\).

Precision metadata include:

- the verified prime-prefix provenance and actual interval count;
- the arithmetic dtype and mantissa bits;
- series-tail estimates and a modeled cumulative-theta roundoff indicator;
- a cancellation indicator comparing absolute local summands with the final
  coefficient, and degrees falling below the roundoff indicator;
- the status of the affine norm subtraction;
- explicit `certified_interval_result=False` and open B-ABS/RH flags.

The roundoff indicators are **not guaranteed error bounds**. They use a
standard arithmetic-error model and do not certify the platform logarithm or
exponential implementations. The optional high-precision experiment reports
measured discrepancies against a different exact finite-sum formula and an
independent direct integration check. Successful finite comparisons establish
neither an interval certificate nor a bound for all large T.
