# ψ − θ correction to the continuous moments, and the R₄ plateau test

**Date:** 2026-11-20
**Scope:** exact identities and floating-point evaluations of them over finite
prime-power prefixes. The decimals below are **not interval certificates**.
Nothing here assumes or proves RH.
**Code:** [`experiments/check_psi_theta_correction.py`](../experiments/check_psi_theta_correction.py) ·
**Data:** [`results/psi_theta_correction_20261120/`](../results/psi_theta_correction_20261120/)

## Summary

- The difference between the ψ residual and the θ residual shifts every
  continuous Legendre moment by an amount that needs only prime powers
  \(p^j\le X\), \(j\ge2\). It is exact and RH-free. At the campaign endpoint
  \(X=5\times10^{10}\), \(c_2[F]-c_2[Q]=-0.201121262828\). Closed-form interval
  integration here matches the campaign's direct summation to 12 digits.
- Its large-window form is \(0.645123988808\ldots/\log X\), and the constant has
  a closed form, \(2\sqrt5\,(C_3-B-\gamma-1)\). The leading term has the wrong
  sign at every reachable endpoint. It is overtaken by the next term until
  \(\log X\approx 533\), because the prime-cube term decays on a log-scale of 6.
- The prime powers do **not** explain the R₄ plateau. The continuous ψ residual
  has 4–5× the R₄ of the θ residual at every window from \(10^3\) to \(10^8\),
  still rising. Within θ, \(T c_2\) has settled, while \(T c_3\) is still
  drifting; the plateau is carried by degree 3 at these windows.

## Definitions

As in [continuous Theorem A](CONTINUOUS_THEOREM_A.md), with \(T=\log X\) and
\(0\le y\le T\):
\[
F(y)=e^{-y/2}\bigl(\psi(e^y)-e^y\bigr),\qquad
Q(y)=e^{-y/2}\bigl(\theta(e^y)-e^y\bigr),\qquad
c_k[R]=\frac1T\int_0^T R(y)\,\varphi_k(y/T)\,dy ,
\]
with \(\varphi_k(u)=\sqrt{2k+1}\,L_k(2u-1)\). Write
\[
g(y)=F(y)-Q(y)=e^{-y/2}\,(\psi-\theta)(e^y)\ \ge 0,
\qquad h=1-g ,
\]
so \(h\) is the function of Theorem A, §2. The correction is
\[
\Delta_k(T)=c_k[F]-c_k[Q]=\langle g_T,\varphi_k\rangle
=-\langle h_T,\varphi_k\rangle\qquad(k\ge1).
\]
Degree-one detrending does not change \(c_k\) for \(k\ge2\).

## 1. Exact evaluation

Since \(\psi-\theta=\sum_{j\ge2}\theta(x^{1/j})\), \(g\) is a step function
times \(e^{-y/2}\) whose jumps sit at prime powers \(p^j\le X\), \(j\ge2\).
Only \(p\le\sqrt X\) enter: 20,639 jumps at \(5\times10^{10}\). On each constant
piece, \(\int e^{-y/2}y^n\,dy\) has a terminating antiderivative, so each
\(\Delta_k\) is a finite sum of elementary terms (mpmath, 40 digits).

| \(X\) | \(T\) | \(\Delta_2\) | \(\Delta_3\) | \(\Delta_4\) |
|---:|---:|---:|---:|---:|
| \(10^{6}\) | 13.816 | −0.256820 | 0.088346 | 0.019091 |
| \(10^{8}\) | 18.421 | −0.244113 | 0.141349 | −0.029261 |
| \(10^{10}\) | 23.026 | −0.212768 | 0.162571 | −0.067427 |
| \(5\times10^{10}\) | 24.635 | **−0.201121262828** | 0.165128 | −0.077683 |
| \(10^{12}\) | 27.631 | −0.180187 | 0.165607 | −0.092788 |

To convert θ-based coefficients to the ψ convention, add \(\Delta_k\).

## 2. Exact three-term structure of \(\Delta_2\)

Because \(\varphi_2(y/T)=\sqrt5\,(1-6y/T+6y^2/T^2)\) and \(\int_0^1\varphi_2=0\),
\[
\boxed{\Delta_2(T)=\frac{\sqrt5}{T}\Bigl[K_0(T)-\frac{6K_1(T)}{T}+\frac{6K_2(T)}{T^2}\Bigr],
\qquad K_m(T)=\int_0^T\bigl(g(y)-1\bigr)\,y^m\,dy .}
\]
This is an identity at every \(T\), not an expansion. At \(X=5\times10^{10}\):

| term | value |
|---|---:|
| \(\sqrt5\,K_0/T\) | +0.017070 |
| \(-6\sqrt5\,K_1/T^2\) | −0.498214 |
| \(+6\sqrt5\,K_2/T^3\) | +0.280023 |
| \(\Delta_2\) | −0.201121262828 |

with \(K_0(T)=0.188062\), \(K_1(T)=22.537\), \(K_2(T)=312.05\).

## 3. The large-window constant in closed form

\(g-1\in L^1([0,\infty))\) unconditionally (classical prime number theorem
error term for \(\theta(\sqrt x)-\sqrt x\); the remaining prime powers are
\(O(e^{-y/6})\)). Hence \(K_m(T)\to K_m(\infty)\) and
\(T\Delta_2(T)\to\sqrt5\,K_0(\infty)\).

**Closed form.** For \(\varepsilon>0\) and \(s=\tfrac12+\varepsilon\),
\[
\int_0^\infty\bigl(g(y)-1\bigr)e^{-\varepsilon y}dy
=\frac1s\sum_{j\ge2}\sum_p\frac{\log p}{p^{js}}-\frac1\varepsilon .
\]
Isolate the \(j=2\) pole with \(\sum_p\log p\,p^{-w}=-\zeta'/\zeta(w)-\sum_{i\ge2}\sum_p\log p\,p^{-iw}\)
and \(-\zeta'/\zeta(w)=1/(w-1)-\gamma+O(w-1)\), then let \(\varepsilon\to0^+\)
(dominated convergence):
\[
\boxed{K_0(\infty)=2\,(C_3-B-\gamma-1)=0.288508218578187\ldots}
\]
\[
B=\sum_p\frac{\log p}{p(p-1)}=0.755366610831688\ldots,\qquad
C_3=\sum_p\frac{\log p}{p^{3/2}-p}=2.476836385022314\ldots
\]
Both sums are evaluated through \(\sum_p\log p\,p^{-w}=-\sum_j\mu(j)\,\zeta'(jw)/\zeta(jw)\).
So \(\sqrt5\,K_0(\infty)=0.645123988808\ldots\), the constant in
\(\Delta_2\sim0.645123988808\ldots/\log X\).

**Independent check.** Cauchy's formula for the derivatives of the same
transform at \(\varepsilon=0\) (a circle of radius 0.05, inside the radius
\(1/6\) set by the prime cubes) reproduces \(K_0(\infty)\) to all printed
digits and gives \(K_1(\infty)=25.6100\), \(K_2(\infty)=409.659\).

**Why the leading term fails.** With these limits, the bracket in §2 changes
sign at \(T\approx16.5\) and \(T\approx516\) (\(X\approx10^{224}\)). The
leading term exceeds the \(1/T^2\) term only for \(T>6K_1/K_0\approx533\).
The cause is the prime-cube term \(e^{-y/2}\theta(x^{1/3})\sim e^{-y/6}\): its
log-moments have scale \(6^{m+1}m!\). It also leaves \(K_0(T)\) short of its
limit at the endpoint, \(0.188\) versus \(0.289\). That deficit, \(-0.100\),
matches the tail \(6X^{-1/6}=0.099\). Finite windows therefore need the exact
correction of §1.

## 4. The R₄ plateau test: θ versus ψ

The continuous, degree-one-detrended ratio
\(R_4=\sum_{k=2}^4c_k^2\big/(N-c_0^2-c_1^2)\) was evaluated for both residuals
by exact integration between jumps, using the local stable form of
[`EXACT_ARITHMETIC_MOMENTS.md`](EXACT_ARITHMETIC_MOMENTS.md) in NumPy
`longdouble`. The θ path reproduces `pbss.arithmetic_moments` at \(X=10^5\)
(maximum coefficient difference \(1.0\times10^{-17}\); identical \(R_4\)).

| \(X\) | \(T\) | \(R_4\), θ | \(T^2R_4\), θ | \(R_4\), ψ | \(T^2R_4\), ψ |
|---:|---:|---:|---:|---:|---:|
| \(10^{3}\) | 6.908 | 0.0295 | 1.41 | 0.2528 | 12.07 |
| \(10^{4}\) | 9.210 | 0.0739 | 6.27 | 0.3557 | 30.17 |
| \(10^{5}\) | 11.513 | 0.1027 | 13.61 | 0.4221 | 55.95 |
| \(10^{6}\) | 13.816 | 0.1121 | 21.39 | 0.4632 | 88.42 |
| \(10^{7}\) | 16.118 | 0.1070 | 27.81 | 0.4862 | 126.30 |
| \(10^{8}\) | 18.421 | 0.0984 | 33.39 | 0.5008 | 169.92 |

**The prime powers are not the cause.** Adding them back makes the ratio
worse at every window. The θ residual is the cleaner object here: in \(F\) the
mean level rises from about −1 at the start of the window to 0, while in \(Q\)
the missing prime squares (\(\theta(x)\approx x-\sqrt x\)) hold it near −1
throughout. Affine detrending removes a constant offset, but not a step.

**Neither residual has reached the \(O(T^{-2})\) regime of Theorem A by
\(10^8\):** \(T^2R_4\) still grows for both.

**Within θ, the drift sits in degree 3:**

| \(X\) | \(T c_2\) | \(T c_3\) | \(T c_4\) | \(R_2\) |
|---:|---:|---:|---:|---:|
| \(10^{3}\) | 0.056 | 0.242 | −0.099 | 0.0013 |
| \(10^{4}\) | 0.462 | 0.204 | −0.285 | 0.0469 |
| \(10^{5}\) | 0.807 | 0.034 | −0.319 | 0.0887 |
| \(10^{6}\) | 1.024 | −0.237 | −0.289 | 0.0989 |
| \(10^{7}\) | 1.104 | −0.531 | −0.177 | 0.0852 |
| \(10^{8}\) | 1.053 | −0.837 | −0.054 | 0.0602 |

\(Tc_2\) has stayed between 1.02 and 1.10 since \(10^6\), and \(R_2\) falls from
\(10^6\) to \(10^8\) by a factor 0.61, close to the \(T^{-2}\) prediction 0.56.
\(Tc_3\) is still moving. Higher-degree coefficients pick up boundary terms
weighted by \(\varphi_k^{(m)}(0)\), which grow roughly like \(k^{2m}\), so a
later settling for degree 3 is plausible. That is an interpretation, not a
derivation.

**Next test.** Extend the continuous θ moments to the campaign endpoint
\(5\times10^{10}\) on SOULKILLER and track \(Tc_3\). If it settles, the plateau
is pre-asymptotic, degree by degree. If it keeps drifting, a further term of
the residual is involved.

These are continuous values. The sampled grand-campaign figures
(\(R_4\approx0.15\)–\(0.19\) through \(10^{10}\), [`STATUS.md`](STATUS.md))
use a different estimator and are not directly comparable.

## 5. Relation to Theorem A and the November 2026 RH claim

In [continuous Theorem A](CONTINUOUS_THEOREM_A.md), the θ step (§2, bound (7))
needs only \(h\in L^1\cap L^2\). That holds unconditionally, as in §3 above;
the RH-strength decay stated in §2 (6) is more than the proof uses. So RH
enters Theorem A only through (EF) and (MS) for \(F\).

On 2026-11-20 Levent Alpöge posted a [note](https://alpo.ge/remark.pdf)
claiming a proof of RH and of GRH for Dirichlet L-functions. The only evidence
offered is a 14-part Lean archive whose top-level theorem is Mathlib's own
`RiemannHypothesis`. No independent build had been reported at the time of
writing. If the claim is confirmed, Theorem A holds unconditionally for the
same continuous θ residual. It would not change any finite-window result
above. The archive works with Mathlib's `vonMangoldt` (ψ-type sums), so
comparisons with it need the \(\Delta_k\) of §1.

## Reproduce

```bash
PYTHONPATH=src python3 experiments/check_psi_theta_correction.py \
    --r4-x 1e3,1e4,1e5,1e6,1e7,1e8
```

Requires `mpmath`. Parts 1–3 take about 3 minutes; the \(10^8\) window of §4
takes about 5 minutes on 2 cores. The run is serial by design: parts 1–3 touch
only primes up to \(\sqrt X\), and §4 is one memory-bound vectorized pass per
window.

## Non-claims

- Not a proof of RH, and not an unconditional proof of Theorem A.
- The plateau mechanism in §4 is open. The table rules out one explanation and
  locates the drift; it does not derive it.
- The decimals are floating-point evaluations of exact identities, not interval
  certificates. The asymptotic statements in §3 follow from the closed forms;
  the crossover values use the computed limits \(K_1(\infty)\), \(K_2(\infty)\).
