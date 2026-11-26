# The 2026-11 RH proof, its simplifications, and what PBSS takes from them

**Date:** 2026-11-26. **Code:** [`src/pbss/rh_moment_asymptotics.py`](../src/pbss/rh_moment_asymptotics.py),
[`experiments/check_rh_moment_asymptotics.py`](../experiments/check_rh_moment_asymptotics.py),
`pbss asymptotics`. **Tests:** [`tests/test_rh_moment_asymptotics.py`](../tests/test_rh_moment_asymptotics.py) (19).
**Artifact:** [`results/rh_moment_asymptotics_20261126/`](../results/rh_moment_asymptotics_20261126/).

**Status in one line.** RH and GRH for Dirichlet L-functions have a Lean proof (Alpöge, with Claude,
2026-11-20). Several groups have machine-checked it, and a shorter edition exists. Human peer review is
still pending. PBSS proves nothing about RH. Its RH-conditional theorems now rest on that formal proof.
This page adds explicit RH asymptotics for the PBSS moments, which explain the R₄ plateau quantitatively.

## 1. What changed since 2026-11-20

| Item | Date | Content |
|------|------|---------|
| [Alpöge's announcement](https://alpo.ge/remark.pdf) | 11-20 | Lean/Mathlib proof of RH and GRH for Dirichlet L-functions, developed with Claude |
| [Full-archive verification](https://mrubinst.github.io/grh-lean-verification/) (Rubinstein) | 11-21 | All 175,187 modules compile; only `propext`, `Classical.choice`, `Quot.sound`; kernel replay passed; nanoda checked 1,484,729 declarations |
| Exposition specialized to ζ (Rubinstein, `zeta-paper_v2`, 265 pp.) | 11-22 | Human-readable route through the proof; see §1.1 |
| [Compact edition](https://github.com/lsernau/grh-lean) (Sernau et al.) | 11-25 | About 76,000 lines in 794 modules (readable edition about 300,000 lines); theorem statements identical to the original |
| [Compact-edition check](https://mrubinst.github.io/grh-lean-compact-check/) | 11-26 | `lake comparator --paranoid` passed on six kernels |
| Further independent checks | by 11-26 | Sutherland ([full](https://math.mit.edu/~drew/grh/), [short](https://math.mit.edu/~drew/shortgrh/)), [Arledge receipts](https://github.com/Infatoshi/rh-verification-receipts), [Naskręcki dashboard](https://nasqret.github.io/riemann-rev/) |

### 1.1 The simplified route (Rubinstein's exposition)

The exposition states the proof as Chebyshev's argument run on a much better "integer":

1. **Zeros to primes, one-sided.** Define the ramp-smoothed prime-power sum and its remainder
   $$S(x)=\sum_{n\le x}\Lambda(n)\Big(\frac{x}{n}-1\Big),\qquad G(x)=S(x)-x\log x+(1+\gamma)x .$$
   If, for every ε > 0, $G(x)\le x^{1/2+\varepsilon}$ throughout windows $[Y,Y^{3/2}]$ for arbitrarily
   large $Y$, then RH holds. Only an **upper** bound is needed, by Landau's theorem on Mellin transforms
   of one sign or Kaczorowski–Pintz.
2. **Chebyshev's principle for lattices.** Divisibility of a nonzero object bounds primes from above;
   the object becomes a linear map between lattices of large rank (slope inequality).
3. **Smooth exponents.** Derivatives of CM functions (the lemniscate sine in the simplest case) have
   divisibility exponents $\approx w\,(x/p^h-1)_+$ at inert primes, so the ramp appears with no floor
   functions.
4. **Size through the linear term.** The map is a chain of $m$ theta lifts. Its norm is known through
   the term linear in $a$.
5. **Denominators** are controlled by a weighted sum of Galois-cohomology lengths (Theorem G).
6. **Rank.** A family of maps fills the target, so the mean slope is controlled.
7. **Balance.** The end-of-chain error $x/m$ and the accumulated error $m$ balance at $m\approx\sqrt x$.
   **This is where the exponent ½ comes from.**
8. **Inert primes to all primes.** Averaging over translates of a residue-class table with two auxiliary
   quadratic characters recovers $G(x)$ (Estimates B and C).

Steps 2–6 give the **main inequality (M)**: on average over translates, the primes are no more numerous
than expected in the cells of the table, up to errors $x/m+xm/k+m+k/n^2$ (times $x^{C\varepsilon}$).
(M) itself is quoted. The exposition reduces it to thirteen numbered hypotheses, none proved outright
there. Steps 1 and 8 are proved in full from cited theorems, and the lemniscate case is worked in full.

PBSS's object is step 1: an explicit-formula statistic of $\psi$ on logarithmic windows. Steps 2–8
are outside PBSS's scope and nothing here touches them.

## 2. The ramp statistic as a PBSS residual

**Explicit formula** (valid for $x>1$, absolutely convergent):
$$G(x)=\log 2\pi-\sum_\rho\frac{x^\rho}{\rho(\rho-1)}-\sum_{m\ge1}\frac{x^{-2m}}{2m(2m+1)} .$$
Under RH, with $A_2=\sum_\rho|\rho|^{-2}=2+\gamma-\log 4\pi=0.0461914179\ldots$ (multiplicity counted):
$$-A_2\sqrt x-(1-\log 2)\,x^{-2}\ \le\ G(x)-\log 2\pi\ \le\ A_2\sqrt x .$$

**PBSS coordinates.** With $g(y)=e^{-y/2}G(e^y)$ and $F(y)=e^{-y/2}(\psi(e^y)-e^y)$,
$$e^{-y}G(e^y)=1+\gamma+\int_0^y F(u)e^{-u/2}\,du,\qquad F=g'-\tfrac12 g .$$
So the residual in Theorem A is the derivative of the proof's statistic, up to a bounded term. Under RH,
$g$ is almost periodic with absolutely summable coefficients $1/(\rho(\rho-1))$.

**Measured** (exact from prime powers, long double):

| Check | Result |
|-------|--------|
| Published values $G(10^4),G(10^5),G(10^6),G(2\cdot10^6)$ = 1.749, 3.671, 11.040, −12.268 | reproduced (max diff 4.6e-4, i.e. rounding) |
| Explicit formula over 200 zeros vs primes, $x=10^3\ldots10^7$ | within the RH zero-tail bound at every point |
| $\max_{x\le10^8}|r(x)|$ with $r=(G-\log2\pi)/\sqrt x$ | 0.0247 at $x=34{,}186{,}499$, i.e. **0.535 A₂** |
| One-sided windows $[Y,Y^{3/2}]$, $Y=10^2\ldots10^5$ | $\max r\le 0.0184$, $\min r\ge-0.0227$; all inside $\pm A_2$ |

These are finite checks of RH consequences, not evidence for RH.

## 3. Exact RH expansion of the PBSS moments

**Prime-sum form** (exact, unconditional). For a residual with jumps $w_n$ at $\log n$:
$$T\,c_k=\sum_{n\le X}w_n\,W_k(\log n)-B_k(T),\quad W_k(y)=\int_y^T e^{-u/2}\varphi_k(u/T)\,du,\quad
B_k(T)=\int_0^T e^{u/2}\varphi_k(u/T)\,du ,$$
with no $B_k$ term for $\psi-\theta$. This agrees with the interval evaluator (θ, $10^{-15}$), the
PR #2 ψ artifact ($2\cdot10^{-14}$), and PR #2's Δₖ at $5\times10^{10}$ and $10^{12}$ ($3\cdot10^{-16}$).

**Mellin constants.** $\Phi(s)=-\zeta'(s)/(s\zeta(s))-1/(s-1)$ is the Mellin transform of $F$ at
$s=\tfrac12+\varepsilon$ and is analytic for $|s-\tfrac12|<\tfrac12$. Put $\mu_j=(-1)^j\Phi^{(j)}(\tfrac12)$:

| $j$ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| $\mu_j$ | −3.3721834192 | −7.2506114432 | −29.334380757 | −176.35494469 | −1411.3630193 |
| $\mu_j/(-2\log2\pi\,2^j j!)$ | 0.9174 | 0.9863 | 0.9976 | 0.9995 | 0.9999 |

$\mu_0=2-2\zeta'(\tfrac12)/\zeta(\tfrac12)$. The second row shows that the pole of $\Phi$ at $s=0$, which
comes from the $-\log 2\pi\,e^{-y/2}$ term (small $x$), dominates every moment beyond the first.

**Theorem (under RH).** For every $T>0$ and $k\ge0$, $T\,c_k[F]=S_k(T)+O_k(T)$ with
$$S_k(T)=\sum_{j\le k}\frac{\varphi_k^{(j)}(0)}{j!\,T^j}\Big[\mu_j-\int_T^\infty y^j e_F(y)\,dy\Big],\qquad
O_k(T)=-\sum_\rho\frac{e^{i\gamma T}}{\rho}\sum_{i\le k}\frac{(-1)^i\varphi_k^{(i)}(1)}{(i\gamma)^{i+1}T^i},$$
$e_F(y)=-\log(2\pi)e^{-y/2}-\tfrac12e^{-y/2}\log(1-e^{-2y})$, and
$$|O_k(T)|\le C_1\sum_{i\le k}\frac{|\varphi_k^{(i)}(1)|}{(\gamma_1T)^i},\qquad C_1=A_2\Big(1+\frac{1}{4\gamma_1^2}\Big)=0.0462492 .$$
*Proof sketch.* Expand $\varphi_k(y/T)$ in powers of $y$. The integrated explicit formula ($g$ above,
absolutely convergent) lets each zero term be integrated exactly. Lower-limit terms sum to the Abel
values $\mu_j$ of $\int_0^\infty y^jF$. Upper-limit terms give $O_k$. The bound uses $|\rho|\ge|\gamma|\ge\gamma_1$
and $\sum_\rho\gamma^{-2}\le C_1$. For θ, $c_k[Q]=c_k[F]-\Delta_k$ exactly, with Δₖ from PR #2 (no RH).

**Measured** ($X=10^3\ldots10^9$, $k=0\ldots4$, 200 zeros): every $|T c_k-S_k|$ is inside the RH envelope.
After the 200-zero sum the remainder is at most $4.8\times10^{-4}$, at most 3.8% of the rigorous zero-tail
bound. The test suite repeats this with the 50-zero table at $10^4\ldots10^7$.

## 4. The R₄ plateau, explained

Two boundary layers at $y=0$ (small $x$) decide the finite-$T$ behaviour.

1. **Numerator.** $S_k(T)\to\varphi_k(0)\mu_0$, but the correction is about $-2.15\,k(k+1)/T$ relative,
   because $\mu_1/\mu_0=2.15$. The window length needed to reach a given fraction of the limit is:

   | $k$ | 50% of limit | 90% of limit |
   |---|---|---|
   | 2 | $T=20.8$ | $T=124.8$ |
   | 3 | $T=39.8$ | $T=247.6$ |
   | 4 | $T=65.2$ | $T=411.4$ |

   For θ, Δₖ adds the slower $e^{-y/6}$ prime-cube layer (PR #2), with its leading term dominant only near
   $\log X\approx 516$.
2. **Normalizer.** Observed: $(N[F]-A_2)\,T=2.504,\,2.510,\,2.513,\,2.511,\,2.511,\,2.511$ for
   $X=10^3\ldots10^8$. So $N[F]\approx A_2+n_1/T$ with $n_1\approx2.511$, and $N$ reaches $A_2$ only when
   $T\gg n_1/A_2\approx54$.

**Prediction vs measurement.** Feeding $S_k+O_k$ (and $-T\Delta_k$ for θ) and the measured $N_{\rm aff}$ into
$R_4=\sum_{k=2}^4c_k^2/N_{\rm aff}$ reproduces the measured R₄ for ψ and θ at every $X$ from $10^3$ to $10^8$,
to a relative $6.4\times10^{-4}$ (θ at $10^8$: 0.09841 predicted and measured).

**Reading.** The plateau is Theorem A's pre-asymptotic regime, fully accounted for by the explicit formula.
It is neither an anomaly nor evidence about zeros. The $T^{-2}$ regime for R₄ starts only at $\log X$ in
the hundreds, which no computation reaches.

## 5. Theorem A with explicit constants

Under RH, $N_{\rm aff}\to A_2^{(m)}=\sum_\gamma m_\gamma^2/|\rho|^2\ge A_2$ (Cramér), and $T c_k[F]\to\varphi_k(0)\mu_0+O_k^\infty(T)$
with $|O_k^\infty|\le C_1\sqrt{2k+1}$. Hence, with $\sum_{k=2}^d(2k+1)=(d+1)^2-4$:
$$\limsup_{T\to\infty}T^2R_d[F]\le\frac{((d+1)^2-4)(|\mu_0|+C_1)^2}{A_2}.$$
If all zeros are simple, $T^2R_d$ eventually oscillates inside
$\big[((d+1)^2-4)(|\mu_0|\mp C_1)^2/A_2\big]$. For θ, replace $\mu_0$ by $\kappa_\theta=\mu_0-K_0(\infty)$ with
$K_0(\infty)=2(C_3-B-\gamma-1)$ from PR #2 (unconditional), because $T\Delta_k\to\varphi_k(0)K_0(\infty)$.

| $d=4$ | $\kappa$ | centre $21\kappa^2/A_2$ | band (simple zeros) |
|---|---|---|---|
| ψ ($F$) | −3.3721834192 | 5169.9 | [5029.0, 5312.7] |
| θ ($Q$, PBSS's residual) | −3.6606916378 | 6092.3 | [5939.4, 6247.3] |

Measured $T^2R_4$ at $X=10^8$ is 169.9 (ψ) and 33.4 (θ), so both are far from these limits, as §4 predicts.

## 6. What this changes in PBSS status

* **Theorem A** ([`CONTINUOUS_THEOREM_A.md`](CONTINUOUS_THEOREM_A.md)) is conditional on RH. That hypothesis
  now has a machine-checked proof, so Theorem A holds relative to it, pending human review, with the constants of §5.
* **B-ABS** ([`ABSOLUTE_MOMENT_CONVERSE.md`](ABSOLUTE_MOMENT_CONVERSE.md)): RH ⇔ subexponential growth of
  $|c_k(T)|$. Relative to the formal proof, the growth statement holds, and in fact $c_k(T)=O_k(1/T)$ by §3.
* The normalized converse is unchanged as a PBSS-internal implication.
* PBSS remains a diagnostic. It did not contribute to, and does not certify, the proof.

## 7. Falsifiable predictions for the campaign endpoint

Continuous θ moments (exact $\Delta_k$; RH expansion with 200 zeros; ± is the rigorous zero-tail bound):

| $X$ | $c_0$ | $c_1$ | $c_2$ | $c_3$ | $c_4$ |
|---|---|---|---|---|---|
| $5\times10^{10}$ | −1.144944 ± 0.000169 | 0.079548 ± 0.000293 | 0.028037 ± 0.000378 | −0.053033 ± 0.000447 | 0.021208 ± 0.000508 |
| $10^{12}$ | −1.130030 ± 0.000151 | 0.086266 ± 0.000261 | 0.016620 ± 0.000337 | −0.049199 ± 0.000399 | 0.030430 ± 0.000453 |

Exact continuous θ moments from the $5\times10^{10}$ prime table must fall inside these bars. Sampled-array
values (fixed grids, smoothing) are different objects and need not.

## 8. Reproduce

```bash
PYTHONPATH=src python3 experiments/check_rh_moment_asymptotics.py   # ~50 s, writes the artifact
PYTHONPATH=src python3 -m pytest -q tests/test_rh_moment_asymptotics.py
PYTHONPATH=src python3 -m pbss asymptotics --x 1e7 --ramp-y 1e4
```

## 9. Non-claims

* PBSS does not prove RH. Every RH-conditional statement here is labelled as such.
* The formal proof is machine-checked. Its human peer review is pending, and this page relies on the
  independent checks listed in §1, not on a review of its mathematics.
* §5's band and §4's $n_1$ are asymptotic or empirical as stated. $n_1$ has no closed form here.
* Nothing here applies to fixed sampled arrays, arbitrary weights, or smoothed residuals.
