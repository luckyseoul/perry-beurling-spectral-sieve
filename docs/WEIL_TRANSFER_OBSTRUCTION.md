# Weil forms do not supply the PBSS projection lower bound

**Date:** 2026-09-05. **Status:** exact transfer obstruction and conditional
finite-dimensional lemma; no proof of arithmetic H*, B-RES, or RH.

The earlier identification of `W := ||P_d q|| / ||q||` as a "Weil-visible channel"
was a proposed dictionary, now retracted in
[`HSTAR_ANTICANCELLATION.md`](HSTAR_ANTICANCELLATION.md).
The cited Weil results do not establish that dictionary or its required lower
bound. The issue persists even before estimating arithmetic remainders.

## What the primary results actually imply

Bombieri's Theorem 1 states a criterion using **every** nonzero smooth compactly
supported Mellin test function. Positivity of the resulting Weil functional is
equivalent to RH. The same paper establishes positivity on sufficiently short
support unconditionally; restricted positivity therefore does not establish the
global criterion. Its finite-truncation eigenvalue result must be read with its
stated hypothesis of finitely many exceptional zeros.
[Bombieri, *Remarks on Weil's quadratic functional*, I, Theorem 1, p. 191 and abstract](https://www.bdim.eu/item?fmt=pdf&id=RLIN_2000_9_11_3_183_0).

The current technical paper linked by Anthropic uses an off-line pair
`{rho, 1 - conjugate(rho)}` as a block of signature `(1, 1)`. Its Lemma 3.1 gives
the **upper** inertia bound

\[
n_+(A^*JA)\le n_+(J).
\]

Proposition 4.1 applies that bound to the evaluation map of a finite family of
test functions. It requires neither injectivity nor a lower singular-value
bound. It supplies no lower bound for a PBSS polynomial projection. This
observation uses the displayed linear algebra, independently of the paper's
zero-proportion conclusions.
[Anthropic technical paper, Lemma 3.1 and Proposition 4.1, pp. 5–6](https://www-cdn.anthropic.com/95c246936988e43127bc6b2ceb7077c1dad2d68e.pdf).

## Exact two-dimensional obstruction

Write the Hermitian contribution of one off-line pair of multiplicity `m > 0`
in evaluation coordinates as

\[
J=m\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\langle J(x,y),(x,y)\rangle=2m\operatorname{Re}(\bar x y).
\]

Its eigenvalues are `m` and `-m`. The map

\[
A:\mathbb C\longrightarrow\mathbb C^2,\qquad Az=(z,iz)
\]

is injective, both evaluation coordinates are nonzero for `z != 0`, and
`||Az||^2 = 2|z|^2`. Nevertheless,

\[
A^*JA=m(i-i)=0.
\]

Thus an indefinite pair can have a completely vanishing compression while its
evaluation map retains nonzero mass. For the perturbed map
`A_eta z = (z, (i+eta)z)` with real `eta > 0`,

\[
A_\eta^*JA_\eta=2m\eta,
\qquad A_\eta^*A_\eta=2+\eta^2.
\]

The compressed form is positive, but its lower bound relative to evaluation
energy tends to zero as `eta -> 0`. Injectivity and positivity at each parameter
therefore do not supply a uniform coercivity constant.

This is a counterexample to the proposed **linear-algebra transfer**, not an
assertion that a particular zeta zero or arithmetic residual realizes these
coordinates. Any arithmetic restriction excluding them must be proved
separately. The existence of a negative eigendirection of `J` does not imply
that the chosen finite test family intersects that direction.

## A recoverable projection lemma with explicit hypotheses

Let `H` be a Hilbert space, `P` an orthogonal projection, and
`B: C^m -> H` an injective map whose range is a specified signal space. Define
the **positive** Gram matrices

\[
G=B^*B>0,\qquad K=B^*PB\ge0.
\]

Then the following statements are equivalent:

1. `P` is injective on `range(B)`.
2. `K` is positive definite.
3. There is a constant `kappa > 0` such that
   `||Ps|| >= kappa ||s||` for every `s` in `range(B)`.

The optimal squared constant is

\[
\kappa^2=\lambda_{\min}(G^{-1/2}KG^{-1/2}).
\]

Indeed, for `s = Bv`, the squared ratio is the generalized Rayleigh quotient
`(v*Kv)/(v*Gv)`. Compactness of the unit sphere gives a positive minimum precisely
when `K` has trivial kernel. Necessarily `m <= rank(P)`; a fixed projector cannot
be injective on signal spaces of unbounded dimension.

For families depending on `T`, positivity of each `K_T` proves only a constant
`kappa_T`. A uniform lower bound requires the additional estimate
`inf_T lambda_min(G_T^(-1/2) K_T G_T^(-1/2)) > 0`. Neither the inertia inequality
nor an empirical collection of positive eigenvalues proves this estimate.

Arithmetic cancellation requires another hypothesis. If `q = s + r`,
`s` belongs to the signal space, and

\[
\|Pr\|\le\eta\|Ps\|\quad(0\le\eta<1),
\qquad \|r\|\le M\|s\|,
\]

then the reverse triangle inequality and the triangle inequality give

\[
\frac{\|Pq\|}{\|q\|}
\ge\frac{1-\eta}{1+M}\frac{\|Ps\|}{\|s\|}
\ge\frac{1-\eta}{1+M}\kappa.
\]

This is a valid conditional anti-cancellation statement. Its assumptions are
about a positive projection Gram matrix and the actual remainder, not about the
signature of a different indefinite form.

## The missing arithmetic implication

To obtain H* by this route one must establish all of the following for the
specified arithmetic residual and test family:

1. An explicit map from arithmetic off-zero contributions to the Hilbert-space
   signal model, with the correct norms and preprocessing.
2. A quantitative lower bound for its projected positive Gram matrix at the
   required scale. A constant bound is stronger than a scale-dependent one.
3. Control of the projected remainder that prevents cancellation of the total
   channel, as well as control of the denominator.

Merely bounding a component `||P q_off|| / ||q_total||` does not bound
`||P q_total|| / ||q_total||`: another component may cancel its projection.
Consequently the informal words "no cancellation" must be replaced by an
inequality involving the total residual or a remainder estimate such as the
one above.

The exact block counterexample identifies why Weil/BGST inertia cannot supply
these implications by itself. It leaves arithmetic B-RES open; it does not
preclude a different theorem using absolute moments, Mellin poles, or a
different family of tests with separately proved estimates.
