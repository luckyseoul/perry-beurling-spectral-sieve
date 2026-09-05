# Fixed-degree visibility of a real off-critical model

**Status: proved model asymptotics, with numerical corroboration.** The result
below concerns one explicitly specified function family. It is not an arithmetic
counterexample, a proof that an off-critical zeta zero exists, or a solution of
B-RES or the Riemann Hypothesis.

It corrects the phase-independent asymptotic previously used to support the
candidate `L_star_local_anticancellation`. That expression is an upper envelope
of the leading term. It is not a uniform lower bound.

## Definitions and result

Fix real constants \(A,t>0\) and an integer \(d\ge0\). In **continuous**
\(L^2([0,1],du)\), set
\[
q_T(u)=e^{ATu}\sin(tTu),\quad
\varphi_k(u)=\sqrt{2k+1}\,P_k(2u-1),\quad
R_d(T)=\frac{\sum_{k=0}^d|\langle q_T,\varphi_k\rangle|^2}
 {\|q_T\|_2^2}.
\]
Let \(\lambda=A+it\), \(r=|\lambda|\), \(\beta=\arg\lambda\),
\(\theta=tT\), and define
\[
 B(\theta)=\Im\frac{e^{i\theta}}\lambda
 =\frac{A\sin\theta-t\cos\theta}{r^2},\qquad
 H(\theta)=\frac1A-\Re\frac{e^{2i\theta}}\lambda.
\]
Then, uniformly in the endpoint phase,
\[
 R_d(T)=\frac{4(d+1)^2B(\theta)^2}{T H(\theta)}+O_{A,t,d}(T^{-2}).
 \tag{1}
\]
In particular,
\[
 R_d(T)\longrightarrow0,\qquad
 \limsup_{T\to\infty}T R_d(T)=\frac{4A(d+1)^2}{A^2+t^2},\qquad
 \liminf_{T\to\infty}T R_d(T)=0.
 \tag{2}
\]
For \(d\ge1\), the sharper minimum scale is
\[
 \boxed{\displaystyle
 \liminf_{T\to\infty}T^3R_d(T)
 =\frac{A d^2(d+1)^2(d+2)^2}{3(A^2+t^2)^2}>0.}
 \tag{3}
\]
Thus the intact model has no positive lower bound for \(R_d\) independent of
\(T\). Its phase troughs decay as \(T^{-3}\), despite the exponential envelope
in the unnormalized function.

These statements keep \(d,A,t\) fixed. They make no assertion uniform as
\(d\to\infty\), \(A\to0\), or \(t\to0\).

## Exact coefficients and norm

For \(0\le j\le k\), write
\[
 D_{k,j}=\varphi_k^{(j)}(1)
 =\sqrt{2k+1}\frac{(k+j)!}{j!(k-j)!},\qquad
 \varphi_k^{(j)}(0)=(-1)^{k+j}D_{k,j}.
\]
Repeated integration by parts terminates after \(k+1\) steps and gives the
**finite exact identity**
\[
 e^{-AT}\langle q_T,\varphi_k\rangle
 =\Im\left[
 e^{i\theta}\sum_{j=0}^k\frac{(-1)^jD_{k,j}}{(T\lambda)^{j+1}}
 -e^{-AT}(-1)^k\sum_{j=0}^k\frac{D_{k,j}}{(T\lambda)^{j+1}}
 \right].
 \tag{4}
\]
Using \(\sin^2x=(1-\cos2x)/2\) also gives
\[
 e^{-2AT}\|q_T\|_2^2
 =\frac{H(\theta)-e^{-2AT}H(0)}{4T},\qquad
 H(0)=\frac{t^2}{A r^2}.
 \tag{5}
\]
The phase-dependent denominator cannot be replaced by \(1/A\): both
\(A T\) and \(tT\) grow, with fixed ratio. The oscillatory term in (5) stays
of the same order as the nonoscillatory term. Its positivity is uniform:
\[
 0<h_-:=\frac1A-\frac1r\le H(\theta)\le\frac1A+\frac1r.
 \tag{6}
\]
An elementary explicit upper bound follows as well. Since
\(|\varphi_k(u)|\le\sqrt{2k+1}\),
\(|\langle q_T,\varphi_k\rangle|\le\sqrt{2k+1}e^{AT}/(AT)\).
For \(T\ge\log(2)/A\), (5), (6), and \(H(0)<2h_-\) imply
\[
 0\le R_d(T)\le\frac{8(d+1)^2}{A^2h_-T}.
 \tag{7}
\]
This already proves vanishing without any phase subsequence or cancellation
operation.

## The leading envelope and phase troughs

The first two terms of (4), uniformly in \(\theta\), are
\[
 e^{-AT}\langle q_T,\varphi_k\rangle
 =\frac{D_{k,0}B(\theta)}T
 -\frac{D_{k,1}C(\theta)}{T^2}
 +O(T^{-3})+O(e^{-AT}/T),
 \quad C(\theta)=\Im\frac{e^{i\theta}}{\lambda^2}.
 \tag{8}
\]
Use \(D_{0,1}=0\) in this formula. Since
\(\sum_{k=0}^dD_{k,0}^2=(d+1)^2\), equations (5), (6), and (8) prove
(1). The identity
\[
 A H(\theta)=r^2B(\theta)^2+\sin^2\theta
\]
shows that its leading coefficient is at most
\(4A(d+1)^2/r^2\), attained when \(\theta=n\pi\). These windows prove
the limsup in (2).

On the different sequence
\[
 T_n=\frac{\beta+n\pi}{t},
\]
we have \(B(tT_n)=0\), \(H(tT_n)=H(0)\), and
\(C(tT_n)=(-1)^{n+1}t/r^3\). Therefore, for \(d\ge1\),
\[
 R_d(T_n)\sim
 \frac{4A d^2(d+1)^2(d+2)^2}{3r^4T_n^3},
 \tag{9}
\]
using
\[
 \sum_{k=0}^d(2k+1)k^2(k+1)^2
 =\frac{d^2(d+1)^2(d+2)^2}{3}.
\]
This proves the liminf in (2). For \(d=0\), the lower endpoint in (4)
supplies the entire coefficient and gives the **exact** expression
\[
 R_0(T_n)=\frac{4A}{r^2T_n(e^{2AT_n}-1)}.
 \tag{10}
\]
There are also exact zeros of the degree-zero coefficient at arbitrarily
large \(T\): its numerator is
\(e^{AT}(A\sin(tT)-t\cos(tT))+t\), whose sign changes on intervals
surrounding successive zeros of its leading sinusoidal factor. Continuity
then gives roots. Hence the degree-zero liminf in (3), with zero right side,
is also zero.

## Proof of the sharp positive trough constant for d at least one

Let \(v=(D_{k,0})_{k=0}^d\), \(w=(D_{k,1})_{k=0}^d\), and put
\[
 K=\|v\|^2=(d+1)^2,\quad
 L=v\cdot w=\frac{d(d+1)^2(d+2)}2,\quad
 M=\|w\|^2=\frac{d^2(d+1)^2(d+2)^2}3.
\]
Equations (5) and (8) yield
\[
 T^3R_d(T)=\frac{4}{H(\theta)+o(1)}
 \left\|T B(\theta)v-C(\theta)w+O(T^{-1})\right\|^2.
 \tag{11}
\]
Any sequence on which the left side stays bounded must have \(TB\) bounded;
otherwise the nonzero vector \(v\) dominates the bounded other terms.
Consequently \(B\to0\), so modulo \(\pi\) its phase approaches \(\beta\).
After multiplication by the irrelevant sign \((-1)^n\) and passing to a
subsequence, every finite limiting value of (11) has the form
\[
 \frac4{H(0)}\left\|xv+\frac{t}{r^3}w\right\|^2,
 \qquad x\in\mathbb R.
\]
Its minimum is
\[
 \frac{4t^2}{H(0)r^6}\left(M-\frac{L^2}K\right)
 =\frac{A d^2(d+1)^2(d+2)^2}{3r^4}.
\]
This proves the lower inequality for the liminf. It is attained by choosing
windows solving
\[
 tT=\beta+n\pi-\frac{tL}{r^2KT},
\]
namely the larger root of
\(tT^2-(\beta+n\pi)T+tL/(r^2K)=0\), for sufficiently large \(n\).
Then \((-1)^nTB\to-tL/(r^3K)\), the minimizing value. This proves (3).

## Detrending, partial cancellation, and fixed endpoint tapers

For any orthogonal projector \(P\), with \(R=\|Pq\|^2/\|q\|^2\),
\[
 R(q-\alpha Pq)
 =\frac{(1-\alpha)^2R}{1-(2\alpha-\alpha^2)R},\qquad0\le\alpha\le1,
 \tag{12}
\]
provided the residual has nonzero norm. Its square root is nonincreasing in
\(\alpha\), but generally is **not** \((1-\alpha)\sqrt R\); the total norm
changes as well.

Continuous orthogonal detrending by \(P_m\), with \(m<d\), gives
\[
 R_d((I-P_m)q)=\frac{R_d(q)-R_m(q)}{1-R_m(q)}.
 \tag{13}
\]
In (1), the leading factor \((d+1)^2\) becomes
\((d+1)^2-(m+1)^2\). Thus ordinary fixed-degree orthogonal detrending
preserves uniform vanishing and the leading-cancellation sequence. If
\(m\ge d\), the remaining degree-\(d\) projection is zero. This requires
matching inner products; it is not an identity for an unrelated sampled
regression.

For a fixed bounded real taper satisfying
\(w(1-v)=c v^p+O(v^{p+1})\) as \(v\downarrow0\), with \(c\ne0\),
\(p\ge0\), and enough endpoint regularity for the stated remainder, Laplace
scaling \(v=T(1-u)\) gives
\[
 \langle wq_T,\varphi_k\rangle
 =\frac{c e^{AT}\Gamma(p+1)\varphi_k(1)}{T^{p+1}}
   \Im\frac{e^{i\theta}}{\lambda^{p+1}}
 +O(e^{AT}T^{-p-2}),
\]
\[
 \|wq_T\|_2^2
 =\frac{c^2e^{2AT}}{T^{2p+1}}\bigl(J_p(\theta)+O(T^{-1})\bigr),
\quad
 J_p(\theta)=\frac{\Gamma(2p+1)}2
 \left[(2A)^{-2p-1}-\Re\frac{e^{2i\theta}}{(2\lambda)^{2p+1}}\right]>0.
\]
The last inequality has a positive lower bound uniform in phase because
\(|\lambda|>A\). Consequently
\[
 R_d(wq_T)
 =\frac{(d+1)^2\Gamma(p+1)^2
       [\Im(e^{i\theta}/\lambda^{p+1})]^2}{T J_p(\theta)}
   +O(T^{-2})=O(T^{-1}).
 \tag{14}
\]
Hann and nontrivial Tukey endpoint tapers have \(p=2\). Such fixed tapers
do not produce a positive limiting visibility for this model. At a zero of
the leading numerator, the stated estimates give \(R_d=O(T^{-3})\).

## Scope for B-RES

The model refutes the assertion that the intact real off-critical mode itself
has \(\liminf R_d>0\), including after the regular tapers described above.
The arithmetic hypothesis in B-RES uses a different denominator,
\(\|q_T^{\rm arith}\|\), and concerns the projection of a total arithmetic
residual in the presence of other contributions. The theorem here does not
identify that denominator or control those other contributions. Therefore it
does **not** refute arithmetic RH or, by itself, settle arithmetic B-RES.
It removes the claimed pure-mode support for a positive, T-independent
self-normalized visibility threshold.

## Implementation and evidence

- `src/pbss/off_critical_exact.py` evaluates continuous moments via a scaled
  modified-Bessel identity and the exact norm. It does not use the sampled
  projector. Results are float64 evaluations, not interval certificates.
- `tests/test_off_critical_exact.py` compares those values with independent
  adaptive integration and verifies projector identities and phase limits.
- `experiments/check_off_critical_exact.py` evaluates (4) and (5) independently
  with optional `mpmath` at 80 decimal digits, checks one integral directly,
  and records the distinct phase subsequences. Its numerical evidence
  corroborates the proof; finite samples do not prove its asymptotic claims.
