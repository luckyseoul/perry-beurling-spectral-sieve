# Proofs of Lemmas M1–M6

All integrals are on \([0,1]\) with measure \(du\). Let \(\{\varphi_k\}\) be the
orthonormal shifted Legendre system of `basis.py`.

---

## Lemma M1 (pure mode energy)

**Statement.** \(R_d(\varphi_m)=1\) if \(m\le d\), else \(0\).

**Proof.** Orthonormality gives \(\langle\varphi_m,\varphi_k\rangle=\delta_{mk}\).
Hence
\[
\|P_d\varphi_m\|^2=\sum_{k=0}^d\delta_{mk}
=
\begin{cases}1&m\le d,\\0&m>d,\end{cases}
\qquad
\|\varphi_m\|^2=1,
\]
so \(R_d(\varphi_m)\) is as claimed. □

---

## Lemma M2 (orthogonal defect formula)

**Statement.** If \(j\le d\), \(\varepsilon\in[0,1]\), \(f\perp V_d\), \(\|f\|_2=1\), and
\(q=\sqrt{1-\varepsilon^2}\,f+\varepsilon\varphi_j\), then \(R_d(q)=\varepsilon^2\).

**Proof.** For \(0\le k\le d\),
\[
\langle q,\varphi_k\rangle
=\sqrt{1-\varepsilon^2}\,\langle f,\varphi_k\rangle+\varepsilon\langle\varphi_j,\varphi_k\rangle
=\varepsilon\,\delta_{jk},
\]
since \(f\perp V_d\). Therefore
\[
\|P_d q\|^2=\sum_{k=0}^d\varepsilon^2\delta_{jk}=\varepsilon^2.
\]
Also
\[
\|q\|^2=(1-\varepsilon^2)\|f\|^2+\varepsilon^2\|\varphi_j\|^2
+2\varepsilon\sqrt{1-\varepsilon^2}\langle f,\varphi_j\rangle
=1-\varepsilon^2+\varepsilon^2=1,
\]
using \(\langle f,\varphi_j\rangle=0\). Hence \(R_d(q)=\varepsilon^2\). □

---

## Lemma M3 (critical-line pure mode decay)

**Statement.** \(R_d(\sin(\omega\,\cdot))=O_d(\omega^{-2})\) as \(\omega\to\infty\).

**Proof.** Fix \(k\le d\). The function \(\varphi_k\) is a polynomial, hence \(C^\infty\).
Integration by parts:
\[
\int_0^1\sin(\omega u)\,\varphi_k(u)\,du
=\Bigl[-\frac{\cos(\omega u)}{\omega}\varphi_k(u)\Bigr]_0^1
+\frac1\omega\int_0^1\cos(\omega u)\,\varphi_k'(u)\,du.
\]
The boundary term is \(O(\omega^{-1})\) and the integral is \(O(\omega^{-1})\)
(again integrate by parts, or bound by \(\|\varphi_k'\|_1\)). Thus
\[
|c_k|=|\langle\sin(\omega\cdot),\varphi_k\rangle|\le\frac{C_k}{\omega}.
\]
Also
\[
\|\sin(\omega\cdot)\|_2^2=\int_0^1\sin^2(\omega u)\,du
=\frac12-\frac{\sin(2\omega)}{4\omega}\to\frac12,
\]
so for large \(\omega\), \(\|\sin\|_2^2\ge1/4\). Therefore
\[
R_d
=\frac{\sum_{k=0}^d|c_k|^2}{\|\sin\|_2^2}
\le 4\sum_{k=0}^d\frac{C_k^2}{\omega^2}
=O_d(\omega^{-2}).
\]
For \(\omega=tT\) with fixed \(t>0\), this is \(O(T^{-2})\). □

**Leading \(k=0\) asymptotics (used in numerics).**
\[
c_0=\int_0^1\sin(\omega u)\,du=\frac{1-\cos\omega}{\omega},
\qquad
R_0=\frac{c_0^2}{\|\sin\|_2^2},
\]
implemented as `predicted_R_d_critical_scaling`.

---

## Lemma M4 (persistent defect blocks vanishing)

**Statement.** Under M2 with fixed \(\varepsilon>0\), \(R_d(q)=\varepsilon^2\) for every
choice of unit \(f\perp V_d\).

**Proof.** Immediate from M2: the right-hand side does not depend on \(f\). □

**Diagnostic corollary.** Along any family \(q_T\), if \(R_d(q_T)\to0\) then the
representation of \(q_T\) cannot keep a fixed low-degree mass \(\varepsilon_0>0\)
orthogonal-defect form. Vanishing of \(R_d\) is necessary for “no persistent
polynomial defect.” Converting that into RH is **not** included in M4.

---

## Lemma M5 (finite critical-line superposition — finite-mode A₀)

**Statement.** Let \(N<\infty\), amplitudes \(a_n\in\mathbb{R}\), ordinates \(t_n>0\),
phases \(\phi_n\in\mathbb{R}\), and
\[
q_T(u)=\sum_{n=1}^{N}a_n\sin(t_n T u+\phi_n)
\]
(or the cosine form with the same frequencies). Then
\[
R_d(q_T)=O_d(T^{-2})\qquad(T\to\infty)
\]
at the **same order** as pure-mode M3 (not a weaker rate). Code:
`bound_R_d_finite_mode_sum`, `finite_cl_superposition`, tests `test_M5_*`.

**Proof.** Fix \(k\le d\). Linearity and the M3 integration-by-parts bound give
\[
\bigl|\langle q_T,\varphi_k\rangle\bigr|
\le\sum_{n=1}^{N}|a_n|\,\frac{C_k}{t_n T}
=\frac{C_k}{T}\sum_{n=1}^{N}\frac{|a_n|}{t_n}.
\]
Hence
\[
\|P_d q_T\|_2^2
=\sum_{k=0}^{d}|\langle q_T,\varphi_k\rangle|^2
\le\frac{(d+1)C_{\max}^2}{T^2}\Bigl(\sum_n\frac{|a_n|}{t_n}\Bigr)^2
=O_{d,N,\{a,t\}}(T^{-2}).
\]
Expanding \(\|q_T\|_2^2\) into diagonal sine terms \(\to\tfrac12\sum a_n^2\) and
cross terms \(O(T^{-1})\), for large \(T\) and \(a\not\equiv0\) we have
\(\|q_T\|_2^2\ge\tfrac14\sum a_n^2>0\). Therefore
\[
R_d(q_T)=\frac{\|P_d q_T\|_2^2}{\|q_T\|_2^2}=O_d(T^{-2}).
\]
□

**Remark.** This is the finite-mode extension of model theorem A₀ needed for
truncated explicit-formula residuals. It does **not** control infinite zero sums
or the arithmetic residual \(q_T=(\theta-x)/\sqrt{x}\); full Theorem A and RH
remain open.

---

## Lemma M6 (admissible weight preserves model-mode decay)

**Statement.** Let \(w\in L^\infty([0,1])\) with \(\|w\|_\infty\le W<\infty\), and suppose
for the family \(q_\omega(u)=\sin(\omega u)\) one has \(\|w q_\omega\|_2^2\ge c_w>0\) for
large \(\omega\) (true for Hanning/Tukey with a positive bulk). Then
\[
R_d(w q_\omega)=O_d(\omega^{-2})\qquad(\omega\to\infty).
\]
The same bound holds for finite critical-line superpositions (weighted finite-mode A₀).

**Proof.** For each \(k\le d\),
\[
|\langle w q_\omega,\varphi_k\rangle|
\le W\,|\langle q_\omega,\varphi_k\rangle|
\le W\frac{C_k}{\omega}
\]
by the M3 integration-by-parts bound on \(\langle q_\omega,\varphi_k\rangle\). Hence
\[
\|P_d(w q_\omega)\|_2^2
=\sum_{k=0}^d|\langle w q_\omega,\varphi_k\rangle|^2
\le W^2\sum_{k=0}^d\frac{C_k^2}{\omega^2}
=O_d(W^2\omega^{-2}).
\]
Dividing by \(\|w q_\omega\|_2^2\ge c_w\) yields \(R_d(w q_\omega)=O_d(\omega^{-2})\).
For a finite sum \(q_T^{(N)}=\sum_{n=1}^N a_n\sin(t_n T u+\phi_n)\), apply the same
estimate modewise and sum (as in M5) to get \(O_d(T^{-2})\). □

**Scope.** Model residuals only. **Not** a proof that \(R_d(w q_T^{\mathrm{arith}})\to0\).

**Code.** `lemmas.bound_R_d_weighted_sine_order`,
`lemmas.bound_R_d_weighted_finite_mode_sum`, weights in `pbss.weights`.


---

## Lemma M7 (\(R_d\) perturbation majorant)

**Statement.** Let \(q_0,r\in L^2([0,1])\) with \(\|q_0\|_2>0\) and
\(\delta:=\|r\|_2/\|q_0\|_2<1\). Set \(q=q_0+r\) and \(R_0=R_d(q_0)\). Then
\[
R_d(q)
\le
\frac{\bigl(\sqrt{R_0}+\delta\bigr)^2}{(1-\delta)^2}.
\]
In particular, if \(R_0\to0\) and \(\delta\to0\), then \(R_d(q)\to0\).

**Proof.** Orthogonal projection \(P_d\) is a contraction on \(L^2\), so
\[
\|P_d q\|_2
\le \|P_d q_0\|_2 + \|P_d r\|_2
\le \|P_d q_0\|_2 + \|r\|_2
= \sqrt{R_0}\,\|q_0\|_2 + \|r\|_2
= \bigl(\sqrt{R_0}+\delta\bigr)\,\|q_0\|_2.
\]
Also \(\|q\|_2\ge\|q_0\|_2-\|r\|_2=(1-\delta)\|q_0\|_2\). Therefore
\[
R_d(q)
=\frac{\|P_d q\|_2^2}{\|q\|_2^2}
\le
\frac{(\sqrt{R_0}+\delta)^2\|q_0\|_2^2}{(1-\delta)^2\|q_0\|_2^2}
=\frac{(\sqrt{R_0}+\delta)^2}{(1-\delta)^2}.
\]
The limit claim is immediate. □

**Code:** `pbss.ab_closure.energy_ratio_perturbation_bound`, discrete check
`verify_m7_on_grid`. **Role:** closes the triangle step in Full Theorem A under
small identification/tail/arithmetic remainders (ANT-1…3).

