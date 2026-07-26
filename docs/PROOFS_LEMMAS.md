# Proofs of Lemmas M1–M4

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
