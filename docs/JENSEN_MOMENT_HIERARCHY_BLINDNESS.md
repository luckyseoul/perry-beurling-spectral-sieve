# Jensen / moment hierarchy blindness for the de Bruijn–Newman constant Λ

**Technical note** · Nicholas D. Perry · August 2026  
**Code:** `pbss.jensen_blindness` · **Tests:** `tests/test_jensen_blindness.py`

---

## Explicit non-claims

- **Not a proof of the Riemann Hypothesis.**  
- **Not a new sharp analytic upper bound on Λ** beyond the optional bookkeeping remark §6.  
- **Not** a claim that Full Theorem B / B-RES is solved, and **not** a reopening of Full A packaging (Full A remains closed *conditionally* under RH + cited ANT in this repository).  
- The false-hyperbolicity ratios at \(t=-0.7\) are **historical project-record numerics** (Jun 2026); they are **not re-run** in the shipped module. The **index argument** is mathematical and checkable via pure helpers.

---

## 1. Setup

### 1.1 de Bruijn–Newman heat flow

Following the standard normalization (e.g. Polymath15),

\[
H_t(z)=\int_0^\infty e^{t u^2}\Phi(u)\cos(zu)\,du,
\]

with \(H_0\) proportional to the Riemann \(\xi\)-function on the critical line. There is a
real constant \(\Lambda\) (de Bruijn–Newman) such that \(H_t\) has only real zeros for all
\(t\ge\Lambda\). Rodgers–Tao proved \(\Lambda\ge 0\). **RH is equivalent to \(\Lambda\le 0\).**

Upper bounds on \(\Lambda\) in the Polymath15 framework use dynamical inequalities and
*local* zero-free / barrier information at large height — not bulk spectral statistics alone.

### 1.2 Moments and central Jensen / Turán certificates

Associated even moments of the (heat-deformed) measure are

\[
m_k(t)=\int_0^\infty u^{2k}\,e^{t u^2}\Phi(u)\,du
\qquad(k=0,1,2,\ldots).
\]

From the Taylor data of an entire function in the Laguerre–Pólya class (or from the
moment sequence), one builds **Jensen polynomials** and **Turán / Laguerre inequalities**
(Csordas–Norfolk–Varga; Griffin–Ono–Rolen–Zagier). A **central** certificate of shift
index \(n\) uses only a finite initial segment of that data.

**Working map used in the Jun 2026 campaign** (and in code
`central_shift_moment_order`):

\[
M(n)=n+4\qquad\Rightarrow\qquad n\le 30\ \Longrightarrow\ M=34.
\]

Even moments through order \(M\) means \(m_k\) for \(2k\le M\), i.e. at most
\(\lfloor M/2\rfloor\) independent even-moment parameters after scaling.

---

## 2. The index argument (checkable)

### 2.1 Lemma (even-moment resolution bound)

**Statement.** An even-moment certificate that uses moments only through maximum order
\(M\ge 0\) has at most

\[
N_{\max}(M)=\Big\lfloor\frac{M}{2}\Big\rfloor
\]

free spectral parameters (zero-pair slots) in the standard Hankel / generating-function
counting. Consequently it cannot resolve a distinguished zero of **ordinal index**
\(N>N_{\max}(M)\).

**Proof.** The even moments \(m_0,\ldots,m_K\) with \(2K\le M\) give \(K+1\) real numbers.
One overall scale is free; at most \(K=\lfloor M/2\rfloor\) independent shape parameters
remain to place or constrain oscillatory / zero features. A zero with ordinal index
\(N>K\) is outside that finite-parameter model. □

**Code:** `max_zero_ordinal_probed_by_even_moment_order(M)`.

### 2.2 Corollary (moment order needed for a given ordinal)

**Statement.** To have \(N_{\max}(M)\ge N\) one needs \(M\ge 2N\).

**Code:** `min_even_moment_order_for_zero_ordinal(N)`.

### 2.3 Application to the accessible hierarchy vs the Lehmer pair

| Object | Value | Source |
|-------|------:|--------|
| Accessible central shifts | \(n\le 30\) | project campaign |
| Accessible moment order | \(M=34\) | \(M(n)=n+4\) |
| Ordinals probed | \(N_{\max}(34)=17\) (bound); campaign text \(\sim 30\) zeros | **Lemma 2.1** / record |
| Approx. height of low zeros | \(\lesssim 101\) (table / record) | project record + `approx_height_of_nth_zero` |
| Binding Lehmer pair | zero **#~6709**, \(\gamma\approx 7005.062866,\ 7005.100565\) | standard tables; project record |
| Moment order needed | \(M\ge 2\cdot 6709=13418\) (\(\sim 13400\)) | **Corollary 2.2** |

Thus the accessible central hierarchy (**\(M=34\)**) and the Λ-binding Lehmer pair
(**\(N\approx 6709\)**) are separated by a **structural** gap of order

\[
M_{\mathrm{needed}}-M_{\mathrm{access}}=13418-34=13384,
\]

not by a few more digits of precision.

**Code entry point:** `index_argument_report()` / `jensen_blindness_report()`.

---

## 3. False hyperbolicity at \(t=-0.7\) (historical numeric)

**From the project record (Jun 2026; not re-run in `pbss.jensen_blindness`):**

> The central Turán / Jensen certificate **falsely certified hyperbolicity** of the
> relevant \(H_t\)-side data **down to heat-time \(t=-0.7\)**, with the minimum Turán
> ratio barely moving: \(\approx 1.04343\) at \(t=0\) to \(\approx 1.04336\) at
> \(t=-0.7\).

Moments were computed to order 34 at high precision; \(m_0(0)=\xi(1/2)/8\) was checked to
~10 digits. Turán and higher Jensen conditions passed at \(t=0\), consistent with the
literature on LP-class inequalities.

**Interpretation.** The certificate remained in the “all-real / hyperbolic” regime deep
into negative heat-time on **central** data, while the mechanism that actually forces
\(\Lambda> -\varepsilon\) is the **local** Lehmer pair at height \(\sim 7005\), invisible
to the central moment budget (Section 2).

Shipped label: `HISTORICAL_FALSE_HYPERBOLICITY` in `jensen_blindness.py` with
`rerun_here=False`.

---

## 4. Model sensitivity demo (re-runnable, not \(H_t\))

To illustrate **scale blindness** without recomputing \(\Phi\), the module evaluates a
product model \(f(x)=\prod_j(1+x^2/\gamma_j^2)\) on the first \(n\) shipped zeta ordinates
and again after adjoining a distant height \(\sim 7005\). The change in the discrete
Laguerre/Turán ratio is \(O(1/h^2)\):

```text
sensitivity = |τ(with distant) − τ(low)|  ≍  C / h²
```

**Code:** `model_false_hyperbolicity_demo`, `laguerre_L1_ratio`.  
This supports the slogan “central certificates are blind at scale \(1/h^2\) to height-\(h\)
features”; it is **not** a substitute for the historical \(H_t\) run.

---

## 5. Conclusion: local-at-height, not bulk/central

1. **Index.** Even-moment certificates of order \(M\) probe at most \(\lfloor M/2\rfloor\)
   zero ordinals.  
2. **Lehmer.** The pair that controls small negative heat-time sits at \(N\sim 6709\),
   requiring \(M\gtrsim 13400\).  
3. **Accessible hierarchy.** Central shifts \(n\le 30\) (\(M=34\)) live far below that
   budget.  
4. **False hyperbolicity (historical).** Central certificates can look “healthy” down to
   \(t=-0.7\) while remaining structurally silent about the binding pair.  

Therefore upper bounds on \(\Lambda\) obtained from **central** Jensen/moment hierarchies
are **local-at-height** phenomena: they cannot, as a matter of parameter count, replace
large-height barrier / Lehmer analysis. This is a **limit theorem for a method class**,
not a value of \(\Lambda\) and not RH.

---

## 6. Optional bookkeeping: \(\Lambda\le 0.20\)

Feeding Platt–Trudgian’s rigorous verification of RH to height \(3\times 10^{12}\)
(arXiv:2004.09765) into the Polymath15 / Ki–Kim–Lee dynamical framework yields the
updated numerical upper bound

\[
\Lambda\le 0.20
\]

(improving the published Polymath figure \(0.22\)). This is **input-height bookkeeping**,
not a new method, and is **not** the main result of this note.

---

## 7. Relation to PBSS Full A/B

| PBSS object | Interaction with this note |
|-------------|----------------------------|
| Full A (`closed_conditional`) | Orthogonal — residual \(R_d\) under RH + cited ANT |
| Full B / **B-RES** | Orthogonal — not solved here |
| Projection diagnostic | Different instrument; same meta-lesson: bulk certificates miss local binding constraints |

---

## 8. How to reproduce the checkable parts

```bash
cd perry-beurling-spectral-sieve
PYTHONPATH=src python3 -c "from pbss.jensen_blindness import jensen_blindness_report; import json; print(json.dumps(jensen_blindness_report(), indent=2, default=str))"
PYTHONPATH=src python3 -m pytest tests/test_jensen_blindness.py -v
```

---

## 9. References (selected)

- Csordas, Norfolk, Varga (1986), Turán inequalities and the zeros of the Riemann \(\xi\)-function.  
- Griffin, Ono, Rolen, Zagier (2019), Jensen polynomials for the Riemann zeta function and related functions.  
- Polymath15 (Tao et al.), effective estimates for the de Bruijn–Newman constant.  
- Rodgers, Tao (2020), \(\Lambda\ge 0\).  
- Platt, Trudgian (arXiv:2004.09765), rigorous RH verification to height \(3\times 10^{12}\).  
- Project record: [`PROJECT_RECORD_2025-11_to_2026-06.md`](PROJECT_RECORD_2025-11_to_2026-06.md) §9 Thread 3.
