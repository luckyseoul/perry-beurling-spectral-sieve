# L*: retracted coercivity proxy and retained model identities

**Status: `retracted_coercivity_proxy` (2026-09-05).** The proposed pure-mode
coercivity interpretation is false. Its phase-independent real-sine
asymptotic and claimed analytic lower envelope are retracted. Arithmetic
B-RES and the Riemann Hypothesis remain open.

The original finite-grid experiments remain historical observations. They do
not override the continuous model theorem in
[OFF_CRITICAL_PHASE_OBSTRUCTION.md](OFF_CRITICAL_PHASE_OBSTRUCTION.md).
No Guinand-Weil quadratic form was evaluated by this proxy, and no equivalence
between that form and polynomial projection energy was established.

## The mathematical correction

Let \(P_d\) be the continuous orthogonal projector onto shifted Legendre
polynomials of degree at most \(d\) on \([0,1]\), and fix \(A,t>0\). For
\[
 q_T(u)=e^{ATu}\sin(tTu),\qquad
 R_d(T)=\frac{\|P_dq_T\|_2^2}{\|q_T\|_2^2},\qquad C_d=\sqrt{R_d},
\]
the corrected model theorem is
\[
 R_d(T)=O(T^{-1})\longrightarrow0.
\]
It holds uniformly in the endpoint phase, without subtracting any projection
or imposing cancellation. Fixed regular endpoint tapers, including Hann and
Tukey, also preserve vanishing for this model.

The old expression, with \(a=AT\) and \(b=tT\),
\[
 \frac{4a(d+1)^2}{a^2+b^2}
 =\frac{4A(d+1)^2}{(A^2+t^2)T},
\]
is the **maximum of the leading term over endpoint phases**. Precisely,
\[
 \limsup T R_d(T)=\frac{4A(d+1)^2}{A^2+t^2},\qquad
 \liminf T R_d(T)=0.
\]
It is not a phase-independent equivalent of \(R_d(T)\). Its square root,
or half that square root, is not a uniform lower bound for \(C_d\).

Set \(\lambda=A+it\), \(\theta=tT\), and
\[
 B(\theta)=\Im(e^{i\theta}/\lambda),\qquad
 H(\theta)=1/A-\Re(e^{2i\theta}/\lambda)>0.
\]
The phase-aware statement is
\[
 R_d(T)=\frac{4(d+1)^2B(\theta)^2}{T H(\theta)}+O(T^{-2}).
\]
For \(d\ge1\), phase troughs have the sharper, positive scaled limit
\[
 \liminf T^3R_d(T)
 =\frac{A d^2(d+1)^2(d+2)^2}{3(A^2+t^2)^2}.
\]
The proof uses exact terminating antiderivatives for the polynomial moments
and an exact oscillatory norm. Both the proof and independent high-precision
corroboration are documented in the linked obstruction note.

The previous argument failed at two points: taking the imaginary part of a
complex endpoint moment can cancel its leading term, and the oscillatory
piece of the squared norm remains of the same order as its nonoscillatory
piece when \(AT\) and \(tT\) grow together.

## What remains valid under partial cancellation

For an orthogonal projector \(P\), write \(q=Pq+(I-P)q\) and
\(R=\|Pq\|^2/\|q\|^2\). If \(q-\alpha Pq\ne0\), then
\[
 C(q-\alpha Pq)
 =\frac{(1-\alpha)\sqrt R}{\sqrt{1-(2\alpha-\alpha^2)R}},
 \qquad 0\le\alpha\le1.
\]
This is nonincreasing in \(\alpha\). It generally differs from
\((1-\alpha)C(q)\), because cancellation also changes the denominator.
At full cancellation, the projection vanishes if the remaining signal is
nonzero. If \(q\) was entirely in the projected space, full cancellation
leaves the zero signal and its ratio is undefined.

These identities apply to the shipped weighted QR projector when every norm
uses the same quadrature. They are Hilbert-space projection facts. Calling
\(C_d\) a “Weil-visible” channel does not turn them into a theorem about a
Weil quadratic form.

## The arithmetic injection bridge is a calibration

The retained `arithmetic_injection_bridge_row` experiment constructs an
arithmetic residual, removes its degree-\(d\) projection, and then inserts a
chosen model bump \(b\). Let the orthogonalized background be \(r\), so
\(P_dr=0\). Linearity gives
\[
 C_d(r+b)=C_d(b)\frac{\|b\|}{\|r+b\|}.
\]
The code now uses matching trapezoid norms throughout this identity. Numerical
agreement checks the projection and normalization. No arithmetic zero forces
the inserted bump, so this experiment does not establish an arithmetic
footprint or prevent cancellation among genuine arithmetic contributions.

## Report and API interpretation

`pbss.hstar_anticancellation` retains historical function names for callers,
with explicit corrected metadata:

| Field or surface | Current meaning |
|---|---|
| `candidate_status` | `retracted_coercivity_proxy` |
| `R_d_leading` / `R_d_phase_leading` | Phase-aware additive asymptotic; not a relative estimate at a phase zero |
| `R_d_limsup_envelope` | Old expression, correctly labeled as the leading phase envelope |
| `R_d_continuous_exact` | Exact continuous identities evaluated in float64, not an interval certificate |
| `rel_err_Rd` | Sampled QR ratio compared with the continuous identity |
| `eps_star_analytic` | Legacy half-envelope reference; `eps_star_analytic_is_lower_bound=False` |
| `numeric_above_half_leading` | Historical finite-sample comparison only |
| `constraint_signal` / `all_constraint_signals` | Retired coercivity flags, now false |
| `sample_separation_signal` | Finite sample comparison between the chosen pure, cancelled, and critical-line signals |
| `eps_star` | Observed finite-grid minima; explicitly not a uniform lower bound |
| `bridge_ok` | Injected-bump projection calibration; not arithmetic B-RES |
| `weil_form_evaluated` | False |

A floor defined as half the observed pure-mode value is self-calibrated and
cannot establish an independent coercivity premise. A minimum over finitely
many sampled windows supplies no positive bound for all later windows.

Continuous asymptotics also require resolution checks before comparing them
with a sampled calculation: the oscillation period and endpoint concentration
width both shrink with \(T\). Keeping a fixed grid while increasing \(T\)
can introduce aliasing and endpoint-quadrature error. The phase-aware rows
report discrepancy from an independent continuous evaluation explicitly.

## Historical artifacts and review

The existing files under `results/hstar_anticancellation/` are preserved as
historical outputs. Their former interpretation as an analytic coercivity
envelope is superseded; their finite numerical values are not a proof of
arithmetic B-RES or RH.

| Historical artifact | Scope retained |
|---|---|
| `grid.json`, `multi_T.json` | Values on the listed grids only |
| `analytic_envelope.json` | Historical comparison with the old phase envelope; lower-bound interpretation withdrawn |
| `partial_cancel.json` | Numerical cancellation curves; corrected normalization identity above |
| `arithmetic_bridge.json` | Injected-bump calibration on arithmetic background |
| `hstar_proxy.png`, `SUMMARY.txt` | Historical presentation of that campaign |
| `referee_lstar.json` | Historical statement-shape review; not a proof and superseded on the disputed analytic claim |

New runs use `results/hstar_corrected/` by default. Rank-4 proportion
feasibility remains STOP and is a separate problem.

## Proof obligations for an arithmetic replacement

| Obligation | Status |
|---|---|
| Orthogonal projection and correctly normalized partial cancellation | Valid model identities |
| Phase-independent pure-real-sine equivalent or half-envelope lower bound | Retracted |
| Fixed-degree pure-real-sine vanishing and sharp phase subsequences | Proved model theorem |
| Identification of polynomial visibility with an actual Weil quadratic form | Not established |
| Off-critical arithmetic contribution forced by a genuine zero, with its true amplitude and denominator | Open |
| Non-cancellation for the total arithmetic residual, with secondary and remainder control | Open |
| Uniformity sufficient to derive arithmetic B-RES and RH | Open |

The new model theorem removes the claimed support for a positive limiting
self-normalized pure-mode threshold. The arithmetic residual has additional
terms and a different normalization question; the model theorem alone neither
proves nor disproves arithmetic B-RES or RH.
