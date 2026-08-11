# Anthropic — related Riemann zeta work (reference)

**Local copies for offline reference.** Third-party material; copyright Anthropic / authors of the linked artifacts. Not part of PBSS, not an endorsement, and **not a claim that PBSS overlaps this result**.

## Primary link

- Blog / announcement: [Learning more about Claude's mathematical capabilities](https://www.anthropic.com/research/riemann-zeta)  
  (Anthropic Research, 10 Aug 2026)

## What they report (one paragraph)

An unreleased research version of Claude improved a **lower bound on the proportion of non-trivial zeros of \(\zeta\) that lie on the critical line**, from the prior literature figure **~41.6% to ~67.2%**. The argument draws on work in the **Levinson–Conrey / zero-proportion** tradition, including results of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh and Bombieri, with a Weil-type quadratic-form rank inequality. Anthropic states the techniques are **not** expected to prove RH. Formalization: Lean repo below.

## Local PDF copies

| File | Description |
|------|-------------|
| [`claude-paper.pdf`](claude-paper.pdf) | Claude’s technical paper |
| [`anthropic-informal-note.pdf`](anthropic-informal-note.pdf) | Anthropic informal note (concise expert statement) |
| [`claude-appendix-methodology.pdf`](claude-appendix-methodology.pdf) | Claude’s methodology appendix |

## Upstream URLs (canonical)

| Artifact | URL |
|----------|-----|
| Blog | https://www.anthropic.com/research/riemann-zeta |
| Paper PDF | https://www-cdn.anthropic.com/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf |
| Informal note | https://www-cdn.anthropic.com/23455459f8832d06bb175cc0f88d019aed962ef8.pdf |
| Methodology appendix | https://www-cdn.anthropic.com/d7f3ecf1d01392d887f8bc974ca187e2a121b1ed.pdf |
| Lean formalization | https://github.com/anthropics/zeta-23-lean |
| Process transcripts (linked from blog) | see Anthropic page “Further reading” |

## Relation to PBSS

| | Anthropic result | PBSS (this repo) |
|--|------------------|------------------|
| Object | Unconditional lower bound on **fraction of zeros on \(\mathrm{Re}s=\tfrac12\)** | Projection diagnostic \(R_d\) on **density residuals** / Beurling systems |
| Main tools | Weil quadratic forms, on/off-line definite subspaces, BGST/Bombieri line | Shifted Legendre projection, model lemmas M1–M6, arithmetic peel |
| RH claim | Explicitly **not** a proof of RH | Explicitly **not** a proof of RH |

**Bottom line for this archive:** keep these files as **context on concurrent AI-assisted zeta research**. PBSS does not implement or reproduce the 41.6% → 67.2% proportion bound.

## Provenance

Downloaded from Anthropic CDN URLs above into this directory for the PBSS research archive (Aug 2026). If PDFs go missing or the CDN moves, prefer the blog page as the index.
