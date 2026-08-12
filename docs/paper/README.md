# PBSS status note (short paper)

- **Source:** `pbss_status_note.tex`
- **PDF:** `pbss_status_note.pdf`
- **JPEGs:** `page1.jpg` … `page5.jpg` (page renders)
- **Figures:** shared from [`../figures/`](../figures/) (not duplicated here)

Build:
```bash
cd docs/paper
tectonic pbss_status_note.tex
pdftoppm -jpeg -r 150 pbss_status_note.pdf page
# rename page-1.jpg → page1.jpg etc. if needed
```

Independent research note — not a proof of RH.

`pbss_status_note.tex` covers lemmas M1–M6+, campaigns, and Theorem-A scaffolding.
**Not a proof of RH.**
