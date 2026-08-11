# PBSS status note (short paper)

- **Source:** `pbss_status_note.tex`
- **PDF:** `pbss_status_note.pdf` (5 pages)
- **JPEGs:** `page1.jpg` … `page5.jpg` (page renders)
- **Figures:** `figures/` (pipeline, status map, grand campaign plots)

Build:
```bash
tectonic -X compile pbss_status_note.tex
pdftoppm -jpeg -r 200 pbss_status_note.pdf page
```

Independent research note — not a proof of RH.

## Status note contents (2026-07-26 refresh)

`pbss_status_note.tex` covers lemmas M1–M6, grand campaign to \(10^{10}\), open-plateau
through \(5\times 10^{10}\), and Theorem-A scaffolding (weight class + remainder).
Rebuild PDF with `pdflatex pbss_status_note.tex` when TeX is available.
**Not a proof of RH.**
