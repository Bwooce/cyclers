# Corpus Document Policy (OCR → digest → index)

**Adopted:** 2026-06-19. **Motivation:** Szebehely 1967 "Theory of Orbits"
— the foundational CR3BP textbook the entire genome rests on — sat in
`cyclers_pdf/papers/` for weeks **undigested**, because every digest wave
swept only newly-*acquired* papers and never the pre-existing corpus. It
was also a 661-page image-only scan needing OCR. This policy closes both
gaps so no document silently goes unprocessed again.

## The rule

A document in the private `cyclers_pdf/papers/` corpus is **"processed"**
only when ALL THREE hold:

### 1. OCR-first (text-searchable) — use a tool, not Claude tokens
On acquisition, probe for a text layer (`pdffonts <f>`; or a `pdftotext`
probe returning empty = image-only).
- **Text-layer PDFs** (most modern arXiv/journal papers): just
  `pdftotext -layout`. No OCR needed.
- **Image-only PDFs** (old scans — Szebehely 1967, 1970s AIAA scans):
  run **`ocrmypdf --skip-text <in> <out>`** (Tesseract under the hood) to
  add a text layer ONCE, then `pdftotext`. This is a cheap deterministic
  CPU step — do NOT vision-read hundreds of page-images through the Read
  tool (the Szebehely digest burned large token cost doing exactly that).
- Tooling: `ocrmypdf` (PyPI, `uv add`) requires the `tesseract-ocr` +
  `ghostscript` system binaries (installed; helper `verify/ocr.py`, #400).
  No document stays a black-box image.

#### Three content classes (the hybrid rule)
OCR is bulk text only; two other content classes need Claude vision.
Cost is bounded by using the cheap OCR text to decide *which* pages to
vision-read — never vision-read a whole document.

| Content | Tool | Note |
|---|---|---|
| Body text | ocrmypdf/Tesseract | cheap, deterministic, greppable |
| Equations / tables (precise values) | Claude vision on the specific page | Tesseract garbles math, subscripts, table structure |
| **Diagrams / figures (semantic)** | Claude vision — **caption-guided** | Tesseract extracts ~nothing from plots/graphs |

- **Diagram rule:** Tesseract reads figure *captions* fine but extracts
  nothing from the figure itself. Many findings live in the diagram —
  Tisserand-Poincaré graphs, B-plane error ellipses, heliocentric
  trajectory plots, family-network figures, orbit-shape plots. Workflow:
  OCR text → read captions (cheap) → identify the load-bearing figures
  the captions reference → Claude-vision-read ONLY those figure pages to
  extract what the diagram shows. Bounded cost (N relevant figures, not
  all), full diagram coverage.
- **Fidelity guard (figure-derived ≠ sourced):** reading a *value off a
  plot* (a curve, a contour, an error ellipse) is **digitization, not
  OCR** — lossy. Vision tells you what a diagram *shows*; it does not
  authoritatively *measure* it. Flag any plot-read value `figure-derived`
  and never put it in a sourced numeric field without the explicit
  digitize rung (cf. the published-values-are-display discipline and the
  digitize step of the never-give-up-reproducing ladder). A published
  table or the original always beats an eyeballed plot.

### 2. Chapter/section-summary digest
A verdict note committed to `docs/notes/YYYY-MM-DD-digest-<slug>.md`.
- **Papers** (journal/conference): full-page read.
- **Books / theses**: chapter-summary scope — TOC + index + the 3–5
  project-relevant chapters deep-read, sample the rest. This is the
  Hintz-2023 / Belbruno-2004 / Parker-2007 / Szebehely-1967 pattern.
  **Never read a whole textbook page-by-page** (that triggers the 32 MB
  tool-result failure that killed an earlier agent).
- Every numeric value / claim carries a section/page citation
  (sourced-only discipline).

### 2b. arXiv-source data-recovery (before declaring "graphical-only")
**Adopted 2026-06-25 (#458; precedent #442).** A paper's rendered PDF
hides machine-readable data the arXiv **e-print source** carries: numeric
values printed full-precision in prose, ancillary `.dat`/`.csv`/`.txt`
data files, `\addplot table` / `\pgfplotstableread` figure data,
TikZ coordinate dumps, full-precision `\def`/`\newcommand` macros, and
commented-out `%` numeric blocks. **Before** declaring any paper value
"graphical-only" / "no IC table" / "not tabulated", or filing a
digitization gap, **PULL THE ARXIV E-PRINT SOURCE** and check it:

```
mkdir -p /tmp/x/<id> && cd /tmp/x/<id>
curl -sL -A "Mozilla/5.0 (research)" "https://arxiv.org/e-print/<id>" -o e.tar.gz
tar xzf e.tar.gz   # may instead be a single .tex.gz, or a bare PDF (no source)
find . -type f \( -name "*.dat" -o -name "*.csv" -o -name "*.txt" -o -name "*.tsv" \)
grep -rnE "begin\{tabular\}|addplot table|pgfplotstableread" *.tex
grep -rnoE "[0-9]+\.[0-9]{4,}" *.tex   # full-precision prose values + macros
```

- **Caveat (no source exists):** `arxiv.org/e-print` sometimes returns a
  bare PDF (author uploaded a PDF, not TeX) — then there is no source to
  mine and the graphical gap stands (e.g. Vasile-Campagnola 2009
  arXiv:1105.1823, #458). Conference/journal-only papers (AAS, IAC,
  MDPI) have **no arXiv ID at all** → no source path; skip.
- **What this recovers (real precedent, #458):** Antoniadou-Libert 2019
  (arXiv:1811.09442) renders all families graphically with "no IC table",
  yet the `SPA.tex` prose carries the v.c.o./bifurcation anchor points to
  full precision (`(e_1,i_1)=(0.0891812,90°)`, `0.000178`, `0.0002327`,
  …) and the per-DS-map constant orbital elements (`a_2/a_1=0.6312` /
  `0.7595` / `0.4806`, ω/Ω/M angles) — anchor data the figures hide.
- **Fidelity win:** a source-printed value is **sourced**, not
  figure-derived — it satisfies §1's digitize rung directly (no lossy
  plot-read), so the recovered value may go straight into a sourced
  numeric field. This is strictly better than digitizing the rendered
  figure. Still confirm topology before adopting.
- **Negatives still harden:** if the source carries only ΔV/TOF/proxy
  tables and no state vectors (Braik-Ross 2026 arXiv:2605.31543;
  Hiraiwa 2026 arXiv:2602.17444 — lobe radii stay histogram-only), the
  graphical-only verdict is now rigorously confirmed, not assumed.

### 3. Registered in the master index
One line in `docs/notes/CORPUS_INDEX.md` per `papers/` file:
`<filename> → <digest-note or mined-by pointer> | <one-line summary> | <status>`
where status ∈ {digested, mined-by-catalogue, mined-by-KNOWN_CORPUS,
undigested}. The index is the discoverability layer **and** the
anti-slip ledger — its absence is precisely what let Szebehely hide.
Updated in the SAME commit as the filing or digest.

## Where this binds

- **New acquisition flow** (see the filing standard in the
  `cyclers_pdf` repo + the project memory): filing a PDF is not "done"
  until its digest note + index line exist. "Filed, digested, indexed."
- **Periodic audit**: sweep `papers/` against the index + `catalogue.yaml`
  (first_published / corroborating_sources DOIs) + `literature_check.py`
  KNOWN_CORPUS (DOI/author) to catch any doc that is undigested AND
  unmined. Re-run after big acquisition waves.

## Honest prioritisation

A document already **mined** by a live catalogue row or KNOWN_CORPUS
anchor (its tables power real rows) is lower-priority for a standalone
digest, but still gets an index line recording *where* it is mined.
**Foundational-theory documents that no row mines** — Szebehely (CR3BP
theory), Gurfil 2007 (methods), the Doedel bifurcation papers — are the
highest-risk for slipping and MUST get the full digest. The genome
depends on them implicitly; the index makes that dependency explicit.

## Status

- Policy + enforcement memory: adopted 2026-06-19.
- Master index `CORPUS_INDEX.md`: to be built by the corpus
  digest-coverage audit (task #397) — sweeps all ~120 `papers/` files,
  classifies each, and seeds the living index.
- OCR-coverage sweep: image-only `papers/` files identified + OCR-status
  recorded as part of the same audit.
