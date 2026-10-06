# Corpus Document Policy (OCR → digest → index)

**Adopted:** 2026-06-19. **Motivation:** Szebehely 1967 "Theory of Orbits"
— the foundational CR3BP textbook the entire genome rests on — sat in
the private paper corpus for weeks **undigested**, because every digest wave
swept only newly-*acquired* papers and never the pre-existing corpus. It
was also a 661-page image-only scan needing OCR. This policy closes both
gaps so no document silently goes unprocessed again.

## The rule

A document in the private paper corpus is **"processed"**
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
- Language packs: `brew install tesseract-lang` (installed 2026-10-06, all languages, about 690 MB).
  Use `-l rus+eng` for KIAM preprints and `-l fra+eng` for French scans. Russian TeX-typeset PDFs
  with Type 3 fonts have garbled text layers: run `ocrmypdf --force-ocr`, file the original PDF
  unchanged, and add the OCR text as a `.txt` sidecar.
- Old solid RAR archives (inside KIAM source zips): `unar` (`brew install unar`); `7z` and `bsdtar`
  cannot open them. Translations and renders: `tectonic` (`brew install tectonic`).

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

#### Numbers from OCR or a text layer (adopted 2026-10-05, `#964`, from the `#960` papercuts)

- **Never trust OCRed 4s.**
  - Tesseract (inside ocrmypdf) reads the typewriter digit 4 as `h`, `L` or `u`. In the Menning 1968
    thesis tables, "4474" came out as "Lush" and "4924" as "hgok".
  - It also drops table rows, which silently shifts any row-by-row alignment.
  - It scrambles double-spaced typewritten body text.
- **Table recipe for scans:**
  1. Render the table pages at 400 dpi in grey:
     `pdftoppm -f P -l P -r 400 -gray -png in.pdf out`.
  2. Run Tesseract with a digit whitelist:
     `tesseract out-P.png stdout --psm 6 -c tessedit_char_whitelist="EV0123456789.- "`.
     Put the column letters the table uses into the whitelist.
  3. Accept only well-formed cells, for example `^0\.\d{3}$` for a 3-decimal column, and compare them
     with an independent copy.
  4. Read every disputed or shifted page by vision.
  5. Record which pages were image-checked.
- **Publisher text layers also corrupt numbers.** Seen in this corpus:
  - lost minus signs (Russell & Strange 2009, Table 6).
  - digit groups split by thin spaces ("0.310 73").
  - "O." printed for "0.".
  - misread labels ("II" for l1, "12" for l2: Henon & Guyot 1970).
- **The rule for text-layer numbers:** any number taken from a text layer or from OCR into a golden test
  value, a catalogue field or a `data/sources` file gets a page-image check. Record the page.
- **The second-witness rule for transcribed tables:** a table transcribed into `data/sources/` needs a
  second witness before it is used as a control or a golden. The witness can be:
  - another copy (a text-layer edition, the thesis behind the paper);
  - an independent OCR pass; or
  - a physics or self-consistency identity between the table's own columns, checked row by row.
  Precedent: the Hollister & Menning 1970 Table 3 YAML, transcribed by eye from an image scan, had 27 wrong
  cells (18 of them turn angles with 3 read as 5). The periapsis identity
  r_p = mu/V^2 (1/sin(theta/2) - 1) against the printed Rmin flags most of them at once. Details:
  `docs/notes/2026-10-05-hollister-menning-1970-table3-recheck.md`.
- **pdftotext stderr:** on many AIAA PDFs, `pdftotext` writes tens of thousands of
  "Syntax Warning: Badly formatted number" lines to stderr. Always use `pdftotext ... 2>/dev/null`; the
  flood truncates tool output and hides real errors.
- **DOIs:** confirm DOIs with `scripts/crossref_check.py`. Crossref's bibliographic search often returns
  the AIAA conference record (10.2514/6.YYYY-NNNN) before the journal version, so a DOI is CONFIRMED only
  when the journal or volume matches too.
- **Index filenames:** every `CORPUS_INDEX.md` row names its file IN FULL, with no "..." abbreviations, so
  index-to-disk checks can be scripted. All abbreviated and stale names were expanded on 2026-10-05
  (`#964`); all 412 PDFs on disk are named in the index.

#### Wanted lists: check against the whole corpus first (adopted 2026-10-05, `#964`)

Before any "wanted papers" or acquisition list goes to the owner, run
`python scripts/check_wanted_vs_corpus.py <list.md>` and resolve every hit by hand. The checker:
- extracts the first 3 pages of every `papers/` PDF, plus any `.txt` companion;
- for each row of the list, reports (1) any DOI in the row found in a corpus text, filename or
  `CORPUS_INDEX.md`; (2) the first author's surname among a corpus filename's author tokens (any
  position) with a filename year within one year, or the same in a `CORPUS_INDEX.md` `[identity: ...]`
  note, with transliteration aliases (Bruno/Brjuno, Henon/Hénon, Olle/Ollé, ...), or the row's whole
  title in an identity note; (3) the opening of the quoted title, spaces removed, found in a corpus text;
  (4) a report, paper or thesis number in the row (AAS, AIAA, JPL TR, NASA TN/TM/CR, MIT TE/RE, NTRS)
  found in a corpus filename or an identity note.

Treat the hits as follows:
- Title hits are mostly citations inside other papers.
- Surname+year and DOI hits need a human look.
- A row whose paper is held is removed. A row whose content is held in another form (a report, preprint or
  thesis version) is marked "acquire only for attribution".

Precedent: the first `#960` wanted list carried Bruno 1978a, b (held as "brjuno-1978") and Henon 2001 (held).
A filename grep on "bruno" missed the first, and nobody searched for the second. The batch-10 list then asked for
Rall 1969, held for months as "hollister-rall-1970-periodic-orbits-NASA-CR.pdf" (the chairman's name and
the report year); the first checker version missed it.

**Identity notes in the index.** When a filename does not show the title-page first author and year (a
thesis filed under its supervisor, a report reissue year, a transliteration), the index row carries
`[identity: Author, I. (Year), "Title", what it is (thesis, report number, NTRS id)]`. Add one too when
the held file is a journal version of a report that is cited by its number (Hollister 1969 JSR = MIT RE-36).
The checker reads these notes: surname and year, the whole title, and report numbers. Prefer a corrective `git mv` in `cyclers_pdf` plus a redirect row in the index when the
filename is actively misleading.

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
  private paper corpus + the project memory): filing a PDF is not "done"
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
