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

### 1. OCR-first (text-searchable)
On acquisition, probe for a text layer (`pdffonts <f>`, or a `pdftotext`
probe). If image-only, it must be OCR'd — either by a digest agent reading
page-images through the Read tool (the Szebehely method: selective, works
on 600+ page scans) or an explicit OCR pass. Text-layer PDFs use
`pdftotext -layout`. No document stays a black-box image.

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
