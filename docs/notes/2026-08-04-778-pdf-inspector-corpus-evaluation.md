# `#778`: `pdf-inspector` corpus evaluation — do NOT use for digest reprocessing

**Task:** `#778` (tooling evaluation only — no code/catalogue changes). User asked to check for
a machine-wide memory documenting a "new PDF tool," and if found, use it to reprocess the whole
`cyclers_pdf/papers/` corpus (251 files) looking for improvements to existing digests.

**Tool found:** `pdf-inspect` (machine-wide CLI wrapping the `pdf-inspector` PyPI package,
github.com/firecrawl/pdf-inspector), installed 2026-08-02 during an unrelated hardware-schematic
session on this Mac. Documented at `~/.claude/memory/pdf_inspector_tool.md` (machine-wide, not a
`cyclers` project memory). That memory already carried one caution (worse than `pdftotext
-layout` on a dense schematic's pin tables) but no evaluation against this project's own corpus
type (academic papers: prose + equations + numeric tables, mixed scan/born-digital).

**Verdict: do NOT adopt for reprocessing/improving this corpus's digests.** Full evidence below;
machine-wide memory updated with the same findings for future sessions/projects.

## Method

Before committing to reprocessing all 251 files (which would touch many existing sourced
digests), ran orientation checks first, per [[feedback_verify_positive_control_source_applicability]]
and [[feedback_check_dont_guess]]-style discipline — verify the tool actually works on THIS
corpus's content type before trusting it, rather than assuming the existing machine-wide memory's
schematic-specific finding does or doesn't generalize.

1. **Corpus-wide mechanical sweep** (`pdf-inspect FILE.pdf --pages 1` over all 251 files in
   `cyclers_pdf/papers/`, capturing stdout byte count + stderr diagnostics + `pdfinfo` Producer;
   pure mechanical, no LLM tokens, ~70s total).
2. **Three targeted head-to-head comparisons** against `pdftotext -layout` on files whose correct
   content was already independently known from prior verified digest/gate work this session:
   - Anderson & Lo 2011 (JAS 58:167) — an old Acrobat-Paper-Capture OCR scan.
   - Vaquero 2013 Table 4.1 — the exact Saturn-Titan IC/eigenvalue table `#765` confirmed to
     near-machine precision this session.
   - Anderson & Lo 2010 (JGCD 33:1899) — a two-column AIAA layout with an interleaved rotated
     "Downloaded by..." watermark, a case where `pdftotext -layout` is known to struggle.

## Findings

**1. Corpus-wide false-confidence empty-output rate: ~10%.**
- 41/251 (16%) produced empty/near-empty (<50 bytes) stdout on page 1.
- 24/251 (9.6%) did so while self-reporting `pdf_type=text_based confidence=1.0` — i.e. the tool
  asserted full confidence and silently returned nothing. This is a **false-confidence failure**,
  strictly worse than the 16/251 that honestly self-flagged as `image_based`/`scanned` (a
  legitimate signal to fall back to OCR).
- No single producer/encoder pattern predicts it — it hit both old Acrobat-scan PDFs (`Adobe
  Acrobat 9.0 Paper Capture Plug-in`) and modern born-digital LaTeX PDFs (`pdfTeX-1.40.24`,
  `pdfTeX-1.40.25`) alike.
- Affected files include foundational papers this project's whole genome depends on: **Howell
  1984** (halo orbits), **Richardson 1980** (analytic construction of periodic orbits about
  collinear points), Russell 2004 dissertation, Vallado 1991.

**2. Fidelity, where it DOES produce output: demonstrated silent corruption of exactly the
content this project's sourced-value discipline cares most about.**
- Vaquero 2013 Table 4.1 (nondimensional/dimensional IC, period, eigenvalue table): pdf-inspect
  **dropped the `×10^exponent` notation** on every numeric row (`1.25869×10⁶` rendered as
  `1.25869 10 6` — ambiguous, silently wrong if read literally) and **dropped an `fi` ligature**
  (`identified` → `identied`) — a systematic per-font glyph-mapping loss, not a one-off typo,
  meaning a `grep` for the correct word would silently miss the occurrence. `pdftotext -layout`
  reproduced this same table with zero corruption (verified against it directly).
- Anderson & Lo 2010 Eq. 1–2 (equations of motion): pdf-inspect **silently deleted the mass-ratio
  symbol μ** from running prose ("the mass of the smaller body ... is , and the larger body ...
  has mass 1") and mangled the equation layout itself beyond use (dropped signs, subscripts,
  coefficients).

**3. One genuine, narrow positive.** On the SAME Anderson & Lo 2010 page — a two-column layout
with a rotated watermark interleaved into the column text — pdf-inspect's position-aware reading
order produced substantially more readable prose than `pdftotext -layout` (which interleaved the
watermark mid-sentence and produced large whitespace-gutter artifacts). `is_complex_layout=True`
was correctly self-reported for this file. This confirms the tool has real value as a fast
prose-navigation aid on complex layouts — just not as a source of truth for any sourced value.

## Verdict

**Do not use `pdf-inspect` to reprocess or "improve" any existing digest in this corpus.** The
~10% false-confidence empty-output rate plus demonstrated silent corruption of exponents,
symbols, and ligatures make it unsafe as a digest source — errors of exactly this kind would be
invisible without a line-by-line cross-check against the original, which defeats the purpose of
using a faster tool at all. This would directly violate
[[feedback_golden_tests_sourced_only]]/[[feedback_published_rounded_values_are_display]]-style
sourced-value discipline if any digest content were updated from its output.

**Narrower, defensible future use:** fast prose-navigation/skimming on genuinely complex
multi-column layouts (`is_complex_layout=True` AND non-trivial stdout length) to locate where a
topic is discussed — never as the source of a sourced numeric/equation/table value, which must
still go through `pdftotext -layout`, `ocrmypdf`, or the existing selective-vision hybrid per
[[feedback_corpus_document_policy]].

**Scope note:** this task evaluated the *instrument*, not the *content* of existing digests — no
digest was read/audited for correctness as part of this task. A digest-content audit (independent
of this tool) would be a separate, dispatchable task if wanted.

**Machine-wide memory updated:** `~/.claude/memory/pdf_inspector_tool.md` (outside this repo,
applies to any project on this Mac) — description field corrected (no longer claims unqualified
"good for prose/text-heavy PDFs"), full findings appended.
