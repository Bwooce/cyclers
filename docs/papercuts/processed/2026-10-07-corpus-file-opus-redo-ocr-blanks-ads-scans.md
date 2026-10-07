# `ocrmypdf --redo-ocr` silently blanks ADS JBIG2 scans

- Date: 2026-10-07
- Agent: corpus-file-opus
- Seen before: no

What happened: on five ADS article scans (jb2pdf.py, one JBIG2 stencil mask per page), `ocrmypdf --redo-ocr`
exited 0 and wrote a good text layer, but the output pages render blank: the stencil images are dropped.
Two digest subagents caught it by rendering pages. NTRS scans were not affected.
Workaround: re-ran with `ocrmypdf --force-ocr` (about 1 MB per paper) and checked every page renders against the original.
Suggested fix: the corpus policy now says to use `--force-ocr` for ADS scans and to render-check one page of
any OCR output before filing (`docs/notes/corpus-document-policy.md`).

Disposition: Promoted: already in docs/notes/corpus-document-policy.md (--force-ocr for ADS scans; render-check one page before filing).
