# Russian-language corpus files need tools this Mac lacks (rus tessdata, solid-RAR unpacker)

- Date: 2026-10-06
- Agent: corpus-file-opus
- Seen before: no

What happened (batch 19):
- The KIAM preprint PDFs (2007-34, 2008-16, 2010-64) have Type 3 Cyrillic fonts, so `pdftotext`
  gives mojibake. The Homebrew tesseract has only eng and osd.
- With `TESSDATA_PREFIX` pointed at a scratch folder holding a downloaded `rus.traineddata`,
  `ocrmypdf` failed with "read_params_file: Can't open hocr" until `configs/` was symlinked in.
- The KIAM 2005-10 zip has solid RAR archives. The installed `7z` reports "Unsupported Method" and
  `bsdtar` reports "RAR solid archive support unavailable".

Workaround:
- Copied eng/osd into the scratch tessdata, linked `/opt/homebrew/share/tessdata/configs`, and fetched
  `rus.traineddata` from the tessdata_best repository.
- `brew install unar` (now installed system-wide).

Suggested fix: `brew install tesseract-lang` (all languages, about 650 MB), or a note in
`docs/notes/corpus-document-policy.md` giving the TESSDATA_PREFIX recipe and `unar` for old RARs.

Disposition: Promoted: TESSDATA_PREFIX + rus.traineddata recipe and `unar` go in docs/notes/corpus-document-policy.md (corpus-file-opus, dispatched 2026-10-07). `unar` is installed.
