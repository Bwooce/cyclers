# pdftotext on AIAA PDFs prints thousands of "Syntax Warning: Badly formatted number" lines

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: pdftotext on the Hollister-Menning 1970 and Campagnola 2019 AIAA PDFs flooded stderr with tens of thousands of warnings. That truncated the useful tool output and hid an error.
Workaround: `2>/dev/null` on every pdftotext call.
Suggested fix: note it in the corpus policy's tool line, or wrap pdftotext in a helper that discards those warnings.
