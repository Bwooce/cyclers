# Tesseract (via ocrmypdf) misreads the typewriter digit 4 as h, L or u, so OCRed tables from 1960s theses are unusable as-is

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: in the Menning 1968 MIT thesis, `ocrmypdf --skip-text` produced table rows like "Lush 0.148" and "hgok" for dates such as 4474 and 4924. It also dropped some rows, which shifted the row-by-row alignment. Double-spaced typewritten body text came out scrambled.
Workaround: rendered the pages at 400 dpi in grey, ran `tesseract --psm 6 -c tessedit_char_whitelist="EV0123456789. "`, and compared only well-formed cells against a trusted table. Shifted or disputed pages were read by vision.
Suggested fix: add the whitelist-and-400-dpi recipe and a "never trust OCRed 4s" warning to docs/notes/corpus-document-policy.md, sec. 1.
