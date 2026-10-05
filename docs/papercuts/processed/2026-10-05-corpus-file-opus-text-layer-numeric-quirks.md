# Publisher text layers silently corrupt table numbers (lost minus signs, split digits, misread labels)

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: Russell & Strange 2009's text layer drops the minus sign on six Table 6 values; the page image shows them. Henon & Guyot 1970 splits digit groups ("0.310 73"), prints "O." for "0.", and turns family labels l1, l'1, l2, h26 into "II", "1'1", "12", "ha6".
Workaround: custom joins in a scratch parser, then a check against page images (row counts plus spot values).
Suggested fix: add to the corpus policy that any number taken from a text layer into a golden or a catalogue field gets a page-image check, as the policy already requires for OCR.

Disposition (main, 2026-10-05): Backlog: #964 (page-image check for any text-layer number used as a golden or catalogue value).
