# Hollister & Menning 1970 Table 3: text-layer swap and transcription recheck (#960)

**Source:** W. M. Hollister & M. D. Menning, "Periodic Swing-By Orbits between Earth and Venus",
J. Spacecraft and Rockets 7(10):1193-1199 (1970), doi 10.2514/3.30134.

**The swap (2026-10-05):** the held image-only scan
`cyclers_pdf/papers/hollister-menning-1970-periodic-swingby-earth-venus-JSR-7-10.pdf` (md5 40dca8db...) was
replaced under the same filename.
- The replacement has a clean typeset copy with a text layer (md5 4f1dbeef3ed43d77bffc2aa6b73ce142;
  cyclers_pdf commit 6364264).
- The old scan stays in git history, and the OCR `.txt` beside it is kept.
- Table 3 is on PDF pp.4-6 = journal pp.1196-1198.

## Method

1. I parsed the new text layer's Table 3 (pdftotext -layout) and diffed it against
   `data/sources/hollister-menning-1970-table3.yaml` (twobody-gen-opus's `#942` transcription, made visually
   from the old scan).
2. I settled every disagreement two ways:
   - by reading the cell on a 400 dpi page image, and
   - where the difference was large, by a physics check: the periapsis from the row's own V_r and theta,
     r_p = mu/V^2 (1/sin(theta/2) - 1), must match the printed Rmin.
3. The old OCR `.txt` was too garbled to serve as a third witness.

## Result: 27 YAML cells differ from the print

Rows are counted 1-26 within each orbit, in YAML order. Format: YAML value -> printed value.

| Orbit, row(s) | Field | YAML -> printed | How confirmed |
|---|---|---|---|
| 3, rows 8-9 | turn | 30.1 -> 50.1 | physics |
| 3, row 10 | turn | 34.3 -> 54.3 | physics |
| 3, rows 18-20 | turn | 37.2 -> 57.2 | physics |
| 7, rows 13-15 | turn | 38.4 -> 58.4 | physics |
| 9, row 18 | turn | 37.1 -> 57.1 | physics |
| 10, rows 16-17 | turn | 54.8 -> 84.8 | physics |
| 15, row 12 | turn | 30.7 -> 50.7 | physics, image |
| 12, row 3 | turn | 17.6 -> 19.9 | physics |
| 4, row 9 | turn | 13.5 -> 13.8 | image |
| 5, row 7 | turn | 15.3 -> 15.5 | image |
| 12, row 13 | turn | 12.0 -> 12.1 | image |
| 14, row 15 | turn | 25.3 -> 25.5 | image |
| 1, rows 6-7 | Rmin | 1.37 -> 1.57 | physics |
| 2, row 15 | Rmin | 3.85 -> 3.83 | image |
| 7, row 17 | Rmin | 5.15 -> 5.45 | physics |
| 8, row 3 | Rmin | 5.31 -> 5.34 | image |
| 14, row 11 | Rmin | 4.39 -> 4.59 | physics |
| 7, rows 23-25 | V_r | 0.173 -> 0.175 | image |

That is 18 turn-angle cells, 6 Rmin cells and 3 V_r cells.

Unchanged:
- The two existing loader fixes match the print: orbit 1 row 12 is planet E, and orbit 6 row 3 is 995.
- Orbit 15 rows 18-26 in the YAML match the image.

**Menning 1968 thesis cross-check** (`docs/notes/2026-10-05-digest-menning-1968-mit-thesis-earth-venus-periodic-orbits.md`):
- It confirms all of the above.
- It fixes one further JSR print error: orbit 13 row 1 V_r is 0.124 in the thesis, against 0.129 in JSR.
- It leaves orbit 5 row 12 Rmin open: 5.08 in the thesis, 5.03 in JSR, about 5.05 by physics.

## Catalogue check (read only)

- Rows `hollister-menning-1970-ev-orbit-01..15` store only the maximum Earth and Venus v_inf
  (max V_r x 29.785 km/s).
- I recomputed all 15 pairs from the corrected table, including the orbit 7 V_r change. All 15 match the
  catalogue. There is no catalogue mismatch.

## Disposition

- I sent the corrections to twobody-gen-opus, owner of the YAML, and to the lead.
- I did not edit the YAML.
