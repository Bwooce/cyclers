# A data table transcribed visually from an image-only scan carried 27 wrong cells into a positive control

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: yes (2026-10-05-twobody-gen-opus-visual-transcription-used-as-control.md, same incident seen from the corpus side; this entry adds the "second witness" rule for the corpus policy)

What happened: data/sources/hollister-menning-1970-table3.yaml, the #942 control, was transcribed from the old image-only scan. The clean text-layer copy showed 18 turn angles (3 read as 5), 6 Rmin values and 3 V_r values wrong. Many failed a basic physics check: the periapsis from V_r and the turn angle against the printed Rmin.
Workaround: diffed against the new text layer, settled each cell on 400 dpi images, and reported the list to the owner. Details: docs/notes/2026-10-05-hollister-menning-1970-table3-recheck.md.
Suggested fix: any table transcribed into data/sources needs a second witness (another copy, OCR, or a physics consistency check) before use as a control. Add this rule to the corpus policy.
