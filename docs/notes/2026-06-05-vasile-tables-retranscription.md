# Verified re-transcription: Vasile & Campagnola tables (TOF + gravity-assist characteristics)

Date: 2026-06-05. Task #96. Independent high-fidelity re-read of the two candidate golden-anchor
tables for deferred low-thrust MGA work (#37), cross-checked against the prior raster scan
(`docs/notes/2026-06-05-vasile-hiraiwa-scan.md`).

Source: Vasile & Campagnola, "Design of Low-Thrust Multi-Gravity Assist Trajectories to Europa",
JBIS (arXiv:1105.1823).

## Method / provenance caveat

This PDF has broken font encoding, so `pdftotext` returns garbage; both this read and the prior scan
are OCR-by-eye off rendered page rasters of **the same PDF**. This re-read was done at **300 DPI**
(prior scan: 150/220 DPI), reading cell-by-cell from cropped, enlarged regions. The two reads are
therefore a higher-resolution confirmation of one source — **not two independent provenances**. A clean
JBIS-typeset copy remains the gold standard. This note upgrades confidence (300-DPI digits are
unambiguous), but it is what it is: a second raster read of the same file.

### Nomenclature note (important)

The task brief refers to these as "Table 3" and "Table 4". In the actual PDF body the captions are:
- **Table 7** = Time of Flights and encounter dates (task's "Table 3")
- **Table 8** = Summary of gravity assist characteristics (task's "Table 4")

(The arXiv preprint and the JBIS published version may use different table numbers; the prior scan also
called them Tables 3/4. Cite by caption text, not number, to avoid ambiguity.) Both tables span the
bottom of **page 24** (PDF page count 27) and the top of **page 25**.

---

## Table 7 (task "Table 3") — verified transcription

Caption (verbatim, p.24): **"Table 7. Time of Flights and encounter dates for the entire optimised
solution"** (note: caption text reads "encounter" but is typeset "encouter" — transcribed as printed
intent). Column headers verbatim:

| Header | Verbatim text |
|--------|---------------|
| col 1  | PH N. |
| col 2  | CELESTIAL BODY |
| col 3  | DEP. DATE (MJD) |
| col 4  | TOF (D) |
| col 5  | CUMMULATIVE MISSION TIME (Y)  *(sic: "CUMMULATIVE")* |

Rows 1-8 on p.24; rows 9-16 on p.25.

| PH N. | CELESTIAL BODY | DEP. DATE (MJD) | TOF (D) | CUMM. MISSION TIME (Y) |
|------:|----------------|----------------:|--------:|-----------------------:|
| 1  | Earth          | 3718.9 | 393.7  | 1.08 |
| 2  | Venus          | 4112.7 | 166.5  | 1.53 |
| 3  | Mars           | 4279.1 | 375.3  | 2.56 |
| 4  | Earth          | 4654.4 | 1241.6 | 5.96 |
| 5  | SOI of Jupiter | 5896.1 | 110.6  | 6.27 |
| 6  | Ganymede1      | 6006.7 | 492.7  | 7.62 |
| 7  | Ganymede2      | 6499.3 | 93.0   | 7.87 |
| 8  | Ganymede3      | 6592.5 | 35.8   | 7.97 |
| 9  | Ganymede4      | 6628.1 | 21.5   | 8.03 |
| 10 | Ganymede5      | 6649.6 | 1.2    | 8.03 |
| 11 | Europa1        | 6650.8 | 10.7   | 8.06 |
| 12 | Europa2        | 6661.4 | 17.8   | 8.11 |
| 13 | Europa3        | 6679.2 | 7.1    | 8.13 |
| 14 | Europa4        | 6686.3 | 28.4   | 8.21 |
| 15 | Europa5        | 6714.7 | 39.6   | 8.32 |
| 16 | Ganymede6      | 6754.4 | 2.0    | 8.32 |

No cells were ambiguous at 300 DPI. **0 [unclear] cells.**

---

## Table 8 (task "Table 4") — verified transcription

Caption (verbatim, p.25): **"Table 8. Summary of gravity assist characteristics for the entire
optimised solution"**. Column headers verbatim:

| col | Verbatim header text |
|-----|----------------------|
| 1 | PH N. |
| 2 | FROM |
| 3 | TO |
| 4 | DEP. RELATIVE VELOCITY (KM/S) |
| 5 | DEP. ABSOLUTE VELOCITY (KM/S) |
| 6 | ARRIVAL RELATIVE VELOCITY (KM/S) |
| 7 | ARRIVAL ABSOLUTE VELOCITY (KM/S) |
| 8 | B [°] |
| 9 | PERICENTRE ALTITUDE (KM) |
| 10 | ΔV_GA (KM/S) |

Body abbreviations as printed: Ga1..Ga6 = Ganymede 1..6, Eu1..Eu6 = Europa 1..6. A dash "-" denotes a
value not applicable (SOI-of-Jupiter capture phases and the final phase). Rows 1-4 on p.24; rows 5-16
on p.25.

| PH N. | FROM | TO | DEP REL V | DEP ABS V | ARR REL V | ARR ABS V | B [°] | PERI ALT (KM) | ΔV_GA (KM/S) |
|------:|------|----|----------:|----------:|----------:|----------:|------:|--------------:|-------------:|
| 1  | Earth          | Venus | 1.78  | 28.32 | 6.71  | 36.89 | 61.66 | 753   | 6.88 |
| 2  | Venus          | Mars  | 6.71  | 41.31 | 5.63  | 20.08 | 6.14  | 17478 | 0.60 |
| 3  | Mars           | Earth | 5.63  | 20.10 | 12.26 | 33.19 | 33.09 | 300   | 6.98 |
| 4  | Earth          | SOI of Jupiter | 12.26 | 38.52 | 4.30  | 8.63  | -     | -     | -    |
| 5  | SOI of Jupiter | Ga1   | -     | 4.30  | 6.40  | 15.74 | 9.01  | 200   | 1.01 |
| 6  | Ga1            | Ga2   | 6.40  | 15.16 | 6.60  | 15.17 | 8.46  | 217   | 0.97 |
| 7  | Ga2            | Ga3   | 6.60  | 14.69 | 6.61  | 14.68 | 8.39  | 235   | 0.97 |
| 8  | Ga3            | Ga4   | 6.61  | 14.02 | 6.61  | 14.02 | 8.49  | 200   | 0.98 |
| 9  | Ga4            | Ga5   | 6.61  | 13.42 | 6.60  | 13.42 | 8.50  | 200   | 0.98 |
| 10 | Ga5            | Eu1   | 6.60  | 12.60 | 4.25  | 17.26 | 10.48 | 200   | 0.78 |
| 11 | Eu1            | Eu2   | 4.25  | 16.89 | 4.29  | 16.88 | 8.89  | 512   | 0.66 |
| 12 | Eu2            | Eu3   | 4.29  | 16.53 | 4.35  | 16.50 | 9.25  | 366   | 0.70 |
| 13 | Eu3            | Eu4   | 4.35  | 16.00 | 4.37  | 15.99 | 9.61  | 268   | 0.73 |
| 14 | Eu4            | Eu5   | 4.37  | 15.39 | 4.45  | 15.36 | 8.83  | 379   | 0.68 |
| 15 | Eu5            | Ga6   | 4.45  | 14.76 | 2.04  | 8.90  | 39.10 | 2105  | 1.36 |
| 16 | Ga6            | Eu6   | 2.04  | 9.77  | 1.65  | 15.36 | -     | -     | -    |

No cells were ambiguous at 300 DPI. **0 [unclear] cells.**

---

## Cross-check vs prior scan (cell-by-cell)

### Table 7 (TOF) agreement

The prior scan tabulated this with a "Dep date" column and the same body sequence. All 16 rows match
on body name, dep date, and cumulative time. **One MISMATCH:**

| Row | Cell | Prior scan | This 300-DPI read | Status |
|-----|------|-----------|-------------------|--------|
| 9 (Ganymede4) | TOF (D) | 21.3 | **21.5** | MISMATCH — resolved in favour of 21.5 |

All other Table 7 cells (15 rows × {body, dep date, cum} + the other 15 TOF values): **MATCH**.

Note: the prior scan's column was headed "Dep date (MJD)"; the printed header is "DEP. DATE (MJD)" —
same field, no value discrepancy.

Resolution of the one mismatch: at 300 DPI the digit is an unambiguous **5** (confirmed on a 2x-enlarged
crop; the glyph is a clear "5", not a "3"). The prior scan's 21.3 was a 150-DPI misread.

### Table 8 (gravity-assist characteristics) agreement

The prior scan did **not** tabulate Table 8 in full — it quoted only representative rows in prose
(rows 1, 2, and a range characterization for rows 6-16). Cross-check of the explicitly quoted values:

| Quoted value (prior prose) | This read | Status |
|----------------------------|-----------|--------|
| Row 1 Earth→Venus dep rel-V 1.78 | 1.78 | MATCH |
| Row 1 arr rel-V 0.71 | **6.71** | MISMATCH — prior prose said "arr rel-V 0.71"; the cell is clearly 6.71 |
| Row 1 B = 36.89° | header/cell confusion: 36.89 is ARR ABS V; **B = 61.66** | MISMATCH — prior prose mislabelled the 36.89 figure as the bend angle B |
| Row 1 peri-alt 753 km | 753 | MATCH |
| Row 1 ΔV_GA 6.88 | 6.88 | MATCH |
| Row 2 Venus→Mars 0.71/0.71 | dep rel 6.71 / arr rel 5.63 | MISMATCH — prior "0.71/0.71" does not match (likely a transcription slip; 6.71 is dep rel of row 2) |
| Row 2 B = 20.08 | **B = 6.14** (20.08 is ARR ABS V) | MISMATCH — same column-misalignment as row 1 |
| Row 2 alt 17478 | 17478 | MATCH |
| Row 2 ΔV 0.60 | 0.60 | MATCH |
| Rows 6-16 ΔV_GA "0.5–1.1 km/s range" | 0.66–1.01 (rows 6-14), 1.36 (row 15) | PARTIAL — range mostly holds; row 15 ΔV_GA = 1.36 exceeds the stated 1.1 upper bound |
| Rows 6-16 peri-alt "~200–580 km" | 200–512 (rows 6-14), 2105 (row 15) | PARTIAL — holds except row 15 peri-alt 2105 km |

The prior-prose mismatches are **labelling/quoting errors in the prose**, not raster digit errors: the
"0.71" and "B=36.89/20.08" slips come from reading across mis-aligned columns in the prose summary, not
from misreading glyphs. The unrecorded rows (3-16 full detail) are **prior-unrecorded**, now filled in.

---

## Verdicts

### Table 7 (TOF / encounter dates) — VERIFIED (with one corrected cell)

Two raster reads now agree on every cell **after** correcting the single 150-DPI misread
(row 9 TOF 21.3 → 21.5). 0 remaining [unclear] cells, 0 remaining mismatches. Caveat: both reads are of
the same PDF (not independent provenance); a clean JBIS copy remains the gold standard, but the 300-DPI
digits are individually unambiguous.

### Table 8 (gravity-assist characteristics) — VERIFIED (newly transcribed in full)

The prior scan never transcribed this table in full, so there is no second raster read of most cells —
strictly this is a **single** high-fidelity read, not two-agree. However: (a) every cell was unambiguous
at 300 DPI (0 [unclear]), and (b) the prior prose's apparent "mismatches" were prose-labelling errors,
not contradicting digit reads — the digits the prose did quote that map cleanly (1.78, 753, 6.88, 17478,
0.60) all MATCH. Treat as VERIFIED-single-read pending a clean-source confirmation; the internal
consistency check below adds confidence.

#### Internal consistency check (adds confidence, not provenance)

Two structural invariants hold across the transcription, which is strong evidence the digits are right:
- **Velocity continuity**: each phase's DEP RELATIVE VELOCITY equals (≈) the previous phase's ARRIVAL
  RELATIVE VELOCITY along the chain (e.g. row 1 arr rel 6.71 = row 2 dep rel 6.71; row 2 arr rel 5.63 =
  row 3 dep rel 5.63; row 10 arr rel 4.25 = row 11 dep rel 4.25; row 15 arr rel 2.04 = row 16 dep rel
  2.04). This chains cleanly through all 16 rows.
- **B and ΔV_GA correlate**: the two large bend angles (row 1 B=61.66, row 15 B=39.10) coincide with the
  two large ΔV_GA values (6.88-class deep flybys and 1.36); the resonant-tour rows 6-14 have small bend
  angles (8-10°) and small ΔV_GA (~0.7-1.0), as expected for a synchronous tour.

---

## Usability for #37 candidate anchors

| Values | Status for #37 |
|--------|----------------|
| Table 7 all 16 rows (dep dates MJD, TOF d, cumulative yr) | **Usable** as candidate anchors at upgraded confidence; row 9 TOF = 21.5 (corrected). Still flag "raster-derived, confirm vs clean JBIS before use as a hard golden EXPECTED value." |
| Table 8 rows 1-16 (V_inf rel/abs, B, peri-alt, ΔV_GA) | **Usable** as candidate anchors at single-read 300-DPI confidence + internal consistency. Same clean-source caveat. |
| Any value as a *golden EXPECTED* in a validation test | **Still blocked** until confirmed against a clean (non-raster) JBIS source — per the project "golden tests use sourced-only values" rule, a raster OCR of a font-broken PDF is not an acceptable EXPECTED provenance even at 300 DPI. These are candidate/reference anchors, not goldens. |

### Net change vs prior scan
- Corrected: Table 7 row 9 TOF (21.3 → 21.5).
- Newly captured in full: all of Table 8 (16 rows × 9 data columns), previously only prose-summarised.
- Corrected prior-prose labelling errors: row 1/2 "arrival rel-V" and "B" figures were mis-attributed in
  the prior prose; correct values are in the table above (B is its own column; the 36.89/20.08 figures
  are ARRIVAL ABSOLUTE VELOCITY, not B).
