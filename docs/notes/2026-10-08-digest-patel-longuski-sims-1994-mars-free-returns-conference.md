# Digest: Patel, Longuski & Sims 1994, "Mars Free Return Trajectories" (AIAA 94-3766): diff against JSR 1998

M. R. Patel, J. M. Longuski & J. A. Sims, "Mars Free Return Trajectories", AIAA/AAS Astrodynamics
Conference, Scottsdale AZ, Aug 1994, paper **AIAA-94-3766-CP** (printed top right of p.1), proceedings
pp. 473-481. **DOI 10.2514/6.1994-3766** (Crossref: title, the three authors, "Astrodynamics Conference",
issued 1994-08-01).

- Source file: upload `c376598e-Mars_Free_Return_Trajectories.pdf`, 9 pp, md5 af2c4391293ab4a205f93ea7432f38c6,
  Ghostscript scan with a poor text layer.
- Filed (OCR copy, was `patel1994-ocr.pdf`) (`ocrmypdf --force-ocr -l eng`, PDF/A-2b, 9 pp), md5
  0730523e58aceb4880d9397edf88e3ca. Render check, pages 1, 5 and 7 at 100 dpi: same pixel size, mean
  grey difference 4.1, 7.8 and 6.7 of 255 (image recompression). A 300-dpi look at Table 2 in the
  copy is sharp and complete.
- Proposed filename: `patel-longuski-sims-1994-mars-free-return-trajectories-aiaa-94-3766.pdf`
  (or `...-aiaa-94-3766-doi-10.2514-6.1994-3766.pdf` to match the JSR file's DOI style).
- Compared against the held JSR 35(3):350-354 (1998) version (doi 10.2514/2.3333) and its digest
  `docs/notes/2026-10-05-digest-patel-longuski-sims-1998-mars-free-return-trajectories.md`. I read
  both on page images.

## 0. Verdict

The conference version has **no table that the journal lacks**. Both versions have the same three
tables (legend, escalator, Henon) with **the same numbers in every cell** (Tables 2 and 3 checked
cell by cell; see `witness-comparison.tsv`). The extra value of the 1994 version is **seven figures
that the journal cut**: Mars-arrival V_inf, outbound and inbound orbit periods, the 2015-2017
launch/arrival plot and Mars flyby altitude. The journal adds context, a nomenclature, two derived
numbers and a new conclusion. It is a "keep and file" item, low priority for gates. It gives no new
positive control beyond the JSR version.

`patel-1994-tables.yaml` holds Tables 2 and 3. They are not missing from the journal. I transcribed
them because the JSR digest left the down-escalator columns of Table 2 unread. All 21 rows are now
transcribed with two witnesses, and the 1998 page agrees.

## 1. Figure map (1994 -> 1998)

| 1994 | Content | 1998 |
|---|---|---|
| Fig. 1 | TOF vs launch date, 1995-2020, V_inf 4-8 km/s, L/D 950101 to 200101 by 15 days, TFMAX 1300 d, ALTMIN -1500 km | Fig. 1 |
| Fig. 2 | Arrival vs launch date, 2000-2002 (L/D 401 to 20401 by 5 days) | Fig. 2 |
| **Fig. 3** | **Arrival vs launch date, 2015-2017 (L/D 150301 to 170301 by 5 days)** | **cut** |
| **Fig. 4** | **V_inf at Mars vs launch date, 2000-2002 (axis 2-18 km/s)** | **cut** |
| Fig. 5 | Conic, opposition-type free return, TOF 1.4 yr, 2000 Nov 7 | Fig. 3 |
| **Fig. 6** | **V_inf at Mars, 2015-2017** | **cut** |
| **Fig. 7** | **Outbound orbit period vs launch date, 2015-2017 (1-5 yr)** | **cut** |
| **Fig. 8** | **Inbound orbit period, 2015-2017** | **cut** |
| **Fig. 9** | **Outbound orbit period, 2000-2002** | **cut** |
| **Fig. 10** | **Inbound orbit period, 2000-2002** | **cut** |
| Fig. 11 | Free-return options E1/E2, M1/M2 | Fig. 4 |
| Fig. 12 | 1.5-yr period, TOF 3.0 yr (E2-M1-E2), 2001 Jan 1 | Fig. 5 |
| Fig. 13 | 3.0-yr period, TOF 3.0 yr, 2001 Mar 22 | Fig. 6 |
| **Fig. 14** | **Mars flyby altitude (Mars radii, 0-100) vs launch date, 2000-2002** | **cut** |
| Fig. 15 | Up escalator conic, 2001 Feb 5 | Fig. 7 |
| Fig. 16 | Down escalator conic, 2006 Mar 12 | Fig. 8 |
| Fig. 17 | Collision-orbit geometry (P, Q, tau) | Fig. 9 |

The cut figures are STOUR scatter plots. They show trends, not data-grade points. Uses:
- Fig. 4/6 show the Mars V_inf "peak" near March 2001 launches (up to about 16-18 km/s on the axis).
- Fig. 7-10 show that most 3-yr TOF returns have outbound and inbound periods near 1.5 yr.
- Fig. 14 is the plot behind the statement "the orbits which have large Mars flyby altitudes are in
  fact the collision orbits predicted by Henon". The JSR keeps the sentence but drops the figure.

## 2. Text only in the 1994 version

1. "Refined searches for the 9 launch opportunities between 2000 and 2018 were conducted" (p.474).
   Detailed plots are given for 2015-2017 (Figs. 6-8). The JSR does not mention the nine refined searches.
2. The 3-yr period orbits near the March 2001 V_inf peak have "very high arrival V_inf's
   (14 < V_inf < 18 km/s)" and large launch energies (p.477, from Fig. 4). The JSR has no such range.
3. The fast trajectories are linked to Figs. 2-4 and the V_inf at Mars plots ("high Mars arrival V_inf's
   (Figure 4)"). The JSR says "high Mars arrival V_inf" without a plot.
4. The conclusion ends with future work: "Similar analyses can be performed ... for Mars Free Returns
   which include Venus as a gravity-assist body". This is the lead-in to Okutsu & Longuski 2002
   (held), which cites the JSR version as its ref. 16.
5. Table 1 legend, Search Event No.: the 1994 text says that event 3 "is for Mars Arrival (for the
   path 3 4 3)". That is an error, because event 3 in path 3 4 3 is the Earth return. The JSR
   corrects it: "the TOF in Fig. 1 corresponds to the third event in the sequence, namely Earth (3)
   arrival."
6. Affiliations: all three authors are at Purdue, and Patel and Sims are graduate students. In 1998,
   Patel is at Lockheed Martin and Sims at JPL. The copyright line is 1994 (the scan prints "@ 1004"
   in the old text layer; the page image reads 1994).

## 3. Text only in the 1998 version

- Nomenclature block. Wider introduction (refs. 1-5, 9: Braun, Striepe, Desai, Tartabini, Walberg).
- The synodic derivation: Mars period about 1 7/8 yr, synodic period about 2 1/7 yr, 1/7 of a circle
  (51.4 deg) per synodic period, 7 synodic periods for about 15 yr. The 1994 version only says
  "approximately every 15 years".
- The escalator numbers 2.12 yr (Earth-1 to Earth-3, up) and 7 months 20 days = 0.63 yr (up vs down
  offset), linked to the 0.6-yr gap in Fig. 1. The 1994 version has neither number.
- Prado & Broucke 1994 (ref. 26) on Henon's problem via Lambert.
- The conclusion paragraph on fast 1.4-yr free returns in 2015 and 2017, "a timely opportunity".
- JPL acknowledgment. References grow from 17 to 26. Every 1994 reference is in the 1998 list. The
  Uranus-Neptune-Pluto paper is cited as IAA-L-0408 (1994) in 1994 and as Acta Astronautica 36(2) 1995
  in 1998.

## 4. Tables

- **Table 1** (legend): same fields. The 1994 wording is longer, and it has the Search Event error in section 2 item 5.
- **Table 2** (escalators, from Byrnes et al. 1993): 21 rows x 4 value columns, 84 cells.
  - Every value is the same in both versions (dates by day, month and year; V_inf and DeltaV to 2 decimals).
  - Differences are only in typography. 1994 uses short months ("Sep", "Jul") and prints "May 1,1997"
    with no space. One cell differs in print: the down-escalator DeltaV for the Maneuver row after
    Mars-8 is blank in 1994 and a dash in 1998.
  - Caption: 1994 "Up/Down Escalator Orbits (From Byrnes et al.[4])"; 1998 "Up/down-escalator orbits
    (from Ref. 10)".
- **Table 3** (Henon): 3 rows x 5 columns. Same 15 numbers. Caption: 1994 "Summary of Henon's Results";
  1998 "Excerpt from Henon's tables". The 1994 version labels the a column "Semi-Major Axis (AU)".
- Two-witness check (`witness-comparison.tsv`, built by `build_witness.py` from `ocr-raw/`):
  - Witness A: my 400-dpi image reading.
  - Witness B: tesseract psm 6. On the numeric columns it used a digit whitelist; on dates and labels it was a full-table pass.
  - B was compared on 87 cells (dash and blank cells excluded) and agrees on all 87. The full pass
    has glyph noise only on the encounter labels ("Mars—2", "Mars-—i0", "Mars~12"). The decided values
    are the same.
  - Column C (the 1998 page image) gives the same value on 98 of 99 cells. The one exception is the
    blank/dash cell above.

## 5. Gate relevance and positive controls

The same as the JSR digest. Nothing new is gated by the 1994 version. The cut Fig. 4/6 (Mars V_inf by
launch date) and Fig. 14 (flyby altitude) can give a qualitative check on a free-return generator's
Mars V_inf peak near March 2001. They are scatter plots, so use them for shape only, not for numbers.

## 6. Citation mining

There are no 1994-only references (see section 3). Nothing new to add to the wanted list. Wolf 1991
(AAS 91-123) and Howell 1987 remain as listed in the JSR digest.

*Filed as `cyclers_pdf/papers/patel-longuski-sims-1994-mars-free-return-trajectories-aiaa-94-3766-doi-10.2514-6.1994-3766.pdf` (the force-OCR copy). Tables 2-3: `data/sources/patel-longuski-sims-1994-tables.yaml`; the witness comparison is filed beside the PDF.*
