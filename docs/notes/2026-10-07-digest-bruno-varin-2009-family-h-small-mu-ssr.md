# Digest: Bruno & Varin 2009, "Family h of periodic solutions of the restricted problem for small mu" (Sol. Syst. Res. 43:2) (#960 batch 31)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, RAS), Solar System Research 43(1):2-25 (2009),
doi 10.1134/S003809460901002X. English translation (Pleiades) of Astronomicheskii Vestnik 43(1):4-27 (2009).
Received 30 April 2008. RFBR grant 08-01-00082.
- Supplied by the owner as `93c99c80-bruno2009.pdf`: 24 pp., publisher text layer, md5 d333e174283aafe953d3fc3bde2c4652.
  Journal page = PDF page + 1.
- **Proposed corpus filename:**
  `bruno-varin-2009-family-h-periodic-solutions-restricted-problem-small-mu-sol-syst-res-43-2-doi-10.1134-s003809460901002x.pdf`
- **How I read it:**
  - I read the whole text layer.
  - **I read every cell of Tables 1-5 on 300 dpi page images**, zoomed to quarter-table crops. These are journal pp.3,
    9, 10, 11, 16 and 21. The text layer agrees with the images cell for cell. I used it only to locate cells and to
    build the YAML rows.
  - On the page images I also read: eq. (1)-(2) (p.2); the T~ definition and eq. (4) (p.6); the mu_J statement (p.13);
    the side values (pp.14, 21-22, 25).
  - I compared the text section by section with our translation of KIAM 67/2005 (`bruno-varin-2005b-...-en-translation.tex`).
  - I re-read the Table 1 cells that disagree with KIAM 34/2007 on the 34/2007 page image (PDF p.23) and on the JAMM
    71:933 page image (p.939).
  - Arithmetic and integration checks are in `checks_bv2009.py`. The output is in `checks_bv2009.out`.
  - The tables are transcribed in `bruno-varin-2009-family-h-small-mu-tables.yaml`.
- Wanted list: **row 17**, the Bruno-Varin KIAM 67 and 64/2005 family h tables.

## 0. Verdict

**This is the journal form of KIAM Preprint 67/2005, and it prints the tables that the 67 TeX source lacks.**
- It covers the same four mu values: 0, mu_J, 0.1 and 0.2.
- It has the same section plan (generating family; Sun-Jupiter; Earth-Moon note; mu = 0.1; mu = 0.2).
- The text is the same text, sentence by sentence, with the differences listed in sec. 1.
- The orbit counts are the same: 16, 35, 28 and 27.
- All 12-digit side values agree with 67, except one dropped digit (sec. 1).

**What it recovers for wanted row 17: the small-mu half, in full.**
- Table 1 (mu = 0, 16 orbits) now has all its columns. KIAM 34/2007 and JAMM lack x1(0), y2(0), the P2 entry velocity
  (v1, v2), Tr_v and w1. They also give T~ to only 2 decimals.
- Table 2 (mu_J, 35 orbits), Table 3 (the Table 1 / Table 2 map), Table 4 (mu = 0.1, 28 orbits) and Table 5 (mu = 0.2,
  27 orbits) are all new to the corpus.
- **Each orbit has a full symmetric initial state (x1(0), 0, 0, y2(0)) and its half-period crossing.** They can be
  integrated directly. I checked 12-14 rows per table by integration (sec. 3).
- **The big-mu half is not here.** Crossref lists the companion: Bruno & Varin, "Family h of periodic solutions of the
  restricted problem for big mu", Sol. Syst. Res. 43(2):158-177 (2009), doi 10.1134/S0038094609020099. It is not held.
  It is very likely the journal form of KIAM 64/2005, which would have Tables 6-8 (mu = 0.3, 0.4, 0.5).
  **Proposal:** row 17 drops 67/2005 and names this paper as the remaining item.

**Three conventions in the tables differ from what the text says.** Anyone who uses the numbers must know them:
1. **Period.** Page 6 says T~ = T/(2 pi). The printed T~ is **(1-mu) T/(2 pi)**, the 67/2005 definition. Integration
   reproduces the printed T~ to 1e-5 only with the (1-mu) factor. At mu = 0.2 the two conventions differ by 25%.
2. **Jacobi constant.** Page 6 says C = -2H. The printed C is **-2H + mu**, which equals the barycentric (Szebehely)
   C + mu(1-mu). This holds on every row of Tables 2, 4 and 5, to the 5-decimal rounding of the states.
   - This matches 64/2005's note that Hénon's C_H = C - 0.25 at mu = 0.5.
3. **mu_J.** Table 2 was computed at **mu = 0.00095**, as the p.13 text says ("on technical grounds"). It was not
   computed at the 0.00095388 of the section heading, nor at the abstract's 10^-3. Both the integrations and the w1
   column reject 0.00095388.

**Misprints (the printed values are kept in the YAML, each with a note):**
- **Table 1, row 7, C:** printed -1.787103; it should be -1.785103.
  - The row's (x1, y2) give -1.78508 and its P2 velocity gives -1.785106.
  - 34/2007 and JAMM print -1.785103.
- **Table 1, row 14, a~(T/2):** printed 5.88769. This disagrees with the row's own e~(T/2), which gives 5.8779.
  34/2007 and JAMM print 5.877.
- **Table 5, row 1, T~:** printed 0.31031; it should be 0.31301 (digits swapped). Integration gives 0.31301.
- **Table 5, row 5, T~:** printed 1.97805; it should be 0.97805. Integration gives 0.97803. The printed value would put
  the row outside the order between rows 4 and 6.
- **Side value, mu = 0.1, orbit 21:** a~(0) is printed -2.18406122002. 67/2005 has -2.184061220002, so the journal
  dropped a 0.

**What it gives the project:**
- **Integrable, checked symmetric periodic orbits of the retrograde family h (Broucke's A1)** at mu = 0.00095, 0.1 and
  0.2, from tiny retrograde circles about P1 out to multi-revolution orbits with T~ up to 5.5. Each table also gives
  the collision orbits and the stability intervals.
- **Positive controls** for any symmetric-family corrector or continuation code at large mu (`#944` X2, `#956` R9). The
  initial states close with corrections of about 1e-5 or less in y2(0).
- **Catalogue implication (PROPOSAL only):** none directly. Family h is a P1-centred retrograde family, not a cycler
  source. If the YAML goes into `data/sources/`, its header must state the two convention traps above, so that no
  code reads T~ or C at face value.

## 1. Identity with KIAM 67/2005

| Part | KIAM 67/2005 (TeX, our translation) | Sol. Syst. Res. 43:2 |
|---|---|---|
| Abstract | mu = 0, 10^-3, 0.1, 0.2 | same, sentence for sentence |
| Introduction | one paragraph ("Within the programme of [1]") | longer: restates the problem, Hamiltonian (1)-(2), symmetry (3), multiplicity, the four coordinate systems, eq. (4), the critical-orbit definition; cites Bruno & Varin 2008 (SSR 42) where 67 cites KIAM 10/2005 |
| T~ definition | (1-mu) T/(2 pi) (sec. 1.1) | "T~ def T/(2 pi)" (p.6) **but the tables use (1-mu) T/(2 pi)** (sec. 3) |
| sec. 1 generating family | 1.1-1.4; Figs. 1.1-1.8 | same text; Figs. 1-8; Table 1 printed |
| Table 1 T/2 columns | a~(T/2), e~(T/2) for orbits 4 and 11, otherwise (a~'_*, e~_*) by (1.1) | a~(T/2), e~(T/2) for all non-collision rows, plus a new w1(T/2) column; (1.1) and the asymptote values of 67 sec. 1.2 are dropped |
| P2 local system | (a~'_*, e~_*), eq. (1.1) and (2.2)-(2.3) | system IV (w1, y2), w1 = mu/(1-x1), eq. (6) |
| sec. 2 Sun-Jupiter | "mu = mu_J def 0.00095388" | heading 0.00095388, but text "mu_J def 0.00095 on technical grounds" |
| side values | x1(T/2) for 5, 6, 19, 20, 33, 34; Tr_v(12), Tr_v(26) | identical (image-checked) |
| sec. 2.5 Earth-Moon | Broucke 1968; Hénon [4, sec. 10.4.7] | same, without the Hénon section number; "as little as orbit 12" where 67 has "slightly beyond orbit 12" (translation) |
| sec. 3 mu = 0.1 | 28 orbits; orbit-10 x1(T/2), e~(0), **a~'_* = -0.0000201814**; a~(0) for 8, 9, 20, 21 | "about 28"; the a~'_* value is dropped; a~(0)(21) has one 0 fewer |
| sec. 4 mu = 0.2 | 27 orbits; 7 side values | "about 27"; side values identical |
| stability intervals | mu_J: 3; mu = 0.1: 3; mu = 0.2: 6 | the same intervals |
| Conclusion | one paragraph (perturbation grows with a~) | absent; replaced by "for a preliminary publication ... see the preprint" |
| References | 6 | 10 (adds Bruno 1994 English, Bruno-Varin 2006 CMDA, 2007 JAMM, 2008 SSR) |

So 67/2005 is the preliminary version and this paper is its journal form, as the paper itself says on p.25.
**The tables recovered are the journal's Tables 1-5.** They are probably the same as 67's absent tables. Two
differences are certain:
- the Table 1 T/2 columns follow a different convention (67's text describes a~'_*, e~_*);
- the w1 column is new.

## 2. The tables (YAML)

`bruno-varin-2009-family-h-small-mu-tables.yaml`:
- one key per table (`table_1` to `table_5`), plus `side_values`;
- every cell is a string, exactly as printed.

**Frame (eq. (1)-(2), p.2, image-checked):**
- Synodic frame with its origin at **P1** (mass 1-mu), and **P2 at (+1, 0)**.
- y1 = x1' - x2 and y2 = x2' + x1 are canonical momenta.
- H = (y1^2+y2^2)/2 + x2 y1 - x1 y2 - 1/r + mu(1/r + x1 - 1/r2). The +mu x1 term is the indirect term.

**I checked the frame two ways:**
- I integrated Table 4, orbit 1 for one time unit in this frame, and also in the standard barycentric CR3BP with
  x_b = x1 - mu. The two agree to 3e-12.
- Every symmetric orbit starts on x2 = y1 = 0, with x1(0) < 0. The orbits are retrograde about P1, on the side away
  from P2.

**Columns:**
- k;
- x1(0), y2(0);
- either v1(T/2), v2(T/2) (Table 1) or x1(T/2), y2(T/2) (Tables 2, 4, 5);
- T~ and C;
- the modified traces Tr~ and Tr~_v;
- a~(0), e~(0), a~(T/2), e~(T/2);
- w1.

**Astronomical coordinates at mu > 0.** The formula is not printed in this paper; it is in Bruno-Varin 2008, which is
not held. I found it from the 12-digit side values: x1 = a~(2-|e~|) and y2|y2| = (1-mu) e~/x1.
- With the (1-mu) factor, all 8 side-value orbits give the printed x1(0), y2(0) and C.
- Without it, y2 is off by 10% at mu = 0.2.

## 3. Checks (`checks_bv2009.py`, output `checks_bv2009.out`)

**Table 1 (mu = 0) against KIAM 34/2007 and JAMM Table 1:**
- C from (x1, y2), C = 3 - V^2 from the P2 entry velocity, and C from (a~(0), e~(0)) all agree with the printed C to
  within 3e-5 (5-decimal rounding) on 15 of 16 rows.
- **Row 7 is off by 2.0e-3**, the misprint in sec. 0.
- C, a~(0) and e~(0) agree with 34/2007 on every row except row 7's C.
- The 5-decimal T~ rounds to 34/2007's 2-decimal values, except that 34/2007 truncates rows 3, 10 and 14
  (1.99740 -> 1.99, 3.99917 -> 3.99, 4.05653 -> 4.05).
- **The T/2 columns, an open question in the JAMM digest sec. 2.1, are now explained.**
  - At T/2 the generating orbit is at P2.
  - If v1 = 0: y2(T/2) = v2 + 1.
  - If v1 ≠ 0, the mu -> 0 limit is a perpendicular crossing at the pericentre of the P2 flyby:
    - w1 = V^2 v1/(V-|v1|);
    - y2(T/2) = 1 + sgn(v2) sqrt(V^2 + 2|w1|).
  - Then e~ = y2|y2| and |a~| = 1/|2-|e~||.
  - **This reproduces every printed w1(T/2) and e~(T/2), and every a~(T/2) except row 14.** That includes rows 7, 8 and
    15, which have |e~| > 2.
  - **This rule is mine.** It reproduces the printed numbers; the paper does not state it.
  - a~(T/2) is printed positive even where |e~| > 2.
  - JAMM/34/2007 row 7 (51.553) is just outside the range that its e~ allows (51.560-51.586). The journal's 51.56234 is
    inside it.

**Tables 2, 4, 5: C on every row, at t = 0 and at T/2:**
- Printed C - (-2H + mu) is at most 5.1e-5, 4.6e-5 and 4.5e-5 at t = 0. (For Table 2 this is at mu = 0.00095.)
- At T/2 the residual is larger only where x1(T/2) is within about 0.02 of P2. There the 5-decimal x1(T/2) cannot fix
  1/r2.

**w1 = mu/(1 - x1(T/2)):**
- Every row agrees within the rounding of x1(T/2).
- Rows with x1(T/2) = 0 give w1 = mu exactly (0.1000 and 0.2000, as printed). These are built-in positive controls.
- At mu = 0.00095388, 17 of the 35 rows of Table 2 fall outside their rounding range. At 0.00095, none does.

**Integrations** (DOP853, rtol 1e-12; y2(0) corrected at fixed x1(0), so that y1 = 0 at the x2 = 0 crossing nearest the
printed T/2):

| Table | rows integrated | correction to y2(0) | (1-mu)T/(2 pi) vs printed T~ | x1, y2 at T/2 vs printed |
|---|---|---|---|---|
| 2 (mu = 0.00095) | 1-4, 10, 11, 13, 14, 22, 24, 25, 27 | <= 6e-6, except rows 2-3 (-5e-5, -6e-5) | within 1e-5 on all 12 | within 3e-5 except rows 2-3 (1.3e-4, 2e-4) |
| 4 (mu = 0.1) | 1-4, 12-14, 20, 21, 24-26 | <= 1e-5 | within 3e-5 on all 12 | within 5e-5, except y2 = 6.2758 near P2 (4e-4) |
| 5 (mu = 0.2) | 1-6, 11-14, 16, 17, 26, 27 | <= 1e-5, except row 6 (-2.8e-4; x1(T/2) = 0.99042, close to P2) | within 1e-5 on rows 2-4 and 11-27; within 4e-4 on row 6; rows 1 and 5 are the misprints | within 6e-4, except row 6 y2(T/2) (-5.474 against -5.454) |

- **The period convention.** T/(2 pi) without the factor misses by exactly 1/(1-mu) on every row: for example, Table 4
  row 1 gives 0.40270 against the printed 0.36243 = 0.9 x 0.40270.
- **mu_J.** With mu = 0.00095388 or 0.001, rows 10, 11 and 22 miss the printed T~ by 1e-4 to 5e-3. With 0.00095 they
  hit it to 1e-5.
- **Table 2, rows 2-3** (T~ ≈ 0.4995, Tr = -2) close less cleanly than any other row. The correction is 10 times the
  rounding of y2(0), and the T/2 point is off by about 1.5e-4. Correcting x1(0) instead does not help. **I do not call
  this a misprint.** The cause is not found. Both rows are near the Tr = -2 orbit, where the family meets family a.
- **Table 5, rows 10-12.** T~ falls from 1.10710 to 1.10211 to 1.10179. Integration reproduces both 1.10211 and
  1.10179. **So the period is really non-monotonic at mu = 0.2.** The paper claims monotonic period only for mu_J
  (and 67 for mu = 0.1).
- Rows with |Tr~| up to 6.7 still closed, once the half-period guess used the (1-mu) convention. My first run guessed
  with T/(2 pi). It locked onto wrong crossings for T~ > 3, and that run is not used.

## 4. Content notes

- **Family h at mu_J** follows the generating family closely; Table 3 maps the critical orbits. Linear stability holds in
  three intervals: start to orbit 3, orbits 11-13 and orbits 25-27. These are the Ir, E+1/2 and E+1/4 parts.
- **mu = 0.1:** stable from the start to orbit 1, from orbit 2 to orbit 3, and from orbit 8 to orbit 10. Families a and
  c meet family h at orbits 1 and 2.
- **mu = 0.2:** stable in 6 intervals: start-1, 2-3, 8-9, 13-14, 16-17 and 26-27.
- **Earth-Moon:** Broucke 1968 (held) computed family h as his A1 "from the beginning to as little as orbit 12 of
  Table 2". No new Earth-Moon numbers are given.
- At mu = 0.1 and 0.2, the departure from the generating family grows with a~(0): the farther from P1 or P2, the larger
  the perturbation.

## 5. Citation mining (10 references)

Checked with `ls cyclers_pdf/papers | grep -i` and a grep of `docs/notes/CORPUS_INDEX.md`, under the spellings
Bruno, Bryuno, Brjuno, Varin, Hénon, Henon, Broucke, Bray and Goudas. Wanted-list rows are numbered as in
`2026-10-05-960-wanted-papers.md` today.

**HELD:**
- Bruno 1994, de Gruyter (`bruno-1994-restricted-3-body-problem-...`).
- Bruno & Varin 2005, KIAM 67 (`bruno-varin-2005b-...`, TeX source and translation; this paper is its journal form).
- Bruno & Varin 2006, CMDA 95:27 (`bruno-varin-2006-...`).
- Bruno & Varin 2007, JAMM 71:933 (`bruno-varin-2007-...-jamm-71-933-...`). The paper cites it as "pp. 1034-106", a
  truncation of the PMM pages 1034-1066.
- Broucke 1968, JPL TR 32-1168 (`broucke-1968-...`).
- Hénon 1997, LNP m52 (`henon-1997-...`).
- Hénon 2001, LNP m65 (`henon-2001-...`).

**Not held, already on the wanted list:**
- Bruno 1996, KIAM Preprint 93: row 55.
- Bray & Goudas 1967, Adv. Astron. Astrophys. 5:71-130: row 56.

**Not held, not on the wanted list (proposals):**
- **Bruno & Varin 2008, "On families of periodic solutions of the restricted three-body problem", Sol. Syst. Res.
  42(2).**
  - The paper prints pp. 163-185. Crossref gives pp. 154-176, doi 10.1134/S003809460802007X (CONFIRMED, Crossref
    query). The printed pages may be the Russian Astron. Vestnik pages.
  - This is very likely the journal form of KIAM 10/2005, which is held and translated.
  - It defines coordinate systems I-IV and formula (19) for w1, which this paper uses. Proposed as an attribution row.
- **Bruno & Varin 2009, "Family h ... for big mu", Sol. Syst. Res. 43(2):158-177, doi 10.1134/S0038094609020099**
  (Crossref CONFIRMED). This is the remaining half of row 17 (sec. 0).
- From the same Crossref query, for row 55's housekeeping:
  - Row 55's "Bruno & Varin (2008), Families c and i ... at mu = 5e-5, Astron. Vestnik 42" is listed by Crossref as
    Sol. Syst. Res. **43**(1):26-40 (2009), doi 10.1134/S0038094609010031. The volume in row 55 may be wrong.
  - Varin, "Closed families of periodic solutions of a restricted three-body problem", Sol. Syst. Res. 43(3):253-276
    (2009), doi 10.1134/S0038094609030071, is probably the journal form of the held KIAM 16/2008.
  - **None of these was opened**, so every link above is from metadata only.

*Filed as `cyclers_pdf/papers/bruno-varin-2009a-family-h-periodic-solutions-restricted-problem-small-mu-sol-syst-res-43-1-2-doi-10.1134-s003809460901002x.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/bruno-varin-2009-family-h-small-mu-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
