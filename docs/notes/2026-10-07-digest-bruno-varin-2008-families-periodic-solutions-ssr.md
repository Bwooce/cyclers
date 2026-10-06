# Digest: Bruno & Varin 2008, "On Families of Periodic Solutions of the Restricted Three-Body Problem" (Solar System Research 42(2)) (#960 batch 34)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, RAS), "On Families of Periodic Solutions of the
Restricted Three-Body Problem", Solar System Research 42(2):154-176 (2008), doi:10.1134/S003809460802007X.
Russian original: Astronomicheskii Vestnik 42(2):163-185. Received 14 December 2006. "Translated by the authors."
RFBR project 08-01-00082.
- File: `0e01a138-bruno2008.pdf`, 23 pp., md5 `4f997bad51c976c74a58116be533f3ee`.
- **Proposed corpus filename:**
  `bruno-varin-2008-families-periodic-solutions-restricted-three-body-problem-generating-families-c-i-ssr-42-154-doi-10.1134-S003809460802007X.pdf`
  (no `bruno-varin-2008-` prefix is in use).
- **How I read it:**
  - I read the whole text layer, and every table on 300 dpi page images (PDF pp. 6, 12 (landscape, rotated), 16),
    in strips of about 12 rows.
  - All 46 rows of Table 2, the 3 rows of Table 1, the 13 rows of Table 3 and the 10 rows of Table 4 were read on
    the image. The text layer agreed with the image except for the plus-minus and minus-plus signs on infinities and
    on 2 (rows 3, 18_1, 25_2), which I took from the image.
  - I compared it with: KIAM 10/2005 (our translation), the held CMDA 95:27 (2006) paper, KIAM 36/2006 Table 1
    (TeX source, exact), KIAM 51/2007 Tables 1-3 (our translation; every differing cell re-read on 250 dpi images
    of the KIAM 51 PDF), and JAMM 71:933 Table 2 (text layer, image-checked in the KIAM 36 digest).
  - Scripts and outputs in this folder: `ssr2008_data.py` (the transcription), `compare_2008.py` ->
    `compare_2008.out`, `check_2008.py` -> `check_2008.out`, `arcs_mu0.py` -> `arcs_mu0.out` (independent mu = 0
    arc solutions), `cr3bp.py` (equations), `make_yaml.py`.
  - YAML: `bruno-varin-2008-ssr-tables.yaml` (Tables 1-4, printed values, 10 notes).

## 0. Verdict

**This is the journal summary of the Bruno-Varin programme for the generating families c and i at mu = 0.**
It is not simply the journal form of KIAM 10/2005. It has three sources:
- Secs. "Introduction" to "Methods of computation" are the general part. Their text is the revised text that was
  already printed in CMDA 95:27-54 (2006, held), which is itself the journal form of KIAM 10/2005. The CMDA paper
  then treats family h; this paper treats families c'' and i instead.
- "Generating family c" is KIAM 51/2007 (family c'', Table 1, Figs. 3-7).
- "Generating family i" is KIAM 36/2006 (pieces K1-K25, the whole-family sequences, Table 4 = KIAM 36 Table 2),
  with the revised Table 2 and new Table 3 of KIAM 51/2007.
- **Identity for the wanted list:** Varin's KIAM 16/2008 cites, as its ref. [7], "Bruno & Varin, The families c and i
  ... at mu = 5e-5, Astronomicheskii Vestnik 2008, vol. 42" with blank issue and pages. The 2009 SSR paper replaces
  that reference with this paper (Astron. Vestn. 42(2):163; SSR 42(2):154). So **wanted row 55's item "Bruno & Varin
  (2008), Families c and i ..., Astron. Vestnik 42" is this paper** (the title changed in print). Proposal: mark it
  received.

**What it gives the project:**
- The best-checked print of the family-i critical-orbit table at mu = 0 (Table 2, 39 orbits). It rounds where the
  preprints truncated. It corrects 13 cells of KIAM 51 that our checks or the KIAM 36 digest had flagged, and it
  fixes the row 25_2 half-period point that the KIAM 36 digest left UNRESOLVED.
- But it **introduces four new errors** (row 1 C, rows 14_2 and 31_3 T~, a decimal slip in row 25_2) and keeps three
  old ones (rows 16, 32_3, 33_3). Sec. 3 lists them all.
- It revises the signs of the infinite traces Tr, Tr_v on 16 rows. Text and table agree with each other. I could not
  test these signs.
- Catalogue implication (PROPOSAL only): none. There are no cyclers here and no catalogue row cites family i or c''
  at mu = 0. If a future control uses family-i generating orbits, use this paper's Table 2 with the corrections in
  sec. 3 (or our YAML with its notes), not KIAM 36, KIAM 51 or JAMM alone.

## 1. Frame and conventions (stated first)

- Planar circular restricted problem. Synodic frame, origin at P1 (mass 1 - mu), P2 at (1, 0), unit angular speed.
- Hamiltonian (eqs. (1)-(2)): H = H0 + mu R, H0 = (y1^2 + y2^2)/2 + x2 y1 - x1 y2 - 1/r, R = 1/r + x1 - 1/r2.
  The y are canonical momenta: y2 = dx2/dt + x1. C = -2H. Symmetry plane Pi: x2 = y1 = 0.
- Four coordinate systems on Pi: I (x1, y2); II (x1, C); III (a~, e~) at P1, eq. (6): e~ = x1 y2 |y2| / (1 - mu),
  a~ = x1/(2 - e~) [in use: |2 - |e~||]; IV (w1, y2) at P2, eq. (8): w1 = mu/(1 - x1).
- T~ = T/(2 pi). Traces are mapped by eq. (10): Tr~ = Tr if |Tr| <= 2, else (1 + log2|Tr|) sgn Tr.
- At mu = 0 the entry velocity (v1, v2) of an arc into P2 is the rotating-frame velocity. 3 - |v|^2 = C.
- Eq. (19) as printed is garbled in the text layer. The form that reproduces Tables 1 and 3 is
  w1 = V^2 |v1| / (V - |v1|) for v1 < 0, and y2 = 1 + sgn(v2) sqrt(2 w1 + V^2), with V^2 = 3 - C.

## 2. Content

- **Introduction and general part** (pp.154-158): the same as CMDA 2006 secs. 1-5, with family lists, Broucke's
  principle and the nine principal families (Strömgren / Broucke / Bruno names). Differences from KIAM 10/2005
  (section by section, from our translation):
  - KIAM 10 uses three coordinate systems (x1, C; a~, e~; z1 = (x1 - 1)/mu, y2). SSR uses four, and system IV is now
    w1 = mu/(1 - x1). CMDA 2006 had a different system IV (hyperbolic a~*, e~*). So SSR's system IV comes from
    KIAM 51.
  - KIAM 10's P2 curve "a~|1 - |e~|| = 1" is printed in SSR as a~|2 - |e~|| = 1, which is correct (x1 = 1).
  - KIAM 10's unbalanced trace formula sign(Tr) log2(|Tr| + 1) becomes SSR eq. (10).
  - KIAM 10's method section (quasi-regularisation, the time change (5.2)-(5.3), Fourier "without saturation") is
    replaced by a shorter one: fixed-step RK5, shooting and Newton, Thiele-Burrau and Levi-Civita near collisions,
    FFT, computer algebra on demand. C is conserved to < 1e-10 (usually < 1e-13).
  - Kept from KIAM 10: mu_M = 0.1215585 (a slip for 0.01215585; CMDA 2006 has the same slip on its first use).
  - KIAM 10 writes mu = m2; SSR writes mu = m2/(m + m2).
- **Family c''** (pp.168-169; Table 1, Figs. 3-7): the part of B1 with a~ <= 1. C falls from 3 to -1. Three critical
  orbits: 1 (start, the point P2), 2 (P1 collision), 3 (end, two-fold on the circular orbit a = 1, e~ = -1).
- **Family i** (pp.169-174; Tables 2-4, Figs. 8-17): pieces K1-K25 and a tail list; the whole family as blocks
  Id, L_p, (2p+1)Ir, M_p; the sequences U and V; Table 4; the w1, y2 view near P2 (Fig. 15) and the zigzag diagrams
  (Figs. 16-17); the arc-family formulas (27)-(31); hypotheses 1-2 on Tr_v; and the Broucke-principle argument
  (32)-(37) that only the part a < a0 of each C_{m,m+1} joins family i.
- **Text changes from KIAM 36** (from our translation): SSR fixes these KIAM 36 slips: K6 ends at E_{3/2}(+1)
  (orbit 11_2); K16 now names 31_3 and 32_3 as the collision orbits and 33_3 as the extremal orbit; K18 says orbit 36
  coincides with orbit 22. SSR keeps "E_{4/5}(+1)" in the tail list (for E_{4/3}(+1)). The trace descriptions of K7,
  K13 and K16 are rewritten (sec. 4). A stray "7.2." heads "Description of the Whole Family".

## 3. Errors and corrections in the tables (checks in `check_2008.out`, `arcs_mu0.out`)

At mu = 0: C = -y2^2 + 2 x1 y2 + 2/|x1| on Pi. SSR rounds (3^(-2/3) = 0.480750 is printed 0.48075, where KIAM
printed 0.48074), so the tolerance is half a unit in the last digit, propagated.

**Checks that pass:** C from (x1, y2)(0) for all 35 testable non-collision rows (row 26_2 has x1 = -0.00094, two
significant digits; not testable). C at T/2 for all even/no-index rows except 1, 16 and 25_2 (below). 3 - |v|^2 = C
for all 13 odd-index rows. Table 3 from the Table 2 velocities by eq. (19) for 11 of 13 rows (not 32_3, 33_3). The
Id and Ir closed forms. Table 4: 20 of 40 cells re-derived from closed forms to 1e-9 (as in the KIAM 36 digest).
Table 1 row 1 T~ = 1.40673 = tau/pi, with tau = 4.41937 the root of tan tau = 3 tau/4 (Hénon 1969 Hill limit).

**Independent arc solutions** (`arcs_mu0.py`): I solved, at mu = 0, for the symmetric Kepler arc that leaves a point
on the x1 axis and reaches P2, at a given C. **Positive control:** at C = 1.793529 the solution has its axis point at
(0.024392, 8.97987), a~, e~ = 0.73776, 1.96694, and half-arc time h/pi = 1.05598. These are row 13_2's T/2 values and
T~(13_2) - T~(6_1) = 2.05598 - 1 to every printed digit.

| row | cell | SSR prints | earlier prints | finding |
|---|---|---|---|---|
| 1 | C | 3.66806 | 3.466806 (KIAM 36, 51); 3.4668 (JAMM) | **SSR misprint** (the digit 4 is lost). Closed form 1/a + 2 sqrt(a), a = 3^(-2/3): 3.466806 |
| 14_2 | T~ | 1.97277 | 1.97370 (36, 51); 1.974 (JAMM) | **SSR misprint, probably.** The C12 half-arc from the printed x1(0) = -1.08920 at C = 1.401876 (solved: x1 = -1.089200) has h/pi = 0.97371, so T~ = 0.97371 + 1 = 1.97371 |
| 31_3 | T~ | 2.97277 | 2.97370 | as 14_2: 2.97371. Both prints keep T~(31_3) = T~(14_2) + 1 |
| 25_2 | x1(T/2), y2(T/2) | 0.01040, 43.8419 | 0.05931, 5.69137 | **SSR corrects the T/2 point** (the KIAM 36 digest's unresolved 25_2 item). Solved B1 arc at C = 1.484506: axis point (0.001040, 43.84190), a~, e~ = 0.71770, 1.99855 = SSR's a~, e~(T/2). But **x1 = 0.01040 is a decimal slip for 0.00104** (0.01040 gives C = -1729) |
| 25_2 | T~ | 3.01130 | 3.08850 (36, 51); 3.088 (JAMM) | **SSR right**: h/pi = 1.01130 for that B1 arc, and T~ = T~(18_1) + h/pi = 2 + 1.01130 |
| 25_2 | a~, e~(T/2) | 0.71770, 1.99855 | 0.75483, 1.92141 | **SSR right** (C from them = 1.48455) |
| 32_3 | v1, v2(T/2) | -0.45439, -1 | the same in KIAM 51; KIAM 36: -0.77237, -0.78096 | **copy error from row 6_1, inherited from KIAM 51.** Table 3 (2.85831, -1.63118) needs v = (-0.77238, -0.78096); so does the solved B1 arc at C = 1.793529 |
| 33_3 | v1, v2(T/2) | -0.03600, -0.35297 | +0.035995, -0.35297 (51); -0.33446, -0.11839 (36) | **wrong in SSR and KIAM 51.** It fits C but not Table 3 (2.07031, -1.06555). The solved B1 arc at C = 2.874117 gives v = (-0.33446, -0.11839), w1, y2 = 2.07031, -1.06555, as KIAM 36 printed. Also T~(33_3) = [T~(12_2) - 1.29405] + 2(1.29405) = 3.57535, as printed |
| 16 | y2(T/2) | +1.14471 | the same in all prints | sign slip, kept (needs -1.14471 for C = -0.436790) |
| 34_3 | v, Table 3 | (-0.73388, -0.61135); (3.02577, -1.63892) | 36: (-0.73387, -0.61134); 51 T3: (-3.02577, -1.63891) | SSR right to the last digit (solved arc: v = (-0.73388, -0.61135), w1 = 3.02578, y2 = -1.63892) |

**Corrections SSR makes to KIAM 51** (all confirmed by the checks above or by the KIAM 36 digest):
- e~(0) of rows 12_2, 19_1, 27_2 and 33_3 (KIAM 51 had already fixed rows 5_1, 7_1, 17_1, 29_1 and 35_1). These are
  the ~1e-4 errors of KIAM 36 sec. 3.1. SSR's values are the closed-form or (x1, y2) values.
- Row 18_1 a~, e~(T/2): 23.93659, -1.95822 (51: 23.87091, -1.95811). Table 3: 2.12073, -1.39937 (51: -2.12065,
  -1.39932). With v1 = -sqrt(3 - C - 1) = -0.717979, eq. (19) gives w1 = 2.12074, so SSR's are the more exact values.
- Row 32_3 a~, e~(T/2): 1.51346, -2.66074 (51: 1.51374, -2.66061), now consistent with SSR Table 3.
- Table 1 (c''): row 1 T~ 1.40673 (51: 3/2), row 1 w1(T/2) 0 (51: 1), row 2 w1(T/2) +2.51823 (51: -2.51823).
- Table 3: all nonzero w1 now positive, as eq. (19) and the figures require (51 printed them negative).
- Row 3 and row 25_2 y2(0) now -+inf; rows 14_2, 26_2 y2(T/2) and rows 3, 18_1, 25_2 e~ now carry plus-minus signs.

**Not changed, still open:** Table 1 row 2 prints x1(0) = 1 with y2(0) = inf and e~(0) = +-2; the t = 0 point is the
P1 collision (x1 = 0), as KIAM 51's translation noted. Table 1 row 2 C = 1.401879 while Table 2 rows 14_2, 26_2 and
31_3 (the same B1 collision arc) print 1.401876 = 1/0.71333.

## 4. Cell-by-cell comparison (`compare_2008.out`)

Table 2 has 12 compared columns x 39 rows = 468 cells per comparison.

| SSR Table 2 against | same | last digit only (rounding vs truncation) | notation (+- / -+) | genuine difference |
|---|---|---|---|---|
| KIAM 51/2007 Table 2 | 324 + 2 same value | 95 (92 SSR larger in magnitude) | 3 | 44 |
| KIAM 36/2006 Table 1 | 293 + 2 | 87 | 9 | 76 + 1 sign |
| JAMM Table 2 (8 columns, 4 decimals) | 98 | 149 rounded or truncated | - | 26 |

- **The 44 genuine differences from KIAM 51** are: the corrections and errors of sec. 3 (row 1 C; 12_2, 19_1, 27_2,
  33_3 e~(0); 14_2, 31_3 T~; 18_1, 25_2, 32_3 T/2 cells; 33_3 v1), the row 3 y2(T/2) sign notation, and **26 trace
  cells on 16 rows** (11_2, 12_2, 13_2, 14_2, 15_2, 23_2, 24_2, 25_2, 26_2, 27_2, 30_3, 31_3, 32_3, 33_3, 34_3, 39_2;
  Tr and/or Tr_v).
- **KIAM 36 adds** the odd-index half-period columns, which KIAM 36 gives in system IV (a~'_*, e~_*) and SSR in
  system III at P2, and the row 28 a~(T/2) sign slip (SSR has +0.82548, correct).
- **JAMM** differs from SSR in the same trace cells and the same corrected cells. JAMM agrees with KIAM 51 throughout.
- **The trace revision.** KIAM 36, KIAM 51 and JAMM agree on the traces. SSR changes them, for example row 11_2
  Tr "2, -inf, +inf" -> "2, -2, +inf" and Tr_v "2, +inf" -> "2, -inf"; rows 13_2, 14_2, 25_2, 26_2 Tr -inf -> +inf;
  rows 31_3, 32_3 Tr +inf -> -inf and Tr_v -inf -> +inf. The SSR text agrees with its table: K7 now has Tr fall
  from 2 to -2 at 11_2, then +inf, -inf at 12_2, +inf again, and a new "additional orbit 14'_2" where it falls to
  -inf. K13 has a similar "26'_2". These are deliberate. I could not test them: the signs come from the authors'
  perturbation rules (Tr_1 signs, hypotheses 1-2), not from the tabulated states.
- **Table 3 against KIAM 51 Table 3:** 5 sign changes (6_1, 30_3, 31_3, 34_3 w1; SSR positive is right), 18_1 and
  32_3 values changed (SSR more exact, see sec. 3), 34_3 y2 last digit, 33_3 w1 2.07031 against -2.07030.
- **Table 4 against KIAM 36 Table 2:** all 40 cells identical.

## 5. Integrations (mu = 0, DOP853, rtol = atol = 1e-12; `check_2008.out` sec. 6)

- Rows 9 (Id), 11_2, 23_2, 39_2 (E(+1) ellipses through P2): from the printed (x1, y2)(0), the x2 = 0 crossing comes
  at T~ = 1.50001, 2.00000, 3.00034, 3.99997 and at (x1, y2) = (0.71138, 1.18563), (0.07969, 4.87724),
  (0.27595, 2.45690), (0.41247, 1.92053). These match the printed T/2 columns within the 5-digit input
  precision (row 23_2: y1 = 1e-3 at the crossing, because its orbit is sensitive).
- Rows 5_1, 7_1, 17_1, 19_1, 29_1, 35_1 (E_N(+-0) rows that end at P2): after pi T~ the state is at (1, 0) within
  1e-4, and the entry velocity matches the printed (0, v2) to 2e-5.

## 6. Citation mining (31 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i` in CORPUS_INDEX.md and the wanted list.

| reference | status |
|---|---|
| Abalakin, Aksenov, Grebennikov, Ryabov 1971 (Nauka handbook) | not held; not wanted |
| Bartlett 1964 | not held; **wanted row 22** |
| Breakwell & Perko 1974, CM 9:437 | HELD (`breakwell-perko-1974-...`) |
| Broucke 1968, JPL TR 32-1168 | HELD (`broucke-1968-...`) |
| Bruno 1978a, b, CM 18:9 and 18:51 | HELD (`brjuno-1978-...`, `brjuno-1978b-...`) |
| Bruno 1978c, KIAM preprint 91 (periodic flybys of the Moon) | not held; the journal form Bruno 1981 is HELD (`bruno-1981-...`) |
| Bruno 1981, CM 24:255 | HELD (`bruno-1981-...`) |
| Bruno 1993, KIAM 66; Bruno 1996, KIAM 93 | not held; **wanted row 55** |
| Bruno 1994, de Gruyter | HELD (`bruno-1994-...`) |
| Bruno 2000, Power Geometry (Elsevier) | not held; not wanted |
| Bruno 2006, Kosm. Issled. 44(3):258 [Cosmic Res. 44:245] | not held; not wanted. Candidate: it is the family-theory source cited for "general properties" |
| Bruno & Varin 2005, KIAM 10 | HELD (`bruno-varin-2005a-...`) |
| Bruno & Varin 2006, KIAM 36; 2007a, KIAM 51 | HELD (`bruno-varin-2006b-...`, `bruno-varin-2007b-...`) |
| Bruno & Varin 2007b, JAMM 71 | HELD (`bruno-varin-2007-...-jamm-...`); printed as "J. Appl. Math. Mech. 71(6):1034-1066", which are the PMM pages; the JAMM pages are 933-960 |
| Hénon 1968, Bull. Astron. 3:377 | HELD (`henon-1968-...`) |
| Hénon 1969, A&A 1:223 | HELD (`henon-1969-...`) |
| Hénon 1973, A&A 28:415 | HELD (`henon-1973-...`) |
| Hénon 1974, A&A 30:317 (vertical stability II, Hill's case) | **not held; not on the wanted list. Candidate** (vertical-trace data on family c of Hill's problem; cited for c') |
| Hénon 1997, 2001 (LNP m52, m65) | HELD (`henon-1997-...`, `henon-2001-...`) |
| Hitzl & Hénon 1977a, b | HELD (`hitzl-henon-1977-...`, `hitzl-henon-1977b-...`) |
| Kotoulas & Voyatzis 2004; Voyatzis & Kotoulas 2005; Voyatzis, Kotoulas & Hadjidemetriou 2005 | not held; **wanted row 56** |
| Perko 1981a, b | HELD (`perko-1981-...`, `perko-1981b-...`) |
| Strömgren 1935, Publ. Copenhagen 100 | not held (Strömgren 1933 Bull. Astron. is HELD, a different paper); not wanted |
| Szebehely 1967 | HELD (`szebehely-1967-...`) |

Proposals for the wanted list: mark row 55's "Bruno & Varin (2008) ... Astron. Vestnik 42" as received (this paper).
Add Hénon 1974 (A&A 30:317-321) as a low-priority row.

*Wanted-list row numbers are the current (batch-30) numbering.*

*Filed as `cyclers_pdf/papers/bruno-varin-2008-families-periodic-solutions-restricted-three-body-problem-sol-syst-res-42-2-154-doi-10.1134-s003809460802007x.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/bruno-varin-2008-ssr-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
