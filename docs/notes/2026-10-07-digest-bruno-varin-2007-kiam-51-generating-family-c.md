# Digest: Bruno & Varin 2007, "Generating family c of periodic solutions of the restricted problem" (KIAM preprint 51/2007) (#960 batch 30)

A. D. Bruno and V. P. Varin (Keldysh Institute, Moscow), "Порождающее семейство c периодических решений
ограниченной задачи" (Generating family c of periodic solutions of the restricted problem), KIAM Preprint No. 51,
2007, 14 pp. (cover says 14 pp.; the PDF has 15 pages with the cover). In Russian. ISSN 2071-2898.
https://library.keldysh.ru/preprint.asp?id=2007-51. Licence: CC BY 4.0 (stated on the cover).
- Supplied as supplied file `prep2007_51.pdf`, md5 b391cf7faf30e5023c37dd31118d5e28. It has a text layer, but the
  layer is in a legacy Cyrillic font encoding and is unreadable. I read everything on the page images.
- Filed as `cyclers_pdf/papers/bruno-varin-2007b-generating-family-c-periodic-solutions-restricted-problem-kiam-preprint-51-2007-russian.pdf`.
  with `-en-translation.tex` / `.pdf` beside it.
- Full English translation (LaTeX, compiled with tectonic, 13 pp.): filed as the `-en-translation.tex` / `.pdf` pair beside the PDF;
  an OCR text sidecar (`-ocr-rus.txt`, force-OCR rus+eng) is filed too because the PDF text layer is unreadable. It uses the same conventions as the held 34/2007 translation.
- **How I read it:**
  - I read all 15 pages at 150 dpi.
  - I read Tables 1-3 (printed pp. 7-10) again on 300 dpi crops, cell by cell. Every number in this digest and in the
    translation comes from those images.
  - Arithmetic checks: `checks-kiam51.py`, output in `checks-kiam51.out` (filed beside the PDF).
- Wanted list: row 59 (part: "Bruno & Varin ... 51/2007 (generating family c)"). This copy closes that part.

## 0. Verdict

**A short (6 pp. of text) follow-up to KIAM 36/2006 and 34/2007. It gives the generating family c at mu = 0
(Table 1, three critical orbits) and re-tabulates the generating family i (Table 2, 39 critical orbits; Table 3, 13
P2-encounter points). It also introduces the local coordinate w1 = mu/(1 - x1) near P2.**
- **Family c at mu = 0 has two parts.** c' is Hill's family c, which leaves L1. c'' is the part of the segment
  family B1 with a~ <= 1. c'' runs from orbit 1 (the point P2, C = 3), through orbit 2 (a P1-collision, C =
  1.401879), to orbit 3 (the retrograde circle a = 1, e~ = -1, C = -1). There it ends as a locally double family
  on family h. C falls monotonically from 3 to -1. Tr = +inf and Tr_v = -inf on the whole of c'', except at orbit 3.
- **What this adds to the held Bruno-Varin work.** The JAMM 2007 paper and 34/2007 describe c only in one sentence (sec. 4.4: "from L1 as Hill's c, then a
  part of B1 ... to a = 1, e~ = -1, where it ends as a double family on h"). This preprint gives the numbers.
  - It is the first held source with family-c critical orbits at mu = 0.
  - It gives the T/2 data (velocity at the P2 encounter, and w1, y2) for the odd-subscript orbits of family i.
    JAMM Table 2 does not have these.
  - It gives the **number of the "Complex families" preprint: KIAM 35/2007, 27 pp.** (ref. [8]). Wanted row 18 says
    "number not found". This closes that gap.
  - No finite-mu data. The Astron. Vestnik 2008 paper (families c and i at mu = 5e-5) is not cited here and is still
    wanted (row 59).
- **Family i: Table 2 agrees with JAMM Table 2 in the t = 0 columns.** I compared spot rows 1, 18_1 and 39_2 with
  the JAMM digest values. I did not do a full cell-by-cell comparison with JAMM.
- **Misprints and inconsistencies found (details in sec. 3):**
  1. Table 1, orbit 1: T~ = 3/2 is printed. Hénon 1969 (held) gives T/2 = 4.41937 at Gamma -> -inf for Hill's
     family c, so T~ = 1.40673. The preprint's own Fig. 4 shows about 1.41.
  2. Table 1, orbit 2: x1(0) = 1 is printed. The orbit's t = 0 point is the P1 collision (x1 = 0).
  3. Table 2, row 16: the sign of y2(T/2) is wrong (+1.14471; it should be -1.14471).
  4. Table 2, row 25_2: the T/2 point gives C = 2.0045, not 1.484505. Unresolved.
  5. Table 2, rows 32_3 and 33_3: the v1, v2 cells do not give the Table 3 (w1, y2) by eq. (1.3). Row 32_3 repeats
     row 6_1's velocities. Unresolved.
  6. The sign of w1 differs between eq. (1.3), Tables 1 and 3 (negative) and Figs. 1, 2 and 6 (positive). The
     magnitudes agree.
  7. Text: "P2 = {x1 = 0}" and "x1 <= 0 / >= 0" should be x1 = 1 and x1 <= 1 / >= 1. The caption of Fig. 6 says
     "family c''" but the figure shows family i.
- **Catalogue implication (PROPOSAL only):** none. No catalogue row cites family c or i as a cycler source. Family c
  (L1 Lyapunov orbits continued to large size) is background for `#944` X2 / `#956` R9 symmetric-family controls.
  Proposal: if a family-c control is ever built at small mu, use Hénon 1969 for the Hill limit and Table 1 here for
  the mu = 0 C range [-1, 3]. Do not use the printed T~ = 3/2.

## 1. Content

### 1.1 Sec. 1: local coordinates near P2 (eqs. 1.1-1.3)
- w1 = mu/(1 - x1), keeping y2. P2 goes to w1 = +-inf, and P1 to w1 = mu. This replaces z1 = (x1 - 1)/mu of
  KIAM 10/2005 eq. (2.3): w1 = -1/z1.
- By Broucke's principle the sign of 1 - x1 does not change along one family's characteristic.
- At mu = 0 the limit point comes from the entry velocity (v1, v2) of a segment solution into P2, with V^2 =
  v1^2 + v2^2 = 3 - C (eq. 1.3):
  w1 = -V^2 v1 / (V - |v1|), y2 = 1 + sqrt(2|w1| + 3 - C) sgn v2.
- At mu = 0, the P2-side T/2 values in Tables 1-2 are found as follows: x1 = 1, y2 from (1.3), then
  e~ = x1 y2^2 sgn y2 and a~ = x1/|2 - |e~||. This is why |e~(T/2)| > 2 can occur (for example -2.483035).
  **This also settles the open question in the JAMM digest (sec. 2.1) for family h, Table 1, rows 7, 8, 14 and 15.**
  Those T/2 cells are P2-encounter values in this convention, not crossings of the plane Pi. For each of those rows,
  a~(T/2)·|2 - |e~(T/2)|| comes out close to 1, as the convention needs (0.98, 0.999, 1.005, 0.999). That is
  consistency only: w1 is not printed for family h, so C cannot be checked there.

### 1.2 Sec. 2: generating family c
- c' (Hill's c): computed by Hénon [9]. Its vertical traces are in Hénon 1974 [10]. Its characteristics in Hill
  coordinates are in [11, 12, 7].
- c'' (part of B1, a~ <= 1, e~ in [1, 2] and [-2, -1]): Figs. 1-4 give y2, C, (a~, e~) and T~ against x1 / w1.
  Fig. 5 shows orbit 2: a symmetric figure-eight-like loop from P2 at (1, 0) out to |x2| ≈ 1.42 and into P1.

### 1.3 Sec. 3: generating family i near P2
- The pieces K3, K9 and K15-K17 begin and end in P2. They are plotted in (w1, y2) in Fig. 6, using C12, C23 and C34
  for K3/K9/K15/K17 and B1 for K16 (= C12 + 2B1).
- K16 runs along c'' (B1). Between orbits 34 and 33 there is a zigzag. So the lower B1 branch carries zigzags, like
  the upper-branch zigzags of 36/2006 item 1.4. Figs. 7-8 show them schematically against n, with the ordinate
  -+1/(3 - C)^(1/3).
- At mu = 0 the i characteristic lies to the left of the B1 characteristics or on them, never to the right. For
  mu > 0 it lies left of family c.

## 2. Tables (image-read at 300 dpi)

- **Table 1 (c'')**:
  - orbit 1: x1 = 1, y2 = 1, C = 3, T~ printed 3/2;
  - orbit 2: v(T/2) = (-0.77337, -1), T~ = 1, C = 1.401879, a~(0) = 0.71333, e~(0) = +-2,
    a~/e~(T/2) = 2.070245 / -2.483035, w1/y2(T/2) = -2.51823 / -1.57577;
  - orbit 3: y2 = -1, v(T/2) = (0, -2), T~ = 1, C = -1. Tr goes from +inf to 2 and Tr_v from -inf to 2.
- **Table 2 (family i)**: 39 rows with x1(0), y2(0), x1|v1(T/2), y2|v2(T/2), T~, C, Tr, Tr_v, a~(0), e~(0), a~(T/2) and
  e~(T/2). These are the same orbits as JAMM Table 2 / 36/2006 Table 1, plus the T/2 columns.
- **Table 3**: (w1, y2) at T/2 for the 13 odd-subscript orbits.
- Note: 14_2, 26_2 and 31_3 share C = 1.401875. 31_3 has the same P2 velocity (-0.77337, -1) and w1 = -2.51824 as
  c'' orbit 2. So family i touches family c'' at its P1-collision orbit. C differs in the 6th decimal (1.401875
  against 1.401879). The Kepler value is 1.401876 (sec. 3).

## 3. Checks (`checks-kiam51.py` / `.out`)

Formulas at mu = 0: C = -y2^2 + 2 x1 y2 + 2/|x1| on Pi; e~ = x1 y2|y2|; a~ = x1/|2 - |e~||. These are the JAMM
digest formulas.
- **Table 2, t = 0:** C, a~ and e~ reproduce from (x1(0), y2(0)) to <= 1.5e-4 in all finite rows. Row 26_2 (x1 =
  -0.00094) is rounding-limited: 2/|x1| alone has an uncertainty of +-11.
- **Table 2, T/2, even subscript:** these reproduce, except:
  - rows 11_2 and 13_2, which are rounding-limited near P1 (1.5e-3 and 7e-3);
  - **row 16**, where only y2 = -1.14471 gives the printed C (-0.436758 against -0.436790) and e~ = -1;
  - **row 25_2**, where (0.05931, 5.69137) gives C = 2.00454 against the printed 1.484505. The printed C equals
    1/a~(0) = 1.484516 for the t = 0 P1 collision. The a~ and e~ at T/2 do agree with (x1, y2). So either the T/2
    pair or the C is wrong. Unresolved.
- **Circular rows (Kepler, independent):** C = 1/a + 2 sgn(e~) sqrt(a). All 14 rows on Id and Ir agree to <= 3e-5.
  T~ is an integer multiple of the circular synodic period 1/|n - 1|: x1 for 1, 2, 8, 9, 10, 20-22 and 36-38; x3 for
  row 4; x5 for row 16; x7 for row 28. Exact values: 3^(-2/3) = 0.480750 is printed as 0.48074, and (5/3)^(-2/3) =
  0.711379 as 0.71137. The values are truncated, not rounded, as the JAMM digest also found.
- **Odd subscripts, P2 encounter:** V^2 = 3 - C holds to <= 3.3e-5 in all 13 rows. Eq. (1.3) reproduces |w1| and y2
  of Table 3 to <= 7e-5 in 11 rows. It fails for 32_3 and 33_3, although both of those Table 3 rows satisfy (1.3)
  inverted, (1 - y2)^2 = 2|w1| + 3 - C, to 1e-5. Table 3 also reproduces a~(T/2) and e~(T/2) from x1 = 1.
- **Table 1, orbit 2 (independent Kepler check):** this is a rectilinear collision ellipse, so C = 1/a. If P3 leaves
  P2 (r = 1), passes apocentre, and hits P1 in time pi (T~ = 1), then a = 0.713330 and C = 1.401876.
  v1 = -sqrt(2 - 1/a) = -0.77338, and v2 = -1 exactly (the frame velocity at P2). Printed: 0.71333, 1.401879
  and -0.77337. They agree to the last digit or within a few units of it.
- **Table 1, orbit 1 (independent):** the c' -> c'' junction is the high-energy limit of Hill's family c. In the
  massless linear Hill limit, a symmetric arc from P2 back to P2 needs tan tau = (3/4) tau, so tau = 4.41937. Hénon 1969
  Table 2 (held; `-tables.txt` line 105) prints T/2 = 4.41937 at Gamma = -inf. That gives **T~ = 1.40673, not the
  printed 3/2**. Fig. 4 here also tops out near 1.41 at x1 = 1.
- **Positive controls:** the same script reproduces the exact circular-orbit C values (for example 3.466806 at
  a = 3^(-2/3)) and the Hénon root. So the formulas are not trivially self-consistent.

## 4. Citation mining (12 references)

Checked with `ls cyclers_pdf/papers | grep -i` and a grep of CORPUS_INDEX.md.
- [1] Bruno 1990 Nauka book: HELD as its English edition (`bruno-1994-restricted-3-body-problem-...`).
- [2] KIAM 10/2005: HELD (`bruno-varin-2005a-...`). [3] 67/2005 and [4] 64/2005: HELD as TeX source and
  translation (`bruno-varin-2005b-...`, `-2005c-...`).
- [5] KIAM 36/2006 (generating family i): was wanted row 59; filed in batch 30 with a full translation
  (`2026-10-07-digest-bruno-varin-2006-kiam-36-generating-family-i.md`).
- [6] Bruno & Petrovich, KIAM 53/2006, Desingularisations: not held. The JAMM digest proposed it for a Tier D row;
  it is not on the wanted list yet.
- [7] KIAM 34/2007: HELD (`bruno-varin-2007-...-kiam-preprint-34-2007-russian...`).
- [8] KIAM **35/2007**, "Complex families ...", 27 pp.: not held. It is wanted row 18. **Update the row with the
  number 35/2007.**
- [9] Hénon 1969 A&A 1:223: HELD (`henon-1969-numerical-exploration-restricted-problem-V-...`).
- [10] Hénon 1974, "Vertical stability of periodic orbits in the restricted problem II. Hill's case", A&A 30:317-321:
  not held, and not on the wanted list (the held items are Hénon 1973 I and 1973b). **New candidate**, low
  priority: it gives the vertical traces on Hill's families, including c.
- [11] Bruno 1994, in Seimenis (ed.), Hamiltonian Mechanics, Plenum, pp. 43-49: not held. Already proposed in the JAMM
  digest (Tier D row); not on the wanted list.
- [12] Bruno KIAM 93/1996: not held. It is on wanted row 59.

Unresolved: items 4 and 5 of the verdict list (Table 2, rows 25_2, 32_3 and 33_3). Checking them would need the
36/2006 Table 1 (not yet filed) or a recomputation of those orbits at mu = 0.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
