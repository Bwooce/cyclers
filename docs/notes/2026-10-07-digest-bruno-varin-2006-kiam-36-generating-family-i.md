# Digest: Bruno & Varin 2006, "The generating family i of periodic solutions of the restricted problem" (KIAM Preprint 36/2006) (#960 batch 30)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, RAS), "Porozhdayushchee semeistvo i periodicheskikh
reshenii ogranichennoi zadachi", Preprint No. 36, Moscow, 2006 (in Russian). UDC 521.1+531.314. RFBR grant 05-01-00050.
The title is confirmed from ATITLE.TEX and from the library page.
- Source: `prep2006_36.zip` from https://keldysh.ru/papers/2006/source/prep2006_36.zip. md5 444f23407847c7534003dd3820470aa8.
  - It holds five TeX files, dated 5-9 June 2006: A.TEX (driver), ATITLE.TEX (title and abstracts), PAR1.TEX (text),
    LIT.TEX (10 references) and TAB1.TEX (Table 1).
  - **There are no figures in the zip.** The text cites Figs. 1-13 (orbits; characteristics in (a~, e~); T~ against a~(0);
    zigzag schemes; characteristics near P2). None of them is in the source.
  - **The encoding is CP866 (DOS Cyrillic), not CP1251.** ATITLE.TEX fails a CP1251 decode. All five files decode
    cleanly as CP866.
- **Licence:** the library page (`k_2006-36.html`) carries a CC BY 4.0 badge (https://creativecommons.org/licenses/by/4.0/).
  So the English translation can be filed and shared as an adaptation, with attribution.
- **Proposed corpus filenames** (`bruno-varin-2006-` is already used by the CMDA 95:27 paper):
  - Filed as `cyclers_pdf/papers/bruno-varin-2006b-generating-family-i-periodic-solutions-restricted-problem-kiam-preprint-36-2006-russian-tex-source.zip`
  - `...-russian-tex-source.utf8.txt` (made here as `prep2006_36-russian-tex-source.utf8.txt`: the five files, CP866 ->
    UTF-8, concatenated in the order A, ATITLE, PAR1, LIT, TAB1)
  - `...-russian-tex-source-en-translation.tex` and `.pdf` (this batch: `translation.tex` / `translation.pdf`, 15 pp.)
- **How I read it:**
  - I read all five TeX files in full after the CP866 -> UTF-8 conversion (kept in `src/`).
  - I translated every sentence, equation, table cell, the references and both abstracts into `translation.tex`.
    It compiles with tectonic. It follows the 67/2005 translation's conventions: banner, translator's note,
    "Tr." footnotes, `%%BODY-BEGIN/END`, and references with Russian titles in brackets. I added a CC BY 4.0
    licence line under the banner. **The 67/2005 template has no licence line**, so this is new.
  - Both tables are in the source. I copied them verbatim from the TeX, and checked the copy with a script
    (see below).
  - I read JAMM 71:933 Table 2 (journal p.946 = PDF p.14) on a 250 dpi page image. All 39 rows of the PDF
    text layer agreed with the image. The check script reads that text layer at run time.
  - Arithmetic checks: `checks_kiam36.py`, output in `checks_kiam36.out` (filed beside the zip).
- Wanted list: **row 59** (current numbering; the JAMM digest's "row 61" is the batch-28 numbering). This preprint is one
  item of that row. Row 59 also holds Bruno 18/1972, 66/1993, 67/1993, 68/1993, 93/1996, Bruno-Varin 51/2007 and the
  2008 Astron. Vestnik paper. It stays open for those.

## 0. Verdict

**This is the source of JAMM 2007 sec. 6.1 and Table 2: the generating family i at mu = 0.**
- It is "Chapter III" of the Bruno-Varin programme. Chapter II was family h (67/2005).
- **Table 1** gives 39 critical orbits. Each row has x1(0), y2(0) and the half-period point or the P2 entry velocity,
  T~, C to 6 decimals, Tr, Tr_v, (a~, e~)(0), and (a~, e~)(T/2) or (a~'_*, e~_*)(T/2).
  - JAMM Table 2 is a 4-decimal copy of this table, without Tr_v and without x1, y2.
  - **So this preprint is the more precise and fuller source**, with two more columns.
- **Table 2** (not in JAMM): C to 9 decimals at the extremal orbits E_{(m+1)/m}(+1), C_{m,m+1}(1),
  E_{(m+1)/m}(+0) and Id((m+1)/m), for m = 1 to 10.
  - I reproduced three of its four columns to 1e-9 from Kepler two-body arithmetic (sec. 2).
- **The whole-family structure (sec. 1.2) is new.** It is the cyclic rule that builds family i from the pieces
  Id, E_N^+-, B_1 and C_{k,k+1} + l B_1, for every p.
  - I read JAMM sec. 6.1 (pp.949-950) and grepped the whole JAMM text layer.
  - JAMM sec. 6.1 condenses this preprint's sec. 1.1 (pieces K_1-K_20), states the cycle idea in one paragraph,
    and shows the zigzag scheme (its Fig. 9). It cites this preprint (its ref. 47, "Section 1.1") for "a description
    of the whole of this family".
  - JAMM has no counterpart of: the L_p, M_p, U, V formulas; the 9-digit Table 2; the sec. 1.5 recipe; Hypotheses
    1-2; or the sec. 1.9 Broucke-principle argument for family i.
- **What it gives the project:**
  - an exact mu = 0 generating skeleton for family i (Strömgren i: direct orbits about P1);
  - 9-digit Jacobi-constant anchors (Table 2) for continuation to small mu;
  - the Broucke-principle argument (sec. 1.9) for which half of each C_{m,m+1} family joins family i.
  - It is the family-i control for `#944` X2 and `#956` R9 symmetric-family work, at mu = 0 only. **There are no
    finite-mu numbers in it.**
- **Errors found** (sec. 3): 9 e~(0) cells too large by about 1e-4; two sign slips (rows 16 and 28); one inconsistent
  half-period entry in row 25_2, which is the same in both prints; a probable copy slip in row 31_3; and a set of text
  slips. They are all small. None changes the family structure.
- **JAMM Table 2 against this table:** there are 2 genuine numeric disagreements, and in both of them the preprint
  cell is the wrong one. All other differences are rounding, truncation, notation or a different coordinate
  system (sec. 4).
- **Catalogue implication (PROPOSAL only):** none. No catalogue row cites family i at mu = 0. If a future row or
  control uses family-i generating orbits, cite this preprint's Table 1 (6 decimals) and not JAMM Table 2. Apply the
  e~(0) corrections of sec. 3.1 first.

## 1. Content (the translation has the full text)

- **Sec. 1.1, Table 1.** Table 1 has 39 critical orbits along the pieces K_1-K_20. The index m on k_m is the number of
  arc-solutions in the orbit.
  - For odd m, the orbit enters P2 at T/2. The half-period columns then hold the entry velocity (v1, v2) and the
    system-IV coordinates (a~'_*, e~_*) of 64/2005 eq. (4'.3).
  - Otherwise they hold the second crossing, (x1, y2)(T/2) and (a~, e~)(T/2).
  - The text then names the pieces up to E_{7/6}^+.
  - Examples:
    - Row 1: Id(3), T~ = 1/2, C = 3.466806.
    - Row 12_2: C_12(1), the cusp of the right characteristic. T~ = 2.28129, C = 2.874117.
    - Row 39_2: E_{5/4}(+1), C = 2.744731, e~(0) = 0.47863.
- **Sec. 1.2, the whole family.** The family is an infinite chain of blocks Id((p+1)/p), L_p, (2p+1)Ir((p+1)/p), M_p,
  for p = 1, 2, ...
  - L_p walks down E_{(p+1)/p}(+1), C_{p-1,p}+B_1, E_{(p-1)/(p-2)}(+1), C_{p-3,p-2}+3B_1, ... and back up along the
    (-1) branches. M_p is the mirror walk through C_{p,p+1}, C_{p-2,p-1}+2B_1, ...
  - In compact form, L_p and M_p are written with the sequences U and V. These are built from gamma_1 and gamma_2
    (the C_12 and C_23 arcs of minimum a) and from beta_1 and beta_2 (the B_1 arcs at the same C).
- **Sec. 1.3, Table 2.** C alternates along U^+ and V^+: minima at E_{(m+1)/m}(+1), maxima at C_{m,m+1}(1). U^- and
  V^- have a single minimum, at an Ir orbit.
  - Table 2 gives the four extremal C values to 9 decimals for m = 1-10. All four columns tend to 3.
- **Sec. 1.4.** The characteristics fold into zigzags along B_1 and the curve P_2** (Figs. 8-9, not held). These come
  apart for mu > 0.
- **Sec. 1.5.** The recipe for the arc families S, from Bruno 1990 Ch. VI:
  - eq. (1.16) e*(N^-1); eq. (1.17) cos(psi/2) = |e*|; eq. (1.18) psi and tau' for a < 1;
  - eqs. (1.19)-(1.20) e~ = c|c|/(a(1 - eps'' e)) and a~ = eps a;
  - the a-ranges and the (theta, k, eps, eps'') values for B_1, C_12, C_23, C_{2l+1,2l+2} and C_{2l,2l+1}.
- **Sec. 1.6.** The initial part was first given in Bruno KIAM 66/1993 sec. 2 without proof. Hénon 1997 Ch. 10
  (his Table 10.8) justified it. This preprint extends the cycles to the whole family.
- **Sec. 1.7.** Period and trace:
  - T = 2 pi/(N - 1) on Id; T = 2 pi k and Tr = 2 on E_{(k+1)/k}; the half-period tau adds over the arcs.
  - The signs of Tr = +-inf come from the sign of Tr_1 on the adjacent E_N piece.
  - Hypothesis 1: sgn Tr_v1 = -sgn Tr_1.
  - Hypothesis 2: Tr_v jumps only where v2 = 0 or -2. That is impossible for |a~| < 1.
- **Sec. 1.9.** Broucke's principle requires v1 <= 0 at entry into P2. Eqs. (1.22)-(1.27) then show that only the part
  a < (m/(m+1))^(2/3) of each C_{m,m+1} joins family i.

## 2. Arithmetic checks (`checks_kiam36.py` -> `checks_kiam36.out`)

At mu = 0 the checks use: C = -y2^2 + 2 x1 y2 + 2/|x1| on the symmetry line; C = 3 - v1^2 - v2^2 for the P2 entry
velocity; and (a~, e~) <-> (x1, y2) by preprint 64/2005 eq. (4'.1)-(4'.2) (64 says (2.2) of [1] is misprinted). The
table **truncates** to 5 decimals: x1(0) of row 1 is 0.48074, but 3^(-2/3) = 0.480750. So the tolerance is one full
unit in the last digit, propagated.

- **C from (x1, y2)(0): 35 of the 36 non-collision rows pass.** Row 26_2 cannot be tested, because x1(0) = -0.00094 has
  only two significant digits. For the three P1-collision rows (3, 18_1, 25_2), C = 1/a holds to 1.1e-5.
- **C = 3 - |v|^2 at P2 entry: 12 of the 13 odd-m rows pass**, to within truncation. Row 31_3 cannot be tested
  (sec. 3.4). This check is new; the JAMM digest could not do it, because JAMM has no velocity columns.
- **C at the T/2 crossing, even m:** all rows pass except row 16 (sign slip, sec. 3.2) and row 25_2 (sec. 3.3).
- **Closed forms:**
  - Id rows 1, 2, 8, 9, 10, 20, 21, 22, 36, 37, 38: C = 1/a + 2 sqrt(a) within 3e-5. T~ = 1/(N - 1) within 1.3e-4.
    Tr = 2 cos(2 pi T~) = +-2, as printed.
  - Ir rows 4, 16, 28: C = 1/a - 2 sqrt(a) within 8e-6.
  - The a values are 2^(-2/3), (2/3)^(2/3), (3/4)^(2/3), (4/5)^(2/3), (3/5)^(2/3), (5/7)^(2/3) and (7/9)^(2/3), as printed.
- **Table 2, Id column:** 1/a + 2 sqrt(a) with a = (m/(m+1))^(2/3) reproduces all 10 values within 1e-9.
- **Table 2, E(+0) column:** the ellipse with a(1+e) = 1 (apocentre at r = 1), with C = 1/a + 2 sqrt(a(1 - e^2)),
  reproduces all 10 values within 1e-9.
- **Table 2, E(+1) column:** I solved for the direct ellipse of E_{(m+1)/m} (a fixed) that passes through P2 in the
  rotating frame, with Kepler's equation and a scan in e. For each m = 2 to 10 the solution with the largest C
  equals the printed E(+1) value to all 9 digits. This is an independent check, and it confirms that E(+1) is the
  type-I collision ellipse.
- **Table 2, C_{m,m+1}(1) column:** not checked. It needs the arc-solution machinery (1.17)-(1.18). The m = 1, 2 values
  agree with Table 1 rows 12_2, 33_3 and 24_2.
- **Cross-table:** 13 Table 2 values agree with the Table 1 C column. Table 1 truncates: for example, 3.057531626 is
  printed as 3.057531.
- **Half-period additivity (sec. 1.7):**
  - T~(32_3) - T~(13_2) = 1.05598 is the B_1 arc at C = 1.793529. So the C_12 arc is 2.05598 - 1.05598 = 1.00000,
    which equals T~(6_1) = 1, the single C_12 orbit at the same C. **Exact agreement.**
  - The same subtraction gives the B_1 arc at four other C values. They are consistent, but there is no second
    witness for them.

## 3. Errors in the preprint (all are noted in the translation as translator's notes; none are corrected in the copied tables)

### 3.1 e~(0) too large by about 1e-4 in 9 rows
- Rows 5_1, 7_1, 12_2, 17_1, 19_1, 27_2, 29_1, 33_3 and 35_1. In each, the 4th decimal is one unit too high.
- The row's own x1(0), y2(0) and C, and (for E_N rows) the closed-form Kepler value, give:

| row | printed e~(0) | from x1, y2 | closed form (E_N rows) |
|---|---|---|---|
| 5_1 / 7_1 | -0.41269 / +0.41269 | 0.41259 | 0.412600 (= 2 - 2^(2/3)) |
| 12_2, 33_3 | 0.41828 | 0.41818 | (C_12(1); no closed form) |
| 17_1 / 19_1 | -1.31047 / +1.31047 | 1.31035 | 1.310372 |
| 27_2 | -1.23588 | -1.23576 | -1.235784; the row's own e~(T/2) = -1.23578 |
| 29_1 / 35_1 | -0.78868 / +0.78868 | 0.78858 | 0.788585 |

- C computed from the printed (a~, e~)(0) misses the printed C by 3.3e-5 to 1.25e-4 in these rows. In all other rows
  it agrees within 3.2e-5, except the collision-line row 26_2 (2.9e-4; see sec. 4).
- The rows with an E(+-1) or E(+-0) collision ellipse that are not affected (11_2, 15_2, 23_2, 39_2) agree with the
  closed form within 1e-5.

### 3.2 Two sign slips
- **Row 16, y2(T/2) = +1.14471.** It should be -1.14471, because this is the 5-fold retrograde Ir orbit. Rows 4 and 28
  are negative. As printed, it gives C = 3.0575, the direct circular value.
- **Row 28, a~(T/2) = -0.82548.** It should be +0.82548. Rows 4 and 16 are positive, the row's own (x1, y2)(T/2) give
  +, and JAMM prints +0.8255.

### 3.3 Row 25_2: the half-period point does not have the row's C (UNRESOLVED)
- (x1, y2)(T/2) = (0.05931, 5.69137) and (a~, e~)(T/2) = (0.75483, 1.92141) agree with each other. Both give
  **C = 2.0000**, not the printed 1.484505.
- The row's t = 0 data (P1 collision, a~ = -0.67362, C = 1/a) and T~ = 3.08850 are consistent.
- JAMM prints the same values (0.7548, 1.9214). So this is not a JAMM copying error.
- For comparison, rows 13_2 and 11_2, the same kind of C + B_1 orbit, check to 5e-5.
- I cannot tell which number is wrong without recomputing the B_1 arc at C = 1.484505.

### 3.4 Row 31_3 half-period cells
- The table prints "0, +inf", which is the P1-collision notation of row 14_2. But m = 3 is odd, so these columns
  should hold a finite P2 entry velocity with 3 - |v|^2 = 1.401875.
- The system-IV columns (-0.62574, 2.63460) are finite. JAMM's (2.0701, -2.4830) are finite and put the T/2 point at
  x1 = 0.9999, at P2, like every other odd-m row (sec. 4).
- Probable copy slip from row 14_2.

### 3.5 Text slips (each footnoted in the translation)
- K_6 ends at "E_{3/2}(-1), i.e. orbit 11_2". Orbit 11_2 is E_{3/2}(+1): its C equals Table 2's E_{3/2}(+1), and K_7
  starts there.
- K_16 says "up to orbit 34_3 ... and at orbit 34_3". The first should be 33_3, by the Table 1 Tr column. K_16 also
  calls 30_3 a collision orbit; the collisions in the table are 31_3 and 32_3.
- K_18: "orbit 36, coinciding with orbit 20" should be orbit 22 (Id(4/3)).
- K_19: (N-1)_{-1} for (N-1)^{-1}.
- In the tail list, E_{4/5}(+1) should be E_{4/3}(+1), by the L_p rule with p = 5.
- U_{2k} and V_{2k} are written for U_{2,k} and V_{2,k}.
- "(10)-(13)" is written for (1.10)-(1.13). "(1.20) and (1.21)" is written for (1.14) and (1.15).
- One interval has unbalanced parentheses: [2l/(2l+1))^{2/3}, ...].
- V_{1,k}^- and V_{2,k}^-: the Ir orbit named as the minimum does not match the Ir orbit at the end of the sequence.
- Table 1, row 24_2: the Tr cell reads `+\in`, which means +inf.
- Reference 8 says "Broucke M.R." (it is R. A. Broucke).

## 4. Comparison with JAMM 71:933 Table 2 (p.946; image-read)

- I compared every JAMM cell with the preprint: T~, C, a~(0), e~(0), and a~(T/2), e~(T/2) for the 26 rows with even m
  or no index. That is 208 cells:

| class | cells | what it means |
|---|---|---|
| identical | 131 | |
| rounded to 4 decimals | 39 | mostly a~ columns, for example -0.62996 -> -0.6300 |
| truncated to 4 decimals | 30 | mostly C and e~, for example C(7_1) 2.872078 -> 2.8720, C(37) 3.021678 -> 3.0216 |
| notation only | 6 | the preprint prints 2 where JAMM prints +-2 or -+2 (rows 3, 14_2, 18_1, 25_2, 26_2) |
| **genuine disagreement** | **2** | below |

- **Row 28, a~(T/2):** the preprint has -0.82548 and JAMM has +0.8255. **JAMM is right** (sec. 3.2).
- **Row 35_1, e~(0):** the preprint has 0.78868 and JAMM has 0.7885. **JAMM is right**: the closed form is 0.788585,
  which truncates to 0.7885. 0.7885 is neither a rounding nor a truncation of the preprint's 0.78868.
  - JAMM does not truncate consistently. It also rounds: the preprint's -0.60676 appears as -0.6068 in row 15_2 and
    as -0.6067 in row 30_3.
  - So the other e~(0) cells of sec. 3.1 cannot show whether JAMM copied the preprint error. For example, JAMM 29_1
    -0.7886, 12_2 and 33_3 0.4182, and 17_1 and 19_1 1.3104 are each a rounding of the correct value and also a
    truncation of the wrong one. Only row 35_1 discriminates.
- **Half-period columns of the 13 odd-m rows are not comparable.** The preprint gives system IV (a~'_*, e~_*). JAMM
  gives system III (a~, e~) at T/2.
  - Check (output sec. 11): in all 13 JAMM rows, x1 = a~|2 - |e~|| at T/2 is 1.0000 +- 0.0002. So JAMM evaluates
    these cells at P2 = (1, 0), the point of entry.
  - JAMM's |a~(T/2)| is not always the |a~(0)| of the arc. For example, row 6_1 has 0.6657 against 0.5576, and row 18_1
    has 23.872. So the velocity (y2) convention JAMM uses at P2 stays open. I did not check further.
- Tr column: it is the same apart from bracket style. JAMM fixes the `+\in` typo of row 24_2. JAMM has no Tr_v column.
- JAMM's 25_2 T/2 entry repeats the preprint's inconsistent values (sec. 3.3).
- The JAMM digest's sec. 2.1 found row 26_2 off by 0.003 in C from (a~, e~). The preprint's 5-decimal e~(0) = -1.99859
  explains this: with it, C(a, e) = 1.40159 against 1.401875. This is still the collision-line sensitivity, not an
  error.

## 5. Citation mining (10 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| ref | work | status |
|---|---|---|
| [1] | Bruno & Varin, KIAM 10/2005, on families of periodic solutions | HELD (`bruno-varin-2005a-...`, TeX source + translation) |
| [2] | Bruno & Varin, KIAM 67/2005, family h for small mu | HELD as TeX source + translation (`bruno-varin-2005b-...`); its Tables 1-8 are absent (wanted row 17) |
| [3] | Bruno & Varin, KIAM 64/2005, family h for big mu | HELD as TeX source + translation (`bruno-varin-2005c-...`); tables absent (row 17) |
| [4] | Bruno, KIAM 66/1993, single periodic solutions, Sun-Jupiter | not held; already on wanted row 59. The first description of the initial part of family i |
| [5] | Hénon 1997, LNP m52 | HELD (`henon-1997-...`). Ch. 10 and Table 10.8 are the initial part of family i |
| [6] | Bruno 1990, Nauka (Russian) | content HELD in another form: the English edition Bruno 1994, de Gruyter (`bruno-1994-...`). The Nauka original is not held; page and figure numbers cited here (Fig. 20, 40, 42, 45, 74, 90; pp. 109, 128-135) may differ in the English edition |
| [7] | Hénon 2001, LNP m65 | HELD (`henon-2001-...`) |
| [8] | Broucke 1968, JPL TR 32-1168 | HELD (`broucke-1968-...`). Broucke's principle is on p.21 |
| [9] | Perko 1981, SIAM J. Appl. Math. 41:181 | HELD (`perko-1981b-...`) |
| [10] | Perko 1981, Celest. Mech. 24:155 | HELD (`perko-1981-...`) |

- No new candidates. The only reference not held, [4], is already on row 59.
- **Proposal for row 59:** mark Bruno-Varin 36/2006 as received (TeX source; full translation). Keep the row open for
  the other preprints. If the filer wants it, add a note that the 2008 Astron. Vestnik paper and 51/2007 (family c)
  are the remaining Bruno-Varin items.

## 6. Files (filed beside the zip; the `src/` conversions are folded into the `.utf8.txt` file)

- `translation.tex`, `translation.pdf`: the full English translation, 15 pp. (tectonic). It has the licence line, the
  translator's notes, both tables, and the references.
- `src/*.TEX`: the CP866 -> UTF-8 conversions of the five source files.
- `checks_kiam36.py`, `checks_kiam36.out`: the checks of secs. 2-4. The script reads JAMM Table 2 from the corpus PDF
  at run time.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
