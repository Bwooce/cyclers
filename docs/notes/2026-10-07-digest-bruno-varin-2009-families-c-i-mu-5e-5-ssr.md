# Digest: Bruno & Varin 2009, "Families c and i of Periodic Solutions of the Restricted Problem for µ = 5×10^-5" (Solar System Research 43:26) (#960 batch 34)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, Moscow), "Families c and i of Periodic
Solutions of the Restricted Problem for µ = 5 × 10^-5", Solar System Research 43(1):26-40 (2009),
doi 10.1134/S0038094609010031. It is the Pleiades English translation of Astronomicheskii Vestnik 43(1):28-43 (the
running head spells the author "Bryuno"). Received 30 April 2008. RFBR grant 08-01-00082.
- Supplied file `bruno2009_1.pdf`, 15 pp., md5 c0de4fb82dc2ac5ce02cf8d4ad6a60a2. Digital, with a clean text layer.
- **Proposed corpus filename:**
  `bruno-varin-2009-families-c-i-periodic-solutions-restricted-problem-mu-5e-5-solar-system-research-43-26-doi-10.1134-S0038094609010031.pdf`
  - The family-h paper of this batch (SSR 43:158) is then `bruno-varin-2009b-...`.
- **Transcription:** `bruno-varin-2009-families-c-i-mu-5e-5-tables.yaml`, for data/sources/. It holds:
  - Table 1 (family c, 8 orbits, 14 columns);
  - Table 2 (family i, 65 orbits, 14 columns);
  - Table 3 (9 crossings of e~ = 1.4, a~ to 12 decimals).
  - 1031 cells; 7 are blank (all in Table 2 row 3).
- **How I read it:**
  - I read the whole text from the text layer, and the tables, the Fig. 6 axis and Fig. 9 on the page images.
  - Tables 1-3 (pp.27, 37-40): I copied the cells mechanically from the text layer (`build_yaml.py`), with asserted
    row counts. I compared every cell with 250 dpi page renders. The text layer agreed with the image everywhere.
  - Checks: `checks_bigmu_ci.py` -> `checks_bigmu_ci.out` (secs. B0-B5, C, D); `explore_hard_rows.py` ->
    `explore_hard_rows.out`; `table3_family_c_e14.py` and `table3_family_i_e14.py` with their `.out` files.
    Integrations use DOP853 with rtol 1e-12 to 1e-13.
- Wanted list: **row 55** (current numbering; the brief's "row 59" is an older number). It lists "Bruno & Varin (2008),
  'Families c and i ... at mu = 5e-5', Astron. Vestnik 42". **This paper is that item.** The correct citation is
  Astron. Vestnik 43(1):28-43 (2009) = SSR 43(1):26-40.

## 0. Verdict

**This is the first finite-µ companion of the held generating-family preprints KIAM 51/2007 (family c, µ = 0) and
36/2006 (family i, µ = 0).** It computes both families at µ = 5×10^-5, close enough to µ = 0 that the orbits look the
same.
- **Table 1** has 8 critical orbits of family c. It runs from L1 (row 1) to the locally double orbit on family h
  (row 8).
- **Table 2** has 65 critical orbits of the initial part of family i, "one cycle ahead" of the µ = 0 table: T~ runs to
  5.0, against 4 in 36/2006.
- **Table 3** has 9 points where the right-hand characteristic of family i (k = 1-8) and of family c (k = 9) cross
  e~ = 1.4. a~ is given to 12 digits. These show the predicted zigzags accumulating on family c.
- **The data are good.**
  - C from the printed (x1(0), y2(0)) agrees in 67 of the 69 rows that do not start at P1 (c6, i3, i25 and i37
    start there). The two exceptions, i24 and i43, are discussed in sec. 2.3.
  - **L1 (Table 1 row 1) reproduces from first principles:**
    - x(L1) = 0.974675 (printed 0.97468) and C(L1) = 3.005756 (exact);
    - from the L1 linearisation: T/2π = 0.47395, Tr~ = 12.046 and Tr~_v = 1.956. The printed values are 0.47395,
      12.045 and 1.956.
  - Nine integrated orbits close at the printed T~, T/2 state and C (sec. 2.3). One of them, c2, needs a 2e-5
    change of y2(0).
  - **Table 3 reproduces to 1e-12** for k = 1, 2, 3, 4, 6 and 9 (sec. 2.4).
- **How it relates to the µ = 0 tables (sec. 3):**
  - Most µ = 0 rows of 36/2006 Table 1 have a recognisable µ = 5e-5 counterpart.
  - The circular and near-circular orbits (Id, Ir, E_N) agree to |ΔC| ≤ 1.3e-3.
  - Orbits built from arcs near P2 (B_1, C_{k,k+1}) move by up to 0.07 in C, and their critical orbits split into
    pairs or quadruples. For example, 12_2 becomes rows 15-18.
  - This paper confirms four corrections made in the held digests:
    - KIAM 51 orbit 1: T~ max = 1.406729, not 3/2 (it is stated in the text);
    - KIAM 51 orbit 2: x1(0) = 0, not 1 (Table 1 row 6 prints 0.00000);
    - KIAM 36 row 16: the sign of y2(T/2) (row 23 here is negative);
    - KIAM 36 row 28: the sign of a~(T/2) (row 44 here is positive).
  - w1 is printed positive throughout, as in the KIAM 51 figures (not its tables).
- **Problems found (kept as printed; noted in the YAML):**
  - Table 2 row 61: the sign of y2(T/2) is probably wrong.
  - Table 2 row 43: the printed state, T~ and C are not mutually consistent at 1e-3 (unresolved).
  - Table 2 row 24: C is off by about 1.6e-4, or the row is very sensitive (unresolved).
  - Σ in Table 2 runs about 0.001 low from row 7 on.
  - About 8 text slips (sec. 4).
- **Catalogue implication (PROPOSAL only):** none. No catalogue row is a family-c or family-i orbit at small µ.
  - Family c here is the planar Lyapunov family of L1 at µ = 5e-5, continued to large size. Family i is the direct
    P1-satellite family.
  - Both are control data for `#944` X2 / `#956` R9 symmetric-family work at small µ, the regime of Saturn-moon or
    Jupiter-moon systems.
  - Proposal for data/sources/: file the YAML as is. Use Table 3 and the circular-type rows of Table 2 as
    high-precision controls.

## 1. Content

- **Introduction.**
  - It continues the programme of BV 2008a (SSR 42:163; not held), which studied the µ = 0 families c and i.
  - New: **uniformization.** The parameter Σ is the total variation of C along the family from a reference orbit M0
    (eq. 1): Σ = ±(Σ_j |C_j - C_(j-1)| + |C - C_k|) over the C extrema between M0 and M. It is used to compare
    different µ.
  - Modified trace: Tr~ = Tr if |Tr| ≤ 2, else (1 + log2|Tr|) sgn Tr. The same holds for Tr_v.
- **Family c (Table 1, Figs. 1-6).**
  - It starts at L1 and ends as a locally double orbit on family h.
  - At µ = 0 it has two parts:
    - c′, which is Hill's family c (Hénon 1969, with vertical traces in Hénon 1974);
    - c″, a part of the segment family B_1 (BV 2008a).
  - **On c′, T~ rises to 1.406729. On c″ it falls from that value.**
  - At µ = 5e-5:
    - rows 1-5 perturb c′ and rows 5-8 perturb c″;
    - Σ = 2.9767666 - C, with M0 = row 5, the maximum of T;
    - row 6 is a collision orbit, the same as orbit 2 of BV 2008a's family c;
    - row 8 is the locally double orbit on h.
  - Figs. 1-2 (in a~, e~ and w1, y2) are close to the right-hand sides of BV 2008a Figs. 5 and 3. T~(Σ) (Figs. 3-4)
    almost coincides with c″ at µ = 0.
  - The traces (Figs. 5-6) perturb the µ = 0 values Tr~ = +inf and Tr~_v = -inf on c″. On row 8, Tr~ = 2 and
    Tr~_v < 2.
  - Figs. 4 and 6 are "fully analogous" to Hénon 1969 (Table 2, Fig. 2) and Hénon 1974 (Figs. 1-2, Table 2),
    "including the number and position of the critical orbits".
- **Family i (Tables 2-3, Figs. 7-17).**
  - These are direct circular orbits of infinitely small radius around P1, continued.
  - Figs. 7-8 show the left-hand and right-hand characteristics in (a~, e~). They are close to the generating ones,
    with bifurcations near e~ = 1, as in Bruno 1990 Ch. VIII sec. 3. Bruno 1993a computed parts of family i at µ_J.
  - Unlike µ = 0, both characteristics are connected curves without self-intersections. The upper and lower parts of
    the right-hand characteristic almost coincide with family c.
  - **Fig. 9 and Table 3:** zigzags of the right-hand characteristic for e~ > 1. The crossings with e~ = 1.4
    accumulate on family c (k = 9).
  - Fig. 10 is the (w1, y2) view near P2. It almost coincides with BV 2008a Fig. 15 (µ = 0).
  - **Figs. 11-13 (C, e~ and T~ against Σ):** the µ = 5e-5 and µ = 0 curves are very close on the initial part, then
    drift apart.
    - BV 2007 (JAMM): for each µ > 0 there is a Σ(µ) beyond which the generated family i differs drastically from the
      generating one.
  - **Traces (Figs. 14-17):**
    - The instability intervals Tr~ < -2 near second-order resonances are very small. Table 2 shows only one orbit
      with Tr = Tr_v = -2 instead of four, at rows 1, 8 and 28.
    - Near rows 18, 43 and 50 (the singular parts) not all ±2 crossings were found.
    - On Σ ∈ [17, 19.5] (eq. 2) the monodromy matrix lost all significant figures; |Tr| and |Tr_v| > 1e12 there.
    - Fig. 15 was used to correct the µ = 0 Fig. 14, and Fig. 17 to correct Fig. 16.
- **Conclusions for µ = 0** (two lists, for Tr and Tr_v; orbit labels are 36/2006 Table 1 numbers):
  - Tr drops from +2 to -2 and then jumps to +inf at the joints between regular and singular parts (11_2, 23_2,
    39_2).
  - On extremal orbits where the right-hand characteristic has a spinode, Tr jumps +inf -> -inf -> +inf (12_2, 24_2,
    34_3).
  - Jumps through ±inf also occur on non-extremal orbits ("additional orbits 14′ and 26′"). These may sit at the point
    of B_1 with the lowest a (curve f, Bruno 1994 Ch. IV sec. 2, Fig. 41).
  - Tr_v drops from 2 to -inf at 5_1, 7_1, 11_2, 15_2, 17_1, 19_1, 23_2, 27_2, 29_1, 35_1 and 39_2. Tr_v jumps
    -inf -> +inf at 12_2, 24_2 and 34_3.
  - **These refine 36/2006 sec. 1.7** (Hypotheses 1-2 and the ±inf sign rules).
- The preliminary publication is KIAM Preprint 22/2008.

## 2. Checks (`checks_bigmu_ci.out`, secs. B0-B5, C, D)

### 2.1 Frame, C, a~, e~, w1 (B0, B1)
- The frame is the same as in the family-h paper: P1 at the origin; y2 = dx2/dt + x1;
  C = x1² - 2µx1 + µ + 2(1-µ)/r1 + 2µ/r2 - v².
- **Positive control (B0):** L1 at µ = 5e-5. x = 0.974675, C = 3.005756 (printed 3.005756).
  - The in-plane linearisation gives T/2π = 0.47395 (printed 0.47395) and Tr = 2cosh(λT) = 2114, so
    Tr~ = 12.046 (printed 12.045).
  - The vertical frequency gives Tr_v = 2cos(ω_z T) = 1.956 (printed 1.956).
  - This is independent of the tables. It fixes the frame, the C convention and the modified-trace definition.
- **C from the t = 0 state:** all 7 testable rows of family c pass (row 6 starts at P1). Family i passes in 60 of the
  62 testable rows; the exceptions are row 24 (-2.3e-4 against a tolerance of 1.5e-4) and row 43 (-8.7e-4). Rows 3,
  25 and 37 start at P1. Row 38 (x1 = -0.00095) passes only because its rounding tolerance is large.
- **C from the T/2 state:** this passes except in rows 6, 24, 26, 43, 45 and 61. Rows 6, 24, 26 and 45 have their T/2
  point within 0.006 of P2. Rows close to P1 or P2 at the 1e-5 level (37-38, 46-47, 50, 59-62; c6) are
  rounding-limited and cannot be judged.
- **a~, e~ and w1:** all reproduce from the printed states, except row 61 (below).
  - w1 = µ/(1 - x1(T/2)) is printed **positive** in every row.
  - This differs from KIAM 51/2007, whose eq. (1.3) and tables carry negative w1 and whose figures positive.

### 2.2 Σ (B2)
- **Family c:** Σ = 2.9767666 - C reproduces all 8 cells. Row 5 is 0, as printed.
- **Family i:**
  - Rows 1-6 agree with Σ = ±(C - C(row 4)): -3.466, -3.174, -1.587, 0.000, 0.364 and 2.869. The cells are
    truncated, not rounded: -3.466518 is printed -3.466 and 2.869764 is printed 2.869.
  - **From row 7 on, the printed Σ is about 0.001 below the smallest value the printed C allows.** Row 7: 3.172
    against at least 3.1732. Row 8: 3.253 against at least 3.2541. Rows 9 and 10: 3.428 and 3.459 against 3.4290 and
    3.4600.
  - An untabulated C extremum can only raise Σ, so this is not an extremum effect.
  - It may come from a Σ summed on a coarse grid that clips a C maximum near rows 6-7. I did not trace it further.
    It is a small offset in an auxiliary column.

### 2.3 Integration (B3, C; `explore_hard_rows.out`)
- I used Newton on (y2(0), T/2) with x1(0) fixed. These rows reproduce T~, x1(T/2), y2(T/2) and C to the printed
  precision (|Δy2(0)| ≤ 5e-6, |ΔC| ≤ 3e-6): family c rows 3, 4 and 7; family i rows 1, 8, 23, 28 and 44.
- **Family c row 2** needs Δy2(0) = -1.85e-5 and gives y2(T/2) 0.954458 against the printed 0.95443. It is the most
  unstable row (|Tr| about 1800), so a 5-decimal x1(0) moves y2(0) by this much. It is not an error.
- **Family i row 8:** T~ = 1.499975 against the printed 1.50000 (2.5e-5). The other half-integer cells, rows 28
  (2.50000) and 64 (4.50000), may also be idealised. Row 28 integrates to 2.500005.
- **Rows 23 and 44** settle two of the 36/2006 sign questions (sec. 3).
- **Row 6:** C, T~ and x1(T/2) reproduce. y2(T/2) = 0.5911-0.5912 against the printed 0.59100 (1e-4). This is minor.
- **Row 24:** both corrections give C = -0.31589 to -0.31590, against the printed -0.315734, and x1(T/2) = 0.99864
  against 0.99849. The row's T/2 point is 0.0014 from P2. Either the C cell is off by 1.6e-4, or the row is too
  sensitive to judge from a 5-decimal start. Unresolved.
- **Row 43:** the printed start is not periodic at the 1e-3 level (dx1/dt = 3.3e-3 at the crossing).
  - Corrections give C = -0.5571 to -0.5573 and T~ = 2.9957, against the printed -0.554105 and 2.99262.
  - C from the printed t = 0 and T/2 states gives -0.55497 and -0.55542.
  - The text names orbit 43 among the singular parts. Unresolved: one or more of x1(0), y2(0), T~ and C is off.
- **Row 45:** too sensitive. The two corrections differ by 5e-3 in C, so no verdict.
- **Rows 59-62:** at T/2 these pass 2e-5 from P2. My integrations did not converge, so no verdict.
- **Row 61, sign of y2(T/2):** printed +1.44785. Rows 59, 60 and 62 are negative. The row's own e~(T/2) = -2.09633
  needs y2(T/2) < 0, because e~ = x1 y2|y2|/(1-µ) and x1 > 0. **Probable sign slip.**

### 2.4 Table 3, independent reproduction (`table3_family_c_e14.out`, `table3_family_i_e14.out`)
- On the line e~ = 1.4 (x1 > 0), x1 = 0.6 a~ and y2 = +(1.4(1-µ)/x1)^(1/2). I integrated and looked for a
  perpendicular y = 0 crossing (dx1/dt = 0).
- **k = 9 (family c):** root a~ = 0.892589912327; printed 0.892589912326. T~ = 1.282842. The T/2 point is
  x1 = 0.99997668, 2.3e-5 from P2 (w1 = 2.14). This lies between Table 1 rows 5 and 6, as it should.
- **k = 1, 2, 3, 4, 6 (family i):** roots agree with the printed a~ to 4e-13 to 9e-13 (the last printed digit). The
  T/2 point is at the 3rd or 4th y = 0 crossing, with T~ = 4.035, 3.090, 2.225, 2.308 and 4.488.
- **k = 5, 7, 8:** my search (the first 24 crossings, t ≤ 70, brackets of ±2e-6 / ±1e-9) found no root. These are
  not called errors: k = 7 and 8 lie within 3e-9 of the family-c root, where the T/2 pass is within 1e-5 of P2.
- **The zigzag accumulation is real.** The spacing to k = 9 shrinks: 1.4e-3, 4.6e-4, 1.0e-4, 3.9e-5, 1.1e-5, 5.2e-6,
  2.8e-9, 1.1e-9.

## 3. Relation to the µ = 0 generating tables (out sec. D)

The paper gives no row map. I paired rows by (x1(0), y2(0), T~, C) and orbit type, against the 36/2006 Table 1
values in the held translation.

| µ = 5e-5 (this paper) | µ = 0 (36/2006 Table 1) | ΔC | note |
|---|---|---|---|
| i1, i2, i3, i4 | 1 Id(3), 2 Id(2), 3 P1-collision, 4 Ir | 2.6e-5, -1.8e-4, 1.3e-4, 3.1e-4 | i2 moves 0.02 in x1: it is the 2:1-resonance critical orbit, not the circle |
| i5, i6 | 5_1, 7_1 (P2 entries) | 0.062, -0.002 | near-P2 pieces move most |
| i7, i8 | 8 Id(2), 9 Id(3/2) | -1.3e-3, 8e-5 | i8 T~ = 1.5 as at µ = 0 |
| i11-i14 | 11_2 | 0.061-0.063 | one µ = 0 orbit splits into 4 critical orbits |
| i15-i18 | 12_2 (spinode) | -8e-4 | splits into 4 (Tr and Tr_v = ±2 pairs) |
| i21-i22 | 15_2 | 0.007 | |
| **i23** | **16 Ir(2)** | 6e-4 | **y2(T/2) = -1.14935: confirms the 36/2006 row 16 sign fix** |
| i24; i39-i42 | 17_1 | 0.035; 3e-4 | i39-42 match the t = 0 state and C of 17_1 closely, but T~ = 2.86 against 2; i24 matches T~ |
| i25 | 18_1 (P1 collision) | -0.013 | |
| i26, i27, i28 | 19_1, 20 Id(2), 21 Id(5/2) | 2e-3, -9e-4, 9e-5 | |
| i31-i34 | 23_2 | 0.056-0.074 | |
| i35-i36 | 24_2 (spinode) | 2.5e-3 | |
| i37 | 25_2 (P1 collision, the row with the inconsistent T/2 point) | 1.8e-3 | T~ 3.01165 against 3.08850; the T/2 point is near P1 (0.00108) against (0.05931, 5.69137). **It does not settle the 25_2 question** |
| i38 | 26_2 | -2.4e-4 | x1(0) -0.00095 against -0.00094 |
| i43 | 27_2 | 4e-4 | i43 itself is inconsistent (sec. 2.3) |
| **i44** | **28 Ir(3)** | 1e-3 | **a~(T/2) = +0.82531: confirms the 36/2006 row 28 sign fix** |
| i45, i46, i47, i48, i49 | 29_1, 33_3, 34_3, 35_1, 36 Id(3) | 0.025, 1.4e-3, 0.049, 2.8e-3, -7e-4 | |
| i9-10, i19-20, i29-30, i50-i65 | none | | i50-65 have T~ = 3.5-5.0, beyond the 39-row µ = 0 table ("one cycle ahead") |

- **Family c against KIAM 51/2007 Table 1:**
  - c5 (Σ = 0, T~ max 1.33779) <-> the c′/c″ junction. At µ = 0, T~ max = 1.406729 = 4.41937/π. That is Hénon 1969's
    T/2 at Γ -> -inf, quoted in this paper's text. **This confirms the KIAM 51 digest's correction of the printed
    T~ = 3/2.**
  - c6 (P1 collision, C 1.397662) <-> orbit 2 (C 1.401879). This paper prints x1(0) = 0.00000, so **the 51/2007
    "x1(0) = 1" is a misprint for 0, as that digest found.** The T/2 values also correspond: y2(T/2) -1.57440 against
    -1.57577; a~, e~(T/2) 2.08847, -2.47881 against 2.070245, -2.483035; |w1| 2.51256 against 2.51823.
  - c8 (C -0.943718, Tr~ 2.000, Tr~_v 1.967) <-> orbit 3 (C = -1, the retrograde circle a = 1, which ends on h).
  - c1-c4 lie on c′ (Hill's c). Rows 2-4 are the vertical critical orbits (Tr~_v = 2, 2, -2). Their µ = 0 analogues
    are in Hénon 1974 (not held).
- **Summary:** at µ = 5e-5 the families stay close to the generating ones where the orbit stays away from P2. Pieces
  built from arcs near P2 shift by up to 0.07 in C, and each µ = 0 critical orbit becomes 2 or 4 critical orbits,
  one for each Tr and Tr_v crossing of ±2. The paper's Figs. 11-13 show the same thing as curves.

## 4. Text slips (all small)

- p.40: "Figure 16 shows the corresponding dependence ... for µ = 5×10^-5" should be Fig. 17.
- p.40: "Base on the results of calculations of the planar trace Tr~_v" should be "vertical trace".
- p.40: "In the interval (5)" should be eq. (2), Σ ∈ [17, 19.5]. There is no eq. (5).
- p.32, Fig. 9 caption: "(see Table 4)" should be Table 3.
- p.31, Fig. 6: the axis label at -0.01 is printed "0.01" (read on the page image).
- p.27: "Bryuo (1994a)".
- References: BV 2007 JAMM is given as "pp. 993-960". It is 933-960.
- References: "Bruno and Varin, 2009b ... for Small µ, Sol. Syst. Res. 43(2)" is the big-µ paper (SSR 43:158, this
  batch); the title is copied from 2009a.
- References: BV 2008a is "pp. 163-185" here and "163-186" in the family-h paper.

## 5. Citation mining (13 references)

Checked with `ls cyclers_pdf/papers | grep -i` (Bruno/Bryuno/Brjuno, Varin, Hénon) and grep of CORPUS_INDEX.md and the
wanted list.

| reference | status |
|---|---|
| Bruno 1994a, de Gruyter book | HELD (`bruno-1994-...`) |
| Bruno 1993a, KIAM 66/1993 (one-fold periodic solutions, Sun-Jupiter) | not held; wanted **row 55** |
| Bruno 1993b, KIAM 67/1993 (two-fold periodic solutions, Sun-Jupiter; instability intervals at µ_J) | not held; wanted **row 55** |
| Bruno 1996, KIAM 93/1996 (zero-fold and retrograde) | not held; wanted **row 55** |
| Bryuno & Varin 2007, JAMM 71(6):933-960 | HELD (`bruno-varin-2007-...-jamm-71-933-...`) |
| Bryuno & Varin 2008a, "On families of periodic solutions ...", SSR 42(2):163-185 | **not held; not on the wanted list.** It is cited here for the µ = 0 families c and i, the w1 formula (8), and Figs. 3, 5, 7 and 15. **New candidate** (also proposed in the family-h digest) |
| Bruno & Varin 2008b, KIAM Preprint 22/2008 (Families c and i ... µ = 5e-5) | not held. It is the preprint form of this paper. Low priority |
| Bruno & Varin 2009a, "Family h ... small µ", SSR 43(1) | **not held; new candidate**, high value (likely the journal form of KIAM 67/2005 with its tables; wanted row 17) |
| Bruno & Varin 2009b, SSR 43(2) (mis-titled; it is the big-µ paper) | received in this batch (SSR 43:158) |
| Bruno 1994b, "Singular perturbations in Hamiltonian mechanics", in Seimenis (ed.), Hamiltonian Mechanics, Plenum | not held; not on the wanted list (already proposed by the KIAM 51 digest) |
| Hénon 1969, A&A 1:223 | HELD (`henon-1969-...`) |
| Hénon 1974, "Vertical stability ... II. Hill's case", A&A 30:317-321 | not held; not on the wanted list. **It is now cited by two held papers (KIAM 51 and this one) for the c′ vertical traces. Propose a low-priority row** |
| (in text) Bruno 1990, Nauka | content HELD as the English Bruno 1994 |

- **Proposals:**
  - Row 55: take out "Bruno & Varin (2008) Families c and i ... Astron. Vestnik 42". It is received as SSR 43(1):26-40
    (2009); the original is Astron. Vestnik 43(1):28-43.
  - Add new rows for SSR 42(2):163 (2008) and SSR 43(1):2 (2009), and a low-priority row for Hénon 1974.

## 6. Files (this folder)

- `bruno-varin-2009-families-c-i-mu-5e-5-tables.yaml`: Tables 1-3.
- `checks_bigmu_ci.py` -> `checks_bigmu_ci.out`: secs. B0-B5, C (family-i rows) and D cover this paper.
- `explore_hard_rows.py` -> `explore_hard_rows.out`: family-i rows 6, 24, 26, 43, 45, 59, 61 and 62.
- `table3_family_c_e14.py` / `.out` and `table3_family_i_e14.py` / `.out`: the Table 3 reproduction.
- `build_yaml.py`: the text-layer copy (shared with the family-h YAML).

*Filed as `cyclers_pdf/papers/bruno-varin-2009b-families-c-i-periodic-solutions-restricted-problem-mu-5e-5-sol-syst-res-43-1-26-doi-10.1134-s0038094609010031.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/bruno-varin-2009-families-c-i-mu-5e-5-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
