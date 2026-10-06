# Digest: Bruno & Varin 2009, "Family h of Periodic Solutions of the Restricted Problem for Big µ" (Solar System Research 43:158) (#960 batch 34)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, Moscow), "Family h of Periodic Solutions of the
Restricted Problem for Big µ", Solar System Research 43(2):158-177 (2009), doi 10.1134/S0038094609020099. It is the
Pleiades English translation of Astronomicheskii Vestnik 43(2):167-186. Received 30 April 2008. RFBR grant 08-01-00082.
- Supplied file `bruno2009_3.pdf`, 20 pp., md5 6f703c020ee599a96729614be7c25a7c. Digital (FrameMaker), with a clean text
  layer.
- **Proposed corpus filename:**
  `bruno-varin-2009b-family-h-periodic-solutions-restricted-problem-big-mu-solar-system-research-43-158-doi-10.1134-S0038094609020099.pdf`
  - The companion c/i paper (SSR 43:26, the same batch) is proposed as `bruno-varin-2009-...`, because it appeared
    first (issue 1).
  - SSR 43(1):2-25 (family h, small µ) is not held. If it is filed later, it needs another letter.
- **Transcription:** `bruno-varin-2009-family-h-big-mu-tables.yaml`, for data/sources/. It holds all 87 rows of
  Tables 1-3 with 14 columns each. It also holds the three in-text "more precise values" mini-tables and the Hénon
  correspondence rows.
- **How I read it:**
  - I read the whole text from the text layer, and the figure captions and key passages on the page images.
  - Tables 1-3 (pp.160-163, 174-175): I copied the cells mechanically from the text layer (`build_yaml.py`). Blank
    cells are handled by row and asserted. Then I compared every cell with 250 dpi page renders, row by row:
    87 rows x 13 value columns = 1131 cells, of which 46 are blank in the print. The text layer agreed with the
    image everywhere. That includes the odd cells: the
    doubled minus "––1.493", the swapped row 31 and the blanks.
  - I also read the mini-tables on pp.158 and 165 and the Hénon correspondence table (p.167) on the images.
  - Checks: `checks_bigmu_ci.py` -> `checks_bigmu_ci.out`; `explore_hard_rows.py` -> `explore_hard_rows.out`;
    `collision_rows.py` -> `collision_rows.out`. All are kept in this folder. Integrations use DOP853 with
    rtol 1e-12 and atol 1e-13.
- Wanted list: **row 17**. It asks for the printed KIAM 67/2005 and 64/2005 preprints, because the TeX sources have no
  tables. This paper supplies the 64/2005 half (µ = 0.3, 0.4, 0.5). See sec. 6.

## 0. Verdict

**This is the journal form of KIAM Preprint 64/2005. Its Tables 1-3 are the missing Tables 6-8 of that preprint:
critical orbits of family h at µ = 0.3 (29 orbits), 0.4 (27) and 0.5 (31).**
- The paper says so itself (p.177): "Originally, family h for µ = 0.3, 0.4, and 0.5 was presented in the preprint
  (Bruno and Varin, 2005)". Its reference is Preprint 64.
- The orbit counts, the mini-tables and the stability-interval lists match the preprint word for word (sec. 3).
- **What it gives the project:**
  - six-decimal C and five-decimal states for 87 critical orbits of the retrograde-satellite family h at large µ;
  - the first held full tables at µ = 0.3 and 0.4 (Hénon & Guyot 1970, held, give 9 of the planar critical orbits
    per µ, and they agree; sec. 4.1), and a second, independent six-digit source for µ = 0.5;
  - a symmetric-family control for `#944` / `#956` at large mass ratio (binary-star-like µ).
- **The data are good.**
  - C from the printed (x1(0), y2(0)) agrees with the printed C in all 87 rows. The worst case is 4.8e-5, within the
    rounding of the 5-decimal state.
  - Ten integrated orbits close and reproduce T~, x1(T/2), y2(T/2) and C to the printed precision (sec. 2.3). One of
    them is row 31, started as printed.
  - At µ = 0.5, all 13 matching Hénon 1973 orbits agree to 8e-6 or better in x0, x1 and C (sec. 4).
  - At µ = 0.3, 0.4 and 0.5, all 27 family-h critical orbits of Hénon & Guyot 1970 (tables pp.369-371) have a row
    here within 1e-4 in x1(0), 2e-5 in x1(T/2) and 2e-5 in C. Their stability type (a = ±1) matches the printed
    Tr~ = ±2.00 in every case (sec. 4.1). This independently confirms that the rows are critical orbits, and not
    only periodic.
- **Frame:** P1 is at the origin and P2 at (1, 0). y2 = dx2/dt + x1 is a canonical momentum.
  C = x1² - 2µx1 + µ + 2(1-µ)/r1 + 2µ/r2 - v², which is the barycentric Jacobi constant plus µ(1-µ).
  T~ = T/(2π). The integrations decide between T/(2π) and (1-µ)T/(2π) clearly (sec. 2.3).
- **Errors and weak cells found (all kept as printed; noted in the YAML):**
  - µ = 0.4 row 16, e~(0): it is printed "––1.493" (a doubled minus). The value is -1.4927.
  - µ = 0.3 rows 6 and 7, e~(T/2): the sign is probably wrong (printed +157.6128 and +204.0113; they should be
    negative).
  - Near-P2 y2(T/2) cells (µ = 0.3 rows 6-8, µ = 0.4 rows 17-19) are accurate to only 2-3 significant figures.
  - In the collision rows, x1(0) is good to about 1e-3 only. The cells end in "00".
  - µ = 0.5 row 31 is listed from its other crossing.
  - Text slips: "Table 1 of Hénon (1965b)" should be Hénon (1973); "Fig. 11 of Hénon (1965b) corresponds to Fig. 23"
    should be Fig. 22; "Figure 22 is similar to Fig. 12" should probably be Fig. 20.
- **Catalogue implication (PROPOSAL only):** none directly. No catalogue row is a family-h orbit at µ ≥ 0.3.
  - If a large-µ symmetric-family control is built, use this YAML for µ = 0.3 and 0.4.
  - For µ = 0.5, use it together with Hénon 1973 Table 1. They agree to 1e-5.
  - Treat the near-P2 y2(T/2) cells and the collision-row x1(0) with the caveats above.
  - Proposal for data/sources/: file the YAML as is.

## 1. Content

- **Introduction (p.158).** The family h begins with retrograde circular orbits of infinitely small radius around
  P1. BV 2009 (SSR 43(1):2-25) covered µ in [0, 0.2]. This paper completes µ = 0.3, 0.4 and 0.5, "following the
  methods and terminology of Bruno and Varin (2008, 2009)".
- **Table columns:** k; x1(0), y2(0); x1(T/2), y2(T/2); T~; C; Tr~, Tr~_v (modified traces); a~(0), e~(0);
  a~(T/2), e~(T/2); w1(T/2).
  - The text also lists "v1(T/2), v2(T/2) (the entrance velocity at P2)", but no table has those columns. The
    P2-entry rows print x1(T/2) = 1 and y2(T/2) = ∞.
- **µ = 0.3 (Table 1, Figs. 1-8).**
  - The orbits and the figures differ only slightly from µ = 0.2.
  - The local maxima of Tr~ now exceed 2. In µ = 0.2 (BV 2009, Fig. 29) they were below 2.
  - It is linearly stable in both directions in 10 intervals: from the start to orbit 1, then (2,3), (6,7), (11,12),
    (13,14), (16,17), (21,22), (23,24), (25,26) and (28,29).
- **µ = 0.4 (Table 2, Figs. 9-15).**
  - It is like µ = 0.3, except for Fig. 4.
  - The oscillations of the characteristics and of the traces become more uniform.
  - 10 stability intervals: (0,1), (2,3), (6,7), (9,10), (11,12), (14,15), (18,19), (21,22), (23,24) and (26,27).
  - The journal has no µ = 0.4 orbit figure.
- **µ = 0.5 (Table 3, Figs. 16-23).**
  - In Fig. 17, the characteristic for x1(0) < 0 crosses the upper dashed line three times. Fig. 19 approaches
    e~(0) = -2 several times for a~(0) < -6.
  - 10 stability intervals: (0,1), (2,3), (6,7), (9,10), (13,14), (16,17), (20,21), (23,24), (27,28) and (30,31).
  - Here vertical stability always holds where planar stability holds.
- **Comparison with Hénon (p.167-169).**
  - Hénon's x0 = x1(0) - 0.5, C_H = C - 0.25, and stability index a = Tr/2.
  - Correspondence (p.167): h1, h2, h3, h4, h5, h_2 6, h_2 7, h_2 8, h_2 9, h_3 10, h_3 11 = k 1, 2, 3, 6, 7, 9, 10,
    13, 14, 16, 17. "The computed values are different in the fourth significant digit." That refers to Hénon 1965b.
  - Hénon 1973 vertical orbits h1v, h2v, h_2 3v, h_2 4v = k 4, 5, 11, 12.
- **Evolution with µ (pp.173-177).**
  - Family h has no self-bifurcations for µ in [0, 0.5]. It is one two-parameter family over the strip T > 0,
    µ ∈ [0, 0.5], and the critical orbits form subfamilies. Some of them with Tr = ±2 are in Hénon & Guyot 1970.
  - Orbits 1 and 2 meet families "a and e" for all µ > 0. The 64/2005 translation flags this as a possible slip for
    "a and c"; the journal repeats it.
  - At µ = 0 the family is built from the pieces Ir, E_N, A_0 and A_1. At µ = 0.5 it can no longer be split into
    pieces that behave differently.
  - **New in the journal:** the further an orbit is from P1 (or P2), the more family h "walks off" from the generating
    family. For regular perturbations this agrees with Tables 1-2 of the Appendix of Bruno 1994.
  - For µ > 0.3, complete linear stability coincides with planar stability.

## 2. Checks (`checks_bigmu_ci.py` -> `checks_bigmu_ci.out`; also `explore_hard_rows.out`, `collision_rows.out`)

### 2.1 Frame and C (out sec. 0, A1)
- **The momentum convention is settled by the data.** With dx2/dt = y2 - x1, the median |dC| is 1.7e-5, 1.2e-5 and
  1.7e-5 at µ = 0.3, 0.4 and 0.5. With dx2/dt = y2 - x1 ± µ it is 1.9 to 3.5.
- **C from (x1(0), y2(0)): all 87 rows pass.** The tolerance is the propagated rounding of one unit in the last
  digit, plus one unit of C. The largest misfit is 4.8e-5 (µ = 0.3 row 28).
- **C from (x1(T/2), y2(T/2)):** this passes in every row whose T/2 point is not close to a primary.
  - It fails in µ = 0.3 rows 6, 7 and 19, µ = 0.4 rows 6, 17, 23 and 24, and µ = 0.5 rows 6, 13 and 14.
  - In all of these the T/2 point lies within 0.0034-0.025 of P2, or within 0.011 of P1 (µ = 0.4 rows 23-24,
    µ = 0.5 rows 13-14).
  - µ = 0.3 row 18 also fails narrowly (4.5e-4 against a tolerance of 2.8e-4). Its T/2 point is not near a primary.
    Integration (`explore_hard_rows.out`) reproduces the printed row: y2(T/2) = -1.74720 to -1.74723 (printed
    -1.74721) and C = -0.49510 to -0.49513 (printed -0.495131). So my linear tolerance was too tight here; the row
    is not in error.
  - Sec. 2.3 shows that the near-primary misses come from the printed y2(T/2), not from C.
- **a~, e~ (eq. 4'.1 of 64/2005) and w1 = µ/(1 - x1(T/2)):** these reproduce from the printed states in every row,
  except as follows:
  - µ = 0.4 row 16, e~(0) "––1.493": the value -1.4927 is right; the print has a doubled minus.
  - µ = 0.3 rows 6 and 7, e~(T/2): the magnitudes agree, but the sign is printed +. The printed y2(T/2) < 0 gives
    e~ < 0, and so does the integration. These are probable sign slips.
  - µ = 0.5 row 26, a~(0) = -27.648 against -27.654 from the state. This is the |e~| -> 2 ill-conditioning
    (e~ = -1.833), not an error.
  - Cells with |e~(T/2)| within 0.01 of 2 cannot be checked from 5-decimal states: µ = 0.3 rows 13-14 and 25-26,
    µ = 0.4 rows 11 and 23-24, and µ = 0.5 rows 13-14 and 27-28. a~ = x1/|2 - |e~|| amplifies the rounding. For
    example, µ = 0.5 row 27 prints a~(T/2) = 7.4695, and its 5-decimal state gives 0.66.
- **The blank a~(0), e~(0) cells are explained.** In all 10 such rows (µ = 0.4 rows 26-27; µ = 0.5 rows 16-18,
  23-25 and 30-31), the printed state gives |e~(0)| > 2: 2.02 to 3.28. So a~ is far outside the figures. The
  authors left the cells empty.
- **The mini-tables (A3):** all 16 precise values round to the table cells. The e~(0) values also agree with e~ from
  the printed (x1, y2) to the state's rounding. For example, -1.3347312 against -1.334759.

### 2.2 T~ normalisation
- The 64/2005 translation flags a conflict. 67/2005 defines T~ = (1-µ)T/(2π), but 64/2005 writes T/(2π).
- **The integrations settle it.**
  - At µ = 0.3 and 0.4 (4 rows, out sec. A5), the residual (y, dx1/dt) at t = πT~, started from the printed state,
    is 1e-5 to 7e-5 with T~ = T/(2π). With T~ = (1-µ)T/(2π) it is 0.6 to 2.5.
  - At µ = 0.5 the test is weaker. (1-µ)T/(2π) maps πT~ onto the full period, where any symmetric orbit returns to
    its start. The decisive evidence there is that the printed T/2 point is reached at t = πT~ in all 6 rows. One
    example is row 31: x1 = -4.915459 at t = πT~.
  - **These tables use T~ = T/(2π).**

### 2.3 Integration (out A5, C; `explore_hard_rows.out`)
- **Newton on (y2(0), T/2), with x1(0) fixed, from the printed state.**
  - µ = 0.5 rows 1, 3, 9, 16, 30 and 31; µ = 0.3 rows 1, 11 and 23; µ = 0.4 row 9.
  - All converge with |Δy2(0)| ≤ 5e-6. They reproduce T~, x1(T/2), y2(T/2) and C to the last printed digit, within
    1-2 units.
  - **Row 31 (µ = 0.5), started from the printed x1(0) = +0.69253, closes** at x1(T/2) = -4.915459 and
    y2(T/2) = 0.556737 (printed -4.91546 and 0.55674).
  - So row 31 is a genuine orbit, listed from its other crossing. Its pair, row 30, starts at -4.91185.
- **Rows that pass near a primary at T/2.** I located the y = 0 crossing nearest πT~ and corrected two ways: with
  x1(0) fixed, and with y2(0) fixed. The spread between the two shows how sensitive the row is to the 5-decimal start.
  - µ = 0.3 row 6: x1(T/2) 0.995298, T~ 1.300890 and C 2.535397 reproduce. **y2(T/2) = -10.282 to -10.327; printed
    -10.52856.**
  - µ = 0.3 row 7: **-12.295 to -12.311; printed -11.97044.**
  - µ = 0.3 row 8: -84.93; printed -81.39340.
  - µ = 0.4 row 17: -4.68007; printed -4.67847. The spread is 1e-5.
  - µ = 0.4 rows 18 and 19: -73.83 and -95.41; printed -73.66824 and -95.08685.
  - In these rows the printed y2(T/2) misses by far more than the spread: 5 times (row 6), 20 times (row 7) and
    over 100 times (µ = 0.4 row 17). In rows 8, 18 and 19 the two corrections agree, but the residual stays at
    2e-4 to 6e-4, so those three values are approximate.
  - The other near-primary rows (µ = 0.4 rows 6 and 7; µ = 0.5 rows 6, 7, 13, 14, 20 and 21) agree with the
    integration within 0.3 percent, which is within or close to their spread.
  - So the C(T/2) failures of sec. 2.1 come from inaccurate near-P2 y2(T/2) cells, not from wrong C values. I give
    these cells 2-3 significant figures.
  - µ = 0.3 row 19 is too sensitive to judge: the two corrections differ by 2e-3 in C.
- **Collision rows (`collision_rows.out`).**
  - All 12 rows with x1(T/2) = 0 or 1 print x1(0) ending in "00", for example -1.82800 and -3.05200.
  - At the printed C, the printed start passes within 1e-5 to 5e-5 of the body. Moving x1(0) by up to 1e-3 brings
    the pass down to a few times 1e-6.
  - **So these x1(0) cells are good to about 1e-3.** They are not 5-decimal values.
  - For µ = 0.5 row 8, Hénon 1973 gives x1 = -1.82707. My collision estimate is -1.827277. The print has -1.82800.

## 3. Identity with KIAM 64/2005 (section by section, against the held translation `bruno-varin-2005c-...-en-translation.tex`)

| KIAM 64/2005 | SSR 43:158 | same / different |
|---|---|---|
| §4′ Additions and corrections: coordinate systems I-V, eqs. (4′.1)-(4′.7), the Tr~ misprint fix, notes on the 67/2005 figures | not in the journal | **dropped**. The journal has a short Introduction and a column list instead |
| §5 µ = 0.3: Table 6, 29 orbits; mini-tables of Tr_v, a~(0), e~(0) | Table 1, 29 orbits; the same three mini-tables | **identical numbers** (16 values checked on the image against the translation) |
| §5 text: "Figure 5.3 shows two new minima of C"; "greatest difference in Tr~, Figs. 5.8 and 4.6" | no "two new minima" sentence; Tr~ maxima above 2 against BV 2009 Fig. 29 | sentence dropped; the rest is the same in substance |
| §5: 10 stability intervals | the same 10 | identical |
| §6 µ = 0.4: Table 7, 27 orbits; Fig. 6.1 (6 orbits); "except for Fig. 5.4"; 10 intervals | Table 2, 27 orbits; **no orbit figure**; "except for Fig. 4"; the same 10 intervals | the same apart from the dropped figure |
| §7 µ = 0.5: Table 8, 31 orbits; mini-table e~(0) for k = 27-28 | Table 3, 31 orbits; the same mini-table | identical |
| §7: Figs. 7.2-7.6 (5 characteristic plots); "Figs. 7.5, 7.6 are similar to 6.5, 6.6" | Figs. 17-20 (4 plots); "Figure 22 is similar to Fig. 12" | the journal has 4 characteristic plots per µ, in (x1, y2), (x1, C), (a~, e~) and (w1, y2). "Figure 22" is probably a slip for Fig. 20: (w1, y2) at µ = 0.5 against Fig. 12, (w1, y2) at µ = 0.4 |
| §7 Hénon: "Fig. 11 of [7] corresponds to Fig. 7.8" (the planar trace) | "Fig. 11 from Hénon (1965b) corresponds to Fig. 23" | **journal slip**: Fig. 23 is the vertical trace; the planar trace is Fig. 22 |
| §7: "differ by about one unit in the fourth digit" | "different in the fourth significant digit" | the same |
| §7: "Table 1 of [8] (Hénon 1973) corresponds to Table 8" | "Table 1 of Hénon (1965b) for family h corresponds to Table 3" | **journal mis-citation**: Hénon 1965b's tables are unnumbered ("the table below"); the numbered Table 1 with family h is Hénon 1973's |
| §8 Evolution: one-to-one on the strip; Hénon-Guyot subfamilies; "a and e"; pieces Ir, E_N, A_0, A_1; "System IV changes a lot for a~′_* < -1"; complete = planar stability for µ > 0.3 | the same, except the System IV sentence is dropped. **Added:** the walk-off paragraph (Bruno 1994 Appendix) and "Originally ... presented in the preprint (Bruno and Varin, 2005)" | the same plus the additions |
| Abstract: "connection with the generating family" | "evolution as µ increased" | reworded |

- **Conclusion: the journal form of KIAM 64/2005.** The tables cannot be compared, because the preprint's tables are
  not in its source.
- One open point: the journal tables have a w1(T/2) column. w1 = µ/(1 - x1) was introduced in KIAM 51/2007, after
  64/2005. The preprint's §4′ defines System IV, (a~′_*, e~_*). So the printed 2005 Tables 6-8 may have had System
  IV columns instead. This cannot be checked without the printed preprint.
- The four values per row that this check uses (x1, y2, T~, C) do not depend on that choice.

## 4. Hénon cross-checks at µ = 0.5 (out A4)

The conversion (stated on p.167): x_H = x1 - 0.5 for both crossings, and C_H = C - 0.25. Hénon puts M1 at x = -1/2 and
M2 at +1/2, with the same sense of rotation. Hénon's y0dot is "+" for family h; Bruno's dx2/dt = y2 - x1 > 0 at
t = 0. They are consistent.

- **Hénon 1973 Table 1 (six digits):**
  - h1, h2, h3, h1v, h2v, h4, h5, h_2 6, h_2 7, h_2 3v, h_2 4v, h_2 8 and h_2 9 (= k 1-7 and 9-14) all agree.
    |Δx0| ≤ 5e-6, |Δx1(T/2)| ≤ 6e-6 and |ΔC| ≤ 8e-6.
  - **These are two independent six-digit computations that agree.**
  - The vertical critical orbits k = 4, 5, 11 and 12 print Tr~_v = -2.00, as h1v, h2v, h_2 3v and h_2 4v need.
- **Hénon's unnamed ejection row** (x1 = +0.5 at T/2, x0 = -2.32707, C_H = 1.7944) = k 8. C agrees (Δ = -2.7e-5,
  within Hénon's 4 decimals). x1(0) differs by 9.3e-4; see sec. 2.3 (collision rows).
- **Hénon 1965b p.1002 (four digits, from the held digest):** for the 11 corresponding orbits, ΔC ≤ 6e-4 and
  Δx0 ≤ 1.9e-3 (h5). This is the "fourth significant digit" difference the paper reports. Hénon 1965b also differs from
  Hénon 1973 by about this much (h1: x0 -1.0875 against -1.088744). The difference is in the 1965 values.
- k = 16 and 17 (h_3 10 and h_3 11) are only in Hénon 1965b: C 0.195859 / 0.194914 against C_H + 0.25 = 0.196 / 0.195.
  They agree.
- **Use:** for µ = 0.5, Hénon 1973 and this Table 3 are interchangeable to 1e-5 on the 13 shared orbits. Table 3 adds
  18 orbits beyond Hénon's (k 15-31), out to x1(0) = -4.91 and T~ = 4.50.

### 4.1 Hénon & Guyot 1970 (held), family-h critical orbits against µ (`hg1970_compare.py` -> `hg1970_compare.out`)
- H&G tabulate, along µ, the critical orbits of the families h1-h5, h2, h3-h4, h26-h29 and h27-h28. Their columns
  are x0, x1 and C (barycentric, P1 at x = -µ, the usual C), and each family carries its type, a = +1 or -1.
- I took the values from the held H&G digest, whose tables were checked against the page images.
- Conversion: x1_BV = x_HG + µ for both crossings, and C_BV = C_HG + µ(1-µ). At µ = 0.5 this is Hénon 1973's
  conversion. The µ = 0.5 rows act as the positive control, and they match.
- **All 27 H&G rows at µ = 0.3, 0.4 and 0.5 have a counterpart here:**
  - |Δx1(0)| ≤ 1e-4 (mostly ≤ 3e-5), |Δx1(T/2)| ≤ 2e-5 and |ΔC| ≤ 2e-5;
  - a = -1 always meets Tr~ = -2.00, and a = +1 always meets +2.00.
  - Counterparts (k at µ = 0.3 / 0.4 / 0.5):
    - h1 = 1/1/1;
    - h5 = 7/7/7;
    - h2 = 2/2/2;
    - h3 = 4/3/3;
    - h4 = 6/6/6;
    - h_2 6 = 11/9/9;
    - h_2 9 = 14/12/14;
    - h_2 7 = 12/10/10;
    - h_2 8 = 13/11/13.
- The vertical critical orbits (Tr~_v = ±2) and the critical orbits beyond k = 14 at µ = 0.3 and 0.4 have no second
  source in the corpus.

## 5. Text and print slips (all small)

- p.167: "Table 1 of Hénon (1965b)" should be Hénon (1973) (sec. 3).
- p.167: "Fig. 11 from Hénon (1965b) corresponds to Fig. 23" should be Fig. 22 (sec. 3).
- p.165: "Figure 22 is similar to Fig. 12" probably should be Fig. 20 (sec. 3).
- p.158: the columns "v1(T/2), v2(T/2)" are named but not printed.
- p.167: "Strömgen" and "Barlett" for Strömgren and Bartlett. p.169: "Figure 5 of from Hénon (1973)".
- Table 2 row 16: "––1.493". Table 1 rows 6-7: the e~(T/2) sign (sec. 2.1).
- The reference list cites Strömgren 1935 (Publ. Copenhagen Obs. 100). The corpus holds Strömgren 1933, Bull.
  Astron. 9:87, which the held digest calls the journal form of Publ. 100.

## 6. Citation mining (11 references)

Checked with `ls cyclers_pdf/papers | grep -i` and grep of CORPUS_INDEX.md and the wanted list.

| reference | status |
|---|---|
| Bruno 1994, de Gruyter (The Restricted 3-Body Problem: Plane Periodic Orbits) | HELD (`bruno-1994-restricted-3-body-problem-...`) |
| Bruno & Varin 2005, KIAM Preprint 64 (this paper's preprint form) | HELD as TeX source and translation (`bruno-varin-2005c-...`); tables absent; wanted row 17 |
| Bruno & Varin 2008, "On families of periodic solutions of the restricted three-body problem", Sol. Syst. Res. 42(2):163-186 | **not held; not on the wanted list.** Its title matches KIAM 10/2005 (held), so it may be that preprint's journal form. The c/i paper (sec. 6 of its digest) cites it for the generating families c and i, for w1 "formula (8)", and for its Figs. 5, 3, 7 and 15. So it may also contain the 36/2006 and 51/2007 material. **New candidate**, medium priority |
| Bruno & Varin 2009, "Family h ... for small µ", Sol. Syst. Res. 43(1):2-25 | **not held; not on the wanted list as such.** It is very probably the journal form of KIAM 67/2005, with the µ = 0.1, 0.2 (and µ_J) tables that row 17 wants ("Table 2 from BV 2009" is cited as the model for Table 1 here, and its Fig. 29 is cited). **High-value new candidate: it would close the rest of row 17** |
| Bartlett 1964, Kong. Dan. Vidensk. Selsk. Mat.-Fys. Skr. 2(7) | not held; wanted **row 22** (no open copy found) |
| Hénon 1965a, Ann. Astrophys. 28:499 | HELD (`henon-1965a-...`) |
| Hénon 1965b, Ann. Astrophys. 28:992 | HELD (`henon-1965b-...`) |
| Hénon 1973, A&A 28:415 | HELD (`henon-1973-vertical-stability-...`) |
| Hénon & Guyot 1970, in Giacaglia (ed.), pp. 349-374 | HELD (`henon-guyot-1970-...`) |
| Strömgren 1935, Publ. Copenhagen Obs. 100 | content HELD in another form: `stromgren-1933-...` (Bull. Astron. 9:87-130), which its digest calls the journal form of Publ. 100. The wanted-list row was removed in batch 27 |

- **Proposals for the wanted list:**
  - Row 17: mark the 64/2005 half as received in journal form (this paper; tables complete). Keep the row open for
    67/2005, and name SSR 43(1):2-25 (2009) as its likely journal form.
  - Add SSR 42(2):163-186 (2008) as a new row. I have no DOI in hand; I did not guess one.

## 7. Files (this folder)

- `bruno-varin-2009-family-h-big-mu-tables.yaml`: Tables 1-3, the mini-tables and the Hénon correspondence.
- `build_yaml.py`: the mechanical text-layer copy with blank-cell assertions (it needs the text dump, which I deleted
  after use).
- `checks_bigmu_ci.py` -> `checks_bigmu_ci.out`: sections 0, A1-A5 and C cover this paper. The script is shared with
  the c/i digest.
- `explore_hard_rows.py` -> `explore_hard_rows.out`: the near-primary rows.
- `collision_rows.py` -> `collision_rows.out`: the collision rows.
- `hg1970_compare.py` -> `hg1970_compare.out`: the comparison with the Hénon & Guyot 1970 critical orbits.

*Filed as `cyclers_pdf/papers/bruno-varin-2009c-family-h-periodic-solutions-restricted-problem-big-mu-sol-syst-res-43-2-158-doi-10.1134-s0038094609020099.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/bruno-varin-2009-family-h-big-mu-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
