# Digest: Bruno & Varin 2009, "Closed Families of Periodic Solutions of a Restricted Three-Body Problem" (Solar System Research 43(3)) (#960 batch 34)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, RAS), "Closed Families of Periodic Solutions of
a Restricted Three-Body Problem", Solar System Research 43(3):253-276 (2009), doi:10.1134/S0038094609030071.
Russian original: Astronomicheskii Vestnik 43(3):265-288 ("Original Russian Text © A.D. Bryuno, V.P. Varin, 2009").
Received 15 May 2008. RFBR project 08-01-00082.
- **Byline:** the printed byline is "A. D. Bruno and V. P. Varin" (running head "BRYUNO, VARIN"), not Varin & Bruno.
  The preprint it comes from, KIAM 16/2008, is by Varin alone.
- File: `c19cfdb1-bruno2009_2.pdf`, 24 pp., md5 `61f12810b25d5accf68041a9c3883a8b`.
- **Proposed corpus filename** (printed author order; `bruno-varin-2009-` is not in use):
  `bruno-varin-2009-closed-families-periodic-solutions-restricted-three-body-problem-ssr-43-253-doi-10.1134-S0038094609030071.pdf`.
  If the filer prefers the order in the task, `varin-bruno-2009-...` with the same tail.
- **How I read it:**
  - I read the whole text layer. I read every table on 300 dpi page images (PDF pp. 2-6), in strips: Tables 1, 2
    and 3.1-3.7, 886 cells. The text layer agreed with the image in every cell.
  - I compared every cell with our KIAM 16/2008 translation (which was transcribed from the preprint's 300 dpi
    images), and the text section by section. I checked two disputed points on the KIAM 16 page images (PDF pp. 4-5).
  - Scripts and outputs in this folder: `ssr2009_data.py` (transcription), `compare_2009.py` -> `compare_2009.out`,
    `check_2009_t1.py` -> `check_2009_t1.out`, `check_2009.py` -> `check_2009.out`, `check_2009_C.py` ->
    `check_2009_C.out`, `summarize_2009.py` -> `summarize_2009.out`, `cr3bp.py`, `make_yaml.py`.
  - YAML: `varin-bruno-2009-closed-families-tables.yaml` (all tables, printed values).

## 0. Verdict

**This is the journal form of Varin's KIAM 16/2008, with a new closing section.**
- All 886 table cells are identical to KIAM 16/2008 (Tables 1, 2, 3.1-3.7). Secs. 1-2 are a close translation of
  the preprint text.
- New in SSR: a section "Other closed families" and a "Conclusions" section, and the figures are re-split
  (16 figures instead of 14).
- **Table 1 is excellent data.** The eight bifurcation orbits (saddle points mu'_k and elliptic points mu''_k of the
  family-i surface, k = 1-4) are given to 8 decimals. From the printed (a~, e~) and mu, our integration closes each
  orbit to |y1(T/2)| <= 5e-9 at the printed T~ (to 1e-7). It gives Tr~ = 2.0000 (as the text requires) and
  reproduces the printed Tr~_v to 1e-6 or better (sec. 3).
- **Tables 3.1-3.7 are coarser.** They give 118 critical orbits of the first cycle at 7 values of mu, to 5 decimals.
  50 of the 118 rows miss the printed C by more than the rounding allows (up to 7e-3). The traces of the
  reproducible rows mostly agree (sec. 3).
- **The printed C is -2H + mu, not -2H.** All 8 rows of Table 1 agree with -2H + mu to 7e-9, with H from SSR 2008
  eq. (2). At mu = 0 this is Bruno's C = -2H. Anyone who uses these C values at mu > 0 must subtract mu. Neither the
  paper nor KIAM 16 says so.
- **What it gives the project:** 8 exactly reproducible, symmetric, planar periodic orbits at 8 values of small mu,
  each at a fold or a saddle of family i, with C, T and both traces. These are good **positive controls** for a
  symmetric-family continuation code at small mu (`#944` X2 or `#956` R9), especially near the 1:1 / 2:1 region of
  family i. mu''_1 = 0.0366863 and mu'_1 = 0.0041313 bracket the Earth-Moon mu (Table 3.6). There are no cyclers.
- Catalogue implication (PROPOSAL only): none for catalogue rows. If a control test is built, use the Table 1 rows
  at their printed mu, with C_test = C_printed - mu.

## 1. Frame and conventions (stated first)

- The same frame as Bruno & Varin 2008 (SSR 42:154): planar CR3BP, synodic, origin at P1 (mass 1 - mu), P2 at
  (1, 0). H = (y1^2 + y2^2)/2 + x2 y1 - x1 y2 - 1/r + mu (1/r + x1 - 1/r2); y are momenta (y2 = dx2/dt + x1).
  The KIAM 16 image (PDF p.4) prints the same R = r^-1 + x1 - r2^-1.
- Symmetric periodic solutions: x2(0) = y1(0) = x2(T/2) = y1(T/2) = 0.
- Tables give the astronomical coordinates of the crossing on the right characteristic of the generating family i:
  x1(0) = a~ |2 - |e~||, y2(0) = sgn(a~ e~) sqrt(|e~| (1 - mu) / |x1|). Table 3.4 row 19 has a~ = -0.63260: a
  left-characteristic crossing. Our integration closes it (raw T~ off by 6e-4, within the 5-digit precision), so the
  minus sign is real.
- T~ = T/(2 pi). Tr~, Tr~_v: (1 + log2|Tr|) sgn Tr outside [-2, 2]. Tr = trace(M) - 2 for the 4x4 plane monodromy
  matrix M (eq. (4): sum m_kk = 4 <=> Tr = 2); Tr_v is the trace of the 2x2 vertical monodromy matrix.
- mu values: mu_J = 9.5388e-4; mu_M = 1.2155e-2 (I used 0.01215585; the difference does not matter at 5 decimals);
  mu'_k, mu''_k from Table 2.

## 2. Content

- **Sec. 1 (self-bifurcations):** the characteristics of a symmetric family, as mu varies, are level lines of a
  surface M in Pi x mu. Generic self-bifurcations are its saddle points; one-orbit "families" are its extremal
  (elliptic) points. An extremum of y2(0) needs m_{4,1} = 0, an extremum of x1(0) needs m_{1,4} = 0; at a saddle
  or elliptic point both hold (eq. (2)), together with the symmetry conditions (3). This 4 x 4 system in
  x1(0), y2(0), T, mu is nondegenerate and is solved directly (following Varin 2000, RCD 5:313). Tr = 2 there.
- **Sec. 2 (formation and evolution):** closed families i_k branch off family i at mu'_k and shrink to one orbit at
  mu''_k; mu'_k < mu''_k, both decrease to 0, and both scale like k^(-8/3) (Table 2). The elliptic point of i_{k+1}
  lies above the saddle point of i_k on the surface, so for any mu < mu'_1 two closed families coexist.
  Figs. 1-16 show the first cycle (i_1) at mu = mu_J, 2e-3, 3e-3, mu'_1, 5e-3, mu_M, 2.3e-2 and mu''_1, the second
  cycle (i_2), and the third and fourth.
- **New "Other closed families" section:** Colombo, Franklin and Munford 1968 computed two closed families (in i_2
  and i_3) at mu_J. Bruno 1993a computed Tr for them. Bruno 1993b, c found closed two-, three- and four-fold
  families branching from resonant orbits. Voyatzis & Kotoulas 2005 computed six families of perturbed Id_p
  fragments, p = -2 to -7, for Sun-Neptune; the p = -7 one is closed and far from its generating family.
- **New "Conclusions":** closed families exist near the intersection of family Id with P2, both inside (a < 1,
  Sun-Jupiter-asteroid) and outside (a > 1, Sun-Neptune-Kuiper belt).

**Text differences from KIAM 16/2008** (our translation, and the preprint images where marked):
- The SSR text says "Tables 1-8 with 8-digit accuracy". The preprint says Table 1 is given to 8 digits. The SSR
  sentence is a translation slip: Tables 3.x have 5 decimals.
- KIAM 16 ref. [7], "families c and i ... at mu = 5e-5, Astron. Vestnik 2008, vol. 42" with blank issue and pages,
  becomes Bruno & Varin 2008 (SSR 42(2):154). Its printed title is "On families of periodic solutions ...".
- KIAM 16 ref. [8] "Colombo, Franklin" becomes "Colombo, Franklin and Munford".
- KIAM 16 Fig. 1 (right characteristic only) becomes SSR Figs. 1 (left) and 2 (right). The axis-name error of
  KIAM 16 Figs. 2 and 6 (our translation's "[sic]") is fixed in SSR Figs. 3 and 8.
- SSR misprints: "SRS families" and "understude" (abstract), "Voyatis and Kotulas", "Water de Gruyte", and two
  different family-h papers both labelled "2009b".
- **Corrections to our KIAM digest** (`2026-10-06-digest-bruno-varin-kiam-preprints-2005-2010.md`, sec. 4), both
  checked on the KIAM 16 images:
  - The digest writes the monodromy conditions as "m24 = m31 = 0". The preprint (PDF p.5) and SSR have
    m_{4,1} = m_{1,4} = 0.
  - The digest lists Table 3.6 at "mu''_1". The preprint caption (OCR line 486) and SSR give Table 3.6 at mu = mu_M.
  - The digest's Tables 1-2 transcription agrees with SSR in every cell (40 + 20).

## 3. Checks

**Table 1 (8 digits; `check_2009_t1.out`).** I integrated the printed state at the printed mu, with no correction
(DOP853, rtol = atol = 1e-12), to the x2 = 0 crossing nearest pi T~:

| orbit | T~ at crossing (printed) | y1(T/2) | -2H + mu (printed C) | Tr~ | Tr~_v (printed) |
|---|---|---|---|---|---|
| mu'_1 | 1.23370550 (1.23370552) | 2e-9 | 3.10346652 (3.10346652) | 2.000001 | -1.84389372 (-1.84389397) |
| mu''_1 | 1.97012501 (1.97012501) | 1e-10 | 2.91324794 (2.91324794) | 2.000000 | 2.00251016 (2.00251016) |
| mu'_2 | 2.26009772 (2.26009772) | 4e-11 | 3.03169188 (3.03169188) | 2.000002 | -1.98293640 (-1.98293642) |
| mu''_2 | 2.97373803 (2.97373803) | 1e-10 | 2.97306701 (2.97306701) | 2.000001 | 2.00068411 (2.00068411) |
| mu'_3 | 3.26879736 (3.26879728) | 5e-9 | 3.01520595 (3.01520595) | 1.999971 | -1.99690262 (-1.99690250) |
| mu''_3 | 3.97406324 (3.97406324) | 1e-9 | 2.98691861 (2.98691862) | 2.000002 | 1.99975054 (1.99975053) |
| mu'_4 | 4.27305761 (4.27305773) | 9e-11 | 3.00890595 (3.00890595) | 2.000050 | -1.99960331 (-1.99960338) |
| mu''_4 | 4.97446585 (4.97446586) | 6e-11 | 2.99230487 (2.99230487) | 1.999996 | 1.99919355 (1.99919355) |

- All eight close. The small T~ and Tr~ misses on the saddle orbits mu'_1, mu'_3, mu'_4 are the expected
  sensitivity of an 8-digit state at a degenerate point.
- The text says only mu''_1 and mu''_2 are spatially unstable (Tr_v > 2). Our Tr~_v values agree.
- Table 3.4 row 3 (mu = mu'_1) is Table 1's mu'_1 orbit rounded to 5 decimals, as it should be.

**Table 2:** all 12 printed ratios agree with mu'_k/mu'_1, mu''_k/mu''_1 and k^(-8/3) rounded to 3 decimals.

**Tables 3.1-3.7 (5 decimals; `check_2009_C.out`, `check_2009.out`, `summarize_2009.out`):**
- **C from (a~, e~):** 68 of 118 rows agree with the printed C within the propagated rounding tolerance (half a unit
  in a~ and e~). **50 rows do not**, by 1e-5 to 7e-3. The largest misses are Table 3.5 rows 22-23 (3e-3, 7e-3),
  Table 3.7 rows 2-3 (5e-3, 2e-3) and Table 3.6 row 2 (3e-3). The misses are spread over all seven tables and are
  not single-digit slips. I think the tabulated critical orbits were located by interpolation along the family, and
  the columns were interpolated separately. That is a guess; the paper does not say how Tables 3.x were made.
- **Integration with correction:** I integrated each row, then corrected y2(0) with x1(0) fixed so that
  y1(T/2) = 0. 96 of 118 rows give a corrected T~ within 1e-3 of the printed one.
  - Of these 96, the traces agree with the printed Tr~ within 0.01 in 68 rows and within 0.05 in 80; Tr~_v within
    0.01 in 88 and within 0.05 in all 96.
  - The corrected C agrees with the printed C within 1e-4 in 83 rows.
- 22 rows cannot be reproduced from 5 digits: 3.3 #1, 8, 17, 18; 3.4 #17; 3.5 #6, 7, 11, 19-26, 30-33; 3.6 #7;
  3.7 #7. Most are highly unstable orbits (Tr~ 15-29, so |Tr| up to about 1e8), where a 5-digit state cannot follow
  the orbit for a full half-period. Table 3.5 row 23 has no x2 = 0 crossing near pi T~.
- Printed values that are not exactly +-2 at a "critical" orbit are kept as printed: Table 3.4 row 6 Tr~ = 2.00001;
  Table 3.5 row 29 Tr~_v = 1.99999.

## 4. Comparison with KIAM 16/2008 (`compare_2009.out`)

- **Tables: 886 cells compared, 0 differences.** That covers Table 1 (40 cells), Table 2 (8 mu values and 12
  ratios) and Tables 3.1-3.7 (826 cells), against our translation of the preprint.
- So every observation in sec. 3 applies to the preprint as well.
- The text differences are in sec. 2.

## 5. Citation mining (15 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i` in CORPUS_INDEX.md and the wanted list.

| reference | status |
|---|---|
| Bruno 1994, de Gruyter | HELD (`bruno-1994-...`) |
| Bruno 1993a, b, c, KIAM 66, 67, 68 (Sun-Jupiter simple, two-fold, multi-fold) | not held; **wanted row 55** |
| Bruno 2006, Kosm. Issled. 44(3):258 [Cosmic Res. 44:245] | not held; not wanted (candidate; also cited by SSR 2008) |
| Bruno & Varin 2007, PMM 71(6):1034 [JAMM 71:933] | HELD (`bruno-varin-2007-...-jamm-...`) |
| Bruno & Varin 2008, SSR 42:154 | this batch (the companion digest) |
| Bruno & Varin 2009a, "Family h ... for small mu", Astron. Vestn. 43(2) [SSR 43(1)] | **not held; not on the wanted list. Candidate**: journal form of KIAM 67/2005, likely with the tables that row 17 says are missing from the TeX source |
| Bruno & Varin 2009b, "Family of h-periodic solutions ... mu = 5e-5", SSR 43(1):26 | **not held; not wanted. Candidate** (family h at the Sun-Neptune mu) |
| Bruno & Varin 2009 ("2009b" again), "Family of h-periodic solutions ... big mu", Astron. Vestn. 43(2):167 [SSR 43(2):158] | **not held; not wanted. Candidate**: journal form of KIAM 64/2005 (row 17) |
| Colombo, Franklin & Munford 1968, AJ 73:111 | not held; **wanted row 56** (listed as Colombo & Franklin; add Munford) |
| Varin 2000, Regul. Chaotic Dyn. 5(3):313 (Beletsky equation degeneracies) | not held; not wanted (method source for the critical-point solver; low priority) |
| Varin 2008, KIAM 16 | HELD (`varin-2008-...`) |
| Voyatzis & Kotoulas 2005, Planet. Space Sci. 53:1189 | not held; **wanted row 56** |

Proposals for the wanted list:
- Add the three 2009 SSR family-h papers to row 17 as the journal forms of KIAM 67 and 64. They are the likely way
  to get the family-h tables that the TeX sources lack. DOIs not checked here.
- Row 56: add Munford as third author of Colombo & Franklin 1968.

*Wanted-list row numbers are the current (batch-30) numbering.*

*Filed as `cyclers_pdf/papers/varin-bruno-2009-closed-families-periodic-solutions-restricted-three-body-problem-sol-syst-res-43-3-253-doi-10.1134-s0038094609030071.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/varin-bruno-2009-closed-families-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
