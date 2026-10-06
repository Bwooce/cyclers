# Digest: Bruno & Varin KIAM preprints 2005-2010 (family h sources; closed families; flat expansions) (#960)

Six Keldysh Institute (KIAM) preprints by A. D. Bruno and V. P. Varin, all in Russian. They were
supplied by the owner in batch 19 (items 42-45, plus three uploads that were not in the brief). The
target was the **family h data at mu = 0.1-0.5** that Bruno & Varin 2006 (CMDA 95:27, held) leave to
"Preprint (Bruno and Varin 2005b) ... and Preprint (2005c)".

## 0. Verdict

**The family h tables at mu > 0 are NOT recovered.**
- Preprints 67/2005 and 64/2005 survive only as LaTeX source zips in the Keldysh library. The library
  offers no PDF.
- In both zips, `A.TEX` inputs only `atitle`, `par1` and `lit`. **The main critical-orbit tables are not
  in the sources:**
  - Tables 1-5: 16, 35, 28 and 27 orbits at mu = 0, mu_J = 0.00095388, 0.1 and 0.2.
  - Tables 6-8: 29, 27 and 31 orbits at mu = 0.3, 0.4 and 0.5.
- **All figures are also missing.** Only the side tabulars of "more precise values" (a few orbits per
  mu) are present. They are transcribed in sec. 2.

Partial substitutes:
- **mu = 0:** Table 1 of preprint 34/2007 lists the same 16 critical orbits of the generating family,
  with a subset of the columns (sec. 3). It was read on the page image.
- **mu = 0.5:** preprint 64 says Hénon 1965b's table (Ann. Astrophys. 28, p.1002, left; wanted rank 18 after batch 19)
  "corresponds to Table 8", with values differing "by about one unit in the fourth digit". It also says
  Hénon 1973 (A&A 28:415) Table 1 corresponds to Table 8 for vertical stability (not held; added to the
  wanted list). The h1-h311 to k correspondence is in sec. 2.
- **mu = 0.1, 0.3, 0.5:** characteristics exist only as Fig. 7 of preprint 34/2007. **That figure is
  not in the PDF either** (the PDF has no figure pages).
- **Practical route:** recompute family h at the wanted mu by continuation from the mu = 0 table, and
  check it against the stability intervals (sec. 2) and the mu = 0.5 Hénon tables. **No numbers in the
  corpus give mu = 0.1-0.4 orbit coordinates beyond the side tabulars.**

**Wanted rank 20 (was 21 before batch 19) stays, marked "source held".** The printed preprints with the tables are still
wanted.

## 1. Filing

| Corpus file (`cyclers_pdf/papers/`) | What it is | md5 |
|---|---|---|
| `bruno-varin-2005a-families-periodic-solutions-restricted-three-body-problem-kiam-preprint-10-2005-russian-tex-source.zip` (+ `.utf8.txt`) | Preprint N10/2005, "О семействах периодических решений ограниченной задачи трех тел" ("On families of periodic solutions of the restricted three-body problem"), 20 pp. The methods paper, cited as [1] by 67 and 64; the Russian precursor of CMDA 95:27. The zip has `FIG1.RAR` and `FIG2.RAR` (solid RAR; `unar` opens them): two schematic figures (orbit parts near collision; a hyperbolic orbit near P2). | be5aa40ca8ddbf5f065ed02fd80e7814 |
| `bruno-varin-2005b-family-h-periodic-solutions-restricted-problem-small-mu-kiam-preprint-67-2005-russian-tex-source.zip` (+ `.utf8.txt`) | Preprint N67/2005, "Семейство h периодических решений ограниченной задачи при малых mu" ("Family h ... for small mu"), 32 pp. Chapter II, secs. 1-4. | 30746fe7a4b08fa253af02be416cce05 |
| `bruno-varin-2005c-family-h-periodic-solutions-restricted-problem-big-mu-kiam-preprint-64-2005-russian-tex-source.zip` (+ `.utf8.txt`) | Preprint N64/2005, "... при больших mu" ("... for big mu"), 31 pp. Sec. 4' (additions and corrections) and secs. 5-8. The zip also has `TITUL.DOC`. | 7bd58b74a55d20f395b0a0dea3074cea |
| `bruno-varin-2007-periodic-solutions-restricted-three-body-problem-small-mu-kiam-preprint-34-2007-russian.pdf` (+ `.txt` OCR) | Preprint N34/2007, "Периодические решения ограниченной задачи трех тел при малых mu", 23 pp. (cited as 22 pp.), CC BY 4.0. | af414a73484458c05b270870e69907a3 |
| `varin-2008-closed-families-periodic-solutions-restricted-three-body-problem-kiam-preprint-16-2008-russian.pdf` (+ `.txt` OCR) | Preprint N16/2008, "Замкнутые семейства периодических решений ограниченной задачи трех тел" ("Closed families of periodic solutions ..."), 28 pp. | e6a7d49fd2fc25e5f1d33f4538860f65 |
| `varin-2010-flat-expansions-ode-solutions-near-singularity-kiam-preprint-64-2010-russian.pdf` (+ `.txt` OCR) | Preprint N64/2010, "Плоские разложения решений ОДУ вблизи особенности" ("Flat expansions of solutions to ODEs at singularities"), 13 pp. | a6772158f8c88535b10ec88f863111d5 |

How the files were made:
- The `.utf8.txt` files are `iconv -f cp866 -t utf-8` of the TeX files, concatenated in `\input`
  order with a `% ===== file:` separator. N10's also carries the two `FIG.TEX` captions.
- The three PDFs have unusable text layers (Type 3 Cyrillic fonts). Their `.txt` sidecars are
  `ocrmypdf --force-ocr -l rus+eng` text. The PDFs are filed unchanged.

**Not filed:**
- `prep2005_48.zip` is byte-identical to `prep2005_67.zip` (zip md5 30746fe7..., all six member files
  match). The Keldysh library's own entry 2005-48 serves the same zip under the same "при малых mu"
  title. It is a library-side duplicate, not a third preprint.
- No TeX render was attempted. With the tables and figures absent, a PDF would add nothing to the UTF-8
  text.

## 2. Preprints 67 and 64: what the sources hold (READ; exact, from TeX)

**Notation table (N10, sec. 1).** The same families under three naming schemes:

| Strömgren | a | b | c | f | g | h | i | l | m |
|---|---|---|---|---|---|---|---|---|---|
| Broucke | I | J1 | G | C | H1, H2 | A1 | BD | E1 | F |
| Bruno | L2 | L3 | L1 | E+1/1 | 2T1 | IR+ | ID1 | ID-1 | IR- |

Family h starts as retrograde circular orbits of infinitely small radius about P1 (the larger mass).

**Columns of the absent Table 1 (67, sec. 1.1).** k; x1(0), y2(0); v1(T/2), v2(T/2) (the entry
velocity into P2); normalised period T~ = (1-mu)T/(2 pi); C; Tr and Tr_v; a~(0), e~(0); and then
either a~(T/2), e~(T/2) or a~'_*(T/2), e~_*(T/2). Tables 2-8 drop v1 and v2 for x1(T/2) and y2(T/2).

**Side tabulars, mu_J = 0.00095388 (67, sec. 2.1):**

| k | 5 | 6 | 19 | 20 | 33 | 34 |
|---|---|---|---|---|---|---|
| x1(T/2) | 0.999999407087 | 0.999999537683 | 0.999999930558 | 0.999999944993 | 0.999999961994 | 0.999999969812 |

Tr_v(12) = 1.99999842294; Tr_v(26) = 1.99999646745.
- Linear stability holds in three intervals only: start to orbit 3, orbits 11-13, and orbits 25-27.
  These are the Ir, E+1/2 and E+1/4 parts of the generating family.
- Sec. 2.5 points to Broucke 1968 (TR 32-1168, held) for mu = mu_M = 0.012155 (family h = Broucke A1).

**mu = 0.1 (67, sec. 3; Table 4, 28 orbits, absent):**
- Orbit 10: x1(T/2) = 0.999999154316, e~(0) = 0.312369180923, a~'_* = -0.0000201814.
- a~(0): orbit 8 -1.290844776096; orbit 9 -1.291462453008; orbit 20 -2.183948330271;
  orbit 21 -2.184061220002.
- Orbit 9: e~(0) = 0.312455778459.

**mu = 0.2 (67, sec. 4; Table 5, 27 orbits, absent):**
- Orbit 20: x1(T/2) = 0.999997451263, a~(0) = -2.06422013438, e~(0) = -0.000009046988.

| k | 7 | 8 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|---|
| a~(0) | -1.115637023969 | -1.116478994229 | -2.064249299804 | -2.068852563918 | -2.069062897608 | -2.087996577682 |
| e~(0) | 0.100598007243 | 0.100565693989 | -0.000005716504 | 0.000027544242 | 0.000027134384 | -0.000187286305 |

- Linear stability (Table 5) holds in 6 intervals: start to 1, (2,3), (8,9), (13,14), (16,17), (26,27).

**mu = 0.3 (64, sec. 5; Table 6, 29 orbits, absent):**

| k | 11 | 12 | 23 | 24 |
|---|---|---|---|---|
| Tr_v | -0.42260717 | -0.42309806 | -0.54642144 | -0.54604180 |

- a~(0): k 21 -2.15095703; k 22 -2.15137645.

| k | 6 | 7 | 8 | 9 | 21 | 22 | 28 | 29 |
|---|---|---|---|---|---|---|---|---|
| e~(0) | 0.0030257707 | 0.0029132519 | 0.0014496498 | 0.0011385457 | -0.1060353340 | -0.1060879730 | -1.3347311917 | -1.3352313333 |

- Stability in both directions holds in 10 intervals: start to 1, (2,3), (6,7), (11,12), (13,14),
  (16,17), (21,22), (23,24), (25,26), (28,29).

**mu = 0.4 (64, sec. 6; Table 7, 27 orbits, absent):** no side values. Stability in 10 intervals: (0,1),
(2,3), (6,7), (9,10), (11,12), (14,15), (18,19), (21,22), (23,24), (26,27).

**mu = 0.5 (64, sec. 7; Table 8, 31 orbits, absent):**
- e~(0): k 27 -1.60792022; k 28 -1.60837961.
- Stability in 10 intervals: (0,1), (2,3), (6,7), (9,10), (13,14), (16,17), (20,21), (23,24), (27,28),
  (30,31). In every planar-stable interval the orbit is also vertically stable.
- Correspondence with Hénon 1965b (stability index a = Tr/2, x0 = x1(0) - 0.5, C_H = C - 0.25):

| Hénon | h1 | h2 | h3 | h4 | h5 | h_2 6 | h_2 7 | h_2 8 | h_2 9 | h_3 10 | h_3 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| k | 1 | 2 | 3 | 6 | 7 | 9 | 10 | 13 | 14 | 16 | 17 |

- Hénon 1973's vertical-critical orbits h1v, h2v, h_2 3v and h_2 4v are k = 4, 5, 11 and 12.
- Cross-reference: Hénon & Guyot 1970 (held) tabulates Tr = +-2 subfamilies on pp.369-371 (64, sec. 8).
  This is an existing comparison point, not a new check.

**Formula notes (64, sec. 4'):**
- Coordinate systems I-V. III is astronomical about P1 (4'.1)-(4'.2). IV is about P2 (4'.3)-(4'.4),
  with a~'_* = a~_* / mu. V is Hill-scaled (4'.5).
- The generating-family limits near P2: e~_* = (V+|v1|)/|v1| sgn(v1 v2), a~'_* = sgn v1 / (3 - C)
  (4'.6)-(4'.7). These are formulas (1.1) of preprint 67, with their derivation.
- Errata to N10: the second expression of (2.2) for a~ is misprinted. The correct trace compression is
  Tr~ = (1 + log2|Tr|) sgn(Tr) for |Tr| > 2.
- 64, sec. 4'.3: Figs. 2.9, 3.8 and 4.7 (six critical orbits each at mu_J, 0.1, 0.2) were left out of
  preprint 67 for lack of space and put in 64. They are also absent from the zip.

**Conclusion of 64, sec. 8:** family h has **no self-bifurcations** for mu in [0, 0.5]. It projects one
to one onto the strip T > 0, mu in [0, 0.5]. Critical orbits k = 1, 2 meet families a and e for all
mu > 0. For mu > 0.3, full linear stability equals planar stability.

## 3. Preprint 34/2007 (READ; OCR text, Table 1 on the page image)

- **Content:** a review of the mu -> 0 theory.
  - Power geometry gives four limit problems: two-body, Hill, Hénon's intermediate problem, and the
    main limit problem.
  - Generating families of the first and second kind.
  - Sec. 5 covers family h.
- **Mu values cited in sec. 1.1:**
  - 3.5e-9: Saturn-Janus-Epimetheus (read on the page image).
  - Saturn-Mimas.
  - 5.178e-5: Sun-Neptune-KBO.
  - 9.538e-4: Sun-Jupiter-asteroid.
  - 1.215e-2: Earth-Moon.
  - Larger mu: 0.4 and 0.5.
- **Sec. 5.2:** family h at mu_J is "practically indistinguishable" from the generating family in
  (a, e). It was computed at mu = 0.1 and 0.2 [45 = preprint 67] and at 0.3, 0.4 and 0.5 [46 = preprint
  64]. Fig. 7 shows mu = 0.1, 0.3 and 0.5. "No new features, no self-bifurcations."
- **Missing parts:**
  - There is no reference list. The text ends "Continued in the preprint 'Сложные семейства
    периодических решений ограниченной задачи' ('Complex families of periodic solutions of the
    restricted problem'); references there." That follow-up preprint is added to the wanted list.
  - Figs. 1-7 are cited but absent: every page is text.
- **Identity, NOT verified:** this may be the first half of Bruno & Varin, PMM 71(6):1034-1066 (2007),
  "... при малом отношении масс". Its translation is J. Appl. Math. Mech. 71:933-960, doi
  10.1016/j.jappmathmech.2007.12.012 (Crossref CONFIRMED). The JAMM paper is not held, so the two were
  not compared. No identity note was added.

**Table 1: generating family h (mu = 0), 16 critical orbits** (p.22, read on the page image):

| k | T/(2 pi) | C | Tr | a~(0) | e~(0) | a~(T/2) | e~(T/2) |
|---|---|---|---|---|---|---|---|
| 1 | 0.50 | -1 | [-2, +inf] | -1 | -1 | 1 | -1 |
| 2 | 1.50 | 2.679465 | +inf | -1.47175 | 0.43384 | | -+inf |
| 3 | 1.99 | 2.970940 | [+inf, -inf] | -1.58720 | 0.63003 | 1.603 | 1.376 |
| 4 | 2.00 | 2.970934 | [-inf, 2] | -1.58740 | 0.62996 | 1.587 | 1.370 |
| 5 | 2.00 | 0.629961 | 2 | -1.58740 | 0 | 1.587 | +-2 |
| 6 | 2.00 | -1.711013 | [2, -inf] | -1.58740 | -0.62996 | 1.587 | -1.370 |
| 7 | 2.19 | -1.785103 | [-inf, +inf] | -1.76225 | -0.53648 | 51.553 | -2.019 |
| 8 | 2.50 | -1.491531 | +inf | -1.96669 | -0.29891 | 0.511 | -3.954 |
| 9 | 3.50 | 2.414539 | +inf | -2.41232 | 0.23485 | | -+inf |
| 10 | 3.99 | 2.929162 | [+inf, -inf] | -2.51977 | 0.39686 | 2.539 | 1.606 |
| 11 | 4.00 | 2.929161 | [-inf, 2] | -2.51984 | 0.39685 | 2.519 | 1.603 |
| 12 | 4.00 | 0.396850 | 2 | -2.51984 | 0 | 2.519 | +-2 |
| 13 | 4.00 | -2.135461 | [2, -inf] | -2.51984 | -0.39685 | 2.519 | -1.603 |
| 14 | 4.05 | -2.141393 | [-inf, +inf] | -2.56117 | -0.38821 | 5.877 | -1.829 |
| 15 | 4.50 | -1.645629 | +inf | -2.82190 | -0.19649 | 0.358 | -4.790 |
| 16 | 5.50 | 2.312067 | +inf | -3.20444 | 0.17058 | | -+inf |

- **Checks against known values (independent of the table):**
  - At the T~ = 2 and 4 resonances, |a~| = T~^(2/3): 2^(2/3) = 1.58740 and 4^(2/3) = 2.51984, as
    printed.
  - The orbit-5 and orbit-12 C values equal 2^(-2/3) = 0.629961 and 4^(-2/3) = 0.396850. This is the
    mu = 0 Jacobi constant of a circular orbit at that radius.
- **Same 16 orbits as the absent Table 1 of preprint 67, but not the same columns:** 67 also has x1(0),
  y2(0), v1/v2(T/2) and Tr_v, and uses a~'_*, e~_* for some T/2 entries.
- Orbit 7's a~(T/2) = 51.553 is as printed. It may be an a~'_* value under the 67 column convention.
  This is not resolved.

## 4. Preprint 16/2008, Varin, "Closed families" (READ; OCR text, tables on the page image)

- **Content:**
  - Natural family i (Strömgren i, Bruno ID1; prograde circular orbits about P1) has a cyclic
    structure at small mu.
  - That structure breaks up through an infinite cascade of self-bifurcations as mu -> 0 (shown in
    ref. [5] = PMM 2007).
  - Closed families i_k exist only for mu in [mu'_k, mu''_k], and shrink to a single orbit at mu''_k.
  - The paper computes the first four (k = 1-4).
  - The method: saddle and elliptic critical points of the characteristic surface M, from monodromy
    conditions m24 = m31 = 0 (1.3) plus symmetry (1.4). At such points Tr = 2 automatically (1.5).
- **Table 1 (p.8, page image): bifurcation orbits, 8 decimals.**

| | a~ | e~ | T~ | C | Tr_v |
|---|---|---|---|---|---|
| mu'_1 | 0.68871719 | 0.72837343 | 1.23370552 | 3.10346652 | -1.84389397 |
| mu''_1 | 0.81395640 | 1.43867005 | 1.97012501 | 2.91324794 | 2.00251016 |
| mu'_2 | 0.80267498 | 0.83798725 | 2.26009772 | 3.03169188 | -1.98293642 |
| mu''_2 | 0.88734729 | 1.24472130 | 2.97373803 | 2.97306701 | 2.00068411 |
| mu'_3 | 0.85530336 | 0.88504109 | 3.26879728 | 3.01520595 | -1.99690250 |
| mu''_3 | 0.91324869 | 1.17586906 | 3.97406324 | 2.98691862 | 1.99975053 |
| mu'_4 | 0.88572210 | 0.91096036 | 4.27305773 | 3.00890595 | -1.99960338 |
| mu''_4 | 0.92909055 | 1.13776964 | 4.97446586 | 2.99230487 | 1.99919355 |

- **Table 2 (p.8): the bifurcation values themselves.**

| k | mu'_k | mu''_k | mu'_k/mu'_1 | mu''_k/mu''_1 | k^(-8/3) |
|---|---|---|---|---|---|
| 1 | 4.13129887e-3 | 3.66863029e-2 | 1 | 1 | 1 |
| 2 | 6.61705554e-4 | 5.27272358e-3 | 0.160 | 0.144 | 0.157 |
| 3 | 2.15292269e-4 | 1.88241384e-3 | 0.052 | 0.051 | 0.053 |
| 4 | 9.54304953e-5 | 8.86552296e-4 | 0.023 | 0.024 | 0.025 |

- **Checks:** 6.61705554e-4 / 4.13129887e-3 = 0.160, and 2^(-8/3) = 0.157, both as printed.
- The text says only the first two elliptic points (mu''_1, mu''_2) are spatially unstable. The table
  agrees: Tr_v = 2.00251016 and 2.00068411 there, just above 2, against 1.99975053 and 1.99919355 at
  mu''_3 and mu''_4.
- **Tables 3.1-3.7:** critical orbits of the first cycle at mu = mu_J, 2e-3, 3e-3, mu'_1, 5e-3, mu''_1
  and 2.3e-2. Figures (pp.15-28) show the evolution of i_1 and i_2 and the traces.
  - **These are on the page images and were not transcribed.**
  - Table 3.1 (mu_J, 12 orbits) is on p.8. Its k = 7, 8 rows have T~ near 1.03 and C 3.17 / 2.88.
- **Relevance:** this is the "family i" control for any `#944` X2 or `#956` R9 symmetric-family work
  near the 1:1 / 2:1 commensurabilities at small mu. It is not a cycler source.

## 5. Preprint 64/2010, Varin, "Flat expansions" (READ; OCR text)

- **Content:** a plane polynomial ODE p(x,y) dy/dx = q(x,y) at a degenerate non-monodromic singular
  point.
  - Solutions that enter the singular point are, to first order, a power series plus an exponential
    addition. They are fixed uniquely as truncated series of any length in flat functions.
  - The series may converge or diverge, but approximate well either way.
  - Four examples; Fig. 1 shows solutions of equation (17).
- **References:** Bailey et al. 2006 (Experimental Mathematics in Action); Varin 2004 (Mat. Sbornik
  195(7), Poincaré maps of polynomial systems); Bruno 1979 (Local method of nonlinear analysis, Nauka);
  Hardy 1951 (Divergent Series).
- **No three-body application.** UDK 521.1 is the only celestial-mechanics tag. It is filed as general
  background for power-series methods near singularities (collision-regularisation context for X2 and
  R4), as the owner asked. No citation is in scope.

## 6. Citation mining (all six)

Held:
- Bruno & Varin 2006 CMDA.
- Bruno 1981 Celest. Mech. 24:255 (= preprint 91/1978 content).
- Broucke 1968 TR 32-1168.
- Hénon 1997 LNP m52 and 2001 LNP m65.
- Hénon & Guyot 1970.
- Hénon 1968 Bull. Astron. 3:377.
- Hitzl & Hénon 1977a,b.
- Perko 1981.
- Breakwell & Perko 1974.
- Szebehely 1967.
- Bruno 1978a,b (Brjuno), as the Celest. Mech. translations.

Already wanted:
- Hénon 1965b (rank 18).
- Bruno 1994 de Gruyter (rank 24; the 1990 Nauka original is the same book).
- Strömgren 1935 (rank 30).
- Bartlett 1964 (rank 31).

Added to the wanted list:
- Hénon 1965a, "Exploration numérique du problème restreint I. Masses égales, orbites périodiques",
  Ann. Astrophys. 28:499-511: the family h characteristic at mu = 0.5 (64's ref. [6]).
- Hénon 1973, "Vertical stability of periodic orbits in the restricted problem", A&A 28:415-426
  (ADS). Its Table 1 corresponds to the absent Table 8. Crossref also lists a short Celest. Mech.
  8:269-272 note (doi 10.1007/BF01231427) with the same title.
- Bruno & Varin 2007, PMM 71(6):1034-1066 / JAMM 71:933-960 (doi 10.1016/j.jappmathmech.2007.12.012),
  plus the follow-up KIAM preprint "Сложные семейства периодических решений ограниченной задачи" (2007,
  number not yet found). These hold the 2007-34 reference list and the cascade proof cited by 2008-16.
- Bruno's KIAM preprints 18/1972, 66/1993 (Sun-Jupiter single periodic solutions) and 93/1996
  (zero-multiple and retrograde periodic solutions; family h at mu_J as IR+J). Also the Bruno-Varin
  KIAM preprints 36/2006 (generating family i) and 51/2007 (generating family c), and Bruno-Varin 2008
  Astron. Vestnik 42 (families c and i at mu = 5e-5, read on the page image). One Tier D row.
- Colombo & Franklin 1968, AJ 73:111-123 (family i closed families, cited by 2008-16 [8]); Kotoulas &
  Voyatzis 2004, CMDA 88:343-363, doi 10.1023/B:CELE.0000023391.85690.31; Bray & Goudas 1967, Adv.
  Astron. Astrophys. 5:71-130. One Tier D row.

Not listed: Abalakin et al. 1971 (handbook), Babenko 2002 (numerical analysis), Szebehely 1982 (Russian
translation, held in English), Varin 2000 RCD (Beletsky equation), and the Varin 2010 references.
