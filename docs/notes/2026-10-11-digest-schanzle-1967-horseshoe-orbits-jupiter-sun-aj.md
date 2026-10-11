# Digest: Schanzle 1967, "Horseshoe-Shaped Orbits in the Jupiter-Sun Restricted Problem" (#960 batch 43)

Allan F. Schanzle (Wolf Research and Development Corp., Bladensburg, Maryland; work done at the F. J. Seiler Research
Laboratory, USAF Academy), The Astronomical Journal 72(2):149-157 (March 1967), doi 10.1086/110211 (Crossref checked
2026-10-11: title, author, AJ 72, page 149, March 1967 all match), ADS 1967AJ.....72..149S (stamped on the page images).
Received 29 August 1966.
- Proposed filename: `cyclers_pdf/papers/schanzle-1967-horseshoe-shaped-orbits-jupiter-sun-restricted-problem-aj-72-149-doi-10.1086-110211-ads-1967AJ-72-149S.pdf`.
  - The file to file is an `ocrmypdf --force-ocr -l eng` copy, md5 fee02673b9400da9a0238b2b6359db14 (re-made 2026-10-11 with ocrmypdf --force-ocr -O 3; render-checked against the ADS original) (4,164,549 bytes).
  - The source is an ADS 200 dpi JPX scan, 9 pp., md5 2135960a6c830510c7982c42cde74593 (654,461 bytes).
  - **Size note:** `--force-ocr` rasterised the JPX pages, so the first OCR copy was 12 times larger than the source (8.0 MB against
    0.65 MB). The pages stay legible. Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.055/0.059,
    p.5 0.031/0.033, p.9 0.042/0.045 (copy slightly darker, no blank pages). The filed copy is the -O 3 re-make (4.2 MB, 600 dpi 16-colour images, render-checked). 8 MB was acceptable for filing; the source
    JPX file could be kept instead if size matters.
  - PDF page n = journal page 148 + n.
- **How it was read:** all nine pages were read on page images (the equations and the conclusions on 300 dpi crops; Tables
  I-III and the figure captions on 400 dpi renders; the scan is 200 dpi, so 400 dpi is an upsample). The OCR text layer was
  used only to find things.
- **Companion `schanzle-1967-tables.yaml`**: Tables I-III and the Figs. 9-22 caption values, 49 rows, format after
  `data/sources/henon-1974-hill-vertical-stability.yaml`. Comparison file `witness-comparison-schanzle.tsv`, 246 cells.
- Checks (all values in these files are ours): `fitmu_schanzle.py/.out`, `integ_schanzle.py/.out` (with
  `integ_schanzle_mu95388.out`), `caption_vs_table1.py/.out`, `schanzle_vs_taylor.py/.out`,
  `crosscheck_schanzle_taylor.py/.out`, `stability_check.py/.out`. The second witness is the tesseract output kept as
  `tess_s_t1.txt`, `tess_s_t2.txt`, `tess_s_t3.txt`, `tess_s_captions.txt`, `tess_s_xd1.txt`, `tess_s_xd2.txt`.

## 0. Verdict

**This is the earliest numerical continuation of Rabe's horseshoe orbit, with ICs to 9 digits for 11 orbits of one family.**
All printed values reproduce in our integration except the collision orbit at x0 = -2.0 (sec. 3).
- It starts from Rabe's orbit (Rabe 1961, parameter d0 = 1.06, x0 = -1.026903545) and traces a class of symmetric horseshoe
  orbits by continuation in x0 up to x0 = -2.0, which "appears to produce a collision" with the Sun.
- It shows that more than 70 periodic orbits start from one x0 (x0 = -1.026; ydot0 varied), so the horseshoe class is
  one of many horseshoe families ("only one of the multiple of periodic horseshoe-shaped families has been traced").
- Rabe's orbit is unstable: cosh(Omega T) = 1.034, whereas Rabe indicated stable. Rabe's d0 = 1.07 orbit is not in this class.
- It compares fourth-order Runge-Kutta (barycentric) with the Steffensen power-series method (Jovicentric, Rabe's own method).

**For the project it is a positive-control source, not a cycler source (PROPOSALS only):**
- 11 rows with 9-digit x0, ydot0, C, half period and the perpendicular-crossing state (Table I), plus Table II (full-period
  return, two methods) and Table III (cosh Omega T for four orbits). These suit a gate for a CR3BP periodic-orbit
  corrector and a monodromy code at Sun-Jupiter mu with a close Sun passage (the orbit at x0 = -1.9 passes 0.12 from the Sun).
  Use for `#945` R2 controls only as co-orbital background; the orbits are not DROs.
- No catalogue row cites these families. No catalogue change is proposed.
- For the held Bruno-Varin 2007 JAMM paper (digest `2026-10-06-digest-bruno-varin-2007-jamm-periodic-solutions-small-mass-ratio.md`):
  its sec. 7 and reference list name Schanzle 1967 and Taylor 1981 as horseshoe sources and the digest says neither is
  held. With this file and Taylor 1981 filed, that gap is closed for these two papers.
- Taylor 1981 did not cite Schanzle. In Taylor's terms Schanzle traced the 28-loop end of h(27,28) (sec. 3.4).

## 1. Content

**Frame and units** (Introduction, p.149, read on the page image).
- Equations of motion "in the barycentric system":
  xddot - 2 ydot = x - [mu (x + mu - 1)/r1^3] - [(1 - mu)(x + mu)/r2^3],
  yddot + 2 xdot = y - (mu y/r1^3) - [(1 - mu) y/r2^3],
  r1^2 = (x + mu - 1)^2 + y^2, r2^2 = (x + mu)^2 + y^2.
- So r1 is the Jupiter distance and r2 the Sun distance. Jupiter (mass mu) is at (1 - mu, 0), the Sun at (-mu, 0). Same frame
  as Taylor 1981.
- "The unit of length is the Jupiter-sun distance, the sum of the masses is the unit of mass, and the unit of time is such
  that the gravitational constant becomes unity." Jupiter's period is 2 pi. Confirmed on p.156: "one unit of time = 1.8879
  earth years" (11.86223 / 2 pi) and "commensurability exists for 1/2 T = 91.1077" (29 pi).
- **mu is not printed anywhere.** Rabe's value 0.000953875 (printed in Taylor 1981, p.289) reproduces every printed
  Jacobi constant of Table I to 2e-8 (ours). A fit of mu per row gives 0.000953844 to 0.000953917 (`fitmu_schanzle.out`).
  The Jacobi constant fit alone cannot tell 0.000953875 from 0.00095388 (both agree to about 1e-8). The crossing closure
  can: at the x0 = -1.1 crossing the x error is 3.3e-8 with 0.000953875 and 1.6e-6 with 0.00095388
  (`integ_schanzle.out`, `integ_schanzle_mu95388.out`), so 0.000953875 is used. This value is ours, inferred.
- **The brief said Table II has Jovicentric rows; it does not need a transform.** Both columns of Table II hold the same
  barycentric abscissa: the P.S. Initial x equals the R.K. value -1.026903545 (a Jovicentric x would be near -2.026), and
  the crossing values agree to about 1e-8. So the power-series results are printed already in barycentric coordinates.
  The Jovicentric system is mentioned in the text (p.149) only as the working frame of the P.S. program.
- Jacobi constant: not printed. C = x^2 + y^2 + 2(1 - mu)/r2 + 2 mu/r1 - xdot^2 - ydot^2 reproduces Table I (ours).

**Method.** Start on the x axis at x0 with xdot0 = 0, y0 = 0, so that symmetry halves the work. For fixed x0 find ydot0
such that xdot_c = 0 at the first crossing of the x axis. Fig. 1 and Fig. 2 show ydot0 against xdot_c for x0 = A and A + e;
corresponding zero crossings define the continuation. Fig. 3 plots ydot0 against xdot_c for x0 = -1.026; "there exist over 70
periodic orbits starting from the same point on the x axis" (read on the image). Dashed parts: orbits that cross the positive
x axis; dash-dot parts: orbits around L3 (Fig. 4). Discontinuities in 0.005 < xdot_c < 0.010 are explained by Figs. 5-8:
the number of crossings changes; orbits terminate in pairs (Fig. 7: orbits 4 and 9, 3 and 10, 2 and 11, 1 and 12).

**Continuation (Figs. 9-22, Table I).** From Rabe's x0 = -1.026903545 down to x0 = -2.0. "The value of ydot0 was considered
final when the resulting xdot_c was less than 5 x 10^-8 or until a change of one unit in the ninth significant digit of
ydot0 no longer reduced the absolute value of xdot_c." Step size so that the Jacobi constant varies by at most 2 units in
the ninth significant figure. The small indentations on the outer curve (short-period terms due to Jupiter) grow into
small loops (x0 = -1.075) which enlarge until "each loop appears to coincide". The Rabe orbit looks smooth in Rabe's own
work only because the power-series step is about ten times longer than the Runge-Kutta step. T/2 over pi is close to 29
(the integer of commensurability) throughout.

**Stability (p.156-157).** Three methods: Hill's equation (Brillouin 1948, De Vogelaere 1950), the Moulton (1914)
variational equations, and Message (1959)'s infinite determinant (found impractical). Table III: cosh(Omega T) =
1.034 (x0 = -1.026903545), 2.194 (-1.03), 52.45 (-1.05), 744.1 (-1.10); -1 < cosh < 1 stable, > 1 unstable. A sequence of
f(t) graphs between x0 = -1.027 and -1.028 shows the number of short-period fluctuations changing from 29 to 30 near
x0 = -1.0274, which explains the sharp peak in the half-period curve of Fig. 23 for 1.026 < |x0| < 1.030.

**Conclusions (p.157).** Brown (1911) and Rabe (1961) conjectured that the family they investigated would pass through the
equilibrium points near Jupiter; the orbits here do not. Rabe's d0 = 1.07 orbit is not in this class. Rabe's d0 = 1.06
orbit is unstable. Only one of the multiple families has been traced.

**Tables, as printed (page image).**

Table I, "Initial and crossing conditions", p.156 (columns x, xdot, ydot, 1/2 T, C; Initial xdot = 0, 1/2 T = 0):

| x0 | Initial ydot | C | Crossing x | Crossing xdot | Crossing ydot | 1/2 T | C at crossing |
|---|---|---|---|---|---|---|---|
| -1.026903545 | +0.039195940 | 3.00149010 | -0.970127615 | -3e-10 | -0.047658759 | 91.7570380 | 3.00149010 |
| -1.1 | +0.181820531 | 2.99587418 | -0.899961908 | -2e-8 | -0.193943363 | 91.7234374 | 2.99587418 |
| -1.2 | +0.369890946 | 2.97044973 | -0.803717732 | -6e-9 | -0.406923716 | 91.5411205 | 2.97044973 |
| -1.3 | +0.552405919 | 2.92380014 | -0.706844760 | +5e-7 | -0.638391873 | 91.4159707 | 2.92380014 |
| -1.4 | +0.731115256 | 2.85444752 | -0.609606594 | -9e-7 | -0.895081160 | 91.3293166 | 2.85444752 |
| -1.5 | +0.907645992 | 2.75985126 | -0.512149311 | +7e-7 | -1.188433481 | 91.2660184 | 2.75985126 |
| -1.6 | +1.083697373 | 2.63588664 | -0.414623774 | +1e-6 | -1.538680234 | 91.2181072 | 2.63588662 |
| -1.7 | +1.261324851 | 2.47577468 | -0.317293190 | -1e-6 | -1.985609228 | 91.1816006 | 2.47577468 |
| -1.8 | +1.443469807 | 2.26771630 | -0.220863716 | +1e-5 | -2.620799381 | 91.1542143 | 2.26771630 |
| -1.9 | +1.639477501 | 1.97492731 | -0.123600857 | +1e-5 | -3.785955508 | 91.1293187 | 1.97492729 |
| -2.0 | +1.998096500 | 1.00776933 | -0.000957480 | +5e-5 | -744.471140 | 91.1341365 | 1.00776931 |

Table II (p.156): R.K. against P.S. for x0 = -1.026903545, -1.1, -1.2 at Initial, Crossing and Final (full period T).
Final R.K. / P.S. x: -1.026903524 / -1.026903532; -1.100000027 / -1.099997564; -1.200035390 / -1.200067861. Final half-T
(elapsed time) 183.514070 / 183.514056, 183.446834 / 183.446186, 183.080098 / 183.080314. Full values in the YAML.

Table III (p.157): x0 -1.026903545, -1.03, -1.05, -1.10 give cosh(Omega T) 1.034, 2.194, 52.45, 744.1.

Figs. 9-22 captions (p.153-155): x0, ydot0, C, T/2 per orbit; Figs. 10, 11, 12 add three orbits not in Table I:
x0 = -1.03 (ydot0 0.04551460, C 3.0014622, T/2 91.97108), -1.05 (0.08506838, 3.0008697, 91.86565), -1.075 (0.13371004,
2.9990076, 91.78101). All 14 captions are in the YAML.

## 2. Witnesses

- Witness C (decides): the 400 dpi render of the page image, read cell by cell.
- Witness A: tesseract 5.5.3, psm 6, on separate crops. Not the OCR text layer. The scan is only 200 dpi, so the small
  figure captions are weaker for tesseract than the tables.
- `witness-comparison-schanzle.tsv`: 246 cells.
  - 163 cells agree exactly.
  - 15 cells: the digits appear inside a glued tesseract token (weaker). 3 cells: a truncated prefix matches (weaker).
  - 23 exponent cells (11 in Table I, 12 in Table II): a tight 2x crop of the xdot columns gives all 23 mantissas and
    signs correctly. Tesseract garbles most exponents, so each exponent has one witness, the image.
  - 34 zero cells are not compared.
  - 8 cells differ: 1 plus sign read as "-+-" in Table II; 7 caption cells (Fig. 9 ydot0, Fig. 13 x0, Fig. 16 ydot0,
    Fig. 17 x0, Fig. 18 x0, Fig. 19 T/2, Fig. 20 T/2) where tesseract garbles the 200 dpi caption.
  - One print slip: Table II P.S. crossing x of the -1.1 block prints "- .899961936" with the leading 0 missing.
- Independent check on the captions: the 14 caption values agree with Table I within caption shortening (7-8 digits; T/2 to
  5 decimals) on all 11 orbits in both (`caption_vs_table1.out`). Two caption T/2 values differ from Table I by 2e-5:
  Fig. 14 (91.54110 against 91.5411205) and Fig. 16 (91.32930 against 91.3293166).

## 3. Cross-checks

**3.1 Jacobi constant and half-period closure** (`integ_schanzle.out`, ours; DOP853 rtol 1e-13 atol 1e-13, mu 0.000953875).
- C from (x0, ydot0): agrees with the printed Initial C to 2e-8 in all 11 rows.
- Integrating to the printed half period, the state against the printed Crossing row (x error, ydot error):
  -1.026903545: 9.5e-9, 1.9e-8; -1.1: 3.3e-8, 5.6e-8; -1.2: 9.8e-7, 1.8e-6; -1.3: 4.0e-7, 8.1e-7; -1.4: 3.9e-7, 9.2e-7;
  -1.5: 3.1e-7, 8.8e-7; -1.6: 4.7e-7, 1.7e-6; -1.7: 1.9e-7, 9.2e-7; -1.8: 2.4e-10, 5.5e-9; -1.9: 2.9e-9, 1.0e-7.
  The residual |y| at that time is below 6e-5 and |xdot| below 1.3e-4 except at -1.9 (3.2e-4). Compare the printed digits:
  9 digits in the table, so the differences up to 2e-6 are a few units in the 6th place and reflect the 9-digit ICs
  amplified by the orbit (Table III).
- **x0 = -2.0 is a collision orbit.** The printed Crossing x is -0.000957480 = -mu, the Sun's x. The printed ydot there is
  -744.471140. Our integration ends 4.6e-4 from that x with ydot -5.7: the orbit passes within about 1e-3 of the Sun and is not
  recoverable from 9-digit ICs without regularisation. We did not reproduce this row. The row is a limit case, not a golden.
- **Table II full period** (x0, final x and ydot):
  -1.026903545: ours -1.026903532 (equals the printed P.S. value, printed R.K. -1.026903524), ydot ours 0.039195915 against
  printed R.K. 0.039195883 and P.S. 0.039195913. -1.1: ours x -1.099997693, ydot 0.181817073 (printed P.S.
  -1.099997564, 0.181816882; printed R.K. -1.100000027, 0.181820927). -1.2: ours x -1.200153838, ydot 0.370107612 (printed
  P.S. -1.200067861, 0.369982831; R.K. -1.200035390, 0.369938990). The 9-digit ICs are amplified by a monodromy eigenvalue
  of about 1.5e3 at -1.1. Our return errors (2.3e-6; 1.5e-4) are the size of the P.S. column errors (2.4e-6; 6.8e-5). We do
  not reproduce the much smaller R.K. drift for -1.1, so it is not evidence of better closure. Schanzle says himself that
  returning to the initial state was not the criterion of periodicity.

**3.2 Table III** (`schanzle_vs_taylor.out`, `stability_check.out`, ours; half-trace of the leading pair of the one-period planar
monodromy matrix, as printed ICs and after a correction of ydot0 and half period at fixed x0):

| x0 | printed cosh(Omega T) | ours as printed | ours corrected |
|---|---|---|---|
| -1.026903545 | 1.034 | 1.0372 | 1.0372 |
| -1.03 | 2.194 | 2.1938 | 2.1943 |
| -1.05 | 52.45 | 50.658 | 50.612 |
| -1.10 | 744.1 | 748.4 | 746.5 |

The -1.03 row agrees to 2e-4. The -1.026903545 row differs by 0.3 percent, and the -1.10 row by 0.6 percent as printed (0.3 percent corrected). The **-1.05 row differs by 3.5
percent and the gap does not close after correction.** We read 52.45 on the page image and the tesseract pass agrees.
It may be a 1967 numerical error; we cannot resolve it. For -1.05 only the 8-digit caption ICs are available.

**3.3 Schanzle against Taylor 1981** (`crosscheck_schanzle_taylor.out`, `schanzle_vs_taylor.out`; ours). The frames are
the same. Details are in the Taylor digest, sec. 3.4. In brief:
- Rabe's orbit (x0 = -1.026903545, T/2 pi = 29.2072) lies between Taylor's asymmetric bifurcation orbit of h(27,28)
  (x -1.02680, T/2 pi 29.1966) and Taylor's orbit 1b (x -1.02745, T/2 pi 29.2386). Interpolation in x gives ydot 0.03919 and
  C 3.001491 against Schanzle's 0.039196 and 3.00149010. Taylor says Rabe's horseshoe is in h(27,28); the numbers agree.
- Taylor's orbit 1a (x -1.06201, ydot 0.10851, C 3.00013, T/2 pi 29.2266) lies on Schanzle's continuation (interpolating the
  captions at x0 = -1.05 and -1.075 gives ydot 0.10844, T/2 pi 29.2288).
- Period test: Schanzle's T/2 over pi runs 29.207 down to 29.007 (x0 = -1.9) and 29.009 (x0 = -2.0). Taylor's rows c to f
  (27 loops) have T/2 pi 28.0 to 28.2. At equal C they are one Jupiter period apart, and their crossing x and ydot differ
  (orbit d: Taylor start x -0.76665, ydot -0.5123; Schanzle at the same C, interpolated over Table I: -0.732, -0.578).
  So Schanzle traced the 28-loop end of h(27,28) (the end that Taylor says "has not been continued as far", with periods
  dropping toward 29 Jupiter periods), and Taylor's c to f are the 27-loop end. Schanzle's loops that "appear to coincide" and
  his Sun-collision end match Taylor's conjecture (a) of loops merging near a period of 29. The table-I interpolation is
  coarse (rows 0.1 apart in x0), but the period gap of 1.0 is far larger than its error.
- Rabe's d0 = 1.07 orbit is not in Schanzle's class (p.157). This is Taylor's main result 14 years later; Taylor does not
  cite Schanzle.

## 4. Citation mining

References (p.157, read on the page image), checked with `ls cyclers_pdf/papers | grep -i` and
`CORPUS_INDEX.md` (and the wanted list):

| Reference | Held? | Wanted list |
|---|---|---|
| Brillouin, L. 1948, Quart. Appl. Math. 6:167 | NOT held (no file; not in index) | no row. Hill-equation numerical method. Low priority. |
| Brown, E. W. 1911, MNRAS 71:438 | NOT held | no row. The conjecture about horseshoe families. Low priority. |
| De Vogelaere, R. 1950, Canadian J. Math. 2:440 | NOT held | no row. Low priority. |
| Message, P. J. 1959, AJ 64:226 | NOT held (the `message` grep hits only a CCSDS file) | no row. Low priority. |
| Moulton, F. R. 1914, MNRAS 75:40 | NOT held (the wanted list has a different Moulton, 1913, Proc. London Math. Soc. 11:367, in row 28) | no row for the 1914 paper. Low priority. |
| Rabe, E. 1961, AJ 66:500 | NOT held (the `rabe` grep hits only Barrabes files; no index row) | no row. **Highest priority**: the source of the horseshoe orbit; also cited by Taylor 1981. |
| Steffensen, J. R. 1956, Kgl. Danske Videnskab. Selskab, Mat.-Fys. Medd. 30, No. 18 | NOT held | no row. Power-series method. Low priority. |

- None of the seven is held. Rabe 1961 is the one worth chasing (it is also the first reference of Taylor 1981 and of the
  Bruno-Varin horseshoe section).
- Companion in the same issue of AJ: Deprit & Henrard 1967, AJ 72:158 (cited by Taylor 1981 for the Hill-equation index
  and the corrector; the acknowledgement thanks Deprit). Not held, no wanted-list row. Proposal: one row.
- **Sequel lead:** D. B. Taylor 1983, Celest. Mech. 29:51-74, doi 10.1007/bf01358598 (see the Taylor digest, sec. 4).
- Wanted-list row 62: Schanzle 10.1086/110211 CONFIRMED (Crossref, this pass). Row 62 still lists Henrard 2002 as partial
  and Kogan, Lidov-Vashkov'yak as not found; Taylor 1981 has no DOI (A&A 1981 is not in Crossref).

*Filing proposal:* `cyclers_pdf/papers/schanzle-1967-horseshoe-shaped-orbits-jupiter-sun-restricted-problem-aj-72-149-doi-10.1086-110211-ads-1967AJ-72-149S.pdf`
(force-OCR copy). Table transcription: `data/sources/schanzle-1967-horseshoe-tables.yaml` (this file is `schanzle-1967-tables.yaml`).
Witness comparison, scripts and outputs beside the PDF as `<pdf stem>-<file name>`.

*Filed as `cyclers_pdf/papers/schanzle-1967-horseshoe-shaped-orbits-jupiter-sun-restricted-problem-aj-72-149-doi-10.1086-110211-ads-1967AJ-72-149S.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/schanzle-1967-tables.yaml`.*
