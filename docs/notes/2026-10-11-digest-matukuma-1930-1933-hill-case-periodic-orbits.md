# Digest: Matukuma 1930, 1932, 1933, "Periodic Orbits in Hill's Case" (three notes) (#960 wanted-list row 54)

T. Matukuma (Astronomical Institute, Tohoku Imperial University), Proceedings of the Imperial Academy (Tokyo):
- 1930: "On the Periodic Orbits in Hill's Case", 6(1):6-8, doi 10.2183/pjab1912.6.6. Rec. and Comm. by K. Hirayama, Jan. 12, 1930.
- 1932: "Periodic Orbits in Hill's Case. Second Paper (Retrograde Variational Orbits)", 8(5):147-150, doi
  10.2183/pjab1912.8.147. Rec. May 12, 1932 (read before the Physico-Mathematical Society of Japan, Oct. 30, 1931).
- 1933: "Periodic Orbits in Hill's Case. Third Paper (Periodic Ejectional Orbit)", 9(8):364-366, doi
  10.2183/pjab1912.9.364. Comm. Oct. 12, 1933 (read Apr. 2, 1933).

Files (J-STAGE `_pdf/-char/en`, free; 300 dpi CCITT page images with a poor J-STAGE text layer):

| Note | Source md5 (bytes, pp.) | OCR copy md5 (bytes) | Proposed corpus filename |
|---|---|---|---|
| 1930 | baa6eb2b2c131423dca8d97b4aa4109c (138,940; 3) | 70e610daed5cfc6a1a4dba13c70af67a (127,501) | `matukuma-1930-periodic-orbits-hills-case-proc-imp-acad-6-6-doi-10.2183-pjab1912.6.6.pdf` |
| 1932 | 8a40fae6b8b980c47858c2a54832e2c2 (204,684; 4) | c65f869b968317feacd4186988af0164 (183,819) | `matukuma-1932-periodic-orbits-hills-case-second-paper-retrograde-variational-orbits-proc-imp-acad-8-147-doi-10.2183-pjab1912.8.147.pdf` |
| 1933 | c457bdeadc7d9239aff5e0dc0a603167 (175,938; 3) | f780defd256b1668d2f955b0be394937 (156,551) | `matukuma-1933-periodic-orbits-hills-case-third-paper-periodic-ejectional-orbit-proc-imp-acad-9-364-doi-10.2183-pjab1912.9.364.pdf` |

- The files to file are the `ocrmypdf --force-ocr -O 3 -l eng` copies. Dark-pixel fractions at 100 dpi, source/OCR:
  1930 p.1 0.0484/0.0485, p.2 0.0215/0.0216, p.3 0.0483/0.0481; 1932 p.1 0.0480/0.0478, p.3 0.0251/0.0249,
  p.4 0.0370/0.0368; 1933 p.1 0.0489/0.0482, p.2 0.0279/0.0278, p.3 0.0596/0.0594. No blank pages. Page 1 of each
  shows the expected title and author on the image.
- PDF page n = journal page 5 + n (1930), 146 + n (1932), 363 + n (1933).
- **How it was read:** every page on page images (110 dpi), tables and equations on 400 dpi crops. The text layer was used
  only to find things.
- **Companion `matukuma-1930-1933-tables.yaml`** (format after `data/sources/taylor-1981-tables.yaml`): every printed
  number of the three notes, with the Hénon 1969 Table 1 entries and our values beside them.
- Checks (all values ours; DOP853, rtol 1e-12, atol 1e-13 or 1e-14; Hill units mu = n' = 1): `work/hill_check.py/.out`,
  `lampos_scan.py/.out`, `lyap_check.py/.out`, `g_and_ejection.py/.out`, `g_cont.py/.out`, `extra_checks.py/.out`.
  Second witness: tesseract on 400 dpi crops, `work/tess_t32.txt`, `tess_c32.txt`, `tess_t33a.txt`, `tess_t33b.txt`.

## 0. Verdict

**Historical hand-computed Hill-problem orbits; Hénon 1969 (held) supersedes all of them. They give no new orbit, but they
are a cheap, independent set of sourced controls for a Hill-problem integrator, and they show two printed slips.**
- 1930: the looped direct variational family B (2C = 2.0, 1.3, 0.8, 0.4; Hénon family g) and the libration family F around
  the Lagrange points (2C = 4.0, 3.7, 3.4; Hénon a, c). Only Jacobi constants and a figure. One analytic formula for the
  small libration ellipse.
- 1932: retrograde variational orbits (Hénon f) by a first-order expansion in lambda = n/n', near the infinite orbit
  lambda = 0. The "present researches" rows are formula values, not integrations. For lambda < 0 they lie near family f in
  size but not in period. **For lambda > 0 they describe no orbit we can find** (sec. 3).
- 1933: Table 1 tracks family B (= g) x0 down to 2C = -0.85; Table 2 integrates ejection orbits and concludes that the
  periodic ejectional orbit exists only as 2C -> -inf. **Hénon's family g (x0 -> 0, T/2 -> 2 pi as Gamma -> -inf) confirms
  this conclusion.** Table 2 reproduces to 1e-3 to 1e-2 (sec. 3).

**For the project (PROPOSALS only):**
- No catalogue row and no route depends on these notes. Hénon 1969 remains the Hill-problem baseline for `#945` R2 and
  `#953` R7.
- Proposal: if a Hill-problem test file is built from Hénon 1969, add the 1933 Table 2 ejection column 2C = -30 as a
  regularisation/near-collision control (it reproduces to 3e-3 from a start at r0 = 1e-9). Its value is that it tests
  launching from the singular point, which Hénon's periodic tables do not.
- Proposal: record the two slips (1930 denominator 39.665; 1932 eq. 8 sign) in the corpus index line, so that nobody uses
  the 1930 ellipse formula as printed.

## 1. Content

**Equations and units.**
- 1930 (p.6, image): Darwin's intrinsic variables, dx/ds = -sin phi, dy/ds = +cos phi,
  dphi/ds = 1/R = (1/V^2){cos(theta - phi)/r^2 - 2V - 3x cos phi}, with Jacobi's integral V^2 = 3x^2 + 2/r - 2C.
  Integrated by Runge-Kutta (Runge and König 1924). Helper: Z. Kikuta.
- 1932 (p.147, eq. 1): d2x/dt2 - 2n' dy/dt - 3n'^2 x = -mu x/r^3, d2y/dt2 + 2n' dx/dt = -mu y/r^3. With mu = n' = 1
  these are Hénon's Hill equations, and 2C = Hénon's Gamma (both put the libration points at 2C = 3^(4/3) = 4.3267).

**1930 (pp.6-8).** Hill computed the variational orbits up to the cusped orbit; beyond it Poincaré proved looped orbits;
Lord Kelvin found one graphically ("Orbit B2"). Matukuma integrates the looped orbits B1-B4 at 2C = 2.0, 1.3, 0.8, 0.4
("family B", between the cusped orbit and an ejection orbit not yet computed) and the oscillating satellites F1-F3 at
2C = 4.0, 3.7, 3.4 ("family F"). Around the libration points (+-1/3^(1/3), 0) the small retrograde ellipses are printed
as xi^2/1 + eta^2/(3.208)^2 = (4.3267 - 2C)/39.665 (p.8, image). The figure (p.7) shows B1-B4, the cusped orbit labelled
"2C = 2.5?778" (third digit hidden by orbit lines) and F1-F3.

**1932 (pp.147-150).** Hill's series x = a0 cos tau + a1 cos 3tau + ..., y = b0 sin tau + ..., tau = (n - n')(t - t0),
m = n'/(n - n'). Jackson (1913) took n negative and large (small retrograde orbits) down to n = -n'. The orbit of
infinite size x = A cos n'(t - t0), y = -2A sin n'(t - t0) (n = 0, m = -1) suggests a family with n near 0. With
lambda = n/n' and the Fourier ansatz (3), eqs. (4)-(5) expand mu/(n'^2 r^3) in Laplace coefficients c_j of argument
1/3 (c0 to c12 printed, p.149). First-order results: B0 = cube root of [(80/27)(mu/n'^2)(c0 + c1)/(-lambda)]
x {1 + beta1 lambda + ...} (6); A0/B0, Ai/B0, Bi/B0 to first order; x(tau = 0)/B0 = 1/2 - 0.3300 lambda,
y(tau = 90)/B0 = 1 + 0.0346 lambda (7); axis ratio "2 - 1.3892 lambda" (8, as printed);
2C = -(n'^2/4) B0^2 {1 - 2.2237 lambda} (9). beta1 is not found; "by some considerations ... nearly equal to 1.5" is
adopted in (6)'. The table (p.150) lists Jackson's four rows and seven "present" rows (lambda = -1/10 to +1/10).

**1933 (pp.364-366).** Table 1: x0 of family B from Hill's cusped orbit (2C = 2.55788, x0 = 0.27180) down to
2C = -0.85 (x0 = 0.015). x0 seems to reach 0 at some 2C0, which would be a periodic ejectional orbit (Fig. 1). Table 2:
ejection orbits from the origin along the x axis with infinite speed, at 2C = 0, -3.0, -10, -30, -10000, giving -x and y
where the normal makes the angle phi = 90, 60, ..., -90 deg with the x axis. x(-90 deg) = 0 would be the periodic orbit; it
never is, and decreases toward 0 as 2C -> -inf (0.2663, 0.186, 0.1155, 0.0633, 0.000618). Conclusion (p.366): the periodic
ejectional orbit "does not exist in the domain of finite values of 2C"; it exists, of infinite size, only in the limit
2C -> -inf. Family B therefore runs from the cusped orbit (2C = 2.55788) to the ejectional orbit (2C = -inf). Conjecture:
beyond it a "Family C" of orbits with a loop near the origin (Fig. 2). Numerical work mostly by Mr. Hitotuyanagi, funded
by the Saitô Gratitude Foundation.

## 2. Witnesses

- Witness C (decides): the page image, 400 dpi crops. Witness A: tesseract 5, psm 6, digit whitelist, on separate crops.
- 1932 table: 6 tesseract misreads, all clear on the image ("0.93219" for 0.33219, "1986" for 1.986, "-1266" for -126.6,
  "275", "61", "014"). The formulas (6)', (7), (9) reproduce every "present" row (sec. 3.2), a third witness.
- 1932 c_j list: all 13 cells agree.
- 1933 Table 1: all 11 rows agree.
- 1933 Table 2: 7 tesseract misreads settled on the image: 0.8235, 1.1391, 0.733, 5.7402, 20.374, 351.28, 0.000618.
  Our integration agrees with each image reading (sec. 3.4).

## 3. Cross-checks (all values ours)

**3.0 Positive controls (Hénon 1969, held tables).** Our symmetric-orbit corrector (unknowns x0 and T/2; conditions
y(T/2) = 0, xdot(T/2) = 0) gives: f at Gamma = 4: xi = -0.204215, T/2 = 0.268370 (Hénon -0.20421, 0.26837); f at -3.5:
-1.870534, 2.846856 (-1.87053, 2.84686); g at Gamma = 2: 0.221684, 2.103488 (0.22168, 2.10349); a at 4.3: 0.664472,
1.518400 (0.66447, 1.51840). Full-period closure 2e-13 to 4e-11.

**3.1 1930 libration ellipse.** Linearising at L (Omega_xx = 9, Omega_yy = -3) gives nu^2 = sqrt(28) - 1, nu = 2.071594,
and the axis ratio (nu^2 + 9)/(2 nu) = 3.2080: **the printed 3.208 is reproduced.** The energy relation is
4.3267 - 2C = k xi_amp^2 with k = ratio^2 nu^2 - 9 = 35.166. Corrected small Lyapunov orbits give
(3^(4/3) - Gamma)/A^2 = 35.150 (Gamma = 4.32), 35.127 (4.31), 35.103 (4.30), with B/A = 3.2076 to 3.2065
(`lyap_check.out`). **The printed denominator 39.665 is not reproduced** (13 percent high; the amplitude it gives is 6
percent small). Observation (INFERRED, not proven): ratio^2 nu^2 - 9/2 = 39.666, so the slip could be a factor 1/2 on the
9 xi^2 term.
- The F orbits: our a-family orbits at Gamma = 4.0, 3.7, 3.4 give xi0/xi1 = 0.580196/0.775219, 0.526820/0.800140,
  0.481099/0.817900, equal to Hénon's Table 1 exact values. The 1930 note prints no x for them.

**3.2 1932 retrograde orbits (family f).**
- Laplace coefficients (`extra_checks.out`): c_j = b_{3/2}^(j)(1/3) exactly as the text says (the (1 - 3/4 cos^2 tau)^(-3/2)
  expansion with the 64/27 factor). Ours c0 = 2.6020650, c1 = 1.2477938, c2 = 0.5124306, ..., c12 = 0.0000182. The
  printed values are truncated, not rounded, from c5 on (2699, 972, 0.00346 against 0.0269964, 0.0097277, 0.0034685), and
  c0 is printed 2.60205 (one unit low in the last place).
- Jackson's rows (Gamma = 1.8851, 1.4297, -0.4188): our f orbits give |xi| = 0.332167, 0.380854, 0.801069 and
  y(90) = 0.347684, 0.405703, 1.061919 (printed 0.33219, 0.38086, 0.79827; 0.34771, 0.40572, 1.05476), equal to Hénon's
  exact values. The half period checks the lambda label: T/2 = 0.523559 against pi/6 = 0.523599 (lambda = -5), 0.628297
  against pi/5 (lambda = -4), but 1.579542 against pi/2 (lambda = -1): the printed Gamma = -0.4188 belongs to
  lambda = -0.989, and its x is 0.35 percent low.
- Matukuma's "present" rows: the first-order formulas (6)', (7), (9) with the real cube root reproduce all six rows
  (2C = -5.19, -27.07, -126.58, -126.78, -27.49, -6.05; x = 2.197, 5.180, 11.246, 11.265, 5.268, 2.605; y and the ratio
  likewise). So these rows are formula values.
- lambda < 0 rows against the exact family f at the printed 2C: |xi| = 2.280 (printed 2.20), 5.206 (5.18), 11.252
  (11.25); y = 4.359 (4.11), 10.370 (10.29), 22.494 (22.48). **The period does not match**: T/2 = 2.9706, 3.1264, 3.1401
  against pi/(1 - lambda) = 2.8560, 3.1105, 3.1385, so the f orbit at 2C = -5.2 has lambda = -0.058, not -0.1. At large
  |2C| x approaches sqrt(-2C) for every large 2:1 epicycle (x = A cos, y = -2A sin, Gamma = -A^2), so x alone does not test the expansion; the period does.
- Eq. (8) sign: (7) gives y/x = 2 + 1.3892 lambda, and the table follows the + sign (lambda = -1/10: 1.86). The printed
  "2 - 1.3892 lambda" is a sign slip (read on the 400 dpi image).
- **lambda > 0 rows.** These imply T/2 = pi/(1 - lambda) > pi. Hénon's family f has T/2 < pi for all Gamma (Table 3: T/2
  rises to pi only as Gamma -> -inf). A scan at Gamma = -6.05 and -27.49 (`lampos_scan.out`: perpendicular starts on the x
  axis, |x0| up to 1.6 sqrt(-Gamma), both senses, each later y = 0 crossing with T/2 in [2.5, 4.5]) finds only family f
  (T/2 = 3.0036, 3.1267) and a near-collision orbit (T/2 = 4.21, 4.38, crossing at x = -0.005, -0.0006), and **no
  symmetric orbit with T/2 near 3.4907 or 3.1733**. So the lambda > 0 rows are formal continuations of the first-order
  formula. This is conditional on the scan (symmetric orbits, first crossings, the x0 window).

**3.3 1933 Table 1 against Hénon 1969 Table 1** (`g_cont.out`: natural-parameter continuation of g in steps of 0.01 from
Hénon's Gamma = 0.5 row; a direct solve from interpolated guesses jumped to other orbits at four values, so those runs were
discarded):

| 2C | 1933 printed | Hénon Table 1 published (ref.) | Hénon exact | ours |
|---|---|---|---|---|
| 2.55788 | 0.27180 | 0.27180 (1, Hill) | 0.271795 | 0.271795 |
| 2.0 | 0.24 | 0.23 (5, Matukuma 1957) | 0.221684 | 0.221684 |
| 1.3 | 0.20 | 0.15778 (5); 0.2 (2, Kelvin) | 0.157138 | 0.157138 |
| 0.8 | 0.15 | 0.13 (5) | 0.113134 | 0.113134 |
| 0.4 | 0.080 | 0.080 (4, Matukuma 1933) | 0.081581 | 0.081581 |
| 0.3 | 0.075 | 0.075 (4) | 0.074512 | 0.074512 |
| 0.0 | 0.060 | 0.05581 (5) | 0.055875 | 0.055875 |
| -0.2 | 0.050 | 0.050 (4) | 0.045770 | 0.045770 |
| -0.38 | 0.040 | 0.040 (4) | 0.038220 | 0.038220 |
| -0.54 | 0.030 | 0.030 (4) | 0.032613 | 0.032613 |
| -0.85 | 0.015 | 0.0150 (5) | 0.024187 | 0.024187 |

- Our values equal Hénon's exact values in all 11 rows (to 6 decimals; closure 7e-9 to 6e-8, 2.3e-6 at -0.85 where the
  orbit passes very close to the origin).
- The 1933 values at 2C = 2.0, 1.3, 0.8 (the 1930 orbits B1-B3) are 8, 27 and 33 percent high. Hénon credits these Gamma
  to the 1957 paper, not to 1933; the 1933 value 0.20 at 1.3 equals Kelvin's 0.2 (Hénon ref. 2).
- Differences of the 1933 values from 2C = 0.4 down, in units of the printed third decimal: 0.4 -1.6, 0.3 +0.5, 0.0 +4.1,
  -0.2 +4.2, -0.38 +1.8, -0.54 -2.6, -0.85 -9.2. Hénon credits 0.4, 0.3, -0.2, -0.38, -0.54 to ref. 4 (1933).

**3.4 1933 Table 2 (ejection orbits).** Start at r0 = 1e-9 on the -x axis with radial speed from the Jacobi integral, then
find the points where the velocity direction is phi + 90 deg (`g_and_ejection.out`, `extra_checks.out`). Largest relative
difference per column (12 cells each): 2C = 0: 1.1e-3; -3.0: 5.4e-3 (3-digit cells); -10: 1.6e-2 (the column is off by
0.07 to 0.6 percent throughout, worst at -90 deg: 0.1155 against 0.1173); -30: 3.2e-3 (5.7403 against 5.7402, 0.0635
against 0.0633); -10000: 8.1e-3 at -x(-90 deg) (0.000623 against 0.000618), other cells within 3.4e-4. The -90 deg point is
reached at t = 4.322, 3.568, 3.306, 3.195, 3.1416 for 2C = 0 to -10000, so the quarter period tends to pi.
- Relation to Hénon's g: g has xi0 -> 0 and T/2 -> 2 pi as Gamma -> -inf (Hénon Table 4). A quarter period of pi from the
  ejection to the top of the orbit gives a half period of 2 pi between two crossings near the origin, consistent with
  Hénon. The 1933 conclusion (periodic ejectional orbit only at 2C = -inf) is confirmed by Hénon's table, not by us; we did
  not test finite-Gamma collision orbits.

## 4. Citation mining

References in the three notes (read on the images), checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Hill, Collected Works vol. 1 (pp.284, 317, 335) | NOT held (no `hill` file) | no row; Hill 1878/1886 are classical |
| Poincaré, Méth. nouv. t. I (1892) p.108 | NOT held | no row (textbook) |
| Lord Kelvin, Collected Papers vol. IV p.516 | NOT held | no row; historical only |
| G. H. Darwin, Scientific Papers vol. IV p.1 (probably the reprint of "Periodic orbits", Acta Math. 21, 1897; INFERRED, not checked) | NOT held | no row. Darwin's Hill-type periodic orbits are a classical numerical set; low priority |
| Runge and König 1924 (numerics textbook) | NOT held | no row |
| J. Jackson 1913, M.N. 74:62 (retrograde satellite orbits) | NOT held | no row. Source of Hénon's ref. 3 values; low priority |
| Matukuma 1957 (Hénon's ref. 5) | NOT held | wanted row 54 says "Matukuma 1957 NOT FOUND" (doi lookup); still open |

- Wanted-list row 54: the three 1930-1933 notes are now in hand; Hénon 1966a and 1966b and Matukuma 1957 remain.
  Proposal: shrink the row accordingly.
- Hénon 1969 digest line "Matukuma ... not held" becomes "1930-1933 held; 1957 not held".

*Filed as `cyclers_pdf/papers/matukuma-1930-periodic-orbits-hills-case-proc-imp-acad-6-6-doi-10.2183-pjab1912.6.6.pdf` (with the other files named above). Tables: `data/sources/matukuma-1930-1933-tables.yaml`.*
