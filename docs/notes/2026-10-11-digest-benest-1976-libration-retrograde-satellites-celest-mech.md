# Digest: Benest 1976, "Libration Effects for Retrograde Satellites in the Restricted Three-Body Problem. I: Circular Plane Hill's Case" (#960, wanted-list row 59)

D. Benest (Observatoire de Nice), Celestial Mechanics 13:203-215 (1976), doi 10.1007/BF01232724, ADS 1976CeMec..13..203B
(the bibcode is stamped on every page image). Received 29 December 1974. Crossref gives the title without the part
subtitle ("I: Circular Plane Hill's Case"); the subtitle is on the title page.
- Proposed filename: `cyclers_pdf/papers/benest-1976-libration-effects-retrograde-satellites-restricted-three-body-problem-I-circular-plane-hill-case-celest-mech-13-203-doi-10.1007-BF01232724-ads-scan.pdf`.
  - The file to file is the `ocrmypdf --force-ocr -O 3 -l eng` copy `1976CeMec13203B-ocr.pdf`, md5
    2d6beee2a0e70e4f2ab7bfa6a27a8f64 (665,114 bytes).
  - The source is the image-only ADS JBIG2 scan `1976CeMec13203B.pdf`, 13 pp., md5 db5e260a97857c480b77fa6c0adcd436
    (283,219 bytes). Fetched from https://articles.adsabs.harvard.edu/pdf/1976CeMec..13..203B (second try; the first
    returned HTTP 504).
  - Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.0643/0.0672, p.7 0.0340/0.0358, p.13 0.0236/0.0247.
    The copy is a little darker; no blank pages.
  - PDF page n = journal page 202 + n.
- **How it was read:** pp.203, 206, 208-213 on 110 dpi page images (the equations, Tables I-III, eqs. 18-22); the rest
  (figure pages, eqs. 7-10, references) on the OCR text layer, which was checked against the images where used.
- **Companion `benest-1976-tables.yaml`**: Tables I, II and III (103 numeric cells). Check script and output:
  `checks/check_benest.py`, `checks/check_benest.out` (all values in those files are ours).

## 0. Verdict

**This is the analytic theory of the slow guiding-centre libration of large, stable, NON-periodic retrograde orbits around
the small primary, in Hill's problem, with a numerical check.** All its tables reproduce (sec. 3 below).
- The fast motion is the 2:1 retrograde "reference ellipse" of period about 2 pi (the epicycle of Hill's linear problem).
  Its centre (X, Y) librates slowly round B2 on a long vertical oval. Averaging the perturbation over the fast period
  gives two integrals: A = const (the ellipse size) and B(Z, W) = (3/4) W^2 + 4 F(Z) (eq. 20), with Z = Y/A,
  W = X A^(1/2), F(Z) = 2 K[4(Z^2 - 3/4)]/pi.
- The minimum B0 = 4F(0) ~ 2.745760 is the simple-periodic orbit of family f. The small-amplitude libration period is
  T ~ 4.049378 A^(3/2) (eq. 22). The paper suggests (p.215) that divergence of B is a stability criterion (case d escapes).

**For the project (PROPOSALS only):**
- It closes the Benest 1976 part of wanted-list row 59 ("Distant-retrograde, quasi-satellite and co-orbital background
  cited by Bruno-Varin JAMM 2007", line 117 of the wanted list; the Hénon 1974 digest calls it row 63 and the Taylor 1981
  digest row 62). It is the first item that row asks for ("Benest 1976 first"). The rest of row 59 stays open (sec. 4).
- `#945` R2 DRO controls: Table III gives a sourced, reproduced control in Hill units. At x0 = 15 the family-f orbit has
  V0 = -30.0022 (ours -30.002222), and three quasi-periodic neighbours have printed B ranges and libration periods.
  These test a DRO code's Hill-frame conventions and its long-time behaviour near a DRO (libration about the orbit,
  not only the Floquet index). Note the scale: x0 = 15 Hill units is a very distant DRO: in these units L1 is at 3^(-1/3) = 0.693, so 15 units is about
  22 Hill radii (our arithmetic).
- It is not a periodic-orbit catalogue. It complements the held Hénon 1969 digest (family f: all stable, Table 3) and the
  held Hénon 1970 digest (non-periodic orbits round f). Benest's Fig. 1 is his own integration (X0 = 15, V0 = -30.1) that
  illustrates the libration Hénon 1970 described; Benest gives the theory of it.
- No catalogue row cites this paper. No catalogue change is proposed.

## 1. Content

**Frame and units** (sec. 2, p.203, read on the image). Circular plane Hill's case, rotating axes, origin at B2:
x' = u, y' = v, u' = 2v + 3x - x/rho^3, v' = -2u - y/rho^3 (eq. 1). These are Hénon's Hill equations (Hénon 1969/1970).
Our check: the family-f orbit at x0 = 15 has Hill Jacobi constant Gamma = 3x^2 + 2/r - v^2 = -225.00000 (ours).

**Osculating elements** (eqs. 3-6, p.206). Far from B2 the solution is x = X + (1/2) A cos(t - phi), y = Y - A sin(t - phi)
(the 2:1 ellipse), with X = 2k3, Y = -3k3 t + k4 drifting. The elements (X, Y, A, phi) are made osculating (eqs. 7-10;
p.207 notes this is variation of constants, after Brouwer and Clemence 1961).

**Averaging** (eq. 11, p.208): X' = dI/dY, Y' = -(3/2) X - dI/dX, A' = -(2/A) dI/dphi, phi' = (2/A) dI/dA, with
I = (1/pi) integral over 0..2 pi of dt/rho. Neglecting X in rho (justified: X is of order A^(-1/2) for large A) gives
I = 4 F(Z)/A (eq. 18), so A' = 0 and Y' = -(3/2) X.

**F(Z)** (eqs. 14-18, pp.208-210). Three cases Z^2 = 3/4, < 3/4, > 3/4 by Gröbner-Hofreiter formulas; all reduce to
F(Z) = 2 K[4(Z^2 - 3/4)]/pi. F is even in Z, F(0) = 0.686440, and F -> infinity as |Z| -> 1 (the orbit then passes
close to B2 and the approximation fails). Table I: F at 24 values of Z from 0 to 1.

**Integral B** (sec. 3, p.211). B = (3/4) A X^2 + A I is an integral of the averaged system; in (Z, W) it is eq. 20. It is
Jacobi's integral in the new variables, averaged over 2 pi. Fig. 3 plots the level curves for the 15 values of Table II.

**Libration period** (sec. 4, eqs. 21-22, pp.212-213). T = (2/sqrt 3) A^(3/2) integral 0..Z_M of [F(Z_M) - F(Z)]^(-1/2) dZ.
For Z_M -> 0: T = 2^(-1/2) pi^(3/2) [K(3/4) - E(3/4)]^(-1/2) A^(3/2) ~ 4.049378 A^(3/2).

**Numerical verification** (sec. 5, Table III, Fig. 4). x0 = 15, four V0 values: a, b, c stable, d escape.
The (W, Z) track oscillates with period 2 pi about the averaged curve. Near the periodic orbit (a) T is close to eq. 22.

Table III as printed (p.213, read on the image):

| Case | V0 | A | B_min | B_bar | B_max | T |
|---|---|---|---|---|---|---|
| periodic orbit (theory) | -30.0022 | 30.0088 | | 2.745760 | | 665.673 |
| a | -30.005 | 30.02 | 2.745789 | 2.746900 | 2.748013 | 665.829 |
| b | -30.01 | 30.04 | 2.748553 | 2.751666 | 2.754779 | 664.487 |
| c | -30.1 | 30.4 | 3.563967 | 3.620006 | 3.676044 | 452.2 |
| d | -30.2 | begin at 30.8 | 6.440 | | up to 271 | |

B_bar is "average value of B_min and B_max" (p.213); ours (2.745789 + 2.748013)/2 = 2.746901, so it is the mid-range.

## 2. Witnesses

- Witness C (decides): the 110 dpi page image, read cell by cell.
- Witness A: the OCR text layer of the -O 3 copy. A scripted search found 98 of the 103 numeric cells of the YAML verbatim
  in the text layer. The other 5 are the V0 cells: their digits are there and the minus prints as a dash and a space
  ("— 30.0022"), confirmed by eye. The text layer adds a full stop after "2.748 553" (Table III b) that is not on the page.
  No cell disagrees.

## 3. Cross-checks (`checks/check_benest.out`, all ours; scipy ellipk/ellipe; DOP853 rtol = atol = 1e-12)

- **Table I:** all 23 finite rows from eq. 18. Max |diff| 9.3e-7 (Z = 0.999: printed 2.423180, ours 2.423181; Z = 0.995:
  1.920888 against 1.920889). The second form for Z^2 < 3/4 gives the same numbers to 1e-15.
- **B0 = 4F(0):** ours 2.745761, printed 2.745760. **Eq. 22 constant:** ours 4.049378, printed 4.049378.
- **The periodic orbit:** solving u = 0 at the first y = 0 return from (15, 0, 0, V0) gives V0 = -30.002222 (printed
  -30.0022), half period 3.140954, x(T/2) = -15.000000. Its osculating A = |6 x0 + 4 V0| = 30.00889 (printed 30.0088,
  so the print truncates). Eq. 22 at our A gives T = 665.676 (printed 665.673).
- **Cases a-d** (integrated from (15, 0, 0, V0); osculating elements from the inversion of eqs. 5-6, X = 4x + 2v,
  A = |(-6x - 4v, -2u)|, Y = y + A sin(t - phi); B with the exact F):

| Case | printed B_min / B_max / T | ours B_min / B_max / T |
|---|---|---|
| a | 2.745789 / 2.748013 / 665.829 | 2.745789 / 2.748013 / 665.4 |
| b | 2.748553 / 2.754779 / 664.487 | 2.748552 / 2.754773 / 664.3 |
| c | 3.563967 / 3.676044 / 452.2 | 3.563981 / 3.676042 / 452.18 |
| d | 6.440 / up to 271 / escape | B from 5.64 up to 283 by t = 400, r to 1729: escape |

  Our T comes from the zero crossings of a 2 pi running mean of Y over 1400 time units (two libration cycles), so it is
  good to a few tenths. The initial A of each case (30.02, 30.04, 30.4, 30.8) equals |6 x0 + 4 V0|, which confirms the
  start state y0 = u0 = 0. Case d: B at t = 0 is 6.441761 (ours; A0 = 30.8, X0 = -0.4, Z0 = 0), close to the printed
  6.440. So the printed 6.440 is probably the starting B (like "begin at 30.8" in the A column), not a minimum over a
  window; our window minimum (5.64) is a different quantity. The last digit differs (6.442 against 6.440). For cases a
  and b, B at t = 0 equals B_max (2.748013, 2.754773); for c it does not (3.6578 against 3.676044).
- **Hénon 1969 comparison** (held tables file, Table 3 family f and Table 6 formula (11)): for Gamma -> -infinity,
  xi0 = -(-Gamma)^(1/2) and T/2 -> pi, index a -> 1 (stable). Our orbit at x0 = 15 has Gamma = -225.00000, so
  (-Gamma)^(1/2) = 15.0000 (the sign is the frame mirror: Hénon starts f on the negative xi axis), and T/2 = 3.140954,
  just below pi = 3.141593, continuing Hénon's trend (T/2 = 3.07495 at Gamma = -10). Benest's periodic orbit is the
  far end of Hénon's family f.
- **Mismatches kept visible:** periodic-orbit A printed 30.0088 against ours 30.00889 (truncation); case d B (6.440
  against our starting 6.442).

## 4. Citation mining

References (p.215, OCR text layer checked on the image), against `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Benest 1971, A&A 13:157 (Sun-Jupiter retrograde satellites) | NOT held (no `benest` file; one unrelated index hit) | no row. The Hénon 1974 digest lists it as low priority. |
| Benest 1974, A&A 32:39 (other mass ratios, up to mu = 0.054) | NOT held | no row. Same. Of the three, the closest to `#945` R2 because it leaves the Hill limit. |
| Brouwer & Clemence 1961, Methods of Celestial Mechanics | NOT held (no file; one index mention, not a held copy) | no row. Textbook. |
| Gröbner & Hofreiter 1965, Integraltafel | NOT held | no row. Table of integrals. Not needed (eq. 18 checks). |
| Hénon 1970, A&A 9:24 | HELD (`henon-1970-...-aa-9-24-ads-1970AA-9-24H.pdf`, digest 2026-10-07) | n/a |
| Jahnke & Emde 1945, Tables of Functions | NOT held | no row. Not needed. |
| Milne-Thomson 1965, in Abramowitz & Stegun | NOT held (no `abramowitz` file) | no row. Not needed. |

- **Sequel lead (Crossref, 2026-10-11, not in the wanted list):** Benest (1979), "Libration of Retrograde Satellite
  Orbits in the Circular Plane Restricted Three-Body Problem", IAU Symp. 81 (Dynamics of the Solar System), 41-44,
  doi 10.1017/s007418090001247x (also 10.1007/978-94-010-9256-2_8). This is likely the promised extension from Hill to the
  circular problem (p.215: "This will be developed in future papers"); we did not read it. Also by Benest: IAU Colloq. 41
  (1978) 285-303, "Existence of retrograde satellites at large distance", doi 10.1017/s0252921100062394. Proposal: one
  low-priority wanted-list row for the two IAU items.
- Wanted-list row 59 after this paper: Benest 1976 in hand. Still open: Kogan 1988?, Lidov & Vashkov'yak 1993/1994,
  Llibre & Ollé 2001 (A&A 378:1087), Barrabés & Mikkola 2005 (A&A 432:1115), Henrard 2002. The two A&A papers were tried
  in this batch and are blocked by a bot check on aanda.org (see `fetch-log.md`).

*Filed as `cyclers_pdf/papers/benest-1976-libration-effects-retrograde-satellites-restricted-three-body-problem-celest-mech-13-203-doi-10.1007-BF01232724-ads-1976CeMec-13-203B.pdf`. Tables: `data/sources/benest-1976-tables.yaml`.*
