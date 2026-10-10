# Digest: Taylor 1981, "Horseshoe Periodic Orbits in the Restricted Problem of Three Bodies for a Sun-Jupiter Mass Ratio" (#960 batch 43)

D. B. Taylor (HM Nautical Almanac Office, Royal Greenwich Observatory), Astronomy & Astrophysics 103:288-294 (1981),
ADS 1981A&A...103..288T (the bibcode is stamped on every page image). Received 19 May, accepted 14 July 1981.
No DOI: Crossref has no A&A record for 1981 (a 1981 filter on ISSN 0004-6361 returns 0 results; checked 2026-10-11).
- Proposed filename: `cyclers_pdf/papers/taylor-1981-horseshoe-periodic-orbits-restricted-problem-three-bodies-sun-jupiter-mass-ratio-aa-103-288-ads-1981AA-103-288T.pdf`.
  - The file to file is an `ocrmypdf --force-ocr -l eng` copy, md5 7cd79bf4dea69ca85b86584db27287ef (698,451 bytes).
  - The source is an image-only ADS JBIG2 scan, 7 pp., md5 44fbbfb5a5059be14b98d16a87d1687a (162,865 bytes).
  - Page images survived the OCR. Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.079/0.090,
    p.4 0.028/0.033, p.7 0.054/0.061. The copy is a little darker, no blank pages.
  - PDF page n = journal page 287 + n.
- **How it was read:** all seven pages were read on page images (150 dpi overview, 300-400 dpi crops for the tables,
  the equations, the Fig. 5 caption, the mu line and the reference list). The OCR text layer was used only to find things.
- **Companion `taylor-1981-tables.yaml`**: Tables 1-3 and the Fig. 5 caption, 18 rows, format after
  `data/sources/henon-1974-hill-vertical-stability.yaml`. Comparison file `witness-comparison-taylor.tsv`, 102 cells.
- Checks (all values in these files are ours): `integ_taylor.py/.out`, `taylor_events.py/.out`, `taylor_correct.py/.out`,
  `stability_check.py/.out`, `crosscheck_schanzle_taylor.py/.out`, `schanzle_vs_taylor.py/.out`. The second witness is the
  tesseract output kept as `tess_t1.txt`, `tess_t2.txt`, `tess_t3.txt`.

## 0. Verdict

**This is the first systematic catalogue of symmetric horseshoe periodic orbits at the Sun-Jupiter mass ratio, with printed
initial conditions to 5 digits and stability indices.** All 18 printed rows reproduce in our integration (sec. 3).
- It finds 7 families of symmetric horseshoe orbits, named h(n,m) (n, m = loops on the two sides of the smooth part). Their
  smooth-horseshoe parts lie at different distances from Jupiter. Rabe's horseshoe is in h(27,28) (p.290, read on the image).
- Main result: the smooth horseshoes do NOT form one family, as Rabe implied. There are many families, each with a smooth
  region, each with a bifurcation to an asymmetric family ha(n,m). Conjectured to be infinitely many (sec. 4 of the paper).
- Data: ICs, Jacobi constant, T/2 pi, horizontal index s and vertical index s_v for 6 + 4 symmetric orbits (Tables 1-2),
  the 7 asymmetric bifurcation orbits (Table 3) and one asymmetric orbit (Fig. 5 caption).

**For the project it is a positive-control source, not a cycler source (PROPOSALS only):**
- It gives published goldens (5-digit x0, ydot0, C, T/2 pi, s, s_v) for planar CR3BP periodic orbits at mu = 0.000953875.
  These suit a monodromy/Floquet code gate (`#945` R2 and the co-orbital controls of wanted-list row 62). The note in
  sec. 3 matters: one must refine the ICs before computing s for the long-period rows. The uncorrected 5-digit ICs give
  wrong s (e.g. -15.9 for Table 2 a, printed -10.2204).
- It is not a DRO paper. The orbits encircle L3, L4 and L5 and pass close to Jupiter; they are retrograde in the
  rotating frame (p.290) but are not distant retrograde orbits about Jupiter.
- No catalogue row cites these families. No catalogue change is proposed.
- Verdict on the family label: the horseshoe class of Schanzle 1967 is the 28-loop end of Taylor's h(27,28) (sec. 3).
  Taylor does not cite Schanzle. Taylor says the 28-loop end "has not been continued as far"; Schanzle had traced it to
  a collision orbit 14 years earlier.

## 1. Content

**Frame and units** (sec. 2, p.289, read on the page image).
- Planar circular restricted problem, rotating frame, origin at the centre of mass. S (Sun, mass 1 - mu) is at (-mu, 0).
  J (Jupiter, mass mu) is at (1 - mu, 0). "mu is taken to be 0.000953875 which corresponds to the value for the mass of
  Jupiter used by Rabe (1961)."
- G = 1, angular rate 1, period of revolution 2 pi. Equations (1): x'' - 2y' = W_x, y'' + 2x' = W_y,
  W = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2, r1 = Sun distance, r2 = Jupiter distance.
- Jacobi constant: C = x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2 - xdot^2 - ydot^2. Larger C = lower energy.
- Tables print x and ydot at a perpendicular x-axis crossing (y = 0, xdot = 0), C, T/2 pi (full period over 2 pi, so
  29.2386 means 29.2386 Jupiter periods), s (horizontal) and s_v (vertical index, "as described by Henon (1973)").
- Our reading of the index: s is the half-trace of the leading eigenvalue pair of the one-period planar monodromy matrix.
  s_v is the Henon 1973 vertical index. |s| < 1, |s_v| < 1 stable.

**Method** (sec. 2). Rabe's recipe: the closest approach to Jupiter lies near the line S to L4, with the speed V in the
rotating frame near zero. Put P at distance lambda from L4 and find lambda and the angle where V and its derivative both
vanish. Values lambda = 0.05 to 0.13 were used. The orbits were then corrected by the Deprit-Henrard (1967) differential
corrector. The step lambda 0.01 gave members of distinct families each time; 0.002 from lambda 0.092 to 0.098 stayed in
one family. Seven families found, "but this was not an exhaustive search".

**Families (sec. 3.1).** h(27,28), h(20,21), h(16,17), h(13,14), h(10,11), h(9,9), h(9,8). Fig. 2 is the schematic:
a region of smooth horseshoes between points B and E (s from stable to mildly unstable), an extremum of C on the family,
two further critical orbits with s = -1 (from which symmetric four-crossing families should bifurcate, after Henon 1965),
one asymmetric bifurcation (s = +1), and loops added on the outer curve when the family is followed away from B-E.
- Family h(27,28) was followed furthest: both ends to a maximum of s. Minimum s in the unstable stretch is -1.0102
  (read as "- 1.0102" on the image) with C changing by under 3.0e-6 ("less then 3.0 10^-6" in the print, read on the image). The orbits of Fig. 3a, b lie between F and E, the
  orbits of Fig. 3c-f between B and A. Max s about 3.1e4 on the way, 1.0154 at orbit f.
- Smooth horseshoes with smooth inside curves do not exist; the inner curve always has a small forced eccentricity.
- Vertical: two critical orbits with s_v = 1 on each family, so 3D families of horseshoes bifurcate there (Henon 1973).
  h(27,28) is vertically stable except in a small range of C.

**Asymmetric families (sec. 3.2).** The bifurcation orbits are in Table 3. ha(27,28) was continued furthest (Fig. 5, an
orbit with 28 loops). The six other families were started only. The mirror images in the line SJ also exist.

**Conjectures (sec. 4).** (a) h(27,28) ends on both sides with orbits where the loops coincide, repeated 28 and 29
times, periods near 28 and 29 Jupiter periods. (b) Infinitely many families, denser toward L3, L4, L5. (c) The symmetric
families bifurcate from the stable part of Stromgren's family a (long-period orbits around L3, Stromgren 1929).
Applications (sec. 5): Saturn's co-orbital satellites S10/S11, ring particles around Uranus (mu about 1e-10).

**Tables, as printed (page image).**

Table 1, orbits of family h(27,28) (Fig. 3), p.289:

| Orbit | x | ydot | C | T/2 pi | s | s_v |
|---|---|---|---|---|---|---|
| a | -1.06201 | 0.10851 | 3.00013 | 29.2266 | 123.25 | -0.0219 |
| b | -1.02745 | 0.04032 | 3.00149 | 29.2386 | 1.1767 | -0.0829 |
| c | -0.97668 | -0.06118 | 2.99893 | 28.2234 | 297.54 | -0.0146 |
| d | -0.76665 | -0.51230 | 2.93589 | 28.1080 | > 10^4 | 0.5083 |
| e | -0.42413 | -1.51123 | 2.61911 | 28.0348 | > 10^4 | 0.6796 |
| f | -0.11239 | -4.00533 | 1.90214 | 28.0064 | 1.0154 | 0.9610 |

Table 2, smooth horseshoes of Fig. 4, p.292: a (-1.06206, 0.09265, 3.00334, 18.1334, -10.2204, 0.4387, h(16,17)); b
(-1.06456, 0.09238, 3.00428, 15.0578, -24.206, 0.7215, h(13,14)); c (-1.07912, 0.11595, 3.00521, 12.3894, 5.9007, -1.0187,
h(10,11)); d (-1.09137, 0.13177, 3.00705, 9.3743, -38.599, -0.8998, h(9,8)).

Table 3, asymmetric bifurcation orbits, p.293 (x, ydot, C, T/2 pi, s_v): h(27,28) -1.02680, 0.03897, 3.001491, 29.1966,
0.1801; h(20,21) -1.03136, 0.04083, 3.002106, 22.2728, -0.3376; h(16,17) -1.05691, 0.08350, 3.003226, 18.3366, -0.7272;
h(13,14) -1.06857, 0.10140, 3.004029, 15.3692, -0.9100; h(10,11) -1.07931, 0.11627, 3.005213, 12.3840, -1.0172; h(9,9)
-1.05828, 0.07617, 3.004842, 10.4102, -0.7694; h(9,8) -1.09147, 0.13205, 3.007025, 9.3878, -0.8200.

Fig. 5 caption (p.293): ha(27,28) asymmetric orbit, x = -1.10662, xdot = 0.03106, ydot = 0.19381, C = 2.99413,
T/2 pi = 29.1734, s = -815.739, s_v = 0.2906 (y = 0 assumed).

Observation (ours): two printed s_v values lie below -1, Table 2 c (-1.0187) and Table 3 h(10,11) (-1.0172), so those
orbits are vertically unstable. The paper says the vertically stable regions shrink as the smooth horseshoes pass closer
to Jupiter, and that all computed orbits of ha(10,11) are vertically unstable.

## 2. Witnesses (Tables 1-3 and the Fig. 5 caption)

- Witness C (decides): the 300-400 dpi page image, read cell by cell.
- Witness A: tesseract 5.5.3, psm 6, on separate 400 dpi crops of each table. Not the OCR text layer.
- `witness-comparison-taylor.tsv`: 102 cells (names not counted).
  - 99 cells agree exactly. 2 cells (s = "> 10^4", rows d and e) are not compared.
  - 1 cell differs: Table 3 h(27,28) s_v, image 0.1801. One tesseract pass gave 0.1501, another pass 0.1801 (a 5/8
    confusion). The image decides: 0.1801.
  - The first values we took from the OCR text layer for the integration run (`integ_taylor.py`) were afterwards
    re-read on the page image and match it on all 18 rows. The OCR text layer shows a spurious "+" before h(10,11) in
    Table 3 that is not on the page.
- Fig. 5 caption: 7 values read on the image (the printed "x =" before 0.03106 is an x-dot with the dot; read on the crop).

## 3. Cross-checks

**3.1 Frame, mu and units confirmed by integration** (`integ_taylor.out`, `taylor_events.out`, `taylor_correct.out`; all
ours; DOP853 rtol 1e-12 to 1e-13, mu = 0.000953875, Sun at (-mu, 0)).
- Jacobi constant from the printed (x, ydot): agrees with the printed C to 2.1e-6 to 4.2e-6 in Tables 1-2 rows a-d, to
  2.0e-6 or better in Table 3 (Table 3 prints C to 6 decimals but x and ydot to 5), and differs by 4.3e-5 (Table 1 e) and
  7.9e-5 (Table 1 f). The e and f differences come from the 5-digit ydot of two very unstable orbits (see below).
- Perpendicularity at the half period. At the y = 0 crossing nearest in time to the printed half period pi T/2 pi (found
  with an event function; offsets 0.000 to 0.076 time units, -0.249 for e) the residual |xdot| is:
  a 2e-6, b 7.5e-5, c 4.1e-3, d 7.7e-3, e 0.18, f 1.7e-2; Table 2 a to d 1.4e-4, 2.5e-4, 2.6e-4, 3.8e-5; Table 3 rows
  1.8e-5 (h(27,28)), 1.8e-4, 7.9e-5, 1.1e-4, 1.1e-3 (h(10,11)), 5.0e-4, 3.9e-4. Compare with the printed digits (5 decimals
  in x and ydot): rows with s above 300 (c to f) lose several digits, as expected. Orbits d to f cross the axis 12 to 32
  times in the window; the nearest-in-time crossing is the one reported, so e's offset of -0.249 shows that choice.
- After a Newton correction of (ydot0, half period) at the printed x0, every one of the 17 symmetric rows closes
  (fsolve residual below 3e-11). The corrections are small: ydot0 by 1.6e-5 or less (e: -1.511214 against -1.51123;
  b +1.4e-5; c -1.2e-5; f +1.1e-5), T/2 pi by 7e-4 or less. So all 17 printed rows are orbits of their families to the
  printed precision.
- Fig. 5 asymmetric orbit: ours C = 2.994126 against 2.99413. One full period (29.1734 x 2 pi) returns within 1.8e-3 in x
  and 2.5e-2 in y with s = -815.739 and 5-digit ICs. Consistent.

**3.2 Stability index** (`taylor_correct.out`, ours; s = half-trace of the corrected orbit's monodromy):

| Row | printed s | ours s |
|---|---|---|
| Table 1 a | 123.25 | 123.2 |
| Table 1 b | 1.1767 | 1.178 |
| Table 1 c | 297.54 | 297.6 |
| Table 1 d, e | > 10^4 | 1.8e4, 1.2e4 |
| Table 1 f | 1.0154 | 1.015 |
| Table 2 a | -10.2204 | -10.22 |
| Table 2 b | -24.206 | -24.25 |
| Table 2 c | 5.9007 | 5.794 (1.8 percent off) |
| Table 2 d | -38.599 | -38.67 |

- Without the correction (5-digit ICs used as they stand), only Table 1 a and b agree (`stability_check.out`: 124.0
  and 1.095). This is a property of the printed precision, not of the definition.
- Table 3 rows must have s = +1 (asymmetric bifurcation, Fig. 2 legend). Ours (x held / C held): h(27,28) 1.002 / 1.010;
  h(20,21) 0.997 / 1.013; h(16,17) 1.013 / 1.038; h(13,14) 1.000 / 1.000; h(10,11) 1.000 / 1.000; h(9,9) 1.24 / 1.23;
  h(9,8) 2.94 / 1.000. Five rows confirm s = 1. h(9,9) stays 1.23 whichever number is held; we did not resolve it.
  h(9,8) is 1.000 only if C (6 decimals) is held, so its printed x is the less precise number.

**3.3 Vertical index.** Not recomputed. The Henon 1973 definition is held (`henon-1973-...aa-28-415...pdf`, digested
2026-10-06). The criterion s_v = +1 at the asymmetric bifurcation is Taylor's statement, not a result of ours.

**3.4 Schanzle 1967 against Taylor 1981** (`crosscheck_schanzle_taylor.out`, `schanzle_vs_taylor.out`; ours). Both use
the same frame and mu (Schanzle does not print mu; the fit is in the Schanzle digest).
- Rabe's orbit (Schanzle: x0 = -1.026903545, ydot0 = 0.039195940, C = 3.00149010, half period 91.7570380, so T/2 pi =
  29.2072) lies between Taylor's asymmetric bifurcation orbit of h(27,28) (x -1.02680, ydot 0.03897, C 3.001491, T/2 pi
  29.1966) and Taylor's orbit 1b (x -1.02745, ydot 0.04032, C 3.00149, T/2 pi 29.2386). Interpolation in x gives ydot
  0.03919 (printed 0.039196), C 3.001491 (3.00149010), T/2 pi 29.2033 (29.2072). So Rabe's orbit is in h(27,28), as Taylor
  says, and lies 1e-4 in x past the asymmetric bifurcation.
- Taylor's orbit 1a (x -1.06201, ydot 0.10851, C 3.00013, T/2 pi 29.2266) lies on Schanzle's continuation: interpolating
  Schanzle's caption orbits x0 = -1.05 and -1.075 (Figs. 11, 12) at x = -1.06201 gives ydot 0.10844, T/2 pi 29.2288.
- The period decides the branch. Schanzle's T/2 over pi falls from 29.207 (x0 = -1.0269) to 29.007 (x0 = -1.9) and 29.009
  (x0 = -2.0). Taylor's rows c to f (the 27-loop end) have T/2 pi 28.0 to 28.2. At equal C, interpolated over Schanzle's
  Table I, Schanzle's orbit has T/2 pi = 29.202, 29.109, 29.034, 29.007 against Taylor c to f 28.2234, 28.1080,
  28.0348, 28.0064: one Jupiter period apart, and different crossing x and ydot (d: -0.766, -0.512 against -0.732,
  -0.578). So Schanzle traced the 28-loop end of h(27,28) (toward Taylor's a and the period near 29 of conjecture (a)), and
  Taylor's c to f are the 27-loop end. The two ends meet at the maximum of C near orbit b. The interpolation over Table I is
  coarse (rows 0.1 apart in x0) but the period gap of 1.0 is far larger than the interpolation error.
- Schanzle p.156 describes the loops growing "until each loop appears to coincide" and ends at a collision orbit
  (x0 = -2.0, "appears to produce a collision" with the Sun). That fits Taylor's conjecture (a) for the 28-loop end.
- Schanzle p.157 states that Rabe's d0 = 1.07 orbit does not belong to the class traced in his Figs. 9-22. This is the
  same point as Taylor's main result, 14 years earlier. Schanzle's cosh(Omega T) = 1.034 for Rabe's orbit (unstable)
  agrees with Taylor's s = +1 at the bifurcation just next to it and s = 1.1767 at orbit b.

**3.5 Mismatches kept visible.**
- Table 1 c, s: printed 297.54, ours 297.6 (agrees). Table 2 c, s: printed 5.9007, ours 5.794 (1.8 percent).
- Table 3 h(9,9) s = 1 not reproduced (1.23). h(16,17) is 1.013 and 1.038 under the two modes.
- Printed reference "Dermott, S. F., Gold, T., Sinclair, A. T.: 1980, Astron. J. 84, 1125". Crossref: Dermott, Gold & Sinclair,
  "The rings of Uranus: nature and origin", AJ 84:1225 (1979), doi 10.1086/112533. The printed year and page differ
  from Crossref. Kept as printed here.

## 4. Citation mining

References (p.294, read on the page image), checked with `ls cyclers_pdf/papers | grep -i` and
`CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Dermott, Gold, Sinclair, AJ 84, 1125 (printed 1980; Crossref says AJ 84:1225, 1979) | NOT held (no `dermott` file; not in index) | no row. Uranus rings; not on the project path. Low priority. |
| Deprit & Henrard 1967, AJ 72:158 | NOT held (no file; index hits only Deprit-Bartholome 1968) | no row. Same issue of AJ as Schanzle 1967. It defines the Hill-equation stability index and the differential corrector. Medium priority as a neighbour of Schanzle. |
| Deprit & Henrard 1970, in "Periodic Orbits, Stability and Resonances" (Giacaglia ed.) | NOT held | no row. Review of periodic motion around L4/L5. Low priority. |
| Henon 1965, Ann. Astrophys. 28:992 | HELD (`henon-1965b-...-ann-astrophys-28-992-...`) | n/a |
| Henon 1973, A&A 28:415 | HELD (`henon-1973-...aa-28-415-ads-1973AA-28-415H.pdf` + `-tables.txt`) | n/a |
| Rabe 1961, AJ 66:500 | NOT held (the `rabe` grep hits only Barrabes files; index has no Rabe row) | no row. The source of the horseshoe orbit and of the starting recipe; high priority for this topic. |
| Stromgren 1929, Publ. Copenhagen Obs. No. 63 | NOT held (the held Stromgren 1933 Bull. Astron. 9:87 is a different paper) | no row. Family a around L3. |
| Szebehely 1967, Theory of Orbits | HELD (`szebehely-1967-theory-of-orbits-...-book.pdf`) | n/a |

- Taylor does not cite Schanzle 1967 (no occurrence in the reference list or the text).
- **New lead (not in CORPUS_INDEX, not in the wanted list):** D. B. Taylor (1983), "Families of asymmetric periodic
  solutions of the restricted problem of three bodies for the Sun-Jupiter mass ratio and their relationship with the
  symmetric families", Celestial Mechanics 29:51-74, doi 10.1007/bf01358598 (Crossref, checked 2026-10-11). It is the
  direct sequel on the ha(n,m) families and is the follow-up to Taylor's conjectures. Proposal: add one wanted-list row.
  The wanted list has no Taylor entry other than row 62, where Taylor 1981 is marked "not found".
- Wanted-list row 62 (Bruno-Varin JAMM 2007 co-orbital background, Schanzle and Taylor named): Taylor 1981 is now in hand
  here. DOI status for the row: Schanzle 10.1086/110211 confirmed (also in the Schanzle digest); Taylor 1981 has no DOI
  (A&A 1981 is not in Crossref).
- Relation to the held Bruno-Varin 2007 JAMM paper (digest `2026-10-06-digest-bruno-varin-2007-jamm-periodic-solutions-small-mass-ratio.md`):
  its sec. 7 (horseshoe and tadpole orbits, pp.953-959) cites Schanzle 1967 and Taylor 1981 as refs in its list [4, 5, 48-51],
  and the 2026-10-06 digest says none of them is held. After filing this paper (and Schanzle 1967) that line is partly
  closed: Taylor 1981 and Schanzle 1967 are held, Llibre-Ollé 2001 and 2004, Barrabes-Mikkola 2005 and Henrard 2002 are not.
  Bruno-Varin treat horseshoes by perturbation theory in mu near mu = 0; Taylor gives the numerical families at mu = 0.00095
  that the theory should reproduce, so it is a numerical check target for their eq. 7.16 (we did not test this).

*Filing proposal:* `cyclers_pdf/papers/taylor-1981-horseshoe-periodic-orbits-restricted-problem-three-bodies-sun-jupiter-mass-ratio-aa-103-288-ads-1981AA-103-288T.pdf`
(force-OCR copy). Table transcription: `data/sources/taylor-1981-horseshoe-tables.yaml` (this file is `taylor-1981-tables.yaml`).
Witness comparison, scripts and outputs beside the PDF as `<pdf stem>-<file name>`.

*Filed as `cyclers_pdf/papers/taylor-1981-horseshoe-periodic-orbits-restricted-problem-three-bodies-sun-jupiter-mass-ratio-aa-103-288-ads-1981AA-103-288T.pdf`. Table transcription: `data/sources/taylor-1981-tables.yaml`.*
