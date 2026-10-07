# Digest: Russell & Brinckerhoff 2008, "Eccentric Orbits around Planetary Moons" (AAS 08-181) (#960)

R. P. Russell and A. T. Brinckerhoff (Georgia Tech; Russell formerly JPL), "Eccentric Orbits around Planetary Moons",
paper **AAS 08-181**, AAS/AIAA Space Flight Mechanics Meeting, Galveston TX, January 2008 (p.1 footnote; "awarded best
paper of conference"). The conference paper has no DOI.
- Journal form: Russell & Brinckerhoff, "Circulating Eccentric Orbits Around Planetary Moons", J. Guid. Control Dyn.
  32(2):424-436 (March 2009), **doi 10.2514/1.38593**. Crossref confirms both authors, volume, issue and title. Crossref
  gives pp. 424-436. The preprint banner (every page) says "pp. 423-435" and asks readers to cite the journal form. **Not
  held**: `ls cyclers_pdf/papers | grep -i "brinckerhoff\|russell"` shows no copy, and CORPUS_INDEX has no row for it. Because the
  journal form is not held, no identity check was possible. This digest covers the preprint only. Wanted-list row:
  `wanted-row.md` (attribution-only, rank 0).
- File: `08Jan_ellipse_AAS-08-181.pdf`, 19 pp., md5 **27bfd2bddc7cdcbd5d3519217e2272a0**. PDF metadata: Word document
  "GalvestonGany29.doc", Acrobat Distiller 8.1.0, created 2008-02-15, modified 2013-10-16 (the banner was added then). The
  text layer is digital and good, so no OCR is needed.
- Proposed filename: `russell-brinckerhoff-2008-eccentric-orbits-around-planetary-moons-AAS-08-181-journal-doi-10.2514-1.38593.pdf`
  (note: the corpus rule is author-year-title-venue-doi; the DOI is that of the journal form, so the name says "journal").
- How I read it:
  - I read the full text layer.
  - I read the page images for Eqs. (1)-(3) (p.3 at 200 dpi), Eqs. (4)-(10) (p.4), Eqs. (11)-(14) (pp.7-8), Eqs. (15)-(18)
    (p.8), Table 2 (p.9, 300 dpi, two crops) and Table 3 (p.11, 300 dpi crop), and Figures 6 and 10 (pp.12, 14).
  - I did not read Table 1 (p.3) or the Appendix (p.18) on a separate image crop; the Table 1 numbers are confirmed by our
    arithmetic (below). The Appendix has one witness (text layer) and a visual scan of the page. Treat the Appendix digits as
    single-witness.
  - Table 2 and Table 3: image read plus text layer (two witnesses, no mismatch).
- Deliverables (scratch folder `dg41-brinckerhoff/`):
  - `russell-brinckerhoff-2008-tables.yaml`: Tables 1-3 and the Appendix gravity field, with a conventions block.
  - `check_russell_brinckerhoff2008.py` and `check_russell_brinckerhoff2008.out`.
  - `wanted-row.md`.

## 0. Verdict

**It is a mission-design review of the doubly averaged third-body problem at moons, plus two long-repeat periodic orbits
at Ganymede. It is not a cycler paper and it has no Europa orbit data. It gives the project sourced formulas and a sourced
region of validity for "circulating eccentric" bound orbits.**
- **What it is.**
  - A review of all motion in Broucke's doubly averaged model (the contour plots of Figs. 2-3), with the full-cycle period
    reduced to a quadrature (Eq. 11).
  - A table (Table 2) of the maximum inclination a "figure-eight" circulating orbit can keep at 25 moons, at a fixed
    moon-to-spacecraft period ratio of 10.
  - Two periodic orbits (12:81 unstable, 9:56 stable) found in the un-averaged Hill model plus a 4x4 Ganymede field, with
    Cartesian initial conditions (Table 3), and a one-year ephemeris check at ten epochs.
- **What it gives the project.**
  1. The doubly averaged equations (Eqs. 4-10), the period quadrature (Eqs. 11-14) and the maximum-inclination relation
     (Eqs. 15-18). All were read on the page image (sec. 1-2).
  2. Table 2: 16 moons with a figure-eight limit, 9 with none. Our recomputation from Eqs. 15-18 reproduces it, with one
     print slip (the Uranus, Neptune and Pluto sub-headers print the Earth mu_p) (sec. 3).
  3. A sourced rule for where averaging holds: spacecraft period at most one tenth of the moon period. Averaged
     predictions "degrade" above that. Un-averaged periodic orbits "cease to exist or find dramatic character changes much
     beyond roughly 2/5 of the Lagrange point distance", which the paper equates to a period ratio of about 7 (sec. 6).
  4. Two Ganymede Cartesian states (Table 3) that are consistent with their printed elements (sec. 4). They need the
     Groove code and a 4x4 field (Appendix) to use as controls.
- **Relation to the Europa CR3BP families (`#956` Europa family identity).** There is no overlap with the 2005 census, so
  this paper does not add a name to the Russell 2005 crosswalk (sec. 7). In short: these orbits are high-altitude, many-
  revolution (56 to 81) orbits in a Hill model, while the 2005 grid caps the search at 16 crossings in a point-mass CR3BP.
- **Gate relevance.**
  - Jovian endgame context: the paper supports the statement that a bound Ganymede orbit at about 0.39 of the Hill-limit collinear-point distance (24 h
    orbit) is outside the formal averaging limit but still survives one year in the ephemeris. It says nothing about
    capture, flybys or leveraging.
  - No cycler and no flyby sequence. The word "circulating" refers to the argument of periapse, not to a cycler.
- **PROPOSALS only (no catalogue or code change):**
  - (a) If a doubly averaged propagator or a Hill-model periodic-orbit corrector is ever added, the Table 2 rows and the two
    Table 3 states are positive-control candidates. The Table 2 T_c column needs the unprinted e0 = 0.001 (sec. 3).
  - (b) Acquire the journal form for attribution and for one identity check of Tables 2-3 (`wanted-row.md`).
  - (c) The Table 2 mu_p slip and the Ganymede altitude inconsistency (sec. 8) should be raised only against the journal
    form, in the respectful-errata framing, not against this preprint.

## 1. Model (pp.2-4, page images)

- The un-averaged model is **Hill's problem** (limit mu_s/mu_p -> 0, "formally (mu_s/mu_p)^(1/3) << 1"), moon-centred, in
  the rotating Hill frame, with the moon rotating synchronously and an n x n spherical-harmonic moon field (Fig. 1). The
  planet is a point mass. The x axis is the moon-fixed meridian and z is the north pole.
- Eq. (1): x'' = 2 N_s v + dGamma/dx, y'' = -2 N_s u + dGamma/dy, z'' = dGamma/dz.
- Eq. (2): Gamma = (1/2) N_s^2 (3x^2 - z^2) + mu/r + U, r = sqrt(x^2 + y^2 + z^2). U is the non-spherical moon potential.
- Eq. (3): C = 2 Gamma - (u^2 + v^2 + w^2), "analogous to the Jacobi constant".
- Units (p.3): N_s normalised to one. Table 1 (Jupiter-Ganymede):
  - mu_Ganymede 9886.99742842995 km^3/s^2.
  - mu_Jupiter 1.26618626797685e8 km^3/s^2.
  - Jupiter-Ganymede distance 1.0704e6 km.
  - Ganymede mean radius 2631.2 km.
  - Hill time unit 98413.2095723724 s, length unit 45749.9268762215 km, mean motion 1.016123754468760e-5 rad/s.
  - Our check: N_s = sqrt((mu_s + mu_p)/a_s^3) and TU = 1/N_s reproduce the printed values. The printed length unit equals
    (mu_s/N_s^2)^(1/3).
- Ephemeris model (p.4, p.15-16): Sun, Jupiter, Saturn and the Galilean moons as point masses (two-body gravity), plus
  oblateness of Jupiter and Ganymede. Ephemerides DE414 and jup230, poles and prime meridians pck00008 and the IAU 2000
  report (ref. 22).
- Gravity field (Appendix, p.18): 4x4 "representative" Ganymede field. Only J2 and C22 were "estimated with reasonable
  confidence based on Galileo flyby data". The rest "are simply representative of an expected field". Values are in
  `russell-brinckerhoff-2008-tables.yaml`.

## 2. The doubly averaged third-body model (pp.4-8, page images)

- **Assumption.** The Hill potential is averaged twice, once over the spacecraft orbit and once over the moon orbit. This
  needs the spacecraft period much smaller than the moon period. "Typically an order of magnitude difference ... justifies
  the use of the averaging approximation [10]" (ref. 10 is Scheeres, Guman & Villac 2001). Point-mass moon only: the
  averaged system has no U term. Semi-major axis is constant.
- **Equations of motion** (Eqs. 4-8, with N_s the moon mean motion and n the spacecraft mean motion):
  - da/dt = 0.
  - de/dt = (15/8)(N_s^2/n) e sqrt(1 - e^2) sin^2 i sin 2w.
  - di/dt = -(15/16)(N_s^2/n) [e^2/sqrt(1 - e^2)] sin 2i sin 2w.
  - dw/dt = (3/8)(N_s^2/n)[1/sqrt(1 - e^2)] [5cos^2 i - 1 + 5 sin^2 i cos 2w + e^2 (1 - 5 cos 2w)].
  - dOmega/dt = -(3/8)(N_s^2/n)[cos i/sqrt(1 - e^2)] (2 + 3e^2 - 5e^2 cos 2w).
  - e, i, w do not depend on Omega, so the system reduces to three variables.
- **Integrals** (Eqs. 9-10, from Broucke 2003, ref. 2):
  - C1 = (1 - e^2) cos^2 i. C1 runs from 0 to 1.
  - C2 = e^2 (2/5 - sin^2 i sin^2 w).
  - The motion stays on a constant (C1, C2) contour in the (e cos w, e sin w) plane. The contour shape is independent of a.
    The traverse rate depends on a through n.
- **Motion types** (Figs. 2-4, text pp.5-7):
  - C1 > 3/5: every contour circulates around the plot centre, so w circulates. Eccentricity grows and shrinks on the
    elongated interior contours.
  - Bifurcation at C1 = 3/5. Two "islands" emerge and move away from the centre as C1 decreases. For e = 0 this is the
    stability boundary i about 39 deg. Our value: arccos sqrt(3/5) = 39.23 deg.
  - C1 < 3/5: contours either **librate** about one island (ovals) or **circulate** around both islands. The circulating
    contour that surrounds both islands, with C2 just above zero, is the **"figure-eight" orbit**. It reaches the highest
    inclination of all circulating orbits (at e = 0) and its maximum e at w = +/-90 deg, where i is about 39 deg.
  - Worked example (p.7): (e0, i0, w0) = (0.001, 56.8 deg, 0 deg) gives C1 about 0.3. It circulates around both islands.
    e ranges about 0 to 0.7 and i about 39 to 56.8 deg (Fig. 4). With w0 = 90 deg the same start librates around the top
    island only, with a similar range of i and e.
- **Frozen orbits.** At the centre of each island "a frozen orbit exists that is stable because the neighboring orbits simply
  librate with near-constant eccentricity and argument of periapse" (p.5). The paper gives no formula and no table of frozen
  (e, i) values. Our derivation from Eqs. (5) and (7) at w = +/-90 deg (de/dt = 0 there; dw/dt = 0 gives
  cos 2i = (1 - 6e^2)/5): i = 39.23 deg at e = 0, 42.36 deg at e = 0.3 and 47.87 deg at e = 0.5. This is ours, not printed.
  It is shown as a check only (`check_russell_brinckerhoff2008.out`). The paper also cites Ely 2005 (ref. 18) for mid-
  inclination frozen orbits at the Moon.
- **Full-cycle period** (Eq. 11): T_c = (16/3)(n/N_s^2) integral from e0 to ef of
  e sqrt(1 - e^2) / sqrt[(2e^2 - 5 C2)(e^2 - 1)(3e^4 + (5C1 + 5C2 - 3)e^2 - 5C2)] de.
  - Lower limit, Eq. (12): e0 = sqrt(5 C2 / 2).
  - Upper limit, Eq. (13): ef = (1/6) sqrt(6 sqrt(25(C1^2 + C2^2 + 2 C1 C2) + 30(C2 - C1) + 9) - 30(C1 + C2) + 18).
  - Eq. (14): the librating lower limit is the same expression with the inner square root negated.
  - The factor 16/3 is for a circulating contour (symmetric in both axes, so the quadrature is multiplied by four). For a
    librating contour, symmetric about the y axis only, the constant is 8/3. T_c is "a function only of n, N_s, C1 and C2".
    It is largest near C2 = 0.
  - The integrand is singular at the limits. The paper says to integrate "slightly within" the bounds, and to avoid C1 = 3/5
    and C2 = 0 exactly.
  - The paper says the quadrature agrees with propagation of Eqs. (5)-(7) ("very accurate").
- **Our checks** (`check_russell_brinckerhoff2008.out`): the integrand, limits and constants as read on the page image give
  Table 2's T_c to within 1 percent in 14 of 16 rows once e0 = 0.001 is chosen (sec. 3). So Eqs. (11)-(13) are read
  correctly. Our value from Eq. (11) for the 12:81 orbit with its printed elements is 68.8 days; the paper quotes 70.3 days
  (not the same input, probably a different start); the Table 3 T_c from the propagated periodic orbit is 77.59 days.

## 3. Table 2 (p.9): figure-eight science orbits at 25 moons

- **Method** (Eqs. 15-18, p.8): for a spacecraft period T = Ts/10, with minimum periapse radius R + 100 km:
  - Eq. (15): a_max = a_s [ (mu_p/mu_s)(Ts/T)^2 ]^(-1/3). This is the two-body semi-major axis whose period is Ts/10.
  - Eq. (16): e_max = 1 - (R + alt_min)/a_max.
  - Eq. (17): i_max = arccos sqrt((1 - e_max^2) 3/5). This equates C1 at w = 0 (e = 0, i = i_max) and at w = 90 deg (e = e_max,
    i about 39 deg, cos^2 = 3/5).
  - Eq. (18): C1 = cos^2 i_max.
  - A negative e_max means no figure-eight orbit exists at this period ratio.
- **Results** (verified on image and text layer; whole table in the yaml). Units km and days:

  | Moon | a_max | e_max | C1 | i_max (deg) | T_c (days) |
  |---|---|---|---|---|---|
  | Moon | 19119 | 0.904 | 0.110 | 70.6 | 832.7 |
  | Io | 3281 | 0.414 | 0.497 | 45.2 | 98.0 |
  | Europa | 4244 | 0.609 | 0.378 | 52.1 | 147.4 |
  | Ganymede | 9856 | 0.723 | 0.286 | 57.6 | 259.6 |
  | Callisto | 15581 | 0.839 | 0.178 | 65.1 | 537.8 |
  | Tethys | 653 | 0.025 | 0.600 | 39.3 | 193.9 |
  | Dione | 1012 | 0.345 | 0.528 | 43.4 | 172.6 |
  | Rhea | 1812 | 0.523 | 0.436 | 48.7 | 210.6 |
  | Titan | 16286 | 0.836 | 0.181 | 64.8 | 515.4 |
  | Hyperion | 691 | 0.663 | 0.336 | 54.5 | 843.4 |
  | Ariel | 1028 | 0.340 | 0.531 | 43.2 | 160.9 |
  | Umbriel | 1366 | 0.499 | 0.451 | 47.8 | 200.4 |
  | Titania | 3234 | 0.725 | 0.285 | 57.8 | 315.5 |
  | Oberon | 4104 | 0.790 | 0.225 | 61.7 | 455.8 |
  | Triton | 4531 | 0.679 | 0.323 | 55.4 | 223.7 |
  | Charon | 1880 | 0.631 | 0.361 | 53.1 | 230.1 |

  - No figure-eight orbit at Ts/T = 10 (N/A): Phobos, Deimos, Amalthea, Thebe, Adrastea, Metis, Mimas, Enceladus, Miranda.
  - Trend (text p.9): maximum inclination grows for larger moons and moons farther from their planet. Moons with larger
    mu_s and larger a_s have higher i_max.
  - Fig. 5 (p.10): i_max and e_max rise as the period ratio falls. Enceladus has a figure-eight orbit at a period ratio of 5, but
    the altitude constraint fails above a ratio of about 7.7. Ganymede, Europa and Titan keep large e across the plotted
    range.
- **Our checks.**
  - Columns a_max, e_max, C1 and i_max reproduce for every Earth, Jupiter and Saturn row (last printed digit) and for the
    Uranus and Neptune rows within 0.1 percent (we use rounded true mu_p).
  - T_c reproduces within 1 percent for 14 of 16 rows when the quadrature starts at initial e0 = 0.001, w0 = 0, that is
    C2 = 0.4 e0^2. The paper does not say which e0 it used for the table. It says "very small but non-zero initial e".
    This is a **hidden input**: T_c grows without limit as C2 goes to zero (e0 = 1e-6 gives about 1.9 times the printed
    value). A reader cannot reproduce the T_c column from the paper alone. Tethys (C1 = 0.600) does not reproduce, because
    it sits on the bifurcation and T_c is very sensitive there. Charon is 1.6 percent off.
  - **Print slip**: the sub-headers for Uranus, Neptune and Pluto each print mu_p = 398479.14 km^3/s^2, which is the Earth
    value. With it, our Ariel a_max is 2507 km against the printed 1028 km. With the true mu_p the rows reproduce. So the
    calculation used the right values, and only the printed header is wrong.
  - Charon: the printed row (1880 km, 0.631) is nearer to Pluto-alone GM (ours 1885 km, 0.632) than to the Pluto-system GM
    (ours 1810 km, 0.617). The paper does not say which was used.
  - Mimas and Enceladus C1 agree to 0.002 (0.172 against 0.174; 0.578 against 0.577).

## 4. Table 3 (p.11): the two Ganymede periodic orbits

- Un-averaged model: Hill's model plus the 4x4 Ganymede field, searched with the JPL "Groove" code (Fortran 90; algorithms
  from refs. 6, 11 and 14). The corrector closes the orbit **in the body-fixed frame** after a stated number of
  revolutions. The search start for 12:81 was {a0 = 12,320 km, i0 = 60 deg, e0 = 0.1, w0 = Omega0 = nu0 = 0}. It needs about 81
  revolutions for one e-w circulation. The text says "The ~80 day period" against "70.3 day period calculated using the
  quadrature".
- State convention: "Initial conditions given in non-rotating frame aligned with the IAU defined Ganymede body-fixed frame
  at epoch" (footnote). Verified by us: elements computed from the printed Cartesian state with the Table 1 mu_s agree with
  all printed element digits (a, e, i, w, Omega, nu) for both orbits, treating the printed velocity as inertial.

  | Property | Units | 12:81 unstable | 9:56 stable |
  |---|---|---|---|
  | x0 | km | -1.10294724E+04 | 1.27215637E+04 |
  | y0 | km | 6.09916977E+02 | 2.74572065E+03 |
  | z0 | km | -9.92043402E-15 | 0.00000000E+00 |
  | u0 | km/s | -2.56020704E-02 | -5.19574390E-01 |
  | v0 | km/s | -4.71573204E-01 | 3.16290562E-01 |
  | w0 | km/s | 8.73907974E-01 | 6.25411458E-01 |
  | a0 | km | 1.23072793E+04 | 1.30393130E+04 |
  | e0 | - | 1.02457203E-01 | 5.05659039E-01 |
  | i0 | deg | 6.16128288E+01 | 5.61929409E+01 |
  | w0 | deg | 2.94588490E-01 | 1.20188898E+02 |
  | Omega0 | deg | 1.76834834E+02 | 1.21794367E+01 |
  | nu0 | deg | -2.94588490E-01 | -1.20188898E+02 |
  | T_c | day | 7.75866851E+01 | 5.70386714E+01 |
  | avg. i | deg | 56.18 | 56.22 |

  - The table prints no Jacobi constant C for either orbit. Figs. 6 and 10 give C only as a plot axis (LU^2/TU^2): about 4.28 to
    4.35 for the 12:81 family and 4.163 to 4.195 for the 9:56 family (read off the plots, approximate).
  - Naming: "12:81" means the spacecraft makes 81 revolutions while Ganymede makes "12+Delta-Omega" revolutions (a node-
    regression term) before the orbit closes in the body-fixed frame. The integer ratio sets the average a (ref. 14).
  - Our arithmetic: Keplerian period at a0 is 0.9986 d (12:81) and 1.0890 d (9:56). Times 81 and 56 revolutions, 80.9 d and 61.0 d.
    These exceed the printed T_c (77.59 d, 57.04 d) by about 4 to 7 percent, as expected for non-Keplerian motion. The 12:81 T_c is
    10.84 Ganymede periods; the 9:56 T_c is 7.97 Ganymede periods.
  - Our C1 and C2 from the printed osculating elements: 12:81 has C1 = 0.224, C2 = +0.0042 (circulating class); 9:56 has
    C1 = 0.230, C2 = -0.030. By the averaged classification (C2 < 0 means the contour lies around one island) the 9:56 start
    would be of librating type. The paper calls both "circulating" or "figure-eight" families. The osculating elements at one
    instant are not averaged, so this is not a contradiction. It is flagged as unresolved.
- **Families** (Figs. 6 and 10, read off the plots; approximate):
  - 12:81: average inclination about 54.1 to 59.05 deg. Node rate -5.1 to -5.65 deg/day. Average e 0.297 to 0.355. Stability
    index b1 about 2 (at the limit), b2 about 2.7 rising to about 13.7 (unstable). Periodicity about 1e-10 (and worse in places).
    Altitude max about 17,000 to 20,000 km; min falls towards zero at the high-inclination end.
  - 9:56: average inclination about 54.6 to 56.8 deg. C about 4.195 to 4.163. Periodicity about 1e-8. The text says all the
    orbits in this family are linearly stable (both indices at most 2).
  - Stability: indices b1, b2 must both be at most two for linear stability. "Most of the orbits of Figure 6 are therefore
    mildly unstable". The periodicity measure: 10^-q is "roughly equivalent to matching q significant digits".
- **Orbit 12:81 over one period** (Figs. 7-9, text p.12): altitude about 1000 to 18,500 km, e from 0.02 to 0.7, i from 45 to 62 deg.
  w circulates once in the ~78 day period. The eccentricity-vector path is a figure eight like the averaged contours.
  The 9:56 orbit has a smaller T_c (about 57 d against about 78 d). The paper links this to a wider "neck" of the figure eight:
  larger neck, shorter period (consistent with the averaged model).
- **Ephemeris propagation** (p.15-16, Figs. 14-16): the 9:56 orbit at epoch 1 January 2028 (JD 2461772.0). Test: it is called
  "long-term stable" if it survives one year at each of ten arbitrary epochs without impact or escape. 9:56 passes. The 12:81 orbit
  also passes, and its minimum altitude dips to about 250 km against about 1000 km in the conservative model. The 9:56
  minimum altitude is about 200 km in both. A spot check with a point-mass Ganymede only gives similar results. No lifetime
  beyond one year is reported, and no failure case is shown.

## 5. Circulating against librating orbits: summary

| Class | Condition (doubly averaged) | w | Close approaches | Inclination |
|---|---|---|---|---|
| Circulating, outer (not around an island) | C1 > 3/5, or C1 < 3/5 outside the islands | circulates | spread over all latitudes and longitudes | maximum at e = 0 |
| Figure-eight (circulating around both islands) | C1 < 3/5, C2 just above zero | circulates | spread over all latitudes and longitudes (the paper's "ball of yarn") | maximum possible among circulating orbits, i about 39 deg to i_max |
| Librating | C1 < 3/5, contour around one island | librates | near the inclination latitude, in one hemisphere | similar range to the figure-eight at the same start |
| Frozen | island centre | fixed at +/-90 deg | fixed | fixed; stable |

- Full-cycle period for the two types: constants 16/3 and 8/3 in Eq. (11) (sec. 2).
- Science relevance (p.7): circulating figure-eight orbits spread close approaches through latitudes and longitudes, reach the
  highest inclinations, have long periods and high altitudes, and "cost less to achieve" than low circular orbits.

## 6. Lifetimes, validity and numbers a reader may want

- **Validity.** Doubly averaged predictions "match closely the full dynamics when the period ratio is at or above about 10";
  for larger ratios they "degrade". Refs. 25 and 26 (Russell 2006; Hénon 1969) show that perturbed Keplerian orbits "cease to
  exist or find dramatic character changes much beyond roughly 2/5 of the Lagrange point distance", "which corresponds to a
  period ratio of about 7". Ganymede at the validity limit: a = 9,856 km (T_s/T = 10); a 100-km-altitude orbiter can keep e up
  to about 0.723 without impact.
- **Ganymede 24-hour orbit** (p.10): an average a of about 12,320 km (our check: 12,319 km, period ratio 7.16), C1 = 0.22, e from
  near circular to about 0.78 and i_max about 61 deg. The paper says this exceeds the formal limit of about 9,856 km.
- **Lifetimes.** The paper reports no impact or escape lifetimes. For circular low orbits it points to Scheeres, Guman & Villac
  2001 (ref. 10), "characteristic instability times", which is the source of its moon list. The only lifetime-type result is the
  ten-epoch one-year ephemeris test of sec. 4.
- **Non-spherical gravity**: "second order" for high-altitude eccentric orbits (p.2). No sensitivity study is given beyond the
  one-line spot check (sec. 4).

## 7. Relation to the CR3BP periodic families in the Russell 2005 Europa digest

Sources: `2026-10-08-digest-russell-2005-global-search-periodic-orbits-near-europa-aas-05-290.md` (sections 1-4), and this paper.

| Item | Russell 2005 (AAS 05-290) | Russell-Brinckerhoff 2008 (AAS 08-181) |
|---|---|---|
| Model | Jupiter-Europa point-mass CR3BP, finite mu 2.528e-5 | Hill's model (mu_s/mu_p -> 0), plus a 4x4 field at Ganymede; ephemeris check |
| Search | brute-force grid, up to N = 16 crossings (32 or 64 per orbit), axisymmetric or doubly symmetric orbits | doubly averaged analysis, then a differential corrector on 56- and 81-revolution repeat orbits |
| Orbit class | N-periodic symmetric orbits, 5 planar families named, 76 printed 3D orbits | long-repeat orbits that close in the body-fixed frame after 56 to 81 revolutions, e and w cycling |
| Europa data | Table 3: 76 rows | none, except one Table 2 row (Europa a_max 4,244 km, e_max 0.609, i_max 52.1, T_c 147.4 d) |
| Named families | Circle-Egg (H1/g1), Egg-Diamond (H2/g2), L1, L2, DRO | none named. Names "12:81" and "9:56" are repeat ratios |
| Reference to 2005 | none | refs. 13 and 25: the 2005 periodic-orbit result is used for (a) the existence of periodic orbits in the un-averaged third-body problem, and (b) the 2/5 Lagrange-point-distance validity bound |

- **No identity link to the 2005 families.** The 2005 grid caps N at 16 and caps the period by that, so an 81-revolution orbit would
  not appear in the 2005 census. In the 2005 grid the closest class by size is the stable 3D direct orbits with x0 of a few thousand
  km (Table 3 rows near x0 = 3.5e3 to 1.1e4 km, pseudo-inclination 50 to 87 deg). They are N-periodic orbits with periods of
  days, not figure-eight repeat orbits.
- **Scale check (ours).** At Europa the Table 2 a_max is 4,244 km. The 2005 Europa L2 is 13,744 km, so a_max is 0.31 of L2, which is
  inside the "2/5 of the Lagrange-point distance" bound. At Ganymede the Hill length unit is 45,750 km (Table 1), so the Hill-limit collinear point is at
  45,750/3^(1/3) = 31,721 km (our arithmetic; the Hill limit, not the finite-mu value). The 24-hour orbit (a = 12,320 km) is then
  0.39 of that distance, at the "2/5" bound, and its period ratio is 7.16, which matches the paper's "about 7". Rough
  comparison only.
- **How the two meet.** The doubly averaged model is the secular limit of the Hill problem at small a/L. The 2005 periodic orbits at
  small x0 and high inclination (the stable direct 3D region) are the periodic-orbit side of the same dynamics. The 2008 paper makes
  this link only in words (refs. 13 and 14).
- Frozen orbits: the 2005 paper has no frozen-orbit concept. This paper's frozen orbits (island centres) are fixed points of the
  averaged system. They would correspond to short-period 3D orbits in the CR3BP. The paper does not make this correspondence
  explicit.

## 8. Print slips and things left open (as printed; not corrected)

- Table 2 sub-headers for Uranus, Neptune and Pluto all print mu_p = 398479.14 km^3/s^2 (Earth's value). The rows use the true values
  (sec. 3). Probable cause: a copy of the first sub-header.
- Table 2 T_c needs an unprinted initial e0 (our fit: 0.001). Tethys T_c is not reproducible.
- p.10 (text) states the 24-hour Ganymede orbit "leads to a 925 x 19,650 km altitude orbit". At a = 12,320 km and e = 0.78 our
  altitudes are about 79 x 19,297 km. The 925 x 19,650 km pair implies a = 12,919 km and e = 0.725 (our arithmetic). The paper's
  conclusion calls the example "200 x 20,000 km" and "9,000 km altitude" for the near-circular phase. These are for osculating
  elements of different orbits (a0 = 12,307 km and 13,039 km in Table 3, and Figs. 7 and 11), and the text does not say which
  orbit each pair belongs to. I leave this open; I do not say it is an error.
- p.9 text says "Table 2 shows that the doubly averaged system is well suited ... when the semi-major axis is less than or equal to
  ~9,856 km". That is a_max at Ts/T = 10, from Table 2. Consistent.
- Eq. (14) and the text: "eccentricity bounds in Eqs. (14) and (13)" for the librating case. The printed pair is (14) and (13). I read it as e0 from
  Eq. (14) and ef from Eq. (13), as the text on p.7 describes.
- p.7: 12:81 start "approximately 81 revolutions" and period "~80 day" against the Table 3 T_c of 77.59 d and the Fig. 7 text "~78 day".
- Citation numbering: refs. 13 and 25 are the same paper (Russell 2006, J. Astronaut. Sci. 54(2):199-226).
- The preprint cites "Groove" as a "prototype" JPL code. It is not public in this form. I did not check for a public release.

## 9. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` (lara, scheeres, paskowitz, johannesen, ely, aiello, kwok, villac, broucke, kozai) and
CORPUS_INDEX / the `#960` wanted list (grep for Broucke, Lara, Scheeres, Guman, Paskowitz, Johannesen, Ely, Aiello, Kwok, Brinckerhoff).

| Ref | Item | Held? |
|---|---|---|
| 2 | Broucke 2003, JGCD 26(1):27-32, "Long-Term Third-Body Effects via Double Averaging", 10.2514/2.5041 | not held. Source of Eqs. (9)-(10) and Fig. 2 |
| 4 | Lara & San-Juan 2005, JGCD 28(2):291-297 | not held |
| 5, 6 | Lara, Russell & Villac AAS 05-384; Lara & Russell AAS 06-168 | not held (AAS 05-384 is the same item as the 2005 digest's ref. 17) |
| 10 | Scheeres, Guman & Villac 2001, JGCD 24(4):778-787, 10.2514/2.4778 | not held. Source of the period-ratio rule and the moon list |
| 11 | Lara & Russell 2007, JGCD 30(1):259-263, 10.2514/1.22493 | not held |
| 12 | Lara, Russell & Villac 2007, JGCD 30(2):409-418, 10.2514/1.22372 | not held |
| 13, 25 | Russell 2006, JAS 54(2):199-226, 10.1007/BF03256483 | not held (the AAS 05-290 preprint is held; see its digest) |
| 14 | Russell & Lara 2007, JGCD 30(4):982-993 "Long-Life Lunar Repeat Ground Track Orbits" | not held; DOI not confirmed on Crossref in one try |
| 15 | Lara, Russell & Villac, Meccanica, 10.1007/s11012-007-9060-z (fast estimation of stable regions) | not held |
| 18 | Ely 2005, JAS 53(3):301-316, "Stable Constellations of Frozen Elliptical Inclined Lunar Orbits", 10.1007/BF03546355 | not held |
| 19 | Russell, "Enceladus Science Orbit Design", JPL IOM 343M-2007-04 (July 2007) | not held; internal memo |
| 26 | Hénon 1969, A&A 1:223 | held (per the 2005 digest) |
| 1, 3, 7-9, 16, 17, 20-24 | Aiello AAS 05-377; Johannesen & D'Amario AAS 99-360; Lara et al. 2005 Chaos; Paskowitz & Scheeres AAS 04-244 and 05-358; Lara 1999 and 2003; textbooks; IAU 2000 report; Kwok AIAA 84-1985 and AAS 91-464 | not held (not individually checked against CORPUS_INDEX beyond the greps above; the `ls` greps found no match) |

- Wanted-list suggestions (proposal): the journal form of this paper (attribution; `wanted-row.md`), and Broucke 2003, Scheeres-Guman-Villac 2001
  and Lara-Russell-Villac 2007 (the three that define the averaged model and the Europa stability regions). Crossref DOIs were
  checked by title, first authors, volume and pages on 2026-10-08.
- Our own corpus check found no row for this paper in CORPUS_INDEX (grep "ellipse", "08-181", "Brinckerhoff").

*To be filed as `cyclers_pdf/papers/<proposed filename above>`. Check script, output and notes are to be filed beside it. Table transcription:
`data/sources/russell-brinckerhoff-2008-tables.yaml` (proposed path).*

*Filed as `cyclers_pdf/papers/russell-brinckerhoff-2008-eccentric-orbits-around-planetary-moons-aas-08-181.pdf`. Tables: data/sources/russell-brinckerhoff-2008-tables.yaml; check script beside the PDF.*
