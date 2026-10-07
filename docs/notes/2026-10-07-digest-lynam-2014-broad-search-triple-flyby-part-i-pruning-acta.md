# Digest: Lynam 2014, "Broad-search algorithms ... Callisto-Ganymede-Io triple flyby sequences from 2024 to 2040, Part I: Heuristic pruning of the search space" (Acta Astronautica) (#960 batch 36, X9 source, #977)

A. E. Lynam (West Virginia University), "Broad-search algorithms for the spacecraft trajectory design of
Callisto-Ganymede-Io triple flyby sequences from 2024 to 2040, Part I: Heuristic pruning of the search space",
Acta Astronautica 94(1):246-252, January 2014. doi 10.1016/j.actaastro.2013.07.018 (Crossref, checked 2026-10-07: title,
author, volume 94, issue 1, pages 246-252 agree). Received 21 March 2013, revised 6 June 2013, accepted 5 July 2013, online
13 July 2013. Open access (CC BY 3.0).
- File given: `5236e3ff-lynam2014.pdf`, 7 pages, text layer (publisher PDF). md5 548a9f72aa82c5791c1cc3e3e5bd624f.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/lynam-2014-broad-search-algorithms-callisto-ganymede-io-triple-flyby-part1-heuristic-pruning-acta-astronautica-94-246-doi-10.1016-j.actaastro.2013.07.018.pdf`
- **Wanted list:** row 67 (low priority). This file completes the Part I / Part II pair.
- **Transcription:** `lynam-2014-part1-tables.yaml`. It holds Tables 1-3, the text counts, and approximate Fig. 2 extents.
- **How I read it.**
  - I read the whole text layer. I then read pp.247 (eqs. 1-7, crop), 248 (eqs. 8-12, Fig. 1, Table 1), 249 (eq. 13, Fig. 2)
    and 250 (eqs. 14-16, Table 2) on 200-dpi renders. I did not render p.251 (Table 3) at 200 dpi in full;
    its text layer sums exactly (see sec. 4), so I trust it. Equations and all numbers in secs. 1-3 were read on the images,
    except where I say "text layer".
  - Page numbers are the printed journal pages (246-252).
  - Arithmetic checks are in `checks_part1.py` and `checks_part1.out`.
  - The paper has **no table of alignment windows, dates or phase-angle values**. It also gives **no coefficients** for
    the pruning polynomials. The window dates are in Part II (held digest, Figs. 4-5).

## 0. Verdict

**A search-pruning method paper. It is not about cyclers, and it prints no window dates.**
- It defines a "Phase Angle Pruning Heuristic". A circular, coplanar, patched-conic model maps all feasible
  Callisto-Ganymede-Io (C-G-Io) triple flybys into four regions of the (dlambda_Ca,Ga, dlambda_Ga,Io) plane. Real
  ephemeris phase angles at each minute are tested against those regions. A minute outside all four is discarded.
- **Result:** of 8,415,358 one-minute epochs in 2024-2040, 82,110 unique epochs survive (0.976 percent). The pruning removes
  99.024 percent. The abstract says "99%".
- **What it gives the project:**
  - one clean statement of the geometric condition for a C-G-Io triple flyby: Callisto within about -4 to +12 deg of
    Ganymede, and Io within about -57 to +95 deg of Ganymede, with the exact region set by flyby type (sec. 2);
  - the phase-angle sign convention and the closed formula (eqs. 15-16) to compute it from any ephemeris;
  - a worked count of how few epochs in 16 years allow a triple flyby (1.5 percent, or 1.0 percent unique).
- **What it does not give:** no window dates, no polynomial coefficients, no mean motions, no mission-level numbers.
  The pruning curves cannot be rebuilt from the paper alone. They must be regenerated with the model of sec. 1 (X9, sec. 7).
- **Catalogue implication (PROPOSAL only):** none. No cycler content.
- **Misprint found:** the text says "55,179 of the 82,110 unique times have two triple flyby solutions" (p.250). Table 2
  gives total 127,289 and unique 82,110. The difference is 45,179. So 55,179 is probably a misprint for 45,179 (sec. 4).

## 1. The patched-conic model (READ, pp.247-248; image-checked)

- **Model name:** "circular, coplanar, ephemeris-free, patched-conic" (CCEFPC), in MATLAB.
  - Moon orbits are circular and coplanar, with radius equal to the semi-major axis.
  - The spacecraft follows four 2-D conics in the moon plane: initial orbit to Callisto's radius, Callisto to Ganymede,
    Ganymede to Io, Io to the final orbit. They are patched by hyperbolic flybys at the moon orbit radii.
  - "Ephemeris free" means the spacecraft meets each moon at its orbit radius. The moon positions in time are ignored.
- **Inputs:** initial a and e (Table 1, sec. 2). The initial orbit is propagated to Callisto's orbit radius.
- **Equations (all read on the image):**
  - (1) cos nu = (a(1 - e^2) - r) / (e r)
  - (2) V = sqrt(2 mu_Jup / r - mu_Jup / a)
  - (3) cos gamma = sqrt( a^2 (1 - e^2) / (r (2a - r)) )
  - (4) V_inf = sqrt( V^2 + V_Ca^2 - 2 V V_Ca cos gamma )
  - (5) cos alpha_in = (V^2 - V_Ca^2 - V_inf^2) / (2 V_inf V_Ca)
  - (6) sin(delta/2) = mu_Ca / ( mu_Ca + (R_Ca + h_p,Ca) V_inf )
  - (7) alpha_out = alpha_in +/- delta (the sign picks energy-reducing or energy-increasing)
  - Here r is the Callisto orbit radius, nu the true anomaly, V the speed and gamma the flight-path angle before the flyby,
    alpha the pump angle (0 to 180 deg), delta the turn angle. I checked eq. (5) by hand: it follows from
    V^2 = V_Ca^2 + V_inf^2 + 2 V_Ca V_inf cos alpha. Eq. (4) is the law of cosines.
  - (8) p = a(1 - e^2)
  - (9) T = sqrt(a^3 / mu_Jup) [ E2 - e sin E2 - (E1 - e sin E1) ] (elliptic transfer)
  - (10) T = sqrt(-a^3 / mu_Jup) [ e sinh H2 - H2 - (e sinh H1 - H1) ] (hyperbolic transfer)
  - E1, E2, H1, H2 are the anomalies at the start and end of each inter-moon transfer. p and T are kept as initial guesses
    for Part II's Lambert solver. Eqs. (6) to (8) are not repeated for Ganymede and Io; the paper says "similar equations".
- **Flyby altitudes (fixed):** Callisto 100 km, Ganymede 1500 km, Io 300 km (p.248, Table 1).
  - Callisto is lowest because it has the least navigation error (ref [21]).
  - Ganymede is high on purpose. It stands in for the energy lost to the 3-D inclination change. The paper says an equatorial
    1500 km flyby has the same energy-change effect as a 300 km flyby at a B-plane angle of 45 deg.
- **Flyby types (Table 1):** Callisto energy-reducing only. Io energy-reducing only. Ganymede energy-increasing or energy-reducing.
  Io before or after perijove. Callisto and Ganymede are always before perijove.
- **Initial orbit range:** a from -1.7e6 km to 1.2e6 km (through infinity at the parabola), e from 0.6 to 1.3. This covers
  ellipses and incoming asymptotes with V_inf under 6 km/s (p.248).
- **Four cases:** (Ganymede reducing or increasing) x (Io before or after perijove). Each has **1780** distinct (a, e) values,
  so 4 x 1780 = 7120 propagations.

## 2. The alignment condition and the pruning booleans (READ, pp.248-249; image-checked)

- **Phase angles at the Ganymede-flyby epoch (eqs. 11, 12):**
  - dlambda_Ca,Ga = (nu1 - nu2) + n_Ca T_Ca,Ga
  - dlambda_Ga,Io = (nu3 - nu4) + n_Io T_Ga,Io
  - dlambda_Ca,Ga is the angle from Ganymede's position to Callisto's position at the time of the Ganymede flyby.
  - dlambda_Ga,Io is the angle from Io's position (at the Ganymede-flyby time) to Ganymede's position.
  - nu1 is the true anomaly just after the Callisto flyby, nu2 just before the Ganymede flyby, nu3 just after the Ganymede
    flyby, nu4 just before the Io flyby. n_Ca and n_Io are the moon mean motions (the paper does not print their values).
  - Reading: Callisto moves n_Ca T_Ca,Ga during the transfer, the spacecraft sweeps nu1 - nu2 (a negative number for
    forward motion, so the printed form is the sum of the spacecraft's swept angle with the sign flipped and the moon's motion).
    I do not re-derive the sign here. The printed definitions match eqs. (15)-(16).
- **Fig. 2 (p.249):** the four sets of 1780 points plot as four roughly quadrilateral regions in the (dlambda_Ga,Io [x],
  dlambda_Ca,Ga [y]) plane. Axes: x from -60 to about +95 deg, y from -4 to 12 deg. Extents read from the image, approximate to
  about 2 deg:
  - Io after perijove (circles = Ganymede energy-reducing, crosses = Ganymede energy-increasing): x about -57 to +48 deg,
    y about -4 to +12 deg. The two sets overlap over a large area, with the crosses spanning the larger range.
  - Io before perijove (dots = reducing, plus signs = increasing): x about +40 to +95 deg, y about -2 to +12 deg. This region is
    much smaller than the after-perijove region.
  - The smallest of the four regions is the energy-increasing, before-perijove set; Table 2 agrees (11,052 minutes).
- **Pruning boolean (eq. 13), for each of the four cases:**
  - B = (dlambda_Ga,Io <= b0 + b1 dlambda_Ca,Ga + ... + bn dlambda_Ca,Ga^n)
    AND (dlambda_Ga,Io >= c0 + c1 dlambda_Ca,Ga + ... + cn dlambda_Ca,Ga^n)
    AND (dlambda_Ca,Ga <= f0 + f1 dlambda_Ga,Io + ... + fn dlambda_Ga,Io^n)
    AND (dlambda_Ca,Ga >= g0 + g1 dlambda_Ga,Io + ... + gn dlambda_Ga,Io^n)
  - bi, ci, fi, gi are polynomial coefficients. The degree n is **not stated**, and **no coefficient is printed**. The four
    curves are the thick black lines of Fig. 2, fitted to the four edges of each quadrilateral. There are 4 sets of 4 curves.
  - If at least one of the four booleans is TRUE, the epoch is feasible. If all four are FALSE, the epoch is discarded.
  - Feasible regions overlap, so one epoch can match two booleans (a Ganymede-reducing and a Ganymede-increasing solution).
    Each match carries its own initial guess to Part II.
- **Tolerance:** there is no explicit tolerance. The bounding curves are the envelope of the propagated points. The condition is
  "inside the polygon". The envelope is exact only to the polynomial fit and to the model error of the CCEFPC.
- **Why the model error is small enough:** the moon orbits are nearly circular and coplanar (all within 1 deg of
  inclination, p.250), and the three flybys come in rapid succession (Discussion, p.251).

## 3. Ephemeris test, discretisation and pruning result (READ, pp.249-250; image-checked)

- **Ephemeris:** `jup230l.bsp` (Jacobson 2003), read by SPICE. Positions of Callisto, Ganymede and Io every **1 minute**
  over eight 2-year intervals between 2024 and 2040. The split is only for MATLAB array limits. One 2-year interval has about
  1.05 million minutes (1,051,920 for 2 x 365.25 d; checked). Per interval: 3 moons x 3 components x 1.05 million = about
  9.45 million numbers (9,467,280 checked) and 2.10 million phase angles (2,103,840 checked).
- **Orbit normal (eq. 14):** n = (r_Ga,1 x r_Ga,2) / ||r_Ga,1 x r_Ga,2||, with r_Ga,2 the Ganymede position **19 minutes**
  after r_Ga,1. It only supplies a positive direction for the angles.
- **Phase angles from ephemeris (eqs. 15, 16):**
  - dlambda_ephem,Ca,Ga = Sgn( n . (r_Ga x r_Ca) ) arccos( (r_Ga . r_Ca) / (|r_Ga| |r_Ca|) )
  - dlambda_ephem,Ga,Io = Sgn( n . (r_Io x r_Ga) ) arccos( (r_Io . r_Ga) / (|r_Io| |r_Ga|) )
  - The printed denominator is written as the norm of the dot product. That is a typesetting slip. The right denominator is the
    product of the two vector lengths, as I wrote above. I read the image; the layout of both fractions is the same.
  - The sign is positive when the second moon in the name is counter-clockwise of the first about n. That gives the angle from
    Ganymede to Callisto for eq. (15) and from Io to Ganymede for eq. (16), as eqs. (11)-(12) require.
- **Test and compression:** apply the four booleans to every minute. Keep the minute in the vector of each case whose boolean is
  TRUE. Aggregate over the eight intervals.
- **Table 2 (p.250; image-read):** see the YAML. Case vectors: 14,747 + 11,052 + 45,512 + 55,978 = 127,289 (1.513 percent of
  8.42e6). Unique 82,110 (0.976 percent). Pruned 8.33e6 (99.024 percent). My recomputation with the abstract's 8,415,358:
  pruned 8,333,248 minutes, 99.024 percent. The per-row percentages all agree to the printed three decimals.
- **Cost saved:** the abstract says a blind search needs 8,415,358 Lambert calls to find 127,289 triple flybys. After pruning,
  82,110 unique epochs (127,289 solutions) go to the Lambert solver. The call count falls by a factor of about 66 (my division),
  for the unique epochs, or by a factor of 102 if each of the 82,110 is one call (8,415,358 / 82,110 = 102.5). The paper does not
  state which. It states only "99%".
- **Run time (Table 3, p.251):** 167.459 s in total. 162.773 s (97.2 percent) is SPICE reading the ephemeris. The model takes
  0.132 s, the interpolation structures 0.239 s, and the phase-angle and reduction step 4.315 s.
  - The author notes that a faster ephemeris reader, compiled code or parallelism could cut this.

## 4. Validation, and what is not validated

- **How the author validates it:** by construction, not by independent test.
  - The regions come from the same simplified model that the heuristic stands for. There is no printed comparison with
    known triple flybys, no false-negative rate and no positive control in Part I.
  - The check on "does pruning lose real solutions" is Part II. Part II runs Lambert pathfinding on the surviving epochs and
    finds triple flybys (the Part II digest reports 11 windows targeted, six feasible with Earth and Mars flybys, and a December
    2029 CGPI case). That shows the surviving epochs contain feasible solutions. It does not show that no feasible epoch was pruned.
  - A missing positive control is a gap for X9. See sec. 7.
- **Arithmetic checks (`checks_part1.out`):**
  - Table 2 rows sum to 127,289. All five percentages recompute to the printed digits. 8,333,248 pruned = 99.024 percent.
  - 16 years x 365.25 d x 1440 min = 8,415,360. The abstract's 8,415,358 differs by 2 minutes (leap-year or end-point choice).
    The text's 8.41 million and the table's 8.42e6 are both roundings (8.415e6 rounds either way).
  - Table 3 rows sum to 167.459 s. The percentages recompute.
  - **Disagreement:** the text says 55,179 of 82,110 unique times have two solutions. Total minus unique is 127,289 - 82,110
    = 45,179. If 55,179 were right, unique would be 72,110. The Table 2 row values and the sum agree with each other, so the
    printed 55,179 is the odd one out. I read "55,179" twice on the p.250 image. Probable misprint for 45,179. Unresolved.
  - The Part II digest describes the same data with the same numbers (127,289 and 82,110) and does not repeat 55,179.

## 5. How Part I feeds Part II

From the held Part II digest (sec. 1) and Part I text (p.249-250, "Part II [23]"):
1. **The list of candidate Ganymede-flyby epochs.** The up to four time vectors (82,110 unique epochs). Each epoch is the
   Ganymede-flyby time. Part II reads the Galilean positions and velocities at those times.
2. **Four interpolation structures** (p_Ca,Ga, p_Ga,Io, T_Ca,Ga, T_Ga,Io) per flyby case, each a 2-D nearest-neighbour map
   (MATLAB `TriScatteredInterp`) from (dlambda_Ca,Ga, dlambda_Ga,Io) to the output. Part II queries them at the epoch's phase
   angles for initial guesses for the Lambert p-iteration and for the transfer times.
   (The paper says "nearest-neighbor" and also names TriScatteredInterp, which is triangulation-based; I record both as printed.)
3. **Double solutions** are handled by giving the energy-reducing and energy-increasing cases separate guess sets. The Lambert
   solver then finds both without any extra logic.
Part II does not need the polynomial coefficients except to build the epoch list.

## 6. Other content

- **Discussion (p.251):** three other orderings (Callisto-Io-Ganymede, Ganymede-Io-Callisto, Io-Ganymede-Callisto) and other
  moon triples (Callisto-Ganymede-Europa, Callisto-Europa-Io) are named as possible, not done. The paper suggests Uranian moons
  and an Earth-Venus-Earth-Mars analogue, with lower accuracy because of Venus and Mars inclinations.
  - This is the author's own pointer toward applying the idea at planets. It is a pointer only; there is no analysis.
- **Introduction:** reports Galileo's Io assist saved 175 m/s of capture delta-V, and that Lynam et al. showed double and triple
  satellite-aided capture save a further 230 m/s and 350 m/s (p.247). I have not checked those figures.
- **Conclusions:** the pruning removes 99 percent of infeasible solutions. The remaining 1 percent can go to a Lambert solver.

## 7. For X9 (the alignment-census pre-screen, #977)

X9 asks which repeat counts k bring a moon pair's phase back to within a corrector's basin. Part I gives the
**geometric membership test**. Lynam 2015 (held; its digest is `2026-10-07-digest-lynam-2015-triple-quadruple-satellite-aided-captures-cmda.md`)
gives the **closed-form drift**. X9 can use both. The rules below are implementable now.

**Rule 1 (phase angle).** At epoch t compute r_Ca, r_Ga, r_Io from the ephemeris. Take n from two Ganymede positions 19 min
apart (any approximately correct orbit normal works). Then
`dl_CaGa = sgn(n . (r_Ga x r_Ca)) * arccos(r_Ga . r_Ca / (|r_Ga||r_Ca|))` and
`dl_GaIo = sgn(n . (r_Io x r_Ga)) * arccos(r_Io . r_Ga / (|r_Io||r_Ga|))`.
Lynam 2015 uses the opposite sign for dl_CaGa (its eq. 1; see that digest). X9 must fix one convention in code and add a unit
test with a hand-made configuration.

**Rule 2 (membership).** Feasible if `dl_GaIo` lies between the lower and upper bounding curves of the case AND `dl_CaGa` lies
between its left and right curves, for any of the four cases. Part I gives no coefficients. X9 must regenerate the region:
1. Build the CCEFPC model of sec. 1 (eqs. 1-10 plus the Ganymede and Io legs, same equations).
2. Sweep (a, e) over a in [-1.7e6, 1.2e6] km (through infinity) and e in [0.6, 1.3], 1780 points per case, with hp 100, 1500, 300 km.
3. Take the convex hull, or a fitted polynomial envelope, of each case's (dl_GaIo, dl_CaGa) points.
The hull is simpler and has no fit error. Part I's own polynomial envelope is slightly looser or tighter by an unknown amount.

**Rule 3 (positive control, mandatory).** Part I prints no control. The ratio check that X9 can use today needs only the
mean motions, which the paper does not print. With standard periods (my values, not from this paper: Io 1.769138 d,
Ganymede 7.154553 d, Callisto 16.689018 d; `checks_part1.out`), the Callisto-Ganymede synodic period is 12.523 d and the
Ganymede-Io synodic period is 2.350 d. Then 4 x 12.523 = 50.09 d and 16 x 2.350 = 37.60 d. These match the 50.12-50.13 d C-G window
spacing and the 37.6 d G-Io clock quoted in the X9 entry (those two numbers come from Lynam 2015; I did not check them
there). The 0.03-0.04 d gap on the first is of the size of the moon-orbit precession, so I do not call it a disagreement. The
real positive control for the membership test: run Rules 1-2 on 2024-2040 and compare the count with Table 2
(127,289 total, 82,110 unique, per-case 14,747 / 11,052 / 45,512 / 55,978). A match within a few percent is the control.
Also check that the epoch lists contain Part II's published windows (CGIP 2 Apr 2026, 14 Jun 2033; CGPI 2 Apr 2026, 2 Dec 2029,
15 Jan 2033, 31 Dec 2034, from the Part II digest Table 1). If X9's membership test prunes any of these, the test is wrong.
Only then trust an all-negative census.

**Rule 4 (drift, from Lynam 2015).** Both phase angles drift linearly at (n_Ca - n_Ga) and (n_Ga - n_Io), with period
12.523 d and 2.350 d. Window spacing for the pair is then k times the synodic period, k integer. For a cycler on one pair
the same arithmetic gives the k that returns the phase to within the corrector's basin. The tolerance is the corrector's
basin, not Part I's polygon.

**What X9 should not copy.** Part I's polygon encodes triple-flyby feasibility for three moons in one pass. A cycler repeats a
two-moon pattern. X9's pair test uses Rule 4, not Rule 2. Rule 2 is for the three-moon C-G-Io chain (candidate C13-type checks), and as
a worked, ephemeris-level control. The CCEFPC model ignores moon eccentricity and inclination. The paper says the effect is small
for the Galilean moons. X9 must not carry the argument over to other systems without checking it.

**Suggested X9 script output.** Per case: epoch count, first and last epoch, longest gap, and window start dates, to compare with
Part II Figs. 4-5.

## 8. Citation mining (26 references, pp.251-252)

Held status checked with `ls cyclers_pdf/papers | grep -i` and `grep -i CORPUS_INDEX.md` on 2026-10-07. Wanted-list row from the
list as it stands (row 67 is this paper).

| ref | work | status |
|---|---|---|
| [19] | Lynam, Kloster & Longuski (2011), CMDA 109:59-84 | HELD (author manuscript, `lynam-kloster-longuski-2011-...-author-manuscript.pdf`) |
| [20] | Lynam & Longuski (2011), JGCD 34:1485-1494 | HELD |
| [22] | Lynam & Longuski (2011), "Laplace-resonant triple-cyclers for missions to Jupiter", Acta Astronaut. 69:158-167 | HELD (`lynam-longuski-2011-laplace-resonant-triple-cyclers-...`) |
| [23] | Lynam (2014) Part II, Acta Astronaut. 94:253-261 | HELD (with the #943 collision check) |
| [21] | Lynam & Longuski (2012), "Preliminary analysis for the navigation of multiple-satellite-aided capture sequences at Jupiter", Acta Astronaut. 70:33-43 | not held; not wanted-listed (checked: no row). Navigation only. |
| [25] | Strange, Russell & Buffington (2007), "Mapping the V-infinity globe", AAS 07-277 | HELD (`strange-russell-buffington-2007-mapping-v-infinity-globe-AAS-07-277.pdf`, from the Part II digest; the Strange name is not in my grep of the Lynam files) |
| [18] | Strange et al. (2012), AIAA 2012-4518 (SEP tours for Jupiter) | not held; no wanted row. Low priority. |
| [17] | Landau, Strange & Lam (2010), AAS/AIAA SFM | not held. Low priority. |
| [16] | Johannesen & D'Amario (1999), AAS 99-330 | not held (the Johannesen file in the corpus is a different paper, per the Part II digest) |
| [9]-[12] | Longman 1968; Longman & Schneider 1970; Cline 1979; Nock & Uphoff 1979 | none held; no wanted rows. Satellite-aided capture history. |
| [13]-[15] | MacDonald & McInnes 2005 (JGCD 28:365); Yam 2008 PhD; Okutsu, Yam & Longuski 2007 (AAS 07-258) | not checked individually beyond a name grep (no hit for these titles); capture and end-of-life escape, not cycler content |
| [24] | Vallado, Fundamentals of Astrodynamics, 3rd ed. | textbook; the corpus holds only Vallado 1991 per the Part II digest |
| [26] | Jacobson (2003), Jup230 ephemeris (ssd.jpl.nasa.gov) | data source, not a paper. **A reproduction of the X9 control needs `jup230l.bsp`.** I did not check whether the project already holds a Jovian-satellite SPK. |
| [1]-[8] | Galileo and Cassini operations papers | none needed |

- New candidates from this paper: none beyond [21] (low priority) and the Jupiter-satellite SPK (data, not a paper; for X9).

## 9. Unresolved

- The polynomial degree n and all coefficients of eq. (13). Not in the paper. May be in the author's thesis (wanted row 13, full
  text) or Part II's appendix (the Part II digest does not mention one). Not checked.
- The 55,179 versus 45,179 count (sec. 4).
- Whether the "factor 100" saving counts epochs or solutions (sec. 3).
- I did not read p.251 (Table 3) and the references at 200 dpi.

*Filed as `cyclers_pdf/papers/lynam-2014-broad-search-algorithms-callisto-ganymede-io-triple-flyby-part1-heuristic-pruning-acta-astronautica-94-246-doi-10.1016-j.actaastro.2013.07.018.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/lynam-2014-part1-tables.yaml`.*

*Wanted-list row numbers in this digest are the pre-batch-36 numbering; the list was renumbered in batch 36.*
