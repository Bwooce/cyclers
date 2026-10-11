# Digest: Bayliss 1970, "Precision Targeting for Multiple Swingby Interplanetary Trajectories", MIT MSL report TE-39 (corpus policy; #960)

S. S. Bayliss, Ph.D. thesis, Department of Aeronautics and Astronautics, MIT, submitted 11 May 1970. Issued as
Measurement Systems Laboratory report TE-39, June 1970 (NASA CR-109936, NTRS 19700024926). The title page says
S.B. MIT 1965, S.M. MIT 1966. The acknowledgements say the work was sponsored in part by NASA grant
NGR-22-009-207 and that computing was done under MIT Computation Center problem M4131.
- Supervisors on the signature page: W. M. Hollister (chair), W. E. Vander Velde, W. McKay, R. H. Miller. The
  report is approved by the Director of MSL (the signature reads W. Markey). The acknowledgements single out Hollister and thank R. Stern.
- Filed as `cyclers_pdf/papers/bayliss-1970-precision-targeting-multiple-swingby-interplanetary-trajectories-mit-msl-te-39-ntrs-19700024926.pdf`
  (the NTRS scan as received, 182 PDF pages, 3.2 MB, md5 f7a1bdd0dfa3f867c784f504e5612278). The OCR copy is
  `cyclers_pdf/papers/bayliss-1970-precision-targeting-multiple-swingby-interplanetary-trajectories-mit-msl-te-39-ntrs-19700024926-ocr.pdf`
  (md5 73c8b5cf5e1156b10b7efde9966631f0, 2.3 MB, `ocrmypdf --force-ocr -l eng`, no recompression needed).
- The scan carries a poor letter-spaced NTRS text layer. The OCR copy replaced it. No DOI (report). Not held before
  this filing: no Bayliss, TE-39 or 19700024926 hit in CORPUS_INDEX.md, the digests or data/sources.
- Tables transcribed to `data/sources/bayliss-1970-te-39-tables.yaml`.
- I read every page of chapters 1, 5, 6, 7 and 8, Appendix C, the end of Appendix D, Appendix E, the contents,
  the figure list and the references as page images. I read the method chapters 2, 3 and 4 (pp.21-76) in the
  OCR text and looked at their figures. I did not follow the derivations in Appendices B and D line by line.

## 0. Verdict

**A method thesis. It builds accurate reference trajectories for multiple-swingby flights. It is not a study of
cycler guidance.** The method turns a simple patched-conic solution into a finite-sphere-of-influence solution,
then adds analytic perturbation offsets, then checks each leg by numerical integration. Three worked examples
test it.
- Example 1 is an Earth-Venus-Mars-Earth free-fall reconnaissance flight (ch. 5). The whole 1.29-yr flight is
  matched to a total correction of 0.2263 m/s.
- Example 2 is a free-fall Grand Tour, launch 1978, 10.77 yr (ch. 6). The total correction is 2.652 m/s.
- Example 3 is a 3.59-yr segment of Earth-Venus periodic orbit I (ch. 7). **No free-fall solution was found.**
  A powered solution needs a total impulse of 220.534 m/s. Its analytic prediction is wrong by a total of
  38.950 m/s against numerical integration.
- This is the only cycler content. It is a real-ephemeris, perturbed-conic test of one segment of Hollister's
  orbit I. It is a candidate negative control for our R1(b) V-tier work, with the caveats in sec. 7.
- **Hickman is not reproduced. See sec. 1.** The wanted-list row for Hickman 1968 stays open.

## 1. The Hickman question (owner: "the content is IN this report")

**Answer: no. TE-39 cites Hickman three times and carries none of Hickman's analysis, numbers or tables. Row 5
of the wanted list must stay open.** The evidence is the full text, read as images.

The reference entry, printed p.187, reference 24:

> 24. Hickman, D.E., "Guidance Requirements for Periodic Orbits", S.M. Thesis in Dept. of Aeronautics and
> Astronautics, M.I.T., July, 1968

The three citations:
1. p.14 (sec. 1.2), a list of mission types with their references: "7) Earth-Venus and Earth-Mars Periodic Orbits
   [22,23,24,25]". Next to it: "The majority of these studies have been concerned primarily with preliminary
   mission planning and guidance requirements studies using simplified models for targeting."
2. p.101 (sec. 7.1): "The existence of such trajectories was first explored by Hollister in [22]. ... This work was
   extended by Menning [23] who found additional Earth-Venus trajectories and by Rall [25] who examined
   Earth-Mars periodic trajectories. The guidance requirements for the Earth-Venus trajectories were examined by
   Hickman in [24]."
3. p.105 (sec. 7.2): the first iteration cycle uses "initial conditions obtained from [22] and [24]".

What that means.
- Bayliss describes Hickman in one sentence. He gives no result of Hickman's: no guidance cost, no dispersion
  figure, no station-keeping number, no table, no figure.
- Bayliss's own contribution to the periodic-orbit case is not a guidance analysis. It is a targeting run on a
  3.59-yr piece of orbit I. He did not set out to "reproduce, extend or supersede" Hickman. His objectives
  (p.13) are reference-trajectory determination, aiming points along the trajectory, and mission analysis.
  He lists "midcourse guidance requirements" and "the effects of guidance and navigation inaccuracies" among
  applications still to be done with his method (sec. 8.2, p.123).
- The one place Hickman's work enters the numbers is the starting point. Sec. 7.2 says the iteration started from
  conditions taken from [22] and [24]. The starting point equals Hollister's published orbit I dates (sec. 4),
  so the data may come from Hickman, from Hollister, or from both. The text does not say which numbers came
  from which.
- The list of tables (p.11) has only C.1 and C.2. There is no table of Hickman's.
- Hickman's thesis is dated July 1968 and sits in the same department as Menning (July 1968) and Rall
  (October 1969). The bibliographic details in the project's wanted row are now confirmed by an independent
  source: S.M. thesis, Department of Aeronautics and Astronautics, MIT, July 1968.
- NTRS search on 2026-10-11 (the public API, title, author and topic queries, 1968 to 1972 filter) returns
  Bayliss and Rall and nothing for Hickman. The thesis was not found on NTRS by these queries.

**Can the row be marked "content held in another form"? No.** The wording is used for rows where a held
document reproduces the content (Hollister 1968 RE-36 is covered by Hollister 1969). Here nothing held reproduces
Hickman's guidance analysis. Annotation for row 5 (I add it with the post-window commits):
"Not reproduced in Bayliss 1970 (TE-39); cited there as ref [24] (confirms: S.M. thesis, Dept. of Aeronautics
and Astronautics, MIT, July 1968; Bayliss p.105 took initial conditions from [22] and [24]). Not on NTRS (searched
2026-10-11). Still wanted. Route: MIT Libraries thesis scan or interlibrary loan."

## 2. Place in the Hollister lineage

| Year | Work | Status here | What it does for the Earth-Venus orbits |
|---|---|---|---|
| 1963 | Hollister, Sc.D. thesis | held, digested | Mars via Venus; no repeating E-V orbit |
| 1968 | Hollister, RE-36 | wanted row 2, content in the 1969 JSR | Defines orbits I, II, III |
| 1968 | Menning, S.M. thesis | held, digested | 15 Earth-Venus orbits; patched conic, zero-radius SOI |
| 1968 | **Hickman, S.M. thesis** | **wanted row 5, not held** | **Guidance requirements of the E-V orbits** |
| 1969 | Hollister, JSR 6(4) | held, digested | Orbits I, II, III in the circular and real models |
| 1969 | Rall, Sc.D. thesis TE-34 | held | Earth-Mars periodic orbits |
| 1970 | **Bayliss, Ph.D. thesis TE-39** | **held now** | Precision reference-trajectory method; orbit I segment as test 3 |
| 1970 | Hollister and Menning, JSR 7(10) | held, digested | 15 E-V orbits, Table 3 |
| 1971 | Rall and Hollister, JSR 8 | held, digested | Earth-Mars orbits |

- The thesis order is Menning and Hickman (July 1968), Rall (Oct 1969), Bayliss (May 1970). Bayliss is the last
  thesis and the only piece in the sequence that leaves the zero-radius patched-conic model.
- Which orbit it targets: orbit I only (Hollister's orbit I, Menning's orbit 1 and 1H). It uses the 16-yr real
  inclined-elliptic version. It takes the first 3.59 yr of the 16-yr period.
- It does not target orbits II or III, any of Menning's orbits 2 to 15, or any Earth-Mars orbit. The word "Mars"
  occurs only in example 1 and in the lineage sentence on p.101.
- Delta-V and dispersion numbers: sec. 4.

## 3. Map of the report (printed page numbers; PDF page = printed page minus 1 to 7)

| Ch. | Pages | Content |
|---|---|---|
| 1 | 13-19 | Objectives, existing targeting methods, synopsis |
| 2 | 21-34 | Simple patched conic model; advanced patched conic model with finite SOI; cost function; accuracy |
| 3 | 35-52 | Perturbed conic analysis: mid-point method, quadrature, velocity offsets |
| 4 | 53-78 | Trajectory segment matching: Lambert partials, gradient, steepest descent, optimum gradient step, conjugate gradient, Newton-Raphson, turn-angle constraint |
| 5 | 79-88 | Example: Earth-Venus-Mars-Earth reconnaissance |
| 6 | 89-99 | Example: Grand Tour, 1978 launch |
| 7 | 101-116 | Example: segment of Earth-Venus periodic orbit I |
| 8 | 117-123 | Conclusions, contributions, recommendations |
| A-E | 125-183 | Notation; Lambert partial matrices; planetary data (ephemeris, SOI); two-body state transition matrix; trajectory data tables |
| refs | 185-188 | 38 references |
| bio | 189 | Biographical note |

**Method in brief (chs. 1 to 4).**
- Step 1, simple patched conic (sec. 2.1): heliocentric conics matched in |v_inf| at each planet, as in Hollister
  and Menning. Dates are varied until the |v_inf| at each flyby match.
- Step 2, advanced patched conic (sec. 2.2): the flight becomes alternating heliocentric and planetocentric conic
  arcs, matched in position and time at the SOI entry and exit points. The SOI is an "enlarged" sphere, the locus
  where the planet's direct acceleration equals the solar gradient acceleration, radius a (mu_p/mu_s)^(1/3). It
  is 1.6 to 2.8 times the Laplace radius (Table C.2). Entry and exit points are given by azimuth, elevation and
  time. An N-encounter flight has 2N-2 matching points. A cost function (sum of squared velocity mismatches) is
  minimised over the point coordinates.
- Step 3, perturbed conic (ch. 3): perturbations of each two-body leg by the other bodies are computed by
  quadrature with the two-body state transition matrix. Perturbations are propagated from the leg mid-point to
  both ends, because propagating from one end diverges (sec. 3.1). Each leg gets velocity offsets at its ends so it
  keeps the same end points and times.
- Step 4 (ch. 4): the cost function includes the offsets and is re-minimised. Methods: steepest descent, a modified
  "optimum gradient step" (OGS), conjugate gradient, generalized Newton-Raphson (GNR). GNR converges
  quadratically but needs a zero minimum. The turn-angle rule keeps flyby periapsis above 1.1 planet radii
  (eqs. 4.87 to 4.90).
- Step 5 (sec. 5.4): each leg is integrated numerically from the analytic start state. The velocity and position
  errors measure the method.
- Error sources (sec. 3.3): the dominant error is the use of the two-body state transition matrix and of
  disturbing accelerations evaluated on the reference arc, plus the linear theory. Appendix D.2 compares
  quadrature with integration for four legs. The differences are about 1 percent or less.

## 4. The Earth-Venus example in detail (chapter 7, Appendix E.4)

**The sequence (p.102).** 1 Earth-to-Venus transfer; 2 direct return to Venus; 3 direct return to Venus; 4
Venus-to-Earth transfer; 5 direct return to Earth; 6 Earth-to-Venus transfer. "The direct return legs are
trajectories which return to the launch planet one planetary year later." It "forms a part of a larger (16 year)
orbit which repeats indefinitely between Earth and Venus". Twelve SOI points, eleven legs, 3.59 yr (1311 d).

**Match to Hollister's orbit I (computed by us from the printed dates).** The SOI midpoints of Bayliss's passes
equal the flyby dates of Hollister's real-model orbit I (Hollister 1969 Table 1; Menning A-2 list for 1H), to
within about one day.

| Flyby | Hollister JD-2440000 | Bayliss SOI midpoint (ours) |
|---|---|---|
| Earth, start | 806 | exit point 811.30 (about 5 d after the flyby, the SOI crossing time) |
| Venus | 971 | 971.44 |
| Venus (return) | 1196 (Menning 1H list) | 1195.79 |
| Venus | 1421 | 1420.09 |
| Earth | 1592 | 1592.92 |
| Earth | 1957 | 1957.47 |
| Venus | 2125 | entry 2122.37 after the correction in sec. 8, about 2.8 d before the flyby |

So the segment is orbit I from its Earth flyby at JD 2440806 to its Venus flyby at JD 2442125, as the owner
expected. The start conditions are Hollister's, or Hickman's, or both (sec. 1).

**Leg data (printed, p.103-105).** Launch C3 18.427 km2/s2. Heliocentric legs have a = 0.864 au (E-V),
0.723 au (V return legs), 0.9999 au (E return), e from 0.076 to 0.177, i from 4.1 to 7.3 deg. The five
planetocentric hyperbolas are listed with periapsis radius and turn angle.

| Pass | Periapsis radius (planet radii) | Turn angle (deg) | v_inf from periapsis speed (ours, km/s) |
|---|---|---|---|
| Venus 2-3 | 2.308 | 28.87 | 4.995 |
| Venus 4-5 | 1.486 | 36.24 | 5.002 |
| Venus 6-7 | 2.361 | 28.54 | 4.989 |
| Earth 8-9 | 2.331 | 35.69 | 4.372 |
| Earth 10-11 | 1.151 | 47.71 | 4.367 |

The Earth 10-11 pass is at 1.151 radii, which is near the 1.1 radius floor of eq. 4.90. Hollister 1969 Table 1 gives, for
these same flybys, v_inf 0.178 EMOS at Venus, 0.154 at Earth and 0.154 at the Earth launch. Bayliss gets 0.168
at Venus, 0.147 at Earth and 0.144 at launch (1 EMOS = 29.78 km/s; launch from C3). All are 5 to 6 percent
lower, including the launch boundary, which carries no impulse. I note the observation only. The cause is not
determined: finite-SOI leg timing, the different ephemeris and the value of the EMOS unit are all candidates.

**The result (p.105-113).**
- The generalized Newton-Raphson step "resulted in a rapid divergence from the minimum" because the cost has no
  zero. The run used the optimum gradient step only, in five cycles, about 450 s on an IBM 360/65.
- The first cycle (finite-SOI two-body model, no perturbations) starts at a cost of 0.61931 km2/s2 and is still at
  0.0006961 km2/s2 after 200 OGS iterations. The later cycles add the perturbation offsets. The cost then
  starts at 0.023412 and ends at 0.011491 at the start of the fifth cycle. The run stopped because the solution
  point, the cost and the weights stopped changing in a useful way. Both are non-zero minima, so the free-fall
  closure fails at the finite-SOI step and the perturbations make the residual about 16 times larger.
- A velocity impulse is needed at each SOI point. In m/s: point 2 0.181; 3 9.935; 4 10.380; 5 69.168; 6 70.561;
  7 0.607; 8 2.240; 9 27.545; 10 27.402; 11 2.515; points 1 and 12 none. **Total 220.534 m/s over 3.59 yr.**
  The large impulses sit at the ends of the three direct-return legs (3-4, 5-6, 9-10).
- Analytic prediction against numerical integration: the per-point error is 0.24 to 8.86 m/s, total 38.950 m/s.
  Final position errors after one integration run from the analytic start velocity are 143 to 8008 km (sec. 7.4).
  The return legs carry offsets of 136 to 181 m/s, "an order of magnitude larger than those for previous
  trajectories".
- Bayliss's reading (sec. 7.5): the solution "may not be assumed to be a global minimum"; the failure of GNR
  "indicates strongly that no zero of the cost function exists locally but this has not been proven"; "some
  experimentation using different sets of initial conditions was conducted but no firm conclusions on the
  existence of multiple solutions were reached"; the direct-return legs "obviously spend a large time in the
  vicinity of the planet" so the two-body transition matrix is a poor model there.

## 5. Tables

All numeric tables are in `data/sources/bayliss-1970-te-39-tables.yaml`. They are:
- Chapters 5, 6 and 7: SOI point dates, azimuth and elevation; leg times, central angles and elements; periapsis
  data; perturbations (end points and mid-points); velocity offsets; comparison with numerical integration;
  periapsis differences; the correction or impulse totals; the cost-function histories.
- Appendix C: Table C.1 (gravitational parameters), Table C.2 (Laplace and enlarged SOI radii), and the
  eight mean-element expansions (Mercury to Neptune, T in centuries from JD 2415020.0).
- Appendix D.2: four quadrature-versus-integration comparisons.
- Appendix E.2, E.3, E.4: planet positions, planet velocities, SOI points in planetocentric coordinates, true
  heliocentric and planetocentric velocities, for the three examples (6, 8 and 12 points).
- The printed errors found by the checks are in the `errata:` block, not corrected in the values.

Selected numbers (printed): the Appendix C enlarged SOI radii in km are Venus 1,458,966.1; Earth 2,157,378.4;
Mars 1,564,377.2; Jupiter 76,638,659.8; Saturn 94,129,489.7; Uranus 101,287,919.7; Neptune 167,883,945.9. Planet
radii used: Venus 6050, Earth 6378.2, Mars 3410, Jupiter 71400, Saturn 60400, Uranus 23500 km (ours, from the
printed periapsis tables).

## 6. Figures

Seventeen figures. I viewed the plates of Figs. 4.1 to 7.3 and B.1 to C.2. Figs. 2.1 to 3.3 are described from their captions and
the text only, because I did not view those plates. I did not digitise any curve.
- Fig. 2.1, 2.2, 2.3 (pp.32-34): sketches of a simple patched-conic flight, the advanced model with SOI circles at
  each planet, and the entry and exit point coordinates (azimuth, elevation, time) on a SOI.
- Fig. 3.1, 3.2, 3.3 (pp.50-52): the two offset calculation methods and the disturbing-acceleration geometry.
- Fig. 4.1 (p.77): contours of cost with a zig-zag steepest-descent path in a narrow ravine. Fig. 4.2 (p.78): the
  turn angle between incoming and outgoing asymptotes inside a SOI.
- Fig. 5.1 (p.88): heliocentric plot of the Earth-Venus-Mars-Earth flight, six numbered points, orbits of Venus,
  Earth and Mars, scale bar 1 au.
- Fig. 6.1, 6.2 (pp.98-99): the Grand Tour in the ecliptic plane to Neptune (10 au bar), and the Earth-to-Saturn
  segment (2 au bar) with points 1 to 5.
- Fig. 7.1, 7.2, 7.3 (pp.114-116): the periodic-orbit segment in three plots, points 1-4, 4-8 and 8-12, with
  the Earth and Venus orbits (1 au bar). The lines wind about Venus's orbit and cross Earth's orbit near the
  flybys. The direct-return legs are the near-circular loops that stay close to Venus's or Earth's own orbit.
- Fig. B.1 (p.155): Lambert problem geometry. Fig. C.1 (p.165): orbital element coordinate system. Fig. C.2
  (p.166): geometry for the SOI perturbing acceleration.

## 7. Relevance to this project

1. **Orbit I with a finite-SOI and perturbed closure test.** This is the only source in the corpus that tries to
   close a real-ephemeris Earth-Venus orbit I segment with finite SOIs and then with planetary perturbations. It
   fails to find a free-fall solution at the finite-SOI step and at the perturbed step. The powered fix is 220.534
   m/s over 3.59 yr, about 61 m/s per year, concentrated on the Venus-return and Earth-return legs. This is useful
   context for any V-tier statement on R1(b) members: the patched-conic closure of orbit I does not carry over
   to the finite-SOI model without small manoeuvres.
2. **Caveats, so the number is not over-read.** The model is a 1970 mean-element ephemeris with a perturbed-conic
   method. The solution is a non-zero local minimum from one starting point, not a proof of non-existence. The
   method's own prediction error on these legs is 38.95 m/s, which is 18 percent of the impulse. The thesis
   itself says the direct-return legs are the weakest part of the method. Our project already treats such a
   result as conditional on method (negative-results rule). It does not show that orbit I is not a cycler.
3. **Positive control for a real-ephemeris check.** Bayliss lists the SOI-point dates, positions and true
   velocities of this segment in Appendix E.4. They are a published data set for a Hollister orbit I segment in a
   real ephemeris. Our own real-ephemeris closure of the same dates can be compared with them. Expected
   difference: finite SOIs and a different ephemeris. The values carry the printed errors in sec. 8.
4. **Not a source for Hickman's guidance numbers.** For the V-tier question on station-keeping cost, Bayliss gives
   only the 220.5 m/s figure above, for a model that is not Hickman's.
5. **Method value.** The advanced patched conic with an enlarged SOI is the model between zero-radius patched
   conics and N-body integration. The thesis argues for the enlarged SOI from numerical tests by the author and by
   Carlson [29] (sec. C.2). Breakwell and Perko [27] are cited in sec. 1.2 as an analytic study of patched-conic
   accuracy.
6. **Examples 1 and 2 are not cycler data.** Example 1 (ch. 5) is a one-shot Earth-Venus-Mars-Earth flight. Example
   2 (ch. 6) is a Grand Tour with a Jupiter swingby at 22.6 planet radii (turn 25.6 deg).

## 8. Notes on the text and printed errors

The three-witness work found six errors in the report. They are recorded in full in the YAML `errata:` block.
- **E.4 point 12 date.** The Julian date of point 12 (sec. 7.1, p.102) is printed 2442123.36905. It should be
  2442122.36905, one day earlier. The leg 11-12 time of flight (159.40312 d) does not equal the printed JD
  difference (160.40312 d). The E.4 Venus position at point 12 matches an independent evaluation of the Appendix C
  expansion at JD 2442122.36905 (108 km, as at the other Venus rows). It is 3.03 million km from the position at the
  printed JD.
- **E.4 point 10, x.** The planetocentric x is printed +1,637,219. It must be negative: point 10 is an SOI entry
  point and the dot product of position and velocity is +1.4e6 as printed, against -7e6 to -1e9 at all other
  entries. The printed azimuth 183.950 deg also needs a negative x.
- **E.2 point 5, z.** The Mars exit point z is printed 1,060,012 km. The radius closure and the printed elevation
  3.885 deg need about 106,000 km.
- **E.3 point 2, v_z.** The true heliocentric v_z of the Jupiter entry is printed +0.060104. Velocity closure
  needs -0.060101. The print shows no minus sign.
- **Calendar dates.** The printed dates "October 6, 1978", "July 13, 1989", "August 11, 1970" and "March 13, 1974"
  do not follow from the printed JDs by the standard conversion. In chapter 6 the printed date is the day after
  the standard date (JD 2443787.0 is 1978 Oct 5 12:00 and JD 2447720.0 is 1989 Jul 12 12:00). In chapter 7 the
  printed dates are about 1 to 4 days before the standard dates (JD 2440811.298 is 1970 Aug 12 19:09; JD
  2442122.369 is 1974 Mar 15 20:51). Treat the JD columns as the authority.
- **Neptune.** The printed Neptune element expansion does not reproduce the printed Neptune position in E.3 (z 0.28
  percent off); the other six expansions I tested do, to 216 km or better. Cause not found.

Per the project's respectful-errata rule, none of this says the work was careless. It was typed from 1970 line
printer output. The closure checks are ours and falsifiable from the printed tables.

## 9. OCR and three-witness record

- Witness A: a fresh `ocrmypdf --force-ocr -l eng --jobs 4` layer (182 pages, about 16 min on a loaded Mac). The
  NTRS text layer was not used.
- Witness B: the page renders read by eye (110 dpi for all pages read; 200 and 260 dpi crops for the cells that
  needed a closer look). The values in the YAML were typed from these.
- Witness C: closure checks by script. They cover the following.
  - Appendix E: 372 numbers against the OCR text (23 not found verbatim; see the list below), then SOI-radius closure (26 points), velocity closure
    (helio minus planet equals planetocentric), azimuth and elevation closure, the sign test of position against
    velocity at entry and exit points, and C3 closure at launch.
  - Planet velocities: evaluating eq. C.4 with the printed expansions reproduces every printed planet-velocity cell
    (all three examples, 26 rows) to 4e-5 km/s for Earth, Venus and Mars, 1.5e-4 for Uranus and Neptune, 1.3e-3 for
    Saturn and 5.8e-3 for Jupiter. E4 point 12 matches only at JD 2442122.36905 (4e-5 km/s), which confirms the
    date erratum a second way.
  - Appendix C: the printed expansions evaluated by us reproduce the printed planet positions at the SOI times to
    216 km or better for Earth, Venus, Mars, Jupiter, Saturn and Uranus. This checks the digits of those
    expansions and of the chapter JDs together. Enlarged SOI radius equals a (mu_p/mu_s)^(1/3) for all planets
    with an astronomical unit of 149599000 km.
  - Chapters 5 to 7: totals (0.2263, 2.652, 220.534, 38.950) equal the sums of their rows; leg times equal the JD
    differences; central angle minus 180 deg equals turn angle for all 10 planetocentric legs; periapsis radius
    equals a (e - 1); heliocentric central angles equal the angle between the SOI-point vectors to 0.13 deg; the
    per-point errors in sec. 7.4 lie between |a - b| and a + b of the neighbouring leg errors.
- Disagreements between A and B were all OCR errors, except where the printed value is itself wrong (sec. 8).
  Examples: OCR read E.2 planet v_x at point 2 as -21.125118 (the render and closure say -21.115118); read
  224.513 as 224.515; dropped the leading 2 of 205.609; read 101.424 as 107.424; read 147,025.929 as 147,028.929.
- The 23 Appendix E values that did not appear verbatim in the OCR text. Planet-velocity cells: E2 points 2, 5, 6;
  E3 points 1, 2, 8; E4 points 2, 4, 5, 7, 10. All were confirmed by the ephemeris velocity check above (the OCR
  had misread the digit or the minus sign; E2 point 2 reads -21.125118 in the OCR and -21.115118 in the render, and
  the closure agrees with the render). Position cells E3 point 2 y, E4 points 3 x and 10 y: punctuation-only OCR
  errors, confirmed by the radius closure and the ephemeris. True heliocentric velocity cells E2 point 5, E3 point
  8, E4 points 9, 11, 12 and planetocentric cells E2 point 4, E4 point 4: re-read on 260 dpi crops. Those crops agree
  with the typed values. Points 8 (E3) and 12 (E4) are arrival points with a dash for the planetocentric velocity, so
  the crop is the only second reading for their heliocentric velocity cells.
- 44 of 613 chapter numbers did not appear verbatim in the OCR text (mostly lost leading zeros, commas read as
  full stops, or garbled JD lines on the speckled p.102). Each was read again from the render. Of those, the JDs are
  cross-checked by the leg times and the ephemeris, the azimuth and elevation by the Appendix E coordinates,
  the central angles by the turn angles, and the periapsis values by a (e - 1). A small group of cells
  has only the render reading: ch. 7 perturbation cells 44.509 and 59081.77 (OCR dropped a leading digit of the
  second), and the periapsis-difference cells 38.27 and -6.991. The per-point error cell 5.532 is covered by the
  triangle check, and 0.397 was confirmed on a 200 dpi crop.

## 10. Citation mining (38 references, printed pp.185-188)

Status key: held, wanted (row in 2026-10-05-960-wanted-papers.md), or neither.

| Ref | Item | Status |
|---|---|---|
| 1 | Minovitch 1964, JPL TR 32-464 | held |
| 2 | Hunter 1964, Astronautics & Aeronautics 2(5) | wanted, row 50 |
| 3 | Battin 1964, Astronautical Guidance (book) | neither; row 37 names Battin 1959 and 1999 |
| 4 | Battin 1959, JAS 26(9), round-trip interplanetary trajectories | wanted, row 37 |
| 5 | Minovitch 1965, JPL TR 32-849 | held |
| 6 | Myers 1969, JSR 6(6), solar probe gravity-assist trajectories | neither |
| 7 | Sturms and Cutting 1965, AIAA 65-90 | wanted, row 61 |
| 8 | Sturms 1967, JPL TR 32-1062 | held |
| 9 | VanderVeen 1969, Bellcomm TM-69-1013-2 | held |
| 10 | Young 1967, MIT EAL PR-3, "The Two-Planet Flyby Problem" | neither |
| 11 | VanderVeen 1968/69, AAS 68-114 and JSR 6(4):383 | JSR held; AAS paper not held (attribution only) |
| 12 | Hollister 1963, Sc.D. thesis | held |
| 13 | Hollister 1964, AIAA 64-647, Mars transfer via Venus | wanted, row 38 (journal form held) |
| 14 | Sohn 1964, JSR 1(5) | held |
| 15 | Crocco 1956, Proc. VII IAC | held |
| 16 | Silver 1967, AIAA 67-613, Grand Tours of the Jovian planets | neither |
| 17 | Friedlander 1966, JSR 3(10), Grand Tour guidance analysis | neither |
| 18 | Deerwester 1969, NASA TN D-5271 | neither |
| 19 | Deerwester 1966, JSR 3(10), Jupiter swingby missions | neither |
| 20 | Kingsland 1968, AIAA 68-1055, Grand Tour trajectory analysis | neither |
| 21 | Long 1969, Astronautics & Aeronautics 7(6) | neither |
| 22 | Hollister 1968, MSL RE-36 | wanted, row 2 (content held in the 1969 JSR) |
| 23 | Menning 1968, S.M. thesis | held |
| 24 | Hickman 1968, S.M. thesis | **wanted, row 5, open** (sec. 1) |
| 25 | Rall 1969, Sc.D. thesis (TE-34) | held |
| 26 | McDonald 1965, MIT EAL RE-19, interplanetary trajectory program for guidance and navigation | neither |
| 27 | Breakwell and Perko 1966, Progress in Astronautics and Aeronautics 17 | content held in the AIAA 65-689 paper (same title) |
| 28 | Slater and Stern 1966, MIT EAL RE-18, simplified midcourse guidance | neither |
| 29 | Carlson 1969, Ph.D. thesis MIT, analytic guidance for many-body trajectories | neither |
| 30 | Zadeh and Desoer 1963, Linear System Theory | neither (textbook) |
| 31 | Danby 1964, AIAA J. 2(1) | wanted, row 34 |
| 32 | Fletcher and Reeves 1964, Computer J. | wanted, row 37 |
| 33 | IBM System/360 Scientific Subroutine Package manual | neither (manual) |
| 34 | Shak, Buehler, Kempthorne 1964, J. SIAM 12(1) | neither |
| 35 | Escobal 1968, Methods of Astrodynamics | neither (textbook; source of the Appendix C expansions) |
| 36 | Stern 1963, Sc.D. thesis MIT, interplanetary midcourse guidance analysis | neither |
| 37 | Goodyear 1966, NASA CR-522 | named as a candidate in two digests, not a wanted row |
| 38 | Hildebrand 1956, Introduction to Numerical Analysis | neither (textbook) |

**Proposals for the wanted list (the lead decides; none is added by this digest).**
1. **Hickman 1968**: keep row 5, add the annotation in sec. 1.
2. **Young 1967, MIT EAL PR-3 (two-planet flyby problem).** Same lab and decade as the E-V orbits. Check NTRS and
   DSpace. Medium priority.
3. **Carlson 1969, Ph.D. thesis MIT** and **Stern 1963 Sc.D. thesis**. They bear on the validity of the
   enlarged-SOI patched conic (Carlson) and on MIT guidance practice (Stern). Medium to low.
4. **McDonald 1965 RE-19** and **Slater and Stern 1966 RE-18**, the MIT numerical tools Bayliss and the earlier
   theses built on. Low.
5. Low priority, attribution only: Silver, Friedlander, Deerwester x2, Kingsland, Long, Myers (all outer-planet or
   solar-probe).
6. Same-lab guidance lineage worth a check: the other MSL TE-series reports. TE-34 is Rall and TE-39 is Bayliss.
   A search for MSL reports between them (guidance of free-fall orbits) is cheap.
