# Digest: Uphoff, Roberts and Friedman 1976, "Orbit Design Concepts for Jupiter Orbiter Missions", JSR 13(6):348-355 (#960)

C. Uphoff (Vector Sciences, Inc., on a JPL contract), P. H. Roberts and L. D. Friedman, Jet Propulsion Laboratory.
Journal of Spacecraft and Rockets 13(6):348-355, June 1976. DOI 10.2514/3.57096. Presented as AIAA Paper 74-781, AIAA
Mechanics and Control of Flight Conference, Anaheim, 5-9 August 1974. Submitted 2 January 1975, revision received
10 November 1975. 8 pdf pages. Wanted-list row "Clipper and Jovian tour background" (second item; Campagnola and
Kawakatsu 2012 is the first). With this paper the row has no unreceived item left.

- **Filed as** `cyclers_pdf/papers/uphoff-roberts-friedman-1976-orbit-design-concepts-jupiter-orbiter-missions-jsr-13-348-doi-10.2514-3.57096.pdf`.
  Table transcription: `data/sources/uphoff-roberts-friedman-1976-tables.yaml`.
  - The filed file is the received file, unchanged, md5 d31e092e08af7c5dc25af7ef6343c13d. It is a 300 dpi bilevel scan
    (CCITT) with a FineReader OCR layer, so no OCR was run here. Pages are printed pp.348-355 (PDF page 1 = p.348).
  - DOI, volume, issue and pages were confirmed from the page headers, the PDF subject field and the Crossref lookup of
    2026-10-11. The reference list is on p.355.
- **How it was read:** all 8 pages on a 130 dpi render, and Figs. 3 and 5 on 300 dpi crops.
  - The "text layer" is an OCR layer. The two text extractions (raw order and layout) come from one OCR pass, so they are
    not independent on OCR misreads. The render is the check that matters.
  - Table 1 (4 by 4) and every number of the text, Table 2 and the figure labels agree between the OCR and the render.
    The OCR layer scrambles the row order of Table 2 and drops minus signs and subscripts in the equations. Those were
    read from the render.
  - One figure label stayed illegible: the leftmost rotated period label on the top axis of Fig. 5.

## 0. Verdict

**The earliest gravity-assist satellite-tour design paper in this corpus, and the source of the terms orbit pumping (credited
to Beckman), orbit cranking and resonance hopping. It contains no periodic or repeating trajectory and no cycler.** It states the method (zero-sphere-of-influence patched
conics with the Jacobi integral conserved), gives one table of ideal transfer parameters between Galilean moon pairs, and
gives one table of a 42-encounter sample tour that is deliberately non-repeating.
- **Concepts (pp.348-351).**
  - Orbit pumping (a term the paper credits to Beckman): satellite gravity changes the jovicentric energy.
  - Orbit cranking: satellite gravity changes the inclination or flight path angle at the encounter with no change of energy.
  - Resonance hopping: pumping with transfers of n times pi (n even), so that the spacecraft returns to the same satellite;
    the hop moves the orbit from one resonance P_sc = (n/m) P_s to another. Each resonance has one jovicentric speed
    V_n,m at the satellite's distance, and the hop must reach the speed of the next resonance.
  - Odd n-pi transfers connect different satellites with high inclination. They need a fixed semilatus rectum
    p = 2 r1 r2 / (r1 + r2) (Table 1).
- **Earlier held work.** Niehoff 1970 and 1971 ("Touring the Galilean Satellites", held) is older. It uses capture orbits whose
  period is commensurate with the moons, with no gravity assists, so repeat encounters come from period matching. Uphoff,
  Roberts and Friedman add the gravity assists that move the orbit between resonances.
- **Repeated moon.** The paper never calls a sequence periodic or says it repeats. The nearest statements are these:
  - "by returning again and again to the same satellite it is possible to change the orbit period by amounts that would
    require in excess of 5 km/s" (p.353, Fig. 5 pump curves);
  - a spacecraft can alternate an inner and an outer satellite to pump the energy down without lowering the periapsis (p.350);
  - "it is equally feasible to transfer from one satellite to another exactly as we fly from one planet to the other"
    (p.353).
  These are tour and endgame statements. None is a closure condition.

## 1. What the paper does

- **Geometry (pp.348-349, Fig. 1).** Zero-sphere-of-influence patched conic. The excess velocity keeps its length across the
  flyby and turns by an angle a, with sin(a/2) fixed by the satellite's gravity and the minimum flyby radius (Eq. 2).
  The largest velocity change is dV = 2 V_inf sin(a/2) (Eq. 3). It cannot exceed the circular speed at the minimum flyby
  radius, which occurs at a = 60 degrees.
- **Jacobi integral (p.349, Eq. 1).** E_b - n' h_b cos i_b = E_a - n' h_a cos i_a with E, h the jovicentric energy and angular
  momentum, n' the satellite mean motion and i the inclination. The paper states that this is the Jacobi constant of the
  circular restricted three-body problem, for a spacecraft far from the satellite before and after an arbitrarily close flyby.
- **Cranking (pp.349-350, Fig. 2, Eqs. 4-5, Fig. 3).** For a spacecraft already in resonance, the largest crank per flyby is
  delta = arccos((cos a - cos^2 beta) / sin^2 beta), the maximum inclination follows from conservation of the angle between
  the jovicentric velocity and the satellite velocity, and the number of flybys to reach it from zero is 90 degrees over
  delta. Fig. 3 shows one-to-one cranking at Callisto. For a 1300 km flyby (0.5 Callisto radius) the text gives 50 degrees
  in 0.75 years, periapsis from 4.2 to 20 R_J, about 20 encounters. Cranking holds the periapsis nearly fixed at low inclination.
- **Perijove (p.350).** Lowering perijove raises the radiation risk. Alternating an inner and an outer satellite pumps the
  energy down without lowering the periapsis.
- **Odd n-pi transfers (pp.350-351, Table 1).** The discrete values of p that permit high-inclination transfers between pairs
  of Galilean moons, for circular coplanar orbits and R_J = 71372 km.
- **A sample high-inclination transfer (p.351).** Ganymede outbound on 24 August 1984 (4300 km altitude, inclination 27
  degrees, 37.4 day orbit, perijove 6.0 R_J) to Europa on 28 September 1984 (35 days). The return is a 540 degree transfer
  from Europa inbound on 2 November to Ganymede outbound on 9 December, with the Europa flyby above 6000 km. Between the
  two the spacecraft hopped to the nearest Europa resonance to wait for Ganymede to move into position.
- **Perturbations (p.351).** Oblateness, solar gravity and non-participating Galilean satellites. Low-order resonance perturbation:
  period about 400 days, amplitude 0.35 R_J, total semi-major-axis variation about 50,000 km in about 200 days.
  Standoff encounters under 100,000 km are avoided.
- **Navigation (pp.351-352, Fig. 4).** A priori satellite position errors of 50-100 km, radio-alone 1-sigma errors of about
  400 km falling to 50-200 km, optical navigation under 50 km. About 5 to 10 m/s of navigation delta-V per orbit and under
  300 m/s for a whole tour. A "multinominal" strategy lets later encounters absorb uncorrected errors.
- **Performance (pp.352-353, Fig. 5).** Launch in 1981-1982 with C3 up to 80 km^2/s^2, arrival V_inf under 6 km/s, July 1984
  to April 1985. A baseline insertion delta-V of 1290 m/s. A powered Ganymede swingby insertion at 13 R_J with 1740 m/s
  saves 450 m/s against direct perijove insertion, and the same delta-V gives an 80 day orbit instead of a 400 day one.
  Pump curves, one-to-one crank curves and the circular-orbit curve for each moon are drawn.
- **Mission design and sample tour (pp.353-355, Table 2, Fig. 6).** The large-initial-orbit strategy, apsidal rotation by
  inbound or outbound encounters (prograde on inbound encounters that lower perijove or outbound ones that raise it), and a
  42-encounter tour with 28 Ganymede, 9 Callisto and 5 Europa passes. It lowers perijove from 14 to 8.4 R_J while the period
  falls from 105 to 14 days, and it includes a bow-shock excursion about two and a half years after insertion.

## 2. Tables (transcribed in `data/sources/uphoff-roberts-friedman-1976-tables.yaml`)

- **Table 1 (p.351), semilatus rectum p in R_J for odd n-pi transfers:**

  | | Io | Europa | Ganymede | Callisto |
  |---|---|---|---|---|
  | Io | 5.9 | 7.3 | 8.5 | 9.6 |
  | Europa | 7.3 | 9.4 | 11.6 | 13.9 |
  | Ganymede | 8.5 | 11.6 | 15.0 | 19.1 |
  | Callisto | 9.6 | 13.9 | 19.1 | 26.4 |

- **Table 2 (p.354), sample satellite tour:** 14 printed rows, from the insertion on 19 July 1984 at Ganymede (1000 km, smear
  velocity 6.8 km/s) to orbit 42 (perijove 11.3 R_J, period 57 days, no encounter printed). Columns: orbit number, perijove
  radius, period, date, satellite, inbound or outbound, flyby altitude, latitude, smear velocity. All 14 rows are in the YAML.
  The "smear velocity" column is not defined in the text.

## 3. Independent checks, all ours

- **Table 1 from its own diagonal.** The diagonal (5.9, 9.4, 15.0, 26.4) is the orbit radius of each moon in R_J. The formula
  p = 2 r1 r2 / (r1 + r2) reproduces all ten distinct off-diagonal values to the printed digit. Io-Callisto is 9.64 from the
  rounded diagonal (printed 9.6); standard orbit radii give 9.7, which is Io's 5.91 rounding. So the table is self-consistent.
- **The sample transfer has one misprinted value.** The text prints a transfer perijove of 5.12 R_J, a semilatus rectum of
  11.396 R_J and a semi-major axis of 44.38 R_J.
  - With a = 44.38 and p = 11.396 the eccentricity is 0.862 and the perijove is 6.12 R_J. With the printed 5.12 and
    11.396 the eccentricity would be 1.226, which is not an ellipse.
  - With 6.12 the true anomalies at Ganymede (15.0 R_J) and Europa (9.4 R_J) are 106.2 and 75.7 degrees. They sum to 181.9
    degrees, so this is an odd n-pi arc with n = 1. The flight time from outbound Ganymede to inbound Europa is 35.2 days,
    the printed 35.
  - So 5.12 is inconsistent with the paper's own p, a and flight time, and 6.12 is consistent with all three. It is probably a
    misprint. The OCR layer and the render both read 5.12. It is recorded as printed.
  - The "45 R_J" of the text is the 37.4 day pre-swingby orbit: a 37.4 day period gives a = 45.17 R_J.
- **Table 2 intervals.** Consecutive printed encounters are 112, 41, 2, 41 days apart in 1984-1985 and 10, 2 days in
  December 1985; Europa 15 and 14 days in July-August 1986; 29 days at orbit 41. The interval agrees with the printed
  period of the later row to within 3 days for orbits 2, 3, 34, 35 and 41 (orbit 1 has 112 days against a printed 104).
  The encounter count is 28 + 9 + 5 = 42, the number of orbits.
- **A probable misprinted date in Table 2.** The note "(Bow shock excursion ~ 12-25-85)" sits between orbits 41 and 42
  (25 November 1986 and later). The text puts the excursion two and a half years after insertion on 19 July 1984, which is
  mid-January 1987. A 25 December 1985 date falls in orbit 20, whose period is 15 days, so it cannot be a bow-shock
  excursion. It is probably 12-25-86. Recorded as printed.
- **Fig. 3 time axis.** 1000 days against 60 revolutions is 16.7 days per revolution, Callisto's period. The
  zero-altitude inclination peak is near 59 degrees at about 18 revolutions, which supports "nearly 60 degrees after 20
  encounters" (p.353). The 1300 km case of the text (50 degrees in 0.75 years, about 20 encounters) is between the
  0 and 2600 km curves and cannot be read off the figure. Twenty Callisto revolutions are 334 days, or 0.91 years, so
  0.75 years (16.4 revolutions) and "about 20 encounters for the entire sequence" fit only if 50 degrees is reached
  before the end of the sequence.
- **Two loose figures.** "14 x 26 R_J (Hohmann transfer between Ganymede and Callisto)" has an ideal perijove of 15.0 R_J
  from the Table 1 diagonal. And a perijove of 8.4 R_J gives an Io distance of 178,000 km against the printed 186,000 km (the
  paper gives no geometry).
- **Figure 5 positions.** The circular-orbit curve passes at periapsis radius about 26.4, 15, 9.4 and 5.9 R_J for Callisto,
  Ganymede, Europa and Io, equal to the Table 1 diagonal.

## 4. Figures (described)

- **Fig. 1 and Fig. 2 (pp.349-350):** the patched-conic vector geometry and the velocity-sphere picture of cranking and
  pumping, with the angles alpha, beta, epsilon and delta.
- **Fig. 3 (p.350):** inclination (solid) and periapsis radius (dashed) against revolutions and days for one-to-one cranking at
  Callisto, for swingby altitudes of 0, 2600 and 5200 km. Higher flybys crank more slowly but reach the same peak near 59 degrees.
- **Fig. 4 (p.352):** navigation concepts: nominal path, actual path, retargeting maneuver and a multi-nominal path.
- **Fig. 5 (p.353):** periapsis radius (0-30 R_J) against insertion delta-V (0-8 km/s) with orbit-period lines, pump curves
  for each moon, one-to-one crank curves, powered-swingby points and the baseline 1290 m/s line.
- **Fig. 6 (p.354):** a Jupiter-centred plan view of the sample tour (after the paper's ref. 14), with orbits 1, 5, 10, 12, 20,
  40 and 42, the orbits of the three moons, the bow shock and the magnetopause.

## 5. Relevance to this project

- **Not a cycler paper. It is the origin of the vocabulary.** Orbit pumping, orbit cranking and resonance hopping are named
  here. Campagnola and Kawakatsu 2012 (digest of the same day) takes the pump and crank angles from this paper as its ref. 6.
  The tour of Table 2 is a one-way, non-repeating, 42-encounter sequence that the authors call "only one of perhaps a hundred".
- **Against the Ganymede-Callisto and Ganymede-Europa two-working-body rows.**
  - Table 1 gives the ideal odd n-pi semilatus rectum for Ganymede-Callisto (19.1 R_J) and Ganymede-Europa (11.6 R_J). It is the
    first printed statement of the closure condition behind a high-inclination moon-to-moon transfer. The Ganymede-Callisto and
    Ganymede-Europa rows are near-planar and use even and odd arcs in the Russell-Strange model. The table is a ready check on any n-pi transfer
    constructor: it is a function of the two orbit radii only.
  - The Ganymede-Europa sample transfer of sec. 1 (35 days, odd n-pi, inclination near 30 degrees) is a published
    Ganymede-Europa-Ganymede-type fragment, with the return by a 540 degree arc. It is a one-off fragment dated August to
    December 1984 and is not a closed sequence.
  - The statement that an inner and an outer moon can alternate to pump the energy down without lowering periapsis (p.350) is
    the earliest published statement in this corpus of a two-moon working-body idea. It is a remark, with no sequence.
  - "Hohmann transfer between Ganymede and Callisto" (a 15 x 26 R_J orbit by our figures) is mentioned as a Ganymede-orbiter option.
- **Against #973 (two-working-body generator at longer periods).** The paper gives no period, no k, no closed sequence and no
  repeated-moon theorem. It does give the pure-resonance condition (P_sc = n/m P_s and the matching V_n,m) that every
  resonant hop in the generator must satisfy. Nothing here collides with a #973 candidate.
- **Resonances, V_inf and periods printed.** Resonance ratios appear only as the general n:m. The sample tour has periods of
  104, 43, 38, 15, 13, 12, 14, 20, 29 and 57 days, and smear velocities of 5.0 to 7.4 km/s (undefined). V_inf is not printed
  as such at any moon. The approach V_inf at Jupiter is under 6 km/s.
- **Literature-gate label the content supports.** No gate anchor for this paper exists in literature_check.py (grep 2026-10-11).
  If one is added: topology {resonant} and {mga-tour}, primary Jupiter, body set {Io, Europa, Ganymede, Callisto}, multi-moon tour.
  Not {repeated-moon} (no repetition statement). Not {pump-tour} in the gate's Titan sense (V-infinity leveraging): the
  paper's pumping and hopping keep V_inf constant at the moon. A candidate two-moon closed repeating sequence is not covered by it.
- **Risk of misuse.** The sample tour is not a cycler and its dates are 1984-1986 mission dates. The paper's perijove 5.12 R_J
  and bow-shock date 12-25-85 should not be reused without the corrections of sec. 3.

## 6. Citation mining (checked against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Held:** Minovitch 1963, JPL TR 32-464 (ref 1). VanderVeen (ref 3, cited as AAS/AIAA Vail, July 1973, "Triple Planet
  Ballistic Flybys of Mars and Venus"): the corpus holds the JSR 6(4):383 form of 1969 under the same title. The 1973 date and venue are
  the paper's citation; the identity of the two is not checked here.
- **Not in the index and not on the wanted list (new candidates), in this order:**
  - **Beckman, J. C. and Smith, D. B. (1973), "The Jupiter Orbiter Satellite Tour Mission", AAS/AIAA Astrodynamics Conference,
    Vail (ref 4).** This is where the term orbit pumping comes from. Medium priority. Conference paper, likely hard to obtain.
  - **"Mariner-Jupiter-Orbiter Study Report", JPL Report 760-120, July 1975 (ref 5).** It holds the complete 42-encounter sequence
    and the long form of the concepts. Medium priority (JPL report; may be on NTRS).
  - **Uphoff, C. (1975), "The Powered Swingby and Its Application to Jupiter Orbit Mission Design", AAS 75-84 (ref 13).**
    Low to medium priority.
  - **Beckman, Hyde and Rasool (1974), Astronautics and Aeronautics 12, pp.24-35 (ref 14)**, the source of Fig. 6. Low priority.
  - Giacaglia 1968 (SAO Special Report 278, ref 6), Russell 1974 (AAS, orbit determination, ref 7), Curkendall 1973, Duxbury
    1972, Kohlhase 1974, Sherman and Jones 1975 (AAS 75-072), Uphoff 1973 (long-term motion of artificial Jovian satellites),
    Beckman and Miner 1975 (AIAA 75-1141), Friedman and Nunamaker 1975 (AIAA 75-1135): navigation, science or background, not
    added.
- **Niehoff, J. C. (1966), "An Analysis of Gravity-Assisted Trajectories to Solar System Targets", AIAA Paper 66-10 (ref 2):** not
  in the index (the held Niehoff items are 1965 IITRI T-12, 1970 and 1971) and not on the wanted list. Low priority.
- **Wanted-list action:** the row "Clipper and Jovian tour background" has both items received and is moved to the Received
  section.

*Filed as `cyclers_pdf/papers/uphoff-roberts-friedman-1976-orbit-design-concepts-jupiter-orbiter-missions-jsr-13-348-doi-10.2514-3.57096.pdf`. Table transcription: `data/sources/uphoff-roberts-friedman-1976-tables.yaml`.*
