# Digest — Landau, Davis, Karimi (2023), "Trajectory Options for a Uranus Orbiter and Probe" (AAS 23-460 preprint)

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, 15 pages; Table 2 and Figure 6 read
from the rendered pages as well as the text layer). The moon-tour flyby list (Figure 6) is raster image
text, not in the text layer; it was transcribed from the rendered page at 220 dpi. **Purpose:** the
Uranus Orbiter and Probe (UOP) moon-tour design read that the six catalogued Uranian two-moon
quasi-cyclers are provisionally labelled "candidate-novel" pending. This paper is the Ref. [13] of
Landau et al. 2025 (the patched-conic tour algorithm reference). This note is facts only; no novelty
verdict is written here.

## Citation
Damon Landau, Alex Davis, Reza Karimi, "TRAJECTORY OPTIONS FOR A URANUS ORBITER AND PROBE".

- Header line as printed (page 1): "(Preprint) AAS 23-460".
- Affiliations as printed: Landau, Mission Formulation Engineer, Mission Systems Concept Development
  Group, Project Systems Engineering & Formulation Section, Jet Propulsion Laboratory, California
  Institute of Technology (JPL); Davis, Mission Formulation Engineer, same group, JPL; Karimi, Mission
  Design Engineer, Outer Planet Mission Analysis Group, Mission Design & Navigation Section, JPL.
- Venue and date: **not stated.** The paper prints no conference name, city or month. The only year
  printed is in the Acknowledgements: "Copyright (c)2023 California Institute of Technology. U.S
  Government sponsorship acknowledged." (Later papers cite it as the AAS/AIAA Astrodynamics Specialist
  Conference, Big Sky, MT, 2023; this paper does not say so.) Pre-decisional notice printed: "It is
  pre-decisional information for planning and discussion purposes only."
- Filed in the private paper corpus as
  `landau-davis-karimi-2023-trajectory-options-uranus-orbiter-probe-AAS-23-460-preprint.pdf`,
  md5 `c6a37303892c5ef41cf620fc4b2c2e0d`.

## What the paper is
Abstract (quoted): "The scientific value of a Uranus Orbiter and Probe mission is driven by the amount of
data relayed from the probe and the diversity of satellite observations. ... Orbiter requirements include
two distinct flybys of Titania, Oberon, Umbriel, Ariel, and Miranda while minimizing mission Delta-V and
flight time. Automated design methods systematically build a Pareto-optimal trade space of trajectories
that satisfy diverse scientific and technical requirements."

Sections as printed: Introduction; Orbit Insertion and Probe Release Methodology (Probe Release on
Approach; Probe Release from Orbit; Probe-Orbiter Targeting Methodology; UOI and relay results); Moon Tour
Trade Space; Conclusions; Acknowledgements; References. Two tables (Table 1 UOI and relay requirements;
Table 2 satellite tour requirements) and six figures. About two-thirds of the paper is probe-relay and
orbit-insertion geometry (Eqs. 1-42); the moon-tour part is the last section (pages 8 and 11-14).

Context printed: a VVEE gravity-assist cruise launching in 2036 arrives "before Uranus' equinox (in Dec.
2049)"; "this cruise arrival ... occurs at 31-Oct-2049 00:00:00 TDB" with V-infinity (Uranus arrival,
ecliptic frame) = [-3.213, 5.941, 0.111] km/s. (The example tour of Figure 6 prints UOI on 03/02/2049, a
different date; the paper does not reconcile the two. The Figure 6 date format is taken as mm/dd/yyyy
because 06/25/2050 and 04/30/2049 occur.)

## The moon tour method (section "Moon Tour Trade Space")
Goal, quoted: "The goal of the tour is to provide a variety of observations of the five major satellites:
Titania, Oberon, Umbriel, Ariel, and Miranda. Specifically, two flybys of each moon are required (at
least 10 total) with different arrival conditions (different regions illuminated) between flybys."

Construction, quoted: "The tour is built by combining plane-change and v-infinity-leveraging maneuvers at
apoapsis to target the next moon flyby. Then a ballistic gravity assist reduces period (and/or inclination
for out-of-plane arrivals [10]) with a flyby at the minimum allowable altitude. The spacecraft trajectory
is propagated from the flyby departure state to the next apoapsis, and the process repeats for subsequent
flybys."

- Starting points: "The tour search begins at the first apoapsis following relay, where there are many
  (order thousands) of options available from each interplanetary cruise trajectory".
- Two phases: "The design is divided into two phases: initial inclination reduction followed by equatorial
  flybys. Inclination reduction is most efficient using flybys of Titania and/or Oberon, where we target
  the B-plane angle beta to maximally crank the orbit down for a given resonance."
- Flyby geometry: a B-plane frame (Eqs. 43-45) built from the arrival V-infinity and the moon velocity,
  "to always reduce orbit period when cos beta > 0 and always reduce inclination when sin beta > 0";
  departure V-infinity from Eqs. 46-47 using the moon's gravitational parameter and minimum flyby radius.
- Resonance, quoted: "When |S . p| > sin delta then an equatorial orbit cannot be achieved with a single
  flyby. Here we target a M : 1 resonance with the moon, and achieve maximum reduction in inclination by
  flying at the minimum allowable altitude" (Eqs. 48-50, with alpha = 2 pi mu / (M P_ga)^(2/3), P_ga the
  moon period). "The resonant period allows repeated flybys of the moon to gradually reduce inclination
  into the orbital plane, at which point targeting beta = 0 will provide the greatest reduction in
  orbital period for subsequent moon encounters." No specific value of M is stated for any flyby.
- The algorithm, as printed (7 steps): (1) collect spacecraft states at the first apoapsis after the
  probe-orbiter relay; (2) for each apoapsis collect next potential moon-encounter states, "The arrival
  time is sampled at 2-deg resolution in mean anomaly and filtered to within +-1/2 moon-period of the
  next spacecraft periapsis passage to avoid large period-change delta-V at apoapsis"; (3) the apoapsis
  position, moon position and time of flight define a Lambert problem (Ref. [12]); "Only solutions with
  <= 1 revolution are considered to favor tours with shorter flight time"; (4) apoapsis maneuver is
  v_apoD - v_apoA, arrival flyby velocity v_infA = v_A - v_moon; (5) prune to a Pareto set (Ref. [13])
  by number of unique flybys ("defined by combination of moon ID and inbound vs outbound v-infinity"),
  total delta-V including UOI, total flight time since UOI, flyby V-infinity magnitude, inclination
  with respect to Uranus' equatorial plane, and "Proximity of (Uranus-centered) periapsis radius to
  120,000 km. This criteria favors tours that intercept Miranda's orbit, which is more difficult to
  target due to its 4.2 deg inclination"; (6) set flyby periapsis radius and beta to minimise
  post-flyby period and inclination per Eqs. 43-50; (7) Kepler's equation to the next apoapsis.
  "Repeat steps 1-7 until all tours satisfy the design requirement of 10 unique flybys of the moons."
- Dynamical model: patched conics for the probe-relay search ("performed in the patched-conic model");
  for the moon tour the paper says "ballistic gravity assist" and uses Lambert, Kepler propagation and
  impulsive maneuvers ("We model all maneuvers as impulsive"). Moon ephemeris or moon orbit model:
  **not stated.** (Note: the later 2025 paper calls this "the patched-conics algorithm detailed in Ref.
  [13]"; this paper does not use the word "patched-conic" for the moon tour specifically.)

### Table 2, exact transcription
Caption as printed: "Table 2: Design requirements and parameters for satellite tour". The text before it
says "The design constraints that define the trade space for the Uranian satellite tour are provided in
Table 1" (as printed; the tour table is Table 2).

| Parameter | Value |
|---|---|
| min. flyby altitude | 100 km |
| max. flyby v-infinity | 5 km/s |
| min. range to Uranus | 103000 km |
| max. total delta-V including UOI | 2.6 km/s |
| max. flight time since | 500 days |
| max. apo-twist^a delta-V | 400 m/s |
| max. single delta-V | 60 m/s |
| max. orbits between flybys | 2 |

Footnote a, as printed: "an 'apo-twist' [11] reduces inclination and raises periapsis at the apoapsis
preceding the first moon flyby". The row "max. flight time since" is printed with that wording (the
reference event is not on the row; the algorithm step 5 and the Figure 5 x-axis say "since UOI").

Figure 5 (page 13): "Trade space of tours (circles) that provide at least 10 distinct moon flybys";
axes "Days since UOI" (about 300 to 500) against "delta-V including UOI (km/s)" (about 1.8 to 2.6),
four colour classes (probe release from orbit or on approach, ring crossing inside or outside). Text:
"Orbit insertion inside the rings provides the lowest UOI delta-V, while release of probe on approach
reduces flight time"; probe release from orbit outside the rings gives "the lowest total tour delta-V at
expense of longer flight time"; release on approach "the lowest flight time options, at expense of
increased delta-V". The number of tours plotted is not stated.

## The exact tour list (Figure 6; the only tour listing in the paper)
Figure 6 caption as printed: "Example tour that provides two flybys (in- and out-bound) each of Oberon,
Titania, Umbriel, Ariel, and Miranda (top) following inclination-reduction phase (bottom)." Text: "An
example tour that begins with UOI outside the rings and probe released from orbit is shown in Fig. 6. The
bottom panel highlights the inclination-reduction sequence driven by flybys of Titania and Oberon."
There is no table of the tour; the labels are in the figure image. Transcribed from the rendered page
(event number, as printed):

| Event | Printed label |
|---|---|
| 1 | UOI, 03/02/2049, delta-V = 1840 m/s |
| 2 | Target, 04/16/2049, delta-V = 182 m/s |
| 3 | Divert, 04/30/2049, delta-V = 22 m/s |
| 4 | End Relay, 05/30/2049 |
| 5 | ApoTwist, 07/14/2049, delta-V = 175 m/s |
| 6 | Oberon, 08/31/2049, V-inf = 3.688 km/s, Alt = 100 km |
| 7 | delta-V = 23 m/s |
| 8 | Titania, 11/07/2049, V-inf = 3.933 km/s, Alt = 100 km |
| 9 | delta-V = 59 m/s |
| 10 | Titania, 12/24/2049, V-inf = 3.765 km/s, Alt = 100 km |
| 11 | delta-V = 15 m/s |
| 12 | Oberon, 01/30/2050, V-inf = 3.349 km/s, Alt = 100 km |
| 13 | delta-V = 19 m/s |
| 14 | Ariel, 03/04/2050, V-inf = 3.990 km/s, Alt = 100 km |
| 15 | delta-V = 15 m/s |
| 16 | Ariel, 03/30/2050, V-inf = 3.897 km/s, Alt = 100 km |
| 17 | delta-V = 8 m/s |
| 18 | Umbriel, 04/23/2050, V-inf = 4.095 km/s, Alt = 100 km |
| 19 | delta-V = 5 m/s |
| 20 | Umbriel, 05/16/2050, V-inf = 4.070 km/s, Alt = 100 km |
| 21 | delta-V = 45 m/s |
| 22 | Miranda, 06/05/2050, V-inf = 2.802 km/s, Alt = 100 km |
| 23 | delta-V = 60 m/s |
| 24 | Miranda, 06/25/2050, V-inf = 3.056 km/s (no Alt line printed for this flyby) |

The ten moon flybys are events 6, 8, 10, 12, 14, 16, 18, 20, 22, 24: Oberon, Titania, Titania, Oberon,
Ariel, Ariel, Umbriel, Umbriel, Miranda, Miranda (two each, as the caption says). Intervals between
consecutive flybys, computed by me from the printed dates (days): Oberon to Titania 68; Titania to
Titania 47; Titania to Oberon 37; Oberon to Ariel 33; Ariel to Ariel 26; Ariel to Umbriel 24; Umbriel to
Umbriel 23; Umbriel to Miranda 20; Miranda to Miranda 20. First flyby to last flyby 298 days. Days from
UOI (03/02/2049) to the first flyby 182, to the last flyby 480. The sum of the printed delta-V values
(UOI plus all listed maneuvers) is 2468 m/s by my addition; the paper prints no tour total for this
example. The paper does not state how many orbits lie between flybys, the resonance of any flyby, or the
turn angle at any flyby.

## Answers to the specific questions
**(3) Alternation of two moons.** No two-moon alternation A, B, A, B, and no A, B, A with a single B,
occurs in the example tour. The sequence is O, T, T, O, A, A, U, U, M, M. The nearest pattern is
Oberon, Titania, Titania, Oberon (events 6, 8, 10, 12): intervals 68, 47 and 37 days; V-infinity 3.688,
3.933, 3.765, 3.349 km/s. Oberon (12) is followed by Ariel (33 days), and Ariel is not followed by
Oberon again. The paper does not describe it as an alternation.

**Repeated flybys of the same moon (consecutive pairs).**
- Titania, events 8 and 10: 47 days; V-inf 3.933 and 3.765 km/s.
- Ariel, events 14 and 16: 26 days; V-inf 3.990 and 3.897 km/s.
- Umbriel, events 18 and 20: 23 days; V-inf 4.095 and 4.070 km/s.
- Miranda, events 22 and 24: 20 days; V-inf 2.802 and 3.056 km/s.
- Oberon: events 6 and 12 are the two Oberon flybys (not consecutive; Titania flybys lie between), 152
  days apart (68 + 47 + 37); V-inf 3.688 and 3.349 km/s.

The resonance for any of these pairs is **not stated**. The only resonance statement is the
generic M : 1 for the inclination-reduction flybys (Titania and/or Oberon) in the section above, with no
value of M given.

**(4) Words.** Text search of the full text layer (figure-image labels not searchable): "cycler" 0,
"cycling" 0, "periodic" 0, "free-return" 0. "repeat" appears 3 times: "the process repeats for subsequent
flybys" (describing the tour-building loop), "The resonant period allows repeated flybys of the moon to
gradually reduce inclination into the orbital plane", and "Repeat steps 1-7 until all tours satisfy the
design requirement of 10 unique flybys of the moons". "resonance" or "resonant" appears 3 times, all in
the inclination-reduction passage quoted above (for a given resonance; target a M : 1 resonance; the
resonant period allows repeated flybys). "period" appears about a dozen times, about capture-orbit
period, relay period, orbit period reduction and moon period. No trajectory is described as returning to
the same two moons repeatedly or indefinitely: every tour is built until "10 unique flybys" are
reached, with a 500-day flight-time limit after UOI. The word "petal" is not used.

**(5) V-infinity at moon flybys.** Lowest printed: 2.802 km/s (Miranda, event 22). Highest printed: 4.095
km/s (Umbriel, event 18). Per moon: Oberon 3.349 to 3.688; Titania 3.765 to 3.933; Ariel 3.897 to 3.990;
Umbriel 4.070 to 4.095; Miranda 2.802 to 3.056. Design cap in Table 2: max flyby V-infinity 5 km/s. (The
Uranus arrival V-infinity vector magnitude from the printed components is about 6.76 km/s by my
arithmetic; that is the cruise arrival, not a moon flyby.)

## What it does NOT contain
- No table of tour flybys; only the single Figure 6 example (10 moon flybys), with dates, V-infinity and
  100 km altitude.
- No cycler, periodic, free-return or repeating two-moon analysis; no statement of any trajectory that
  returns to the same moons indefinitely.
- No moon-pair resonance values, no turn angles, no per-leg time-of-flight statement in text, no
  statement of how many revolutions lie between flybys beyond the limit "max. orbits between flybys 2".
- No count of tours found, no Pareto-set size, no Lambert-solution statistics.
- No moon ephemeris or moon-orbit model for the tour.
- No venue or date for this paper itself (see Citation).
- No flyby of Oberon, Umbriel or Ariel in an A, B, A pattern.

## Relevant bibliography (as printed in this paper)
- [3] M. Hofstadter, A. Simon, K. Reh, and J. Elliott, "Ice Giants Pre-Decadal Survey Mission Study
  Report," Tech. Rep. JPL D-100520, NASA, June 2017.
- [4] M. Hofstadter, A. Simon, S. Atreya, D. Banfield, J. J. Fortney, A. Hayes, M. Hedman, G. Hospodarsky,
  K. Mandt, A. Masters, et al., "Uranus and Neptune missions: A study in advance of the next Planetary
  Science Decadal Survey," Planetary and Space Science, Vol. 177, 2019, p. 104680.
- [5] J. McAdams, C. Scott, Y. Guo, J. Dankanich, and R. Russell, "Conceptual mission design of a polar
  Uranus orbiter and satellite tour," Spaceflight Mechanics, Vol. 140, 2011. (Cited in the text for
  inclined tours where "flyby opportunities are more constrained".)
- [6] D. Landau, S. Campagnola, and E. Pellegrini, "Star Searches for Patched-Conic Trajectories,"
  Journal of the Astronautical Sciences, 2022, DOI :10.1007/s40295-022-00350-y. (Used for the
  interplanetary cruise search.)
- [10] N. Strage, D. Landau, and J. Longuski, "Design of initial inclination reduction sequence for
  Uranian gravity-assist tours," AAS/AIAA Spaceflight Mechanics Meeting, American Astronautical
  Society, 2014, pp. 1469-1485. AAS 13-801.
- [11] R. W. Luidens and B. A. Miller, "Efficient planetary parking orbits with examples for Mars," Tech.
  Rep. NASA TN D-3220, National Aeronautics and Space Administration, January 1966. (The apo-twist.)
- [12] R. Gooding, "A procedure for the solution of Lambert's orbital boundary-value problem," Celestial
  Mechanics and Dynamical Astronomy, Vol. 48, No. 2, 1990, pp. 145-165.
- [13] Y. Cao, "Pareto Set," https://www.mathworks.com/matlabcentral/fileexchange/15181-pareto-set, 2020.
  [Online; accessed October 31, 2020].
- [14] C.-Y. Wu and R. P. Russell, "Reachable Set of Low-Delta-v Trajectories Following a Gravity-Assist
  Flyby," Journal of Spacecraft and Rockets, Vol. 60, No. 2, 2023, pp. 616-633, 10.2514/1.A35464.
  (Cited for the apoapsis position and velocity from the eccentricity vector.)
- Also about mission design, listed for completeness: [1] National Academies of Sciences, Engineering,
  and Medicine, Origins, Worlds, and Life: A Decadal Strategy for Planetary Science and Astrobiology
  2023-2032, 2022, 10.17226/26522; [2] A. Simon, F. Nimmo, and R. Anderson, "Uranus Orbiter and Probe,"
  tech. rep., NASA, June 2021.
- Omitted as not about tour design: [7] and [8] Odell and Gooding (Kepler's equation), [9] Battin.

## Differences from Landau et al. 2025
(Based on the sibling digest, `2026-10-03-digest-landau-2025-uranus-cruise-tour-design.md`.)

In this 2023 paper and not in the 2025 digest:
- The tour algorithm itself: 7-step apoapsis-to-moon Lambert construction (at most 1 revolution, 2-deg
  mean-anomaly sampling, +-half moon-period filter), Pareto pruning axes, the 120,000 km periapsis
  proximity criterion for reaching Miranda, the inclination-reduction phase using Titania and Oberon at
  an M : 1 resonance, apo-twist.
- Tour constraints in Table 2: min flyby altitude 100 km, max flyby V-infinity 5 km/s, min range to
  Uranus 103000 km, max total delta-V 2.6 km/s, max 500 days flight time since UOI, max apo-twist 400
  m/s, max single maneuver 60 m/s, max 2 orbits between flybys.
- One dated 10-flyby example tour (Figure 6) with per-flyby V-infinity (2.802 to 4.095 km/s), 100 km
  altitude and the inter-flyby maneuvers (5 to 60 m/s) and a UOI-to-end timeline (UOI 03/02/2049 to last
  flyby 06/25/2050).
- Probe-relay and UOI design equations (about two-thirds of the paper).

In the 2025 paper and not in this paper:
- The 26-flyby example moon tour (Table 3) with latitude, longitude, solar longitude, solar phase and
  V-infinity per flyby, at 50 km altitude and V-infinity 2.6 to 4.3 km/s, including the stretch
  Umbriel, Oberon, Umbriel, Oberon (rows 23-26) and the runs of three to four Titania, Ariel flybys; this
  paper has no alternation of that kind.
- Moon-tour science drivers (Table 1 of the 2025 paper: distinct encounters, three flybys per moon),
  the magneto-tour, the Juno-like close passes, and the cruise, mass and mission-duration trade space.
- The 2025 paper's statement of the algorithm is two sentences deferring to this paper (its Ref. [13]);
  it names a 50 km minimum altitude where this paper's Table 2 and Figure 6 give 100 km.
- Different mission epoch: 2025 paper first Titania flyby "6.2 Sep 2045"; this paper's example tour flies
  2049-2050.

## Coordinator's verdict for the six Uranian quasi_cycler rows (added 2026-10-03, `#869`)

A judgement by the coordinating session, which read the tour sections and Figure 6 independently
of the digest above (the two transcriptions of Figure 6 agree); everything above is the factual
digest.

**Nothing here bears on the six rows.**

- The paper has no tour table. Its one example tour (Figure 6) is Oberon, Titania, Titania, Oberon,
  Ariel, Ariel, Umbriel, Umbriel, Miranda, Miranda: no two moons alternate A, B, A. This is further
  from the catalogued geometry than the later papers' stretches (Umbriel-Oberon in Landau et al.
  2025 Table 3; Ariel-Oberon in AAS 25-668 Table 8).
- Every flyby is preceded by an apoapsis manoeuvre (5-60 m/s in the example), so the tour is
  powered between encounters; the six rows are ballistic closures.
- Flyby V-infinity is 2.8-4.1 km/s (design cap 5 km/s) against 0.89-2.16 km/s for the six rows.
- The only repetition in the method is an M:1 same-moon resonance used to bring inclination down,
  and tours stop at ten unique flybys within 500 days. Cycler, periodic and free-return do not
  occur.

**Consequence:** none for the catalogue. The rows' wording ("no published repeating or periodic
two-moon trajectory at Uranus, and none in this V-infinity regime") stands, and their
`corroborating_sources` are unchanged. With this paper the Uranus Orbiter and Probe era tour
papers (this one, Landau et al. 2025, Ellison et al. AAS 25-668, Simon et al. 2026) have all been
read in full. Two OLDER Uranian tour papers that they cite are still not held: McAdams, Scott,
Guo, Dankanich & Russell, AAS 11-188 (this paper's Ref. [5]) and Strange, Landau & Longuski,
AAS 13-801 (Ref. [10], the inclination-reduction sequence). This paper is not added as a separate literature-gate
anchor: the Landau et al. 2025 anchor (same group, same algorithm, and the only one of the two
with a tour table) already covers the class as `mga-tour`.
