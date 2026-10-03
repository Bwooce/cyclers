# Digest — McAdams, Scott, Guo, Dankanich, Russell (2011), "Conceptual Mission Design of a Polar Uranus Orbiter and Satellite Tour" (AAS 11-188)

**Digested:** 2026-10-03 (text-layer PDF, 14 pages, read in full; page 1, page 8 and page 9 (Table 4, Figure 7)
viewed as rendered images and Table 4 checked against the text layer). **Purpose:** the Uranian moon-tour
literature check for the six catalogued Uranian two-moon quasi-cyclers. This paper is Ref. [5] of Landau,
Davis & Karimi 2023 (AAS 23-460) and Ref. [4] of Strange, Landau & Longuski 2013 (AAS 13-801). This note
is facts only; no novelty verdict is written here. Statements marked "by my arithmetic" are computed from
printed values; everything else is quoted or paraphrased from the paper.

## Citation
James McAdams, Christopher Scott, Yanping Guo, John Dankanich and Ryan Russell, "CONCEPTUAL MISSION DESIGN
OF A POLAR URANUS ORBITER AND SATELLITE TOUR".

- Header line as printed (page 1): "AAS 11-188". Printed page numbers run 1257 to 1270.
- Affiliations as printed (page-1 footnotes): McAdams, Mission Design Lead Engineer, The Johns Hopkins
  University Applied Physics Laboratory (JHU/APL), 11100 Johns Hopkins Rd, Laurel MD, 20723; Scott,
  Mission Design Analyst, JHU/APL; Guo, Mission Design Lead Engineer, JHU/APL; Dankanich, Mission Analyst,
  Gray Research, Inc, 21000 Brookpark Rd. M/S 77-4, Cleveland, OH, 44135; Russell, Assistant Professor,
  Georgia Institute of Technology, Atlanta, Georgia, 30332. (Dankanich and Russell carry the same footnote
  marker on the rendered page; the two footnotes are listed in that order.)
- Venue, place, date, volume: **not stated.** The paper prints no conference name, city, month or year of
  publication. Years that appear refer to the mission (launch 2020-2023). Other papers cite it as
  "Spaceflight Mechanics, Vol. 140, 2011" (Landau et al. 2023) and "Space Flight Mechanics Conference, Feb.
  2011. AAS Paper 11-188" (Strange et al. 2013); this paper does not say so itself.
- Sponsorship printed (Acknowledgments): "NASA sponsorship for the Ice Giant Orbiter/Probe Mission Decadal
  study under contract NNN06AA01C, task NNN08AA03T with The Johns Hopkins University Applied Physics
  Laboratory (JHU/APL), where Helmut Seifert provided management oversight."
- Filed in the private paper corpus as
  `mcadams-scott-guo-dankanich-russell-2011-conceptual-mission-design-polar-uranus-orbiter-satellite-tour-AAS-11-188.pdf`,
  md5 `3b5f664fbefb602f0cb5c80afa8ec0d4`.

## What the paper is
Abstract (quoted): "In response to NASA's planetary science decadal survey, this paper outlines the
conceptual mission design of a Uranus orbiter. In the baseline design the spacecraft launches during a
21-day launch period in 2020, followed by a 13-year cruise with solar electric propulsion and a single
Earth flyby. Repeatable launch opportunities are available from 2021-2023. An atmospheric probe is
released 29 days prior to Uranus orbit insertion. After completion of the probe descent phase the
spacecraft inserts into a highly inclined elliptical orbit for 431 days, followed by the satellite tour
with targeted flybys of five satellites."

A whole-mission concept study "managed by The Johns Hopkins University Applied Physics Laboratory, with key
direction from a science steering committee". Sections as printed: Introduction; Journey to Uranus
(Propulsion Trades, Launch Opportunities, Sensitivity to Earth Flyby Altitude and Cruise Phase Duration);
Arrival at Uranus (Probe Release and Atmospheric Entry, Orbit Insertion); Orbit Phase and Satellite Tour
(Primary Science Orbit, Satellite Tour Options); Conclusion; References. Five tables and fourteen figures.
Only the last subsection, "Satellite Tour Options" (printed pages 1264-1269), concerns the moon tour; the
rest is low-thrust cruise, probe, orbit insertion and primary-science-orbit design.

Mission context printed: a 13-year SEP cruise (20 kW, two NEXT thrusters, Atlas 551, one Earth flyby at
1,000 km) arriving "mid-2033"; probe separation 29 days before entry; UOI "June 28, 2033", 1661 m/s,
into a "1.3 RU (33,425-km altitude) periapsis by 21.0-day Uranus orbit with inclination of 97.7 deg"; a
"20.5-orbit, 431-day baseline orbit phase" after UOI "ends with the start of the optional satellite tour".
The tour is described as "optional" and "secondary" in the text.

## The tour design method (subsection "Satellite Tour Options")
- Model, quoted: "The preliminary tour design is based on the zero-radius sphere of influence patched-conic
  model with satellite ephemeris locations provided by 'ura083.bsp'." Later: "Integrating the Uranus
  centric orbit using the patched-conic satellite ephemeris model".
- Method, quoted: "The tour was designed in part using the graphical methods based on V globe maps that
  reduce feasible options for post flyby orbits onto maps with contours of desired quantities." Figure 7 is
  a globe map of the final (Oberon) flyby, titled in the rendered image "v-infinity globe map, ring
  crossing countours (km), v-infinity = 4.6438 (km/s)", with "line is reachable with min alt flyby, x
  before flyby, o after flyby" and a 2:1 resonance band marked.
- Flyby altitude limit: "the circle represents the possible V vector locations after a ballistic 50 km
  flyby. Locations outside the circle require flybys with altitudes less than 50 km and are therefore
  unreachable." Minimum targeted flyby altitude is therefore 50 km. The paper prints no per-flyby altitude.
- Moon-to-moon transfers (manoeuvred, not ballistic): "Maneuvers at flybys occur at the sphere of
  influence. Large maneuvers target the next moon." "Consecutive maneuvers indicate a two-impulse sequence
  that guarantees phasing in about two spacecraft orbits. The first maneuver is targets the ring plane
  crossing, while the second maneuver is at the ring plane and adjusts the period to ensure proper
  phasing. The two-impulse solutions are less favorable for Titania and Oberon because their periods are
  much larger than the other moons ... Therefore, a single-impulse solution is used at Titania and
  Oberon." Between repeat flybys of the same moon: "In most cases, nearby satellite:spacecraft resonances
  were targeted for orbits between repeat flybys of the same satellite (16:1 for Miranda, 10:1 for Ariel,
  6:1 for Umbriel, 3:1 for Titania, and 2:1 for Oberon). Many small maneuvers not shown in the satellite
  tour summary (see Table 4) are required in the design because precise resonances are absent in a full
  ephemeris model."
- Resonance and free-return passage, quoted in full: "Successive encounters with a single body are achieved
  by targeting the resonant bands while encounters with new moons are enabled by targeting intersections of
  the 50 km flyby circle with the ring plane crossing corresponding to the targeted moons' orbital radius.
  (Note that a later maneuver and multiple revolutions then lead to correct phasing.) Successive resonant
  free-returns can be used to reduce delta-v requirements to reach the next moon, although the efficiency
  is low due to the low mass of all the moons. Instead, in this tour we favored short flight times and used
  only one resonant return for each moon in order to slightly reduce delta-v but more importantly achieve a
  second flyby to enhance the science return."
- Constraint structure, quoted: "The near-polar inclined tour is highly constrained because 1) very high
  excess velocity limits flyby capabilities, and 2) the only potential for moon encounters is at the node
  crossings of the spacecraft." "The tour is not optimized end to end. Each leg is a global minimum for its
  associated single-impulse maneuver, or an approximate global minimum for the two-impulse maneuver."
- Maneuver margin: "An additional 5 m/s/flyby is allocated (see Table 5) to accommodate for navigation
  errors and model fidelity errors (based on numbers from Cassini tour design)."
- Flyby speeds and turning, quoted: "With such highly elliptical orbits, spacecraft-satellite encounter
  velocities vary from 10.9 km/s at Miranda to 8.8 km/s at Ariel to 7.3 km/s at Umbriel to 5.6 km/s at
  Titania to 4.7 km/s at Oberon. For Oberon and Titania, these encounter speeds allow for only ~1 deg of
  turning due to a hyperbolic flyby while the other moons provide ~0.25 deg or less turning." (The paper
  says "encounter velocities" and does not use the symbol V-infinity for these.)
- Tour size, quoted: "Lasting 424 days, the baseline tour has ten targeted flybys passing by Miranda (body
  ID 705), Ariel (ID 701), Umbriel (ID 702), Titania (ID 703), and Oberon (ID 704) twice each and four
  additional close untargeted flybys with Umbriel." "Out of 13 such non-targeted encounters with Miranda,
  Ariel, and Umbriel, the four closest encounters are with Umbriel and occur between 1,400 and 3,700 km."
  "The total delta-V for the tour is 619 m/s."
- Alternative tour, quoted: "Results of an alternative satellite tour with 5 Miranda, 1 Ariel, 7 Umbriel,
  4 Titania, and 5 Oberon flybys are not presented here because the tour flight time was more than 50%
  longer than the 424-day duration for the baseline satellite tour. An alternative, Galileo-style Uranus
  satellite tour was offered by Heaton and Longuski [9]."

## The exact tour list (Table 4, the only tour listing; Table 5 is the delta-V budget)
Caption as printed: "Table 4. Satellite flyby and maneuver dates for the baseline satellite tour." Columns:
Julian Date, Body ID, Event, delta-v (m/s). Transcribed from the rendered page 1265 and checked against the
text layer. Blank cells are blank in the print. Rows with a Julian Date and no body are marked "-".
The Table 4 rows from 2464547.64 to 2464557.96 are set in italics in the print (the text says "the
italicized labels in Table 4 indicate a very long flight time after the maneuver that achieves the Titania
transfer").

| Julian Date | Body ID | Event | delta-v (m/s) |
|---|---|---|---|
| 2464214.50 | - | ref_orbit | 0 |
| 2464228.63 | - | maneuver | 88.4 |
| 2464239.46 | - | maneuver | 8.239 |
| 2464261.51 | 705 | flyby(16:1) | 4.979 |
| 2464271.33 | - | maneuver | 1.401 |
| 2464284.13 | (blank) | (blank) | (blank) |
| 2464295.36 | - | maneuver | 79.107 |
| 2464307.06 | - | maneuver | 8.932 |
| 2464330.76 | 701 | flyby(10:1) | 0 |
| (blank) | - | maneuver | 0.984 |
| 2464355.96 | 701 | flyby(10:1) | 0 |
| 2464368.42 | - | maneuver | 99.205 |
| (blank) | - | maneuver | 2.687 |
| 2464407.76 | 702 | flyby(6:1) | 0 |
| 2464419.13 | - | maneuver | 0.956 |
| 2464432.63 | 702 | flyby(6:1) | 0 |
| 2464547.64 (italic) | - | maneuver (italic) | 181.205 (italic) |
| 2464557.96 (italic) | 703 | flyby(3:1) (italic) | 0 (italic) |
| 2464569.19 | - | maneuver | 2.094 |
| 2464584.08 | 703 | flyby(3:1) | 0 |
| 2464598.08 | - | maneuver | 139.335 |
| 2464611.57 | 704 | flyby(2:1) | 0 |
| 2464622.64 | - | maneuver | 1.548 |

Observations on the table as printed (by my arithmetic from the printed values):
- The table has no altitude column and no V-infinity column. The "resonance" is printed in the Event cell.
- The table lists one Miranda flyby (705) and one Oberon flyby (704), while the text says each of the five
  moons is targeted twice. The row at 2464284.13 has a date and nothing else; the paper does not say what it
  is. The 16:1 label on the Miranda row is the resonance quoted in the text "between repeat flybys".
- Calendar dates of the Julian Dates, by my conversion (the paper prints no calendar dates for the tour):
  2464214.50 = 2034-09-09; first listed flyby (Miranda) 2034-10-26; Ariel 2035-01-03 and 2035-01-28;
  Umbriel 2035-03-21 and 2035-04-15; Titania 2035-08-18 and 2035-09-13; Oberon 2035-10-11; last row
  2035-10-22.
- Printed delta-v values sum to 619.072 m/s, matching the stated 619 m/s.
- Intervals between listed flybys, from the printed Julian Dates, in days: Miranda to Ariel 69.25; Ariel to
  Ariel 25.20; Ariel to Umbriel 51.80; Umbriel to Umbriel 24.87; Umbriel to Titania 125.33; Titania to
  Titania 26.12; Titania to Oberon 27.49. First listed flyby (Miranda) to last listed flyby (Oberon) 350.06
  days; ref_orbit row to last row 408.14 days (the text says the tour lasts 424 days; the paper does not
  reconcile the two numbers).

### Table 5, mission delta-V budget (as printed; satellite-tour rows and total)
| Event | delta-V (m/s) | Propulsion type (1 = hybrid, 2 = biprop) | Comment |
|---|---|---|---|
| Miranda Targeting | 103.0 | 2 | |
| Miranda Statistical | 10.0 | 1 | 5 m/s per flyby (2x) |
| Ariel Targeting | 89.0 | 2 | |
| Ariel Statistical | 10.0 | 1 | 5 m/s per flyby (2x) |
| Umbriel Targeting | 102.8 | 2 | |
| Umbriel Statistical | 10.0 | 1 | 5 m/s per flyby (2x) |
| Titania Targeting | 183.3 | 2 | |
| Titania Statistical | 10.0 | 1 | 5 m/s per flyby (2x) |
| Oberon Targeting | 140.9 | 2 | |
| Oberon Statistical | 10.0 | 1 | 5 m/s per flyby (2x) |
| Total (whole mission) | 2501.1 | | |

(Other Table 5 rows: interplanetary statistical 30.0, orbit deflection 30.0, UOI B-plane targeting 10.0,
orbit insertion 1661.0, cleanup 25.0, orbit maintenance 20.0 "~1 m/s per orbit", periapsis reduction 56.0.)
The targeting rows for the five moons sum to 619.0 m/s by my arithmetic, matching the stated tour 619 m/s;
the statistical rows are separate (50 m/s).

## Answers to the specific questions
**(4) Alternation of two different moons.** None. The listed moon order is Miranda, Ariel, Ariel, Umbriel,
Umbriel, Titania, Titania, Oberon (Table 4), or Miranda, Miranda, Ariel, Ariel, Umbriel, Umbriel, Titania,
Titania, Oberon, Oberon by the text ("twice each"). Each moon is visited, then the next moon in order of
increasing orbital radius; no moon recurs after a different moon. The four extra untargeted Umbriel
passes are described only as "four additional close untargeted flybys with Umbriel"; no dates, order or
distances for them are given beyond "between 1,400 and 3,700 km" for the four closest of 13 non-targeted
encounters.

**Repeated flybys of the same moon and the resonance stated** (consecutive pairs in Table 4):
- Miranda 16:1 (second flyby not itemised in Table 4).
- Ariel 10:1: JD 2464330.76 and 2464355.96, 25.20 days.
- Umbriel 6:1: JD 2464407.76 and 2464432.63, 24.87 days.
- Titania 3:1: JD 2464557.96 and 2464584.08, 26.12 days.
- Oberon 2:1 (second flyby not itemised in Table 4).
The resonance is stated as the integer ratio of satellite to spacecraft revolutions as printed (for
example "16:1 for Miranda"); the paper does not define the ratio direction in words.

**(5) Word search** (full text layer, lines joined, ligatures normalised, case-insensitive; figure-image
labels not searched): "cycler" 0; "cycling" 0; "periodic" 0; "repeating" 0; "petal" 0; "free-return" 2;
"repeat" 4; "resonant" 3; "resonance" 4 ("resonances" counted); "free return" (no hyphen) 0. The sentences:
- "free-return" (1): "Successive resonant free-returns can be used to reduce delta-v requirements to reach
  the next moon, although the efficiency is low due to the low mass of all the moons. Instead, in this tour
  we favored short flight times and used only one resonant return for each moon ..." (2): the title of
  Ref. 8, Russell and Ocampo, "Geometric Analysis of Free-Return Trajectories Following a Gravity-Assisted
  Flyby", Journal of Spacecraft and Rockets, 2005.
- "repeat" (4): "Repeatable launch opportunities are available from 2021-2023." (abstract); "launch
  opportunity repeatability" (introduction); "20 primary science orbits providing repeated coverage of the
  same region at ~2-hour intervals for cloud tracking" (an arrival requirement, about Uranus orbits, not
  moons); "resonances were targeted for orbits between repeat flybys of the same satellite".
- "resonant" (3): "targeting the resonant bands"; "Successive resonant free-returns"; "only one resonant
  return for each moon".
- "resonance" (4): "nearby satellite:spacecraft resonances were targeted"; "precise resonances are absent
  in a full ephemeris model"; "return to Oberon on a 2:1 resonance"; "the 2:1 resonance band".
No trajectory is described as returning to the same two moons repeatedly or indefinitely. The only
repetition is one resonant return to the same moon (two flybys per moon), and the tour ends at Oberon after
424 days.

**(6) V-infinity, duration, delta-v.** The paper gives no per-flyby V-infinity. It gives per-moon
"spacecraft-satellite encounter velocities": Miranda 10.9, Ariel 8.8, Umbriel 7.3, Titania 5.6, Oberon 4.7
km/s (lowest 4.7 at Oberon, highest 10.9 at Miranda). Figure 7 (final Oberon flyby globe map) is titled
"v-infinity = 4.6438 (km/s)". Tour duration 424 days; tour delta-v 619 m/s (targeting; plus 50 m/s
statistical allowance in Table 5); whole-mission delta-v 2501.1 m/s (Table 5). Moon-flyby altitude: 50 km
minimum for the ballistic flyby circle; no per-flyby altitude given.

## What it does NOT contain
- No per-flyby V-infinity, altitude, turn angle, B-plane angle or calendar date for the tour; Table 4 has
  only Julian Dates, body ID, event with resonance, and delta-v.
- No itemisation in Table 4 of the second Miranda flyby, the second Oberon flyby, or the four untargeted
  Umbriel flybys, and no explanation of the row with a date only.
- No two-moon alternation, no cycler or periodic trajectory, no trajectory returning to the same moons
  indefinitely; the words "cycler", "cycling", "periodic" and "petal" do not occur.
- No ballistic moon-to-moon transfers: every transfer to a new moon uses large manoeuvres ("Large
  maneuvers target the next moon").
- No equatorial or Galileo-style tour (the paper points to Heaton and Longuski [9] for that); the tour is
  near-polar with node-crossing encounters only.
- No statement of venue, place or publication date for the paper itself.
- Nothing on Hohmann-type or V-infinity-leveraging sequences, Tisserand graphs, or multi-body dynamics.

## Check of the literature-gate anchor sentence against this paper
Anchor text: "Sims, Finlayson, Rinderle, Vavrina & Kawalkowski et al., 'Conceptual mission design of a polar
Uranus orbiter and satellite tour' (2014). Baseline 424-day, 619 m/s tour with two targeted flybys of each
major moon. One-shot insertion tour, NOT cycler."

| Element | Verdict | Paper text |
|---|---|---|
| Authors Sims, Finlayson, Rinderle, Vavrina, Kawalkowski | Contradicted | The authors printed on page 1 are McAdams, Scott, Guo, Dankanich and Russell. The five names are the authors of Reference 1 only: "Sims, J. A., Finlayson, P. A., Rinderle, E. A., Vavrina, M. A., and Kawalkowski, T. D., 'Implementation of a Low-Thrust Trajectory Optimization Algorithm for Preliminary Design,' AIAA 2006-6746, AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, August 21-24, 2006." It is cited as "Mission Analysis Low-Thrust Optimization (MALTO) tool.1" in the propulsion trades. The paper prints this surname as "Kawalkowski"; the spelling "Kowalkowski" does not occur in it. |
| Year 2014 | Not stated | The paper prints no publication year. The year 2014 appears nowhere. (Strange et al. 2013 prints this paper as "Feb. 2011. AAS Paper 11-188"; Landau et al. 2023 prints "Vol. 140, 2011".) |
| Title | Confirmed | "CONCEPTUAL MISSION DESIGN OF A POLAR URANUS ORBITER AND SATELLITE TOUR" (the same words). |
| 424-day | Confirmed | "Lasting 424 days, the baseline tour has ten targeted flybys ..." and "a baseline 424-day, 619-m/s satellite tour design" (Conclusion). |
| 619 m/s | Confirmed | "The total delta-V for the tour is 619 m/s."; Table 4 values sum to 619.07 m/s by my arithmetic. |
| Two targeted flybys of each major moon | Confirmed | "ten targeted flybys passing by Miranda ..., Ariel ..., Umbriel ..., Titania ..., and Oberon ... twice each"; Conclusion: "two targeted flybys each of the five largest moons of Uranus". (Table 4 itself itemises only one Miranda and one Oberon flyby; see above. The paper also lists four additional untargeted Umbriel flybys.) |
| One-shot insertion tour | Not stated | The phrase does not occur. The paper describes a baseline tour that follows a 431-day primary science orbit after UOI, and calls the tour "optional" and "secondary". |
| Not a cycler | Not stated | The word "cycler" does not occur. What the paper does say is in section (5): a 424-day tour with one resonant return per moon, with "Successive resonant free-returns" mentioned and not used. |

Names Sims, Finlayson, Rinderle, Vavrina and Kawalkowski: each appears exactly once in this paper, all in
Reference 1, as co-authors of the MALTO low-thrust optimisation paper cited for the SEP cruise trades. None
appears as an author of this paper, in the acknowledgments, or in the body text. The spelling
"Kowalkowski" does not appear (the paper has "Kawalkowski").

## Relevant bibliography (as printed in this paper)
- [1] J. A. Sims, P. A. Finlayson, E. A. Rinderle, M. A. Vavrina, and T. D. Kawalkowski, "Implementation of a
  Low-Thrust Trajectory Optimization Algorithm for Preliminary Design," AIAA 2006-6746, AIAA/AAS
  Astrodynamics Specialist Conference, Keystone, CO, August 21-24, 2006. (MALTO; cruise only.)
- [2] D. Landau, T. Lam, and N. Strange, "Broad Search and Optimization of Solar Electric Propulsion
  Trajectories to Uranus and Neptune," AAS 09-428. (Cruise.)
- [7] N. J. Strange, R. P. Russell, B. Buffington, "Mapping the V-infinity Globe," Paper AAS 07-277,
  AAS/AIAA Astrodynamics Specialist Conference and Exhibit, Mackinac Island, MI, Aug 2007.
- [8] R. P. Russell, C. A. Ocampo, "Geometric Analysis of Free-Return Trajectories Following a
  Gravity-Assisted Flyby," Journal of Spacecraft and Rockets, Vol. 42, No. 1, pp. 138-151, 2005.
- [9] A. Heaton and J. Longuski, "Feasibility of a Galileo-Style Tour of the Uranian Satellites," Journal of
  Spacecraft and Rockets, Vol. 40, No. 4, July-August 2003, pp. 591-596.
- Where [7] and [8] are cited in the body is not visible in the text layer (superscripts dropped); citation
  numbers 1 to 6 and 9 are visible in the text layer. The body describes "graphical
  methods based on V globe maps" and "free-returns" in the passages where [7] and [8] would apply.
- Also listed, not about tour design: [3] Titan Saturn System Mission Final Report (2009); [4] Lindal 1992
  (Neptune atmosphere); [5] Kazeminejad et al. 2004 (Huygens); [6] Tauber et al. 1994 (probe study).
