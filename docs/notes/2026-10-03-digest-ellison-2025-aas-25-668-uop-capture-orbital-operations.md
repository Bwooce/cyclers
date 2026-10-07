# Digest — Ellison et al. (2025), "Uranus Orbiter and Probe: System Capture and Orbital Operations" (AAS 25-668)

> **Notice (2026-10-07, `#911`):** the six Uranian (1,1) quasi-cycler rows, including `#312`'s
> `umbriel-oberon-1-1-uranian-quasi-cycler-2026`, were WITHDRAWN from the catalogue on 2026-10-04
> (`#888`: a demanded-turn check found they are not ballistic trajectories; upheld by `#937`), and
> `umbriel-1-2-torus-homoclinic-uranus-2026` was WITHDRAWN on 2026-10-03 (`#882`: its connection is not a
> trajectory of its own model). `europa-3-4-crnbp-torus-jupiter-2026` is `known-class-member`. Any
> statement below that calls these rows novel, validated or catalogued is superseded. This note is a
> dated record and is otherwise not rewritten.

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, 23 pages). **Purpose:** the
Uranus Orbiter and Probe (UOP) era mission-design read that the six catalogued Uranian two-moon
quasi-cyclers were provisionally labelled "candidate-novel" pending. This note is facts only; the
novelty verdict is written elsewhere.

## Citation
Donald H. Ellison, Noble A. Hatten, Jacob A. Englander, Zachary R. Putnam, Jackson L. Shannon,
Kyle M. Hughes, Alec J. Mudek, Taabish Z. Rashied, "URANUS ORBITER AND PROBE: SYSTEM CAPTURE AND
ORBITAL OPERATIONS", paper **AAS 25-668** (printed top right of page 1).

- Affiliations as printed: Ellison, Hatten, Englander, Putnam, Shannon: Mission Design Engineer,
  Space Astrodynamics and Controls Group, Johns Hopkins Applied Physics Lab, Laurel, MD, USA, AAS
  Member. Hughes, Mudek, Rashied: Navigation and Mission Design Branch, NASA's Goddard Space Flight
  Center, Greenbelt, MD.
- Conference name: **not stated** on the paper. Date: **not stated**. (Only the sibling paper
  AAS 25-669, reference [31], is printed with "2025 AAS/AIAA Astrodynamics Specialists Conference,
  Boston, MA, 2025"; the paper does not say AAS 25-668 appeared there.)
- Filed in the private paper corpus as
  `ellison-hatten-englander-putnam-shannon-hughes-mudek-rashied-2025-uranus-orbiter-probe-system-capture-orbital-operations-AAS-25-668.pdf`,
  md5 `0bd37dca2812815b7c62b4bb534bde02`.

## What the paper is
Abstract (quoted): "This paper presents a Uranus system capture trade study, which poses the
end-to-end capture, probe delivery, and tour initialization sequence as a multi-agent optimization
problem. Additionally, the accommodations required to establish a successful science tour concept
of operations using an RTG-powered spacecraft at 20 astronomical units (au) are discussed."
Interplanetary cruise is stated to be out of scope ("beyond the scope of this work").

Mission phases listed in Sec. 3: (1) system approach and capture (orbit insertion); (2) atmospheric
probe trajectory targeting; (3) post-separation orbiter divert burn; (4) probe primary mission and
orbiter flyover; (5) orbit raise / first satellite encounter targeting; (6) orbit equatorialization;
(7) primary system science tour.

Sections: 1 Introduction; 2 Accessing the Outer Solar System (Hohmann TOF, capture delta-v,
rings, JWST ring imagery, 97.77 deg obliquity); 3 UOP Mission Components (3.1 System Entry: Orbit
Insertion and Avoiding the Rings; 3.2 Probe Targeting, Divert, and Flyover; 3.3 Equatorialization;
3.4 System Science Tour); 4 Practical Design of System Capture and Orbital Operations (4.1 System
Approach and Capture; 4.2 Orbiter-Probe Multi-Vehicle Optimization Problem; 4.3 Probe Descent
Modeling; 4.4 Orbiter Trajectory to First Moon Encounter; 4.5 Equatorialization; 4.6 Science
Tour); 5 Example Mission Design; 6 Conclusion. Eight tables and 18 figures.

Example trajectory (Sec. 5, stated to be "an example, intended as a demonstration"; "They do not
necessarily represent the decisions that would be made in practice"): 13.75-year, 2-Earth-gravity-
assist solar-electric cruise, launch 14 Oct. 2035, patched-conic Uranus arrival 14 Jul. 2049
(Table 3); capture rp 50559 km, Bθ 210 deg, 130-day capture orbit, INC 120.6 deg (Uranus TODBCI),
patched-conic capture delta-v 1296 m/s, finite-burn 1331 m/s over 1.62 h with a Leros 4-ET (Tables
4-5); probe release to atmospheric entry interface (AEI) 20 days, probe-targeting delta-v 57.1 m/s,
divert 85.6 m/s, Titania-targeting delta-v 240.3 m/s, maximum Titania arrival V-infinity 4.00 km/s
(Table 6).

## Equatorialization (Sec. 3.3 and 4.5)
- Need: "For most capture scenarios at Uranus, including the UOP baseline system science tour, the
  orbiter must equatorialize its initially inclined orbit in order to access the ring/moon plane."
  Footnote: an exception is a polar tour, such as the APL 2010 Decadal design [3]. Capture orbits
  are "highly inclined due to the obliquity of Uranus's equator and the large angle between the
  spacecraft's arrival velocity vector and Uranus's velocity vector."
- Mechanism: "Post-capture equatorialization at Uranus is achieved via a combination of maneuvers
  and flybys of Titania, the most effective gravity-assist body in the Uranus system [23]."
- Convention footnote: equatorialization is called "reducing" INC "even though equatorialization
  actually increases INC to 180 deg due to the retrograde motion of Uranus's moons."
- Apo-twist option: "If delta-ECC is small, then one option is to use an apo-twist maneuver [24] to
  perform all or nearly all of the equatorialization in one fell swoop. This option can require
  significant delta-v, but saves time by avoiding repeated Titania flybys." A small delta-ECC "is
  only achievable when the spacecraft arrives at Uranus near equinox and the arrival asymptote is
  favorable." (delta-ECC is the periapse latitude of the post-capture orbit, Sec. 3.2.)
- Otherwise: "an apo-twist does not produce a simple INC change; instead, both INC and right
  ascension of the ascending node (RAAN) are changed." "A straightforward approach is to use an
  apo-twist to alter the orbit's line of nodes so that it is approximately orthogonal to the
  orbital eccentricity vector. With this geometry, each flyby of Titania—which necessarily occurs
  at the nodal crossing—produces a relatively 'pure' INC change." Depicted in Figure 10, "used in
  the 2021 UOP Decadal design" (Figure 11, "Inclination reduction ... using repeated Titania
  flybys").
- Cost: "the greater the INC change that must be performed, the more expensive it is to do, in
  terms of delta-v, TOF, or both."
- Related relocation: "Reorienting the capture orbit periapse to the day side, via petal rotation
  using repeated flybys of Titania, also requires significant time, as shown in designs by Landau
  et al. [13] and Ellison and Scott [25] ... this relocation phase takes 2-to-3 years."
- Tools (Sec. 4.5): "Design of this mission phase is generally accomplished via a low-fidelity,
  multiple-flyby, path-planning, grid-search tool, such as those described in [13] and [25]. These
  tools can be used to find trajectories that simultaneously reduce the orbit period, in addition
  to equatorializing, in order to decrease the time between moon encounters ..." Period reduction
  "cannot be performed ad infinitum for operational reasons". Then "trajectory optimization
  software can be used to convert zero-sphere-of-influence patched-conic flybys into fully
  propagated trajectories with high-fidelity dynamics."
- In the example (Sec. 5): "Equatorialization is performed with the first 10 Titania flybys (above
  the horizontal line in Table 8), which take place over 1.48 years and use 4.95 m/s of
  deterministic delta-v." (See the table note below on the count.)

## Science tour (Sec. 3.4, 4.6 and 5)
**Requirements (Table 1, 2021 UOP baseline and threshold):**

| Category | Baseline | Threshold |
|---|---|---|
| orbital tour duration | 4 years | 2 years |
| major satellite flybys | 3 targeted, 2 non-targeted | 2 targeted, 1 non-targeted |
| minor satellite flybys | targeted and non-targeted | non-targeted only |
| satellite ground tracks | polar and low-inclination passes | polar only |
| Uranus orbits | close (1.1 Uranus radii) polar and low-inclination dayside passes | polar passes only |

Footnote: "The major Uranian satellites are Miranda, Ariel, Umbriel, Titania, and Oberon." Stated
challenges: total mission duration, uncertainty about the inner zeta ring, and the capture-orbit
periapse being on the night side.

**Operations constraints (Sec. 4.6, Table 7):** minimum TOF between encounters 30 days; minimum TOF
between maneuver and encounter 7 days; minimum flyby altitude 50 km; minimum 3 encounters with each
major moon. Rationale quoted: "If the orbit period is too short (and satellite re-encounter cadence
too quick), then the spacecraft may not be able to downlink all of the data". The 2021 Decadal
"assessed orbital periods in the 30-35 day range for a 3-RTG system"; Europa Clipper is cited at
"flybys every ~18 days".

**Design method (Sec. 4.6 and 5):** "Preliminary science tour design is performed using the same
tools referenced in section 4.5. The primary differences are that (1) inclination reduction is no
longer required, and (2) Titania is no longer the only flyby body under consideration, so the
search space is larger." "the primary science tour was designed using a C++ implementation of the
Star algorithm [30]. The software computes conic transfers (Lambert, resonant, powered leveraging,
etc.) between body states on a fixed-time grid. It is similar to the software that was employed to
design the Europa Clipper pump down trajectory[25], except that it does not perform flyby
C3-matching, which enables the fixed-time grid and leg evaluation in an arbitrary order, leading to
considerable computational time savings."

Search statistics (Sec. 5): "The search assembled a patched-conic leg graph of 30 million legs from
which it identified 79,272 valid tours. Of these tours, 12 encountered the five major satellites at
least 3 times each. The search took ~75 seconds executing with 20 threads using an Intel Core Ultra
9 processor." "The algorithm searched for tours with up to 18 encounters (including Titania-1)."
"In-depth detail on tour selection criteria is beyond the scope of this paper, but tradable
criteria include TOF, delta-v, and the suitability of each flyby to obtain measurements needed by
the various scientific disciplines" (e.g. location of moon within its orbit, flyover groundtrack,
lighting). Also: "the sequence of flybys in Table 8 represents only one possible sequence of many
that can be generated by tools like [30] and [25]."

**Table 8, exact transcription** (columns as printed; dates are M/D/YYYY as printed; "delta-v" rows
are maneuvers). The horizontal rules in the printed table fall after the 9/29/2051 row (end of the
equatorialization part) and after the 6/7/2053 row (table end); the rule position was confirmed from
the rendered page, not the text layer.

| Date | Event | Body | V-inf (km/s) | delta-v (m/s) | Flyby altitude (km) |
|---|---|---|---|---|---|
| 4/7/2050 | flyby | Titania | 4.000 | N/A | 50 |
| 6/3/2050 | delta-v | N/A | N/A | 1.86 | N/A |
| 7/4/2050 | flyby | Titania | 3.985 | N/A | 50 |
| 8/22/2050 | delta-v | N/A | N/A | 1.51 | N/A |
| 9/20/2050 | flyby | Titania | 3.980 | N/A | 50 |
| 11/1/2050 | delta-v | N/A | N/A | 0.91 | N/A |
| 11/29/2050 | flyby | Titania | 3.969 | N/A | 50 |
| 1/6/2051 | delta-v | N/A | N/A | 0.29 | N/A |
| 1/29/2051 | flyby | Titania | 3.973 | N/A | 50 |
| 2/25/2051 | delta-v | N/A | N/A | 0.18 | N/A |
| 3/22/2051 | flyby | Titania | 3.970 | N/A | 50 |
| 4/7/2051 | delta-v | N/A | N/A | 0.01 | N/A |
| 5/4/2051 | flyby | Titania | 3.967 | N/A | 50 |
| 5/27/2051 | delta-v | N/A | N/A | 0.11 | N/A |
| 6/17/2051 | flyby | Titania | 3.970 | N/A | 50 |
| 7/4/2051 | delta-v | N/A | N/A | 0.05 | N/A |
| 7/22/2051 | flyby | Titania | 3.972 | N/A | 50 |
| 8/6/2051 | delta-v | N/A | N/A | 0.02 | N/A |
| 8/25/2051 | flyby | Titania | 3.968 | N/A | 50 |
| 9/29/2051 | flyby | Titania | 3.966 | N/A | 50 |
| *(horizontal rule)* | | | | | |
| 10/6/2051 | delta-v | N/A | N/A | 11.15 | N/A |
| 11/16/2051 | flyby | Umbriel | 4.396 | N/A | 552 |
| 12/30/2051 | flyby | Umbriel | 4.392 | N/A | 62 |
| 2/3/2052 | flyby | Umbriel | 4.395 | N/A | 1721 |
| 3/10/2052 | flyby | Titania | 3.880 | N/A | 784 |
| 4/17/2052 | flyby | Oberon | 3.478 | N/A | 765 |
| 5/23/2052 | flyby | Ariel | 4.384 | N/A | 791 |
| 6/23/2052 | flyby | Oberon | 3.448 | N/A | 5164 |
| 7/28/2052 | flyby | Ariel | 4.325 | N/A | 116 |
| 9/3/2052 | flyby | Titania | 3.922 | N/A | 592 |
| 10/6/2052 | flyby | Ariel | 4.487 | N/A | 227 |
| 11/12/2052 | flyby | Ariel | 4.488 | N/A | 3662 |
| 12/21/2052 | flyby | Ariel | 4.487 | N/A | 3671 |
| 1/27/2053 | flyby | Ariel | 4.488 | N/A | 194 |
| 3/2/2053 | flyby | Oberon | 3.479 | N/A | 130 |
| 4/3/2053 | flyby | Miranda | 3.370 | N/A | 50 |
| 4/17/2053 | delta-v | N/A | N/A | 1.10 | N/A |
| 5/6/2053 | flyby | Miranda | 3.342 | N/A | 50 |
| 5/21/2053 | delta-v | N/A | N/A | 1.11 | N/A |
| 6/7/2053 | flyby | Miranda | 3.311 | N/A | N/A |

(Note: the last row's flyby altitude is printed "N/A", as printed.)

**Stated totals (Sec. 5):** equatorialization (first 10 Titania flybys, "above the horizontal line")
1.48 years, 4.95 m/s deterministic delta-v; tour "another 1.69 years and 13.37 m/s of deterministic
delta-v"; "Altogether, the moon-encounter phases take 3.17 years, use 18.32 m/s of deterministic
delta-v, and satisfy the Baseline requirement of 3 targeted encounters of each major satellite
(Table 1)."

**Checks against the table (arithmetic by the digester, not statements of the paper):**
- Eleven Titania flybys (4/7/2050 to 9/29/2051) sit above the horizontal rule, not ten as the text
  says. 4/7/2050 to 9/29/2051 is 540 days (1.48 yr); 9/29/2051 to 6/7/2053 is 1.69 yr; 4/7/2050 to
  6/7/2053 is 3.17 yr. These durations match the stated totals.
- Summing the printed delta-v rows gives 4.94 m/s above the rule (stated 4.95), 13.36 m/s below
  (stated 13.37), 18.30 m/s total (stated 18.32); differences are consistent with rounding of the
  printed rows.

## Observations (facts only)
All intervals below are calendar days between the printed dates of consecutive flybys.

**(a) Runs of consecutive flybys of the same moon (Table 8)**
- Titania, 11 flybys: 4/7/2050, 7/4/2050, 9/20/2050, 11/29/2050, 1/29/2051, 3/22/2051, 5/4/2051,
  6/17/2051, 7/22/2051, 8/25/2051, 9/29/2051. Intervals: 88, 78, 70, 61, 52, 43, 44, 35, 34, 35
  days. This is the equatorialization phase; all at 50 km altitude.
- Umbriel, 3 flybys: 11/16/2051, 12/30/2051, 2/3/2052. Intervals: 44, 35 days.
- Ariel, 4 flybys: 10/6/2052, 11/12/2052, 12/21/2052, 1/27/2053. Intervals: 37, 39, 37 days.
- Miranda, 3 flybys: 4/3/2053, 5/6/2053, 6/7/2053. Intervals: 33, 32 days.
- No run of two or more consecutive flybys exists for Oberon (3 flybys, none adjacent) or, in the
  tour part, Titania (3/10/2052 and 9/3/2052, separated by other moons).

**(b) Alternation of two different moons (A, B, A, B)**
- One place: **Oberon 4/17/2052, Ariel 5/23/2052, Oberon 6/23/2052, Ariel 7/28/2052** (four
  consecutive flybys). Intervals: 36 (O to A), 31 (A to O), 35 (O to A) days. V-infinity:
  Oberon 3.478, Ariel 4.384, Oberon 3.448, Ariel 4.325 km/s. Altitudes 765, 791, 5164, 116 km. It
  is preceded by Titania 3/10/2052 (38 days before Oberon 4/17) and followed by Titania 9/3/2052
  (37 days after Ariel 7/28), so the alternation lasts exactly four flybys.
- A shorter A, B, A pattern also appears at Ariel 7/28/2052, Titania 9/3/2052, Ariel 10/6/2052
  (intervals 37, 33 days; V-infinity 4.325, 3.922, 4.487 km/s); it is not extended to a fourth
  flyby because Ariel repeats at 11/12/2052.
- No other alternation of two moons occurs in the table. The tour contains no Umbriel-Oberon,
  Titania-Oberon, Umbriel-Titania, Ariel-Titania (beyond the A, B, A above) or Ariel-Umbriel
  alternation of length four.

**(c) V-infinity range at each moon (km/s, from Table 8)**
- Titania: 3.880 to 4.000 (equatorialization 3.966 to 4.000; tour 3.880 and 3.922).
- Umbriel: 4.392 to 4.396.
- Oberon: 3.448 to 3.479.
- Ariel: 4.325 to 4.488.
- Miranda: 3.311 to 3.370.
- Overall 3.311 to 4.488 km/s; no flyby in Table 8 is below 3.3 km/s.

**(d) Use of the words "cycler", "periodic", "free-return", "resonant", "repeating"**
- "cycler", "periodic", "free-return", "repeating": not used anywhere in the paper.
- "resonant": once, Sec. 4.6: "The software computes conic transfers (Lambert, resonant, powered
  leveraging, etc.) between body states on a fixed-time grid." It is a description of the Star
  software's leg types, not of a trajectory in Table 8.
- Related words: "repeated Titania flybys" (Sec. 3.3, twice; Sec. 4.5; Figs. 10 and 11 captions,
  "repeated satellite flybys", "repeated Titania flybys"); "re-encounter cadence" (Sec. 4.6, three
  times); "encounter cadence" (Sec. 3.4); "pump down trajectory" (Europa Clipper, Sec. 4.6).

## What it does NOT contain
- No cycler, periodic-orbit or free-return analysis; no use of those words about any Uranian
  trajectory.
- No two-moon shuttle analysis, no enumeration of moon-pair resonances, no statement of which
  moon pairs support repeated ballistic transfers, and no V-infinity-versus-bend feasibility table
  for the moons.
- No ballistic-only requirement: the example tour uses deterministic delta-v (18.32 m/s total, with
  an 11.15 m/s maneuver on 10/6/2051); tour selection criteria are explicitly "beyond the scope".
- No ephemeris or dynamical-model description for the tour beyond "patched-conic" and "zero-sphere-
  of-influence" flybys converted to high fidelity; the paper does not state the moon ephemeris used.
- No printed conference name or date, no tour figures of the moon-to-moon sequence (the figures
  cover capture, probe delivery, entry corridor and the first-year event sequence), and no list of
  the 12 tours that met the three-encounters-per-moon test beyond the single Table 8 example.
- No interplanetary cruise design (deferred to Englander et al., AAS 25-669, reference [31]).

## Bibliography (references concerning Uranus mission design, satellite tours, tour software)
Reference numbers as printed.

- [1] National Academies of Sciences, Engineering, and Medicine, "Origins, Worlds, and Life: A
  Decadal Strategy for Planetary Science and Astrobiology 2023-2032," 2023. Report of the Committee
  on a Decadal Strategy for Planetary Science and Astrobiology 2023-2032, 10.17226/26522.
- [2] W. B. Hubbard, M. Marley, H. B. Hammel, A. A. Simon-Miller, K. Khurana, B. Hesman, J. Clarke,
  L. Dudzinski, K. Lindstrom, E. Turtle, H. Seifert, D. Eng, R. Gold, E. Adams, D. Royster,
  S. Oleson, M. Mcguire, D. Grantier, Y. Guo, J. Dankanich, C. Scott, R. Russell, J. Drexler,
  L. Wofarth, S. Whitley, J. Drexler, M. Jones, A. Tenteris, M. K. Lockwood, D. Powell, D. Way,
  B. Sequeira, J. Warner, R. Vaughn, W.-J. Shyong, M. Martini, R. Reinders, J. White, S. Cooper,
  D. Napolillo, J. Gyekenyesi, M. Briere, E. Abel, T. Colozza, S. Williams, M. Fraeman,
  G. L. Williams, J. Fincannon, K. Bury, S. Bushman, and J. Fittje, "Ice Giants Decadal Study,"
  June 2010.
- [3] J. McAdams, C. Scott, Y. Guo, J. Dankanich, and R. Russell, "AAS 11-188 Conceptual mission
  design of a polar Uranus orbiter and satellite tour," Advances in the Astronautical Sciences,
  Vol. 140, January 2011.
- [4] J. Dankanich and J. McAdams, "AAS 11-189 Interplanetary Electric Propulsion Uranus Mission
  Trades Supporting the Decadal Survey," Advances in the Astronautical Sciences, Vol. 140,
  January 2011.
- [5] M. Hofstadter, A. Simon, S. Atreya, D. Banfield, J. Fortney, A. Hayes, M. Hedman,
  G. Hospodarsky, A. Masters, K. Mandt, M. Showalter, K. Soderlund, D. Turrini, E. Turtle,
  J. Elliott, K. Reh, P. Agrawal, T. Anderson, D. Atkinson, N. Arora, C. Borden, B. Martin,
  J. Cutts, H. Hwang, M. Le, Y. Lee, A. Petropoulos, S. Saikia, T. Spilker, and W. Smythe, "Ice
  Giants Pre-Decadal Survey Mission Study Report," June 2017.
- [6] S. Bayon et al., "ESA M* Ice Giant Concurrent Design Facility (CDF) study 1," tech. rep., ESA
  Concurrent Design Facility, 2019.
- [7] A. Simon, F. Nimmo, R. Anderson, I. Cohen, R. Gold, A. Azarbarzin, D. Cattopadhyay,
  M. Harrow, J. Arrieta, M. Ozimek, C. Scott, A. Calloway, A. Berman, P. McCauley, H. Hwang,
  D. Prabhu, G. Allen, J. Monk, J. Thornton, and M. Steerman, "Planetary Mission Concept Study for
  the 2023-2032 Decadal Survey: Uranus Orbiter Probe (UOP) Mission Concept Design Study Final
  Report," June 2021.
- [8] D. Ellison, D. Landau, J. Englander, K. Hughes, and N. Hatten, "Uranus Flagship Mission
  Design," Uranus Flagship Workshop: Investigating New Paradigms for Outer Planet Exploration,
  2024.
- [9] D. Ellison and K. Hughes, "Uranus Flagship Mission Design: Challenges and Opportunities,"
  Outer Planets Assessment Group, 2024.
- [10] D. Landau, "Uranus Cruise and Tour Design Impacts on Science Return, Cost, and Risk," Outer
  Planets Assessment Group, 2024.
- [11] R. Restrepo, "Maximizing Science Return for a Uranus Flagship Mission Using Aerocapture,"
  Outer Planets Assessment Group, 2024.
- [12] I. J. Cohen, E. J. Smith, G. B. Clark, D. L. Turner, D. H. Ellison, B. Clare, L. H. Regoli,
  P. Kollman, D. T. Gallagher, G. A. Holtzman, J. J. Likar, T. Morizono, M. Shannon, and
  K. S. Vodusek, "Plasma Environment, Radiation, Structure, and Evolution of the Uranian System
  (PERSEUS): A Dedicated Orbiter Mission Concept to Study Space Physics at Uranus," Space Science
  Reviews, Vol. 219, October 2023, 10.1007/s11214-023-01013-6.
- [13] D. Landau, R. Persinger, R. Karimi, M. Hofstadter, J. Castillo-Rogez, K. Mitchell,
  J. Elliott, S. Weinstein-Weiss, and C. Raymond, "Uranus Cruise and Tour Design Impacts on
  Science, Cost, and Risk," IEEE Aerospace Conference, January 2025.
- [14] R. R. Karimi, D. F. Landau, S. Weinstein-Weiss, and J. Elliott, "Destination Uranus:
  Interplanetary and Capture Trajectory Architecture Analysis for a Flagship-Class Orbiter and
  Probe," AIAA SciTech Forum, Orlando, FL, 2025.
- [15] R. Restrepo, D. Mages, M. Smith, R. Deshmukh, S. Dutta, and L. Benhacine, Mission Design and
  Navigation Solutions for Uranus Aerocapture. AIAA, 2024, 10.2514/6.2024-0715.
- [23] N. J. Strange, D. F. Landau, and J. M. Longuski, "AAS 13-801 Design of Initial Inclination
  Reduction Sequence for Uranian Gravity-Assist Tours," AIAA/AAS Astrodynamics Specialist
  Conference and Exhibit, 08 2013.
- [24] R. W. Luidens and B. A. Miller, "NASA TN D-3220: Efficient Planetary Parking Orbits with
  Examples for Mars," tech. rep., Lewis Research Center, 1965.
- [25] D. Ellison and C. Scott, "Europa Clipper Mission Analysis: A Flexible Graph-Based Search
  Automaton For Gravity-Assist Trajectory Design," Astrodynamics Specialist Conference, 2024.
- [29] R. Anderson, "Uranus Orbiter and Probe (UOP): Mission Architecting Challenges," Uranus
  Flagship Workshop: Investigating New Paradigms for Outer Planet Exploration, 2024.
- [30] D. Landau, S. Campagnola, and E. Pellegrini, "Star Searches for Patched-Conic Trajectories,"
  Journal of Astronautical Sciences, 2022.
- [31] J. A. Englander, D. H. Ellison, N. Hatten, K. M. Hughes, A. J. Mudek, and T. Z. Rashied,
  "Robust Access to Uranus via Solar Electric Propulsion," 2025 AAS/AIAA Astrodynamics Specialists
  Conference, Boston, MA, 2025. Paper AAS 25-669.

Not Uranus-mission or tour-design references (omitted, atmosphere/rings/Juno/Cassini sources):
[16] Kowalkowski et al. Juno launch period; [17] Cassini trajectory web page; [18] Elliot et al.
1977; [19] Smith et al. 1986; [20] Showalter and Lissauer 2006; [21] Hedman et al. 2023;
[22] JWST news release; [26] Bishop et al. 1990; [27] Herbert et al. 1987; [28] Lindal et al. 1987.

Where the paper cites them: [3] and [4] in the introduction (conceptual heritage); [13], [14],
[15] and Englander et al. for interplanetary cruise (Sec. 1); [13] and [25] for equatorialization
and petal-rotation tools and timelines (Sec. 3.4, 4.5); [23] for Titania as the most effective
gravity-assist body (Sec. 3.3); [30] for the Star algorithm (Sec. 4.6); [31] for the cruise
trajectory (Sec. 5).

## Coordinator's verdict for the six Uranian quasi_cycler rows (added 2026-10-03, `#869`)

This section is a judgement by the coordinating session, written after reading the tour material
directly; everything above it is the factual digest.

**Question:** does this paper publish a repeating ballistic two-moon geometry at Uranus, which
would bear on the `candidate-novel` label of the six catalogued rows?

**Answer: no periodic or cycler trajectory is published here, but one passage is closer than
anything else in the corpus and must be cited.**

- The paper designs a one-off tour with a patched-conic broad search (the "Star" algorithm). It
  never uses the words cycler, periodic, free-return or repeating about a trajectory, and it
  analyses no trajectory as a repeating pattern.
- Inside the example tour (Table 8) there is one stretch of four alternating flybys of two moons:
  Oberon, Ariel, Oberon, Ariel on 2052-04-17, 05-23, 06-23 and 07-28, at intervals of 36, 31 and
  35 days. The V-infinity nearly repeats from one lap to the next (Oberon 3.478 then 3.448 km/s;
  Ariel 4.384 then 4.325 km/s). That is, in effect, two laps of an Ariel-Oberon alternation,
  bracketed by Titania flybys, and it is not discussed as such by the authors.
- The catalogued `ariel-oberon-1-1-uranian-quasi-cycler-2026` is the same moon pair but a
  different geometry: legs of 7.75 days each (against 31-36 days), V-infinity 1.52 km/s at Ariel
  and 1.83 km/s at Oberon (against 4.3-4.4 and 3.4-3.5 km/s), and it is a closure that repeats
  over epoch windows through 2083 rather than two laps inside a longer sequence.
- No other moon pair alternates A, B, A, B. There is one A, B, A (Ariel, Titania, Ariel).
- All tour V-infinity values are 3.3-4.5 km/s; the six catalogued rows sit at 0.89-2.16 km/s.

**Consequence proposed (owner's call, not applied here):** keep `candidate-novel` on all six
rows, and add this paper to the Ariel-Oberon row's `corroborating_sources` as the nearest
published geometry, with the differences above stated. The statement "no published two-moon
alternating geometry exists at Uranus" should not be made any more; "no published repeating or
periodic two-moon trajectory, and none in this V-infinity regime" is what the corpus supports.
The Landau et al. 2025 IEEE Aerospace paper ([13] here, DOI 10.1109/AERO63441.2025.11068400) and
the Simon et al. 2026 concept update (DOI 10.3847/PSJ/ae680c) are still unread, so `#869` is not
closed.
