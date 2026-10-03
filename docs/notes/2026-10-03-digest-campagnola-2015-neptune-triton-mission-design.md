# Digest — Campagnola, Boutonnet, Martens, Masters (2015), "Mission Design for the Exploration of Neptune and Triton" (IEEE AESM)

**Digested:** 2026-10-03 (12-page PDF; read in full from the rendered pages). **Text-layer caveat:** the PDF's
text layer is font-encoded (body text comes out as shifted characters with no spaces, and numerals are lost), so
`pdftotext` is not usable as the primary source. Every number below, including Tables 1-8, was read from the
rendered pages; Table 6 was read from a 220 dpi crop. The term counts below come from a shift-decoded copy of the
text layer and are approximate (they cannot see text inside figures; figure labels are noted separately).
**Purpose:** record exactly what in-system trajectory content at Neptune is published before the project works on
Neptune-Triton resonant-orbit objects. Facts only; no novelty verdict is written here.

## Citation
Stefano Campagnola, Arnaud Boutonnet, Waldemar Martens, Adam Masters, "Mission Design for the Exploration of
Neptune and Triton", IEEE Aerospace and Electronic Systems Magazine, July 2015, pages 6-17 (as printed in the
running footers: "IEEE A&E SYSTEMS MAGAZINE", "JULY 2015").

- DOI as printed (page 6): "DOI No. 10.1109/MAES.2015.140119". Also printed: "0885/8985/15 $26.00 (c) 2015
  IEEE"; "Review handled by M. Jah."
- Dates as printed: "Manuscript received June 23, 2014, revised May 4, 2015, and ready for publication May
  6, 2015."
- Volume and issue numbers: not printed on any page read (the task record gives 30(7); the pages print only the
  month and year).
- Affiliations as printed: Campagnola, ISAS/JAXA, Kanagawa, Japan; Boutonnet and Martens, ESOC/ESA, Darmstadt,
  Germany; Masters, Imperial College London, UK. (The authors' address block prints the first author's name as
  "S. Campagnoia"; the title block prints "Stefano Campagnola".)
- Filed in the private paper corpus as
  `campagnola-boutonnet-martens-masters-2015-mission-design-exploration-neptune-triton-ieee-aes-magazine-30-7-doi-10.1109-MAES.2015.140119.pdf`,
  md5 `b4c750065eb317bd85fe1b2624f94331`.

## What the paper is
An ESA L-class (L2/L3) mission concept to Neptune and Triton, written from the trajectory-design part of a 2013
white paper. Introduction, quoted: "We present here the trajectory design for the mission concept and a
high-level discussion on its feasibility, which was included in the white paper. We also present details on the
example Triton tour, a new search for interplanetary transfers, and a preliminary analysis of the gravity losses
at Neptune orbit insertion (NOI), which suggests the use of a large chemical propulsion system." And: "The last
section presents an example 2-year moon tour at Neptune, which demonstrates that all scientific questions can be
effectively addressed using Triton as a tour engine. The tour includes 55 flybys and covers a range of Neptune
orbits and Triton flyby geometries. Another example of a Triton tour can be found in [3]."

Sections as printed: Introduction; A Mission to Neptune and Triton; Payload Options (Tables 1, 2); Enabling
Technologies (tracking network, RTGs, solar electric propulsion, enhancing technologies); Gravity Losses at NOI;
Mission Summary (Table 3); Interplanetary Transfer to Neptune; Chemical Propulsion Option (Table 4, Figures
1-4); SEP Option (Table 5, Figures 5-6); Neptune Orbital Tour (Table 6, Figures 7-10); Appendix (Ariane 5 ECA
performance, Table 7; gravity losses at NOI, Table 8); References [1]-[22]. About two-thirds of the paper is
mission architecture and the interplanetary leg; the Neptune-system trajectory is the "Neptune Orbital Tour"
section (page 14-15) and Table 6 (page 13).

## Neptune-system trajectory content (the part that matters here)

### Which bodies are used
Triton only. The tour section states: "Although not included here, close flybys at other moons and
opportunities for Neptune or Triton occultation can be added." Proteus, Nereid and Larissa are not mentioned in
the text; "Nereid" (an orbit label) and "Larissa" and "Proteus" (small labels next to Neptune) appear only inside
the raster Figure 7 as orbit/body labels. No Neptunian moon other than Triton is flown by in the tour.

### Model
"The tour is computed using linked conic approximation and Triton ephemerides." Neptune-centred conic arcs
between Triton gravity assists (GAs). The interplanetary part uses the same "linked-conic approximation" with
tabular planetary ephemerides and a numerically integrated check, via the in-house global search tool SOURCE
(refs [7], [8]). The tour is designed on the v-infinity sphere (refs [18], [19], [20]); the text calls it "a
powerful graphical tool in which a complicated trajectory is represented by a discrete sequence of dots". A
Tisserand graph is not used in the body text; the word occurs once, in the title of reference [21]. No
three-body (CR3BP) model is used for the tour. NOI is assumed impulsive (with a gravity-loss analysis in the
appendix).

### Tour structure
Starting conditions, quoted: "The spacecraft flies inside the inner rings and executes NOI; at the first
apocenter, a PRM raises the pericenter outside the rings and targets the first Triton flyby T1." The first orbit
period is "100 days, since longer periods do not significantly reduce the NOI Delta-v"; the pericenter altitude
"should be as low as possible, and it is assumed to be at an altitude of 3,000 km (following previous NASA
mission concepts [1])". T1 is placed "close to the line of nodes between the orbital plane of Triton and
Neptune's equator" and "close to the Sun-Neptune direction".

"Starting with T1, a first sequence of flybys (called the crank-over-the-top, or COT, sequence [16]) increases
the inclination of the spacecraft orbit over Neptune's equatorial plane." High v-infinity at T1 "maximizes the
inclination achievable with the COT, while a low v-infinity maximizes the bending angle of each flyby and hence
the number of flybys and the total tour duration; for this example tour, we choose a v-infinity of ~3 km."
(The paper prints "~3 km" in the section text and "v-infinity ~ 3 km/s" is what Table 6 supports; the unit is
printed without "/s" in that sentence.)

The rest of the tour "has no deterministic orbital maneuvers and is split into three phases":

| Phase | Flybys | Table 6 colour | Text description (quoted or closely paraphrased) |
|---|---|---|---|
| I | 20 (T1-T20) | pink | "the first COT sequence with 20 Triton flybys (altitudes between ~150 and ~1,000 km) at the same Triton orbital location. The initial flybys have high altitudes to cope with the larger uncertainties in Triton's ephemerides. Neptune orbits in this phase range in inclination (115-160 deg), pericenter radius (75,000-250,000 km), and solar local time (2-11 p.m. at the apocenter)." Ground tracks "concentrated on the sunlit Neptune-facing hemisphere" (Figure 8). |
| II | 14 (T21-T34) | green | "The 14 flybys rotate the line of apsides anticlockwise using the 'petal strategy' technique, which alternates long and short nonresonant transfers [16], [17]. In this phase, the apocenter varies between 800,000 and 1,300,000 km in altitude and between 2 p.m. and 5 a.m. solar local time. The flyby ground tracks are equatorial, centered alternately at 0 deg and 180 deg Triton longitude, with minimum altitudes between 300 and 2,000 km. Flybys occur all along Triton's orbit at intervals of ~30 deg. The last flyby of this phase lies on the opposite side of Triton's orbit compared to the first flyby, close to the orbital node to maximize the achievable inclination of the following phase." |
| III | 21 (T35-T55) | yellow | "a second COT sequence with 21 Triton flybys. Neptune orbits are again varied in inclination (115-160 deg) and solar local time (7 a.m.-5 p.m. at the apocenter)." Ground tracks "concentrated in the sunlit anti-Neptune-facing hemisphere" (Figure 9). |

20 + 14 + 21 = 55, matching the abstract-level statement of 55 flybys. The paper prints the Phase I and III
inclination range as "115-160 deg" (the reference plane is not stated). Phase boundaries by
colour are those of Table 6 and Figure 7 (Phase I magenta, II green, III yellow).

### Resonances
No resonance is given as an m:n pair with Triton anywhere in the paper. Resonances enter only through the v-
infinity sphere description, quoted: "Figure 10 shows a two-dimensional map of the v-infinity sphere, with
'pump angle' and 'crank angle' used as spherical coordinates ... The figure shows three sequences of dots,
which are the three parts of the tour, and several contour lines that identify special orbits, in particular
resonant orbits (vertical lines, one per resonant period), fixed-inclination orbits (contours, where the labels
indicate the inclination in degrees), and orbits that collide with Neptune or its rings (shaded areas). Finally,
some closed curves are added to show the change in pump and crank angles provided by the GAs." Figure 10's
panel label reads "resonant orbits" above the vertical lines and "minimum-altitude flybys" for a shaded region;
the numerical labels on the vertical lines are not explained in the text. The Phase II description says the
transfers are "nonresonant" (alternating long and short). Which Phase I or III flyby pairs are resonant: not
stated.

### Table 6, exact transcription (page 13, "Tour Events")
Header block as printed: "Event | Epoch | Delta-v (km/s)": NOI, "2044 DEC 20", 2.45; PRM, "2044 FEB 12", 0.29.
Then columns "Event | Epoch | v-infinity-in (km/s) | v-infinity-out (km/s) | h (km)". The v-infinity columns
each have three sub-columns; the paper does not label the components or the frame (the text says "the incoming
and outgoing v-infinity vectors"; the first column is the largest component, near 3 km/s, and the three
components together are consistent with the ~3 km/s stated). h is the flyby altitude. Dates are printed as
"YYYY MON DD". Values are copied as printed (including "-0.00" and "-0.00" signs).

Note on the NOI/PRM rows, as printed and not reconciled in the paper: NOI is dated 2044 DEC 20 and the PRM 2044
FEB 12 (earlier), while Table 5 prints the Neptune arrival as 2043 Dec. 20 and the first flyby T1 as 2044 APR 02.

| Event | Epoch | v-inf in (a) | (b) | (c) | v-inf out (a) | (b) | (c) | h (km) |
|---|---|---|---|---|---|---|---|---|
| T1 | 2044 APR 02 | 2.80 | 0.87 | -0.50 | 2.91 | 0.52 | -0.39 | 1000 |
| T2 | 2044 MAY 08 | 2.91 | 0.52 | -0.39 | 2.97 | 0.13 | -0.12 | 500 |
| T3 | 2044 MAY 25 | 2.97 | 0.13 | -0.12 | 2.95 | -0.20 | 0.31 | 250 |
| T4 | 2044 JUN 06 | 2.95 | -0.20 | 0.31 | 2.84 | -0.49 | 0.75 | 250 |
| T5 | 2044 JUN 24 | 2.84 | -0.49 | 0.75 | 2.68 | -1.01 | 0.80 | 250 |
| T6 | 2044 JUN 30 | 2.68 | -1.01 | 0.80 | 2.48 | -1.01 | 1.30 | 250 |
| T7 | 2044 JUL 05 | 2.48 | -1.01 | 1.30 | 2.18 | -1.01 | 1.76 | 250 |
| T8 | 2044 JUL 11 | 2.18 | -1.01 | 1.76 | 1.90 | -1.47 | 1.76 | 250 |
| T9 | 2044 JUL 29 | 1.90 | -1.47 | 1.76 | 1.57 | -1.90 | 1.67 | 250 |
| T10 | 2044 AUG 16 | 1.57 | -1.90 | 1.67 | 1.22 | -2.30 | 1.45 | 250 |
| T11 | 2044 AUG 21 | 1.22 | -2.30 | 1.45 | 0.76 | -2.30 | 1.73 | 150 |
| T12 | 2044 AUG 27 | 0.76 | -2.30 | 1.73 | 0.54 | -2.00 | 2.14 | 250 |
| T13 | 2044 SEP 20 | 0.54 | -2.00 | 2.14 | 0.97 | -1.69 | 2.25 | 250 |
| T14 | 2044 OCT 02 | 0.97 | -1.69 | 2.25 | 1.46 | -1.47 | 2.14 | 250 |
| T15 | 2044 OCT 19 | 1.46 | -1.47 | 2.14 | 1.73 | -1.01 | 2.20 | 250 |
| T16 | 2044 OCT 25 | 1.73 | -1.01 | 2.20 | 2.13 | -1.01 | 1.82 | 250 |
| T17 | 2044 OCT 31 | 2.13 | -1.01 | 1.82 | 2.44 | -1.01 | 1.38 | 250 |
| T18 | 2044 NOV 06 | 2.44 | -1.01 | 1.38 | 2.66 | -1.01 | 0.88 | 250 |
| T19 | 2044 NOV 12 | 2.66 | -1.01 | 0.88 | 2.78 | -1.01 | 0.35 | 250 |
| T20 | 2044 NOV 18 | 2.78 | -1.01 | 0.35 | 2.92 | -0.49 | 0.28 | 250 |
| T21 | 2044 DEC 05 | 2.92 | -0.49 | 0.28 | 2.96 | -0.35 | -0.00 | 1510 |
| T22 | 2044 DEC 14 | -2.96 | -0.35 | 0.00 | -2.98 | -0.09 | -0.00 | 2149 |
| T23 | 2044 DEC 28 | 2.98 | -0.09 | 0.00 | 2.96 | -0.35 | -0.00 | 2149 |
| T24 | 2045 JAN 06 | -2.96 | -0.35 | 0.00 | -2.98 | -0.09 | -0.00 | 2149 |
| T25 | 2045 JAN 21 | 2.98 | -0.09 | 0.00 | 2.96 | -0.35 | -0.00 | 2149 |
| T26 | 2045 JAN 29 | -2.96 | -0.35 | 0.00 | -2.98 | -0.09 | -0.00 | 2149 |
| T27 | 2045 FEB 13 | 2.98 | -0.09 | 0.00 | 2.96 | -0.35 | -0.00 | 2149 |
| T28 | 2045 FEB 21 | -2.96 | -0.35 | 0.00 | -2.98 | -0.09 | -0.00 | 2149 |
| T29 | 2045 MAR 08 | 2.98 | -0.09 | 0.00 | 2.96 | -0.35 | -0.00 | 2149 |
| T30 | 2045 MAR 16 | -2.96 | -0.35 | 0.00 | -2.97 | 0.18 | -0.00 | 293 |
| T31 | 2045 APR 06 | 2.97 | 0.18 | -0.00 | 2.96 | -0.35 | 0.00 | 294 |
| T32 | 2045 APR 14 | -2.96 | -0.35 | -0.00 | -2.97 | 0.18 | -0.00 | 293 |
| T33 | 2045 MAY 04 | 2.97 | 0.18 | 0.00 | 2.96 | -0.35 | 0.00 | 294 |
| T34 | 2045 MAY 13 | -2.96 | -0.35 | -0.00 | -2.97 | 0.18 | -0.00 | 293 |
| T35 | 2045 JUN 02 | 2.97 | 0.18 | 0.00 | 2.95 | 0.20 | 0.39 | 250 |
| T36 | 2045 JUN 14 | 2.95 | -0.20 | -0.39 | 2.82 | -0.49 | -0.83 | 250 |
| T37 | 2045 JUL 01 | 2.82 | -0.49 | -0.83 | 2.66 | -1.01 | -0.87 | 250 |
| T38 | 2045 JUL 07 | 2.66 | -1.01 | -0.87 | 2,44 | -1.01 | -1.37 | 250 |
| T39 | 2045 JUL 13 | 2,44 | -1.01 | -1.37 | 2.13 | -1.01 | -1.82 | 250 |
| T40 | 2045 JUL 19 | 2.13 | -1.01 | -1.82 | 1.85 | -1.47 | -1.81 | 250 |
| T41 | 2045 AUG 06 | 1.85 | -1.47 | -1.81 | 1.52 | -1.90 | -1.71 | 250 |
| T42 | 2045 AUG 23 | 1.52 | -1.90 | -1.71 | 1.18 | -2.30 | -1.48 | 150 |
| T43 | 2045 AUG 29 | 1.18 | -2.30 | -1.48 | 0.71 | -2.30 | -1.75 | 250 |
| T44 | 2045 SEP 04 | 0.71 | -2.30 | -1.75 | 0.18 | -2.30 | -1.88 | 250 |
| T45 | 2045 SEP 10 | 0.18 | -2.30 | -1.88 | 0.09 | -1.90 | -2.29 | 150 |
| T46 | 2045 SEP 28 | 0.09 | -1.90 | -2.29 | -0.07 | -1.47 | -2.59 | 250 |
| T47 | 2045 OCT 15 | -0.07 | -1.47 | -2.59 | -0.27 | -1.01 | -2.79 | 250 |
| T48 | 2045 OCT 21 | -0.27 | -1.01 | -2.79 | -0.80 | -1.01 | -2.68 | 250 |
| T49 | 2045 OCT 27 | -0.80 | -1.01 | -2.68 | -1.31 | -1.01 | -2.48 | 250 |
| T50 | 2045 NOV 02 | -1.31 | -1.01 | -2.48 | -1.76 | -1.01 | -2.18 | 250 |
| T51 | 2045 NOV 08 | -1.76 | -1.01 | -2.18 | -2.15 | -1.01 | -1.80 | 250 |
| T52 | 2045 NOV 14 | -2.15 | -1.01 | -1.80 | -2.46 | -1.01 | -1.35 | 250 |
| T53 | 2045 NOV 19 | -2.46 | -1.01 | -1.35 | -2.67 | -1.01 | -0.85 | 250 |
| T54 | 2045 NOV 25 | -2.67 | -1.01 | -0.85 | -2.82 | -0.49 | -0.80 | 250 |
| T55 | 2045 DEC 13 | -2.82 | -0.49 | -0.80 | -2.95 | -0.20 | -0.36 | 250 |

Printed oddities, copied as they are: T38 prints "2,44" (T39 in) with a comma in place of the decimal point, in
one cell each for T38 out and T39 in. Row-to-row chaining: in Phases I and III each row's "out" vector equals
the next row's "in" vector. In the Phase II rows and at the phase boundaries it does not match literally: the
first component changes sign between a row's "out" and the next row's "in" (for example T21 out 2.96, -0.35,
-0.00 against T22 in -2.96, -0.35, 0.00; T34 out -2.97, 0.18, -0.00 against T35 in 2.97, 0.18, 0.00), and T35
out (2.95, 0.20, 0.39) against T36 in (2.95, -0.20, -0.39) differs in the sign of the second and third
components. The paper does not comment on these sign changes. The Phase II rows alternate between two vector
pairs (2.96 with -0.35, and 2.98 with -0.09, as printed) with flyby altitude 2149 km for T22-T29 and 293-294 km
for T30-T34. Epoch spacing between consecutive Phase II flybys ranges from 8 to 21 days (T21 2044 DEC 05 to T34
2045 MAY 13).

### Flyby altitudes, delta-v, durations (tour)
- Altitudes: T1 1000 km, T2 500 km, most later flybys 250 km (the paper: "the minimum altitude is 250 km. At
  250 km, the dynamic pressure (and hence heat loads and controllability) is similar to that experienced by
  Cassini during the Titan flybys"); 150 km at some flybys (the text prints "The flybys T10, T45, and T45 are
  examples of lower-altitude flybys (150 km)"; Table 6 prints 150 km at T11, T42 and T45). Phase II altitudes up
  to 2149 km (Table 6). The text explains the 150 km choice: "where the dynamic pressure is similar to that of a
  750- or 800-km altitude flyby at Titan with v-infinity ~ 6 km/s (Cassini's lowest-altitude flyby was at 880
  km, T70). The altitude of these flybys can be increased with a marginal penalty in the tour duration."
- Navigation delta-v: "For the tour, Delta-v ~ 250 m/s is allocated for navigation: ~3 m/s per flyby with the
  exception of the first 10 flybys, where ~10 m/s are allocated for the initial large uncertainties on Triton
  ephemerides. For the same reason, the altitude of the first flyby is kept at 1,000 km and the altitude of the
  second flyby is 500 km."
- Deterministic maneuvers: NOI Delta-v 2.45 km/s and PRM 0.29 km/s (Table 6); "The rest of the tour has no
  deterministic orbital maneuvers." Table 3: "Propulsion: Chemical, total Delta-v of 3 km/s (a, b)" with
  footnotes "(a) Assuming an impulsive NOI. (b) Including deterministic Delta-v as tour navigation."
- Duration: Table 3 says "Two-year orbital tour (covering all solar local times and a range of orbital
  inclinations)"; "55 Triton flybys (providing global surface coverage)". Table 6 epochs run from T1 2044 APR 02
  to T55 2045 DEC 13 (about 20 months by the printed dates).
- v-infinity at Triton: about 3 km/s chosen; Table 6 vectors from roughly 3 km/s down to components near zero
  and back; see the transcription above.

### Triton orbit insertion or endgame
Not designed. Text, quoted: "Although not included here, a Triton orbit phase could also be considered at the
end of the Neptune tour, or replacing part of it, at the additional cost of Delta-v ~ 300 m/s. Over the course
of 3-4 months, a sequence of high-altitude Triton flybys and orbital maneuvers would reduce the spacecraft
velocity relative to Triton, eventually leading to a gravitational capture orbit that is stabilized by a Triton
orbit insertion maneuver [21]. This type of transfer was used for the lunar orbit insertion of the Small
Missions for Advanced Research in Technology 1 spacecraft and is planned for JUICE Ganymede orbit insertion to
reduce the required propellant mass." Reference [21] is Campagnola et al., "Tisserand-leveraging transfers",
JGCD 2014.

## Interplanetary and mission-architecture content (for completeness)
- Table 3 (mission summary): launch 19 Dec. 2028 on Ariane 5 ECA; two Earth GAs and one Jupiter GA; transfer
  time 15 years; SEP module (ejected before the Jupiter flyby), propellant mass 695 kg; 10 European RTGs, 500 W
  at Neptune; mass at launch 6,467 kg; SEP module wet mass 1,516 kg; dry mass at Neptune 1,921 kg; payload mass
  70 kg.
- Table 5 (SEP trajectory): Launch 2028 Dec. 19, 6,467 kg, v-inf 0.76 km/s; Earth GA 1 2030 Feb. 21, 6,102 kg,
  4.83; Earth GA 2 2032 Apr. 13, 5,817 kg, 10.46; Jupiter GA 2033 Sep. 03, 5,772 kg, 11.73; Neptune arrival 2043
  Dec. 20, 5,772 kg (includes the SEP module), 10.56 km/s; post-NOI 2,267 kg (impulsive NOI). Low-thrust arcs
  are modelled as small impulsive maneuvers; three ion engines, total thrust 0.465 N max at 1 AU.
- Chemical option (SOURCE global search, launches 2025-2041): best Earth/Venus/Mars/Jupiter GA sequences listed
  in the text; Table 4 gives eight candidates (M2026, T2026, T2028, C2037, M2037, C2039, T2039, T2040) with
  launch dates, launch v-inf 2.94-4.09 km/s, arrival v-inf 10.4-11.9 km/s, NOI 2400-3060 m/s, durations 16.3-20.5
  years, final masses 1445-2089 kg. The text concludes chemical options "would not deliver sufficient mass into
  Neptune orbit unless the transfer time is allowed to increase to 16-17 years".
- Table 7: Ariane 5 ECA escape mass vs v-infinity (0.5 km/s: 6,686 kg to 5.5 km/s: 2,441 kg). Table 8: NOI
  gravity losses vs arrival v-infinity (10-12 km/s) and engine thrust (450 and 900 N): Delta-v 2.63-4.27 km/s,
  gravity losses 11.3-22.5 percent in mass and 18.9-38.2 percent in Delta-v; capture into a 100-day orbit with
  3,000 km pericenter altitude.
- Table 1 payload (about 70 kg total) and Table 2 science-theme matrix: payload only.

## Term search (shift-decoded text layer; approximate; figures not searched)
cycler 0; cycling 0; periodic 0; repeating 0; free-return 0; resonant 7 (includes "resonant" in the SOURCE
sentence "transfer types, including resonant transfers and v-infinity-leveraging transfers", "a 1:1 Earth
resonant orbit" for the Ariane 5 launch, "nonresonant transfers" in the petal strategy sentence, "resonant orbits
(vertical lines, one per resonant period)", "connect multiple flybys with resonant orbits", and the title of ref
[19], "3D resonant hopping strategies"); resonance 0; heteroclinic 0; homoclinic 0; quasi-periodic 0; torus 0;
four-body 0; Umbriel 0; Titania 0; Oberon 0; Ariel 0; Miranda 0; Triton about 46 (text); Proteus 0 and Nereid 0
in the text (both appear as labels inside Figure 7, with "Larissa"); Tisserand 1 (title of [21]); leveraging 3.

Sentences containing the requested terms:
- cycler, cycling, repeating, free-return: none. The words do not occur in the paper.
- Nearest to a repeating pattern, quoted: "Phase I consists of the first COT sequence with 20 Triton flybys
  (altitudes between ~150 and ~1,000 km) at the same Triton orbital location." and "Flybys occur all along
  Triton's orbit at intervals of ~30 deg." These are sequences of flybys of one moon (Triton) that are not
  stated to be periodic or to repeat.

## Does the paper contain a trajectory that returns to the same two moons, or a periodic or quasi-periodic orbit that encounters two moons
No. The only moon flown by in the Neptune system is Triton, in a 55-flyby tour with no other moon encounter. No
periodic, quasi-periodic or repeating orbit is constructed. The only repetition is Triton flybys repeated by a
single spacecraft on a non-periodic (deterministically designed, ephemeris-dependent) tour. There is no
two-moon trajectory at Neptune in the paper.

## What it does NOT contain
- Any flyby of Proteus, Nereid, Larissa or any Neptunian moon other than Triton.
- Any Uranus, Uranian-moon or two-moon content (Uranus appears only as a planet in a heliocentric figure
  label).
- Any m:n resonance with Triton stated numerically, any resonant periodic-orbit family, any CR3BP, three-body or
  four-body model of the Neptune-system tour, any invariant-manifold or Poincare-map computation.
- A Tisserand-graph construction (only the reference title).
- A designed Triton orbit insertion or Triton capture endgame (only described as an option at about 300 m/s).
- A computed ballistic-capture, periodic or free-return trajectory, or any labelled component frame for the
  Table 6 v-infinity vectors.
- A statement of how the three-component Table 6 vectors are expressed (frame not stated).

## Relevant bibliography (reference list as printed, in full)
[1] Marley, M., Dudzinski, L., Spilker, T. L., and Moeller, R. Planetary science decadal survey: JPL rapid
mission architecture Neptune-Triton-KBO study. NASA, Washington, DC, Final Rep., 2010.
[2] Christophe, B., Spilker, L., and Anderson, J. (2012) OSS (outer solar system): a fundamental and planetary
physics mission to Neptune, Triton and the Kuiper belt, Experimental Astronomy. 34 pp. 203.
[3] Ingersoll, A., and Spilker, T. R. A Neptune orbiter with probes mission with aerocapture orbit insertion. In
Progress in Astronautics and Aeronautics, NASA Space Science Vision Missions, 224. Reston, VA: AIAA, 2008,
81-113.
[4] Pessina, S. M., Campagnola, S., and Vasile, M. Preliminary analysis of interplanetary trajectories with
aerogravity and gravity assist maneuvers. In 54th International Astronautical Congress, Bremen, Germany, Oct.
2003.
[5] Yam, C. H., and McConaghy, T. T. Design of low-thrust gravity-assist trajectories to the outer planets. In
55th International Astronautical Congress, Vancouver, BC, Canada, 4-8 Oct. 2004.
[6] Landau, D. F., Lam, T., and Strange, N. J. Broad search and optimization of solar electric propulsion
trajectories to Uranus and Neptune. Advances in the Astronautical Sciences, Vol. 153, 3 (2009), 2093-2112.
[7] Boutonnet, A., Martens, W., and Schoenmaekers, J. SOURCE: A Matlab-oriented tool for interplanetary
trajectory global optimization-fundamentals (part I). AAS/AIAA SFMM, Kauai, HI, Ser. Paper AAS 13-300.
[8] Martens, W., Boutonnet, W. A., and Schoenmaekers, J. SOURCE: A Matlab-oriented tool for interplanetary
trajectory global optimization-applications (part II). AAS/AIAA SFMM, Kauai, HI, Ser. Paper AAS 13-301.
[9] Masters, A., et al. (Dec. 2014) Neptune and Triton: essential pieces of the solar system puzzle, Planetary
and Space Science, 104 pp. 108-121.
[10] Arridge, C. S., et al. Uranus Pathfinder: exploring the origins and evolution of ice giant planets.
Experimental Astronomy, Vol. 33, 2-3 (Apr. 2012), 753-791.
[11] Ambrosi, R., Williams, H., Samara-Ratna, P., Jorden, A., Slade, R. M., Jaegle, M., Koenig, J., Bannister,
N., Deacon, T., Stuttard, T., Crawford, E., and Vernon, D. Thermoelectric converter system for small-scale
RTGs. ESA TRP Tech Rep., TECS-RTG-TR-D8-001-UL, 2012.
[12] Dougherty, M. JUICE: Exploring the emergence of habitable worlds around gas giants. Assessment Study
Report, Tech Rep. ESA SRE(2011)18, 2011.
[13] Janhunen, P., Quarta, A., and Mengali, G. (2013) Electric solar wind sail mass budget model, Geoscientific
Instrumentation, Methods and Data Systems Discussions, 2(2) pp. 429-455.
[14] Shimazaki, K., Takahashi, M., Imaizumi, M., and Al. E. First approach to lightweight solar panel using space
solar sheet. In Proceedings of 8th European Space Power Conference, Constance, Germany, 14-19 Sep. 2008.
[15] Khan, M. Summary of generic results on Jupiter orbit insertion. ESA/ESOC, MAS-TN-058, Sept. 2006.
[16] Buffington, B. B., Strange, N. J., and Campagnola, S. Global moon coverage via hyperbolic flybys. In 23rd
International Symposium on Space Flight Dynamics, Pasadena, CA, 2012.
[17] Anderson, R. L., Campagnola, S., and Buffington, B. B. Analysis of petal rotation trajectory
characteristics. In Astrodynamics Specialist Conference, San Diego, CA, 2014.
[18] Strange, N. J., Russell, R. P., and Buffington, B. B. Mapping the V-infinity globe. In AAS/AIAA
Astrodynamics Specialist Conference and Exhibit, Mackinac Island, MI, 2007.
[19] Campagnola, S., and Kawakatsu, Y. 3D resonant hopping strategies and the Jupiter magnetospheric orbiter.
Journal of Guidance, Control, and Dynamics, Vol. 35, 1 (2012), 340-344.
[20] Kawakatsu, Y. V-infinity direction diagram and its application to swingby design. In Proceedings of the
21st International Symposium on Space Flight Dynamics, Vol. 1, Toulouse, France, 2009.
[21] Campagnola, S., Boutonnet, A., Schoenmaekers, J., Grebow, D. J., Petropoulos, A. E., and Russell, R. P.
Tisserand-leveraging transfers. Journal of Guidance, Control, and Dynamics, Vol. 37, 4 (Mar. 2014), 1202-1210.
[22] Lux, A. Feasibility mission analysis: trajectory and performance study-Laplace mission on Ariane-5-ECA.
Report AE-NT-1-H-049-AE, Issue 1, Mar. 2006.

Transcription notes: reference [8] author list printed as "Martens, W., Boutonnet, W. A., and Schoenmaekers, J.";
the PDF's reference text is image-read, so minor punctuation may differ from the print. Of these, the in-system
Neptune trajectory works are [3] (another Triton tour), [16] and [17] (hyperbolic-flyby coverage, petal
strategy), [18]-[21] (v-infinity sphere and Tisserand-leveraging methods, applicable to moon tours), and [2];
Uranus: [6] (SEP to Uranus and Neptune) and [10].

## Suggested topology label
**pump-tour.** Justification from the text: the Triton tour is built from COT sequences designed on the
v-infinity sphere with "pump angle" and "crank angle" coordinates and uses the petal strategy of alternating long
and short transfers, which are v-infinity leveraging and flyby-sequence methods; it is a one-off non-repeating 55-
flyby sequence of a single moon, so it is not a repeating-encounter cycler (a plain "mga-tour" label would also
describe its non-repeating nature).

## Provenance of this digest
All 12 pages were read from the rendered PDF; Table 6 was transcribed from a 220 dpi crop in two halves. The
text layer was shift-decoded only for the term counts. Not checked: any other version of the paper, and the
volume and issue numbers (not printed on the pages).
