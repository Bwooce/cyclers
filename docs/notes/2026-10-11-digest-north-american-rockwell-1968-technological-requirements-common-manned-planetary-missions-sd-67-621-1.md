# Digest: North American Rockwell 1968, "Technological Requirements Common to Manned Planetary Missions, Technical Summary", SD 67-621-1 (corpus policy; map digest)

North American Rockwell Corporation, Space Division, Downey, California. Report SD 67-621-1, "Technological
Requirements Common to Manned Planetary Missions, Technical Summary", contract NAS2-3918, January 1968,
NTRS 19680009382. The title page names R. D. Meston (Project Engineer) as preparer and A. Codik (Project
Manager) as approver. The abstract sheet gives the publication date as JAN68.
- Filed as `cyclers_pdf/papers/north-american-rockwell-1968-technological-requirements-common-manned-planetary-missions-technical-summary-sd-67-621-1-ntrs-19680009382.pdf`
  (the NTRS scan as received, 120 PDF pages, md5 2776c48ecd05d0b7eff94d09b895948e). The OCR copy is the same
  name with `-ocr.pdf`.
- The scan has a text layer, but it is poor. It is readable on the text pages and useless on the rotated table
  pages (tables 1-16). The OCR copy was made with `ocrmypdf --force-ocr --rotate-pages --oversample 250`
  and recompressed with Ghostscript (161 MB down to 23 MB, text kept).
- This is volume 1 of 5. The detailed volumes are SD 67-621-2 (Appendix A, Mission Requirements), -3
  (Environments), -4 (Subsystem Synthesis and Parametric Analysis) and -5 (System Synthesis and Parametric
  Analysis). The project holds none of them.
- Tables transcribed to `data/sources/north-american-rockwell-1968-sd-67-621-1-tables.yaml` (Tables 1-18 and
  21-30; Tables 19 and 20 hold words only).
- Corpus check before filing: grep of CORPUS_INDEX.md, the wanted list and the digests for 19680009382,
  "SD 67-621", "Technological Requirements Common" and NAS2-3918 found nothing.

## 0. Verdict

**A systems study of one-shot manned planetary round trips. No cycler, no repeating trajectory, no free
return as such. It carries a useful set of Venus-swingby Earth-Mars round-trip baselines in Tables 13-16,
the same E-V-M-E and E-M-V-E class as Hollister and Prussing 1965 and Young 1966.**
- Mission modes studied: direct stopover, Venus swingby (outbound and inbound), and flyby. Flybys of Mars and
  Venus are not studied.
- Targets: Mercury, Venus, Mars, Jupiter, Ganymede, and the asteroids Vesta and Ceres. Opportunities cover
  1980 to 2000.
- Trajectory content is in Tables 1-16 only. The remaining 14 tables and almost all of the text are about
  modules, subsystems and masses.
- It has no reference list. The only outside sources named are unpublished NASA data and one contemporary
  contract study (sec. 7).
- It carries its own errata sheet on the last page (sec. 2).

## 1. Map of the report (printed page numbers; PDF page = printed page + 7 up to p.81, then + 6 and later + 5)

| Section (printed p.) | Content | Trajectory classes and key statements |
|---|---|---|
| Introduction (1-3) | Aim: find requirements common to near-term (Mars and Venus stopover, early to mid 1980s) and advanced missions | Constraints: high-thrust propulsion only; circular parking orbits in the nine-month base study, elliptical orbits in a three-month amendment; no abort or launch-window analysis; crew 3 to 20 |
| Mission requirements: performance (5-37) | Mission-opportunity and baseline selection, Tables 1-16, Figs. 1-13 | Opposition-class stopovers only. Conjunction-class (about 1000 days) is named and not studied (p.6, p.12). Two-plane transfers for heliocentric angles near 180 deg. Swingby modes for Mercury and Mars. Entry-speed cap 19.8 km/s (65,000 ft/s) |
| Aerobraking technology (38-40) | Heatshield weight fractions, Figs. 14-15 | Mars entry 6.1 to 9.2 km/s: heatshield 6.8 to 14.6 percent of entry weight. Venus 9.2 to 15.2 km/s: 12.7 to 40 percent |
| Planetary excursion module (41-45) | Descent and ascent delta-V, Table 17, Figs. 16-18; Ganymede orbit-insertion guidance | Hohmann-type descent and ascent for elliptical parking orbits |
| Environments (47-60) | Meteoroid, thermal, radiation; Tables 18-19, Figs. 19-20 | Two-plane transfer over the asteroid belt keeps shield weight reasonable (p.113) |
| Subsystem synthesis (61-81) | Life support, communications, propulsion, electrical power; Tables 20-26, Figs. 21-26 | Not trajectory content |
| System synthesis (83-110) | Weight scaling, module masses, common modules; Tables 27-30, Figs. 27-37 | Midcourse allowance 60 m/s per direct return leg and 120 m/s per swingby return leg (p.93) |
| Conclusions (111-114) | Commonality at module and subsystem level | Mars missions limited to the Venus-swingby mode would need only the Apollo-class reentry module (p.104, p.111) |
| Errata (last page) | Five corrections | See sec. 2 |

Mission modes by table:

| Mode | Tables | Targets and notes |
|---|---|---|
| Flyby (one pass, no stopover) | 1 | Vesta, Ceres, Jupiter; two opportunities each |
| Direct stopover | 2-10 | Mercury; Venus aerobraker and retrobraker; Mars aerobraker and retrobraker; Vesta; Ceres; Jupiter; Ganymede. Three opportunity years each, stay times 0, 30 and 60 days |
| Venus swingby, outbound to Mercury | 11 | Earth, Venus, Mercury, then direct return; the swingby adds delta-V (printed column dVswing) |
| Venus swingby, inbound from Mercury | 12 | Earth, Mercury, then Venus swingby on the way home |
| Venus swingby, outbound to Mars | 13 (aerobraker), 14 (retrobraker) | Earth, Venus swingby, Mars stopover, direct return |
| Venus swingby, inbound from Mars | 15 (aerobraker), 16 (retrobraker) | Earth, Mars stopover, Venus swingby, Earth |

Trajectory classes asked about:
- One-shot manned round trips: all of them.
- Conjunction and opposition class: named; only the opposition class was studied.
- Venus swingby: Tables 11-16, with the Venus passage distance given in Venus radii for the Mars cases.
- Free return, periodic or repeating trajectory: **not treated.** The word "repeat" appears once, for the
  statement that performance requirements roughly repeat each synodic cycle (p.6).

## 2. Tables

All in `data/sources/north-american-rockwell-1968-sd-67-621-1-tables.yaml`.
- **Tables 1-10** (printed pp.14-22): the direct and flyby baselines. Per opportunity year there are three
  stay times. Columns carry Julian dates, trip and stay times, delta-V increments, entry speeds, and totals.
- **Tables 11-12** (pp.23-24): Mercury swingby baselines, each in two sub-tables (outbound and return legs).
- **Tables 13-16** (pp.25-28): Mars Venus-swingby baselines. Columns: Earth departure JD and delta-V1, leg
  times, Venus swingby JD and passage altitude hp in Venus radii, Mars arrival and departure, Earth arrival
  JD and entry speed, total time and total delta-V. Velocities are in km/s.
- **Tables 17-30**: planetary excursion module delta-V (17), meteoroid flux models (18), life-support scaling
  equations (21), communications parameters (22), liquid propellants (23), power-subsystem matrices
  (24-26), design summary (27), module masses (28), common-propulsion grids (29-30). Tables 19 and 20 hold
  words only.
- **Errata sheet** (last page, PDF 120): page 2 constraint sentence; page 56 t in weeks; page 107 delete the
  Ganymede exclusion; Table 29 EOE flyby entries on three lines; Table 30 add an F under the 100,000 kg
  chemical module for the Vesta, Ceres and Jupiter flybys and set the Vesta flyby nuclear module mass to
  75,000 kg, not 1,200,000 kg. The yaml keeps the values as printed and carries the errata separately.

Printed rows that fail arithmetic closure (dates adding up, delta-V columns adding up) are listed in the yaml
under `closure_exceptions`. They are not claimed as errors. Rounding, midcourse allowances and trip times
rounded to 5 days could explain several of them. The notable ones:

| Row | Printed | Closure |
|---|---|---|
| Table 3, 1990, 60-day stay | Arrive Venus 7990 | 7830 + 100 = 7930 |
| Table 7, 1985, 6070 | dvtp 9.82; arrival 6410 | dv1 + dv2 + dv3 = 9.92; 6070 + 330 = 6400 |
| Table 8, 1980, 4550 | Trip 300, arrival 4880, total time 745 | Closes if the trip time is 330 |
| Table 8, 1989, 7805, 60-day stay | Arrive Mars 8140 | 7805 + 310 = 8115 |
| Table 11 outbound, 1988 | dvttp 14.50 | 3.96 + 2.32 + 0.0 + 7.27 = 13.55; the mission total 21.46 matches 14.50 |

## 3. Figures (37 plates; described from the list of illustrations, not digitised)

1-5. Trans-planet, total trans-planet, total trans-Earth velocity requirements and contour overlays for the
   1990 Mars opposition, as functions of arrival or launch date (JD - 2400000) and transfer time 100 to 500
   days. These are the plots from which the baseline missions were picked.
6-13. Incremental velocity requirements for elliptical planetary orbits: Mercury retrobraker; Venus aerobraker
   escape and retrobraker; Mars aerobraker escape and retrobraker; Jupiter orbit insertion and escape (1990
   mission); Ganymede retrobraker.
14-15. Mars and Venus aerodynamic-braking heatshield weight fraction against entry velocity.
16-18. Planetary excursion module characteristic velocity against parking-orbit eccentricity (Mercury, Mars,
   Ganymede).
19-20. Mission-module radiation shielding requirements and added shielding (1985 to 1989 mission).
21-26. Life-support mass comparison, interplanetary communications comparison, candidate power systems, power
   system weights, radiator area.
27-37. Reentry-module, mission-module and excursion-module masses; trans-Earth mass; nuclear and chemical
   propulsion-module masses.
- The contour plots (Figs. 1-5) were seen on the page renders only as an overlay of date axes and numbered
  contour lines. I did not digitise any of them. Their trip-time legend (1 to K, 100 to 500 days in steps of
  20) is printed on each plate.

## 4. Relevance to this project

- **Cycler, periodic, repeating, free-return: none.** No catalogue cycler row gets a source, a number or a
  check from this report.
- **What bears on the Earth-Mars Venus-swingby lineage.** Tables 13-16 are Venus-swingby Earth-Mars round trips
  (E-V-M-E outbound, E-M-V-E inbound) of the Hollister and Prussing 1965 type, with powered terms listed and a
  Venus passage altitude. They are an independent 1968 industry set from NASA Mission Analysis Division data.
  - Outbound opportunities in the table: first-row Earth departure JD 6143 (Table 13), 8509 and 10835, which
    are 2366 and 2326 days apart (about 6.4 years).
  - Inbound opportunities: 4940, 7354 and 9660 days (Table 15), 2414 and 2306 days apart.
  - The 6.3 to 6.6 year spacing is our own subtraction of printed JDs. It is the same 6.4-year Earth-Venus-Mars
    repeat stated by Young 1966 (digested 2026-10-06).
  - Mission durations are 545 to 692 days with Venus passage altitudes of 0.11 to 1.55 Venus radii as printed.
    Total delta-V in the nine-row sets: 5.87 to 8.10 km/s aerobraker (Tables 13 and 15), 9.62 to 13.58 km/s retrobraker (Tables 14 and 16).
- **Cross-reference: NASA SP-35 Part 6 (Douglas, 1968), digest 2026-10-11 by corpus-titus-sonnet.** Same
  opportunities, different source. Part 6 (NTRS 19680020626) treats Mars stopover missions using Venus
  swingbys, with contour charts for 16 opportunities from 1975 to 1999 and a stated 2338-day (about 6.4-year)
  repeat, and its trajectory data were verified by the NASA/OART Mission Analysis Division. This report names
  the same division as the source of its opportunities (p.5). Tables 13-16 here give point baselines
  (dates, delta-V, Venus passage altitude) for the opportunities that Part 6 charts. A date-by-date match
  between the two was not made. The 2338-day figure is consistent with the 2326 to 2366 day outbound
  spacing found above.
- **Hollister lineage.** The report does not cite Hollister or Minovitch. It is a mission-system study built on
  NASA data. It is a possible V0-level data set for the Hollister-Prussing E-V-M family. It is not a
  replacement for the primary papers and has no constraint check of its own.
- **Earth-Venus rows.** The Venus direct stopovers (Tables 3 and 4: 100 to 200 day outbound legs, 270 to 345
  day returns, delta-V 7.44 to 11.97 km/s) are conventional round trips. They do not touch any Earth-Venus
  cycler row.
- **Mercury via Venus (Tables 11-12).** These are Venus-swingby Mercury missions of the Cutting and Sturms
  type, related to the VanderVeen 1969 Bellcomm memo (digested 2026-10-06) and the catalogue rows citing it.
  I did not compare dates.
- **Value as a negative result.** The study ranks the Venus-swingby mode as the more attractive Mars mode on
  entry speed and propulsion (p.104). It never considers a repeating trajectory, which supports the project's
  view that cycler ideas were outside the 1967 mission-study toolbox.
- **Data-quality warning for any future use.** Several printed rows do not close (sec. 2). Do not use a row as
  a golden value without its closure check.

## 5. Notes on the text

- The text gives the baseline selection rule: minimise the sum of the incremental velocities for the chosen
  stay time, and shows that minimum mass in Earth orbit gives essentially the same mission (p.12). So the
  trajectories are delta-V optimal, not mass optimal.
- Stay times of 0, 30 and 60 days were evaluated for every case except Mercury (p.13).
- The Jupiter and Ganymede cases were picked to minimise entry speed, not the delta-V sum (p.13).
- The report carries no calendar dates. All dates are JD minus 2440000 as printed.

## 6. OCR and three-witness record

- Witnesses: (1) the NTRS text layer; (2) a new tesseract layer (250 dpi, `ocrmypdf --force-ocr
  --rotate-pages`); (3) Apple Vision on 200 dpi renders rotated upright, for the table pages; (4) every number
  read by eye on the 200 dpi renders, with 300 dpi crops where the scan is damaged.
- The NTRS layer is useless on pages 21-35 (tables 1-16). On those pages the two usable witnesses are the new
  OCR and Vision. A text-presence test (does the printed number appear anywhere on the page in the witness)
  gave, for the 1934 numeric entries of Tables 1-16:

  | Found in | Entries |
  |---|---|
  | Both the new OCR and Vision | 1714 |
  | Only the new OCR | 115 |
  | Only Vision | 71 |
  | Neither (plus the six target names) | 28 numbers |

  This test is weak: it checks presence on the page, not position. The decisive check was arithmetic closure
  and the eye reading. The 28 numbers found in neither witness are in damaged tables 3, 4, 6, 7, 8, 9, 13
  and 16 and were each read on the render.
- Tables 17-18 and 21-28: about 410 numeric entries, of which 160 appear in both text layers, 202 in one, and
  51 in neither (some of those 51 are page numbers and exponents, not table entries). Tables 22 and 23 carry
  most of the misses because of their dense or rotated layout. All were read on the render.
- Damaged scan regions: Table 7 (last column and the entry speed in the 1987 rows), Table 8 (the last
  column, 1980 and 1991 rows, and several time cells). The yaml notes list the affected cells. Values there
  were fixed by closure, and in each case agree with the visible digits.
- Table 1 closed on every row (arrival = departure + first leg; total time = sum of the legs).
- Not three-witnessed: the text of the report beyond the tables was read from the NTRS layer and the new OCR
  layer only. No figure was digitised.

## 7. Citation mining

The report has no reference list or bibliography. Outside sources are named only in the text:

| Item as printed | Status | Why it matters |
|---|---|---|
| "NASA/OART Mission Analysis Division" data on mission opportunities (p.5) | Not a citable item; no report number | The source of the baseline opportunities in Tables 1-16 |
| "Technology Requirements for Atmosphere Braking to Orbit About Mars and Venus", NASA contract NAS2-4135, due January 1968 (p.41) | Neither (not on the wanted list or in the index) | Aerobraker companion study for Mars and Venus; no trajectory content expected |
| SD 67-621-2, "Appendix A - Mission Requirements" (same series) | Neither | The detailed mission and trajectory requirements behind Tables 1-16; most likely to hold per-leg data |
| SD 67-621-3, -4, -5 (Environments; Subsystem Synthesis; System Synthesis) | Neither | Subsystem and mass detail; low value for the project |

### Proposed wanted rows (the lead decides which go on the list)

| Priority | Item | Identifier |
|---|---|---|
| 1 | North American Rockwell, "Technological Requirements Common to Manned Planetary Missions", Appendix A: Mission Requirements, SD 67-621-2, contract NAS2-3918 (1968) | Report number as printed; NTRS not checked |
| 2 | North American Rockwell, "Technology Requirements for Atmosphere Braking to Orbit About Mars and Venus", contract NAS2-4135 (1968) | Contract number as printed; report number not given |
| 3 | SD 67-621-3, -4, -5 | Report numbers as printed; low priority |

Catalogue implications: none.
