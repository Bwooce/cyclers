# Digest: Ehricke 1966, "Interplanetary Maneuvers in Manned Helionautical Missions" (AIAA chapter form of preprint 65-695; #960 wanted-list row 42)

K. A. Ehricke (Director, Advanced Studies, General Dynamics/Convair, San Diego; footnote: now Assistant Director, Astrionics Division, Autonetics, North American Aviation),
chapter in *Methods in Astrodynamics and Celestial Mechanics*, Progress in Astronautics and Aeronautics (R. L. Duncombe and V. G. Szebehely, eds.), AIAA, 1966, pp. 325-352,
DOI 10.2514/5.9781600864919.0325.0352. Crossref (2026-10-11): type book-chapter, container "Methods in Astrodynamics and Celestial Mechanics", issued 1966, pages 325-352; the book DOI is 10.2514/4.864919
(editors Duncombe and Szebehely; ISBN 9781563478666 and 9781600864919). The series volume number is not printed on the chapter pages and Crossref returns none; it is not stated here.
First presented as AIAA/ION Preprint 65-695 (the wanted-list identifier 10.2514/6.1965-695), Astrodynamics Specialist Conference, Monterey, California, 16-17 September 1965. Acknowledges B. Brown, P. Horio and W. Riddell.
Filed as `cyclers_pdf/papers/ehricke-1966-interplanetary-maneuvers-manned-helionautical-missions-aiaa-progress-chapter-doi-10.2514-5.9781600864919.0325.0352.pdf`
(md5 `7a623a631a5259b69fc26eae0f9406eb`, 30 pp., an AIAA ARC download with a text layer; pp.29-30 are blank; PDF page n = journal page n + 324; original upload not modified).
The year in the filename is the volume's (1966); the preprint is 1965. Copy checked: CORPUS_INDEX.md and `docs/notes/*.md` were grepped for "Ehricke", "65-695" and "AZM" on 2026-10-11 (no held copy; the only Ehricke hits are citations).

Status: READ in full (text layer; Table 1 on a 120 dpi page render and two 220 dpi crops; Eqs. 1 to 3 on crops; Figs. 1 to 9 viewed at 75 dpi, Fig. 3 at 300 dpi). Numbers: `data/sources/ehricke-1966-helionautical-maneuvers-tables.yaml`.
Check script and output: `cyclers_pdf/papers/...-check.py` and `...-check.out` beside the PDF. COMPUTED = our arithmetic 2026-10-11. INFERRED = our reading. Page numbers are journal pages.
Related digests, not repeated: `2026-10-05-digest-hollister-prussing-1965-optimum-transfer-mars-via-venus.md` and its RE-16 digest (the Venus-flyby optimum-impulse study of the same conference, same month),
`2026-10-07-digest-nasa-tm-x-53049-1964-symposium-manned-planetary-missions.md`.

## 0. Verdict

A mission-analysis study of what a mid-course heliocentric rocket burn does for a one-shot manned Mars round trip. Nothing repeats. The paper has one table (Table 1), nine figures and three equations, and no reference list.
It names three classes of heliocentric maneuver (apsidal, off-apsidal retro or acceleration, planet-approach), shows that the perihelion brake (PB) on the Mars-to-Earth leg is the strongest, and weighs the result in orbital departure weight (ODW)
for chemical, nuclear and solar-heat-exchange propulsion. It mentions Venus only for combined use of a heliocentric burn with a Venus powered fly-by (PFB). It contains no cycler, repeating-trajectory or periodic-orbit statement.
Nothing here bears on a literal collision with the project's cycler gates.

## 1. Content (READ)

- **Abstract and introduction (pp.325, 328-329).** Manned missions allow two classes of maneuver: controlled perturbation by planetary gravity fields (not treated) and heliocentric maneuvers (this paper).
  Three groups: apsidal maneuvers at or near perihelion or aphelion (perihelion brake PB, aphelion brake AB; Fig. 1a), planet-approach maneuvers 10 to 30 days before entering the planet's sphere of influence (Fig. 1b, the heliocentric
  Earth approach retro HEAR), and retro or acceleration maneuvers in intermediate parts of the transfer orbit (off-perihelion). Objectives: lower mission energy, a wider choice of propulsion (low thrust, high specific impulse), wider departure windows.
  "Exchange ratio" = delta-V saved at the planetocentric maneuver per unit heliocentric delta-V: PB 1.5 to over 2; retro or acceleration 1 to about 1.6 (best 35 to 65 degrees from an apsis); planet-approach less than 1 to about 1.2.
  Because the heliocentric burn can be long, the same total delta-V can be placed at the highest specific impulse available, giving multi-elliptic transfer paths.
- **Why PB on the return from an outer planet (p.329).** The heavy outbound vehicle must take the cheap outbound corridor, which forces a poorer return corridor with perihelion at or below 0.65 AU and a steep crossing of Earth's orbit (high hyperbolic excess).
  Braking near perihelion lowers the aphelion distance and eccentricity.
- **Eqs. 1 to 3 (pp.329-330).** (1) da/dV = (R_P/2)[V/(K/R_P)](n+1)^2; (2) dR_A/dV = R_P(n+1)[2V/(K/R_P) + 2/(V(e + cos eta))]; (3) de/dV = 2(e + cos eta)/V; n = R_A/R_P, eta true anomaly, K the
  gravitational parameter. Fig. 2 plots the three partials against true anomaly for n = 3, 6, 9 and R_P = 0.4 AU. **Checks (ours):** Eq. 1 equals 2 a^2 V/K and matches the Fig. 2 left panel at eta = 0; Eq. 3 is the Gauss equation for
  tangential acceleration; Eq. 2 as printed does NOT reproduce the exact perihelion result (4 a^2 V/K) and our derivation R_P(n+1)[nV/(K/R_P) + (e + cos eta)/V]; Fig. 2 middle panel agrees with neither (data file, `equations`).
  Do not use Eq. 2 or the Fig. 2 middle panel as a numeric control.
- **HEAR (pp.330, 335-336, Fig. 6).** Reducing the hyperbolic excess by delta v_inf at the planet-approach point costs more than the equivalent periapsis reduction by the factor v_p/v_inf; the advantage is time: a low-thrust,
  higher-I_sp drive can be used. Results typical for several Mars years (Fig. 6a, 1979/80; 6b): exchange ratios 0.74 and 0.86 for HEAR started 10 d and 180 d before arrival.
- **Search method (p.331).** IBM 7094-class machine (printed "IBM 709^"), three-dimensional two-point heliocentric transfer orbits by the Gauss-Lambert method: first a mono-elliptic path of given flight time T2 (Mars to Earth) or T1,
  then the heliocentric maneuver changes the path.
- **Mono-elliptic Mars reference missions (pp.331-332, Fig. 3).** Capture orbits n = 1, r* = 1.3 at Mars; Earth departure r* = 1.1, n = 1. Reference missions 1973 to 1990 chosen at low-ODW windows for solid-core nuclear and chemical vehicles with Earth entry at no more
  than 50,000 fps (15.3 km/s). The energy variation among missions is mainly the Earth return velocity. Fig. 3 labels (data file): nine Earth departure dates and mission periods of 420 to 450 d (for example 8-16-73, 180 + 30 + 220 = 430 d).
- **PB results (pp.332-334, Figs. 4a to 4d, 5, Table 1).** The highest exchange ratio occurs when PB is at T2 - 50 d and the flight time is extended by 30 to 50 d (PB-to-Earth time 80 to 100 d); PB velocity 11,000 to 15,000 fps against 13,000 to 29,000 fps
  for a geocentric retro maneuver for the same entry speed. Exchange ratios mostly 1.5 to 2.0, except 1984 to 1988. For 1986 the mono-elliptic T2 = 240 d would give ratios below one (excessive plane change); T2 = 220 d is the compromise
  (T2 = 200 d would need a 22,810 fps Mars departure maneuver). PB takes place shortly after perihelion passage at small true anomaly in the seventies; the PB distance R_PB is lower in the seventies than in the eighties.
  **Table 1** gives the change of orbital parameters for the nine reference return transfers, two flight-time extensions each (data file).
- **ODW effects (pp.336-341, Figs. 7 to 9).** Propulsion letters C (chemical, 450 s), S (solar heat exchange, 700 s), G (graphite solid core, 800 s), P (nuclear pulse, for Jupiter missions). Payload per maneuver and entry speeds in the data file (`text_numbers`).
  A chemical vehicle with PB by the solar drive (C2S) beats a nuclear-departure vehicle with mono-elliptic return (CG) throughout the seventies (Fig. 7); PB smooths the ODW variation in unfavourable years for a nuclear vehicle but is of
  little use in the eighties. For missions ending in an elliptic Earth capture orbit with solid-core engines (Fig. 8), PB lowers ODW in the early eighties and raises it slightly in 1986 (new engines each maneuver, so a weight penalty).
  HEAR (Fig. 9): reduces ODW below mono-elliptic for a chemical vehicle but is mostly inferior to PB, comparable only if the HEAR starts 30 d before arrival.
- **Combination with a planetary flyby (p.338).** A solar-heat-exchange heliocentric burn before or after a Venus encounter on the Mars return could reduce or remove the Venus-centred powered fly-by (PFB) delta-V while keeping a wide Mars departure window;
  fusing the two (modifying the Venus encounter hyperbola with the solar drive) is suggested for study. No numbers.
- **Summary and conclusions (pp.341-342).** Heliocentric maneuvers broaden the propulsion choice, often reduce mission velocity (best at perihelion on the outer-planet return), can reduce ODW or widen windows,
  and are worth including in trade studies; their effect is usually but not always advantageous. No reference list follows.

## 2. Table 1 (transcribed; three witnesses)

18 rows (9 Mars departure dates, delta T = 30 and 50 d), 14 columns. The data file holds every cell. Two internal consistency checks (COMPUTED, `-check.out`):
(a) column 3 = column 2 minus 50 d in all 9 rows; (b) the true-anomaly identity column 14 - column 7 = column 12 holds exactly in 11 of 18 rows and within 0.2 degree in all 18 (rounding or typesetting scatter in the original; the render shows each cell as transcribed).
The text layer has damage in several cells; the render decides.

## 3. Figures (one line each)

1a, 1b PB and HEAR sketches (heliocentric orbits of Earth, Mars and the transfer with the burn point; viewed). 2 partials of a, R_A and e against true anomaly, three ellipses, R_P = 0.4 AU (viewed, approximate). 3 mono-elliptic
reference-mission bar chart (labels transcribed; bar heights about 7 to 72 thousand fps, not read to a stated precision). 4a to 4d delta-V and entry velocity against added flight time for the nine launch opportunities (viewed; curves not read).
5 characteristic PB data, per mission year (horizontal bars; viewed). 6a, 6b HEAR velocity change against Earth entry velocity for several initiation times (viewed). 7 relative orbital departure weight by year for several propulsion combinations
(chemical C4 about 4 to 9 in relative units near 1975-77; viewed). 8 ODW for solid-core capture missions, with and without PB, five years (viewed). 9 ODW comparison, four years, chemical and solar propulsion (viewed).

## 4. Relevance to this project

- **Heliocentric-maneuver classes.** The three classes (apsidal PB and AB, off-apsidal retro or acceleration, planet-approach HEAR) are impulsive deep-space burns applied to a single Earth-Mars-Earth round trip. They are
  not gravity-assist or resonance constructions. The project's cycler definitions (repeating ballistic orbit with encounters) do not include them; they are relevant only as the pre-1966 literature on adding a deep-space delta-V to
  a one-shot Mars round trip, comparable to the deep-space-maneuver (DSM) content of later Earth-Mars studies.
- **Round-trip or repeating-trajectory statements.** "Round trip" here is one Earth departure, a Mars stopover (capture time 30 to 50 d) and one return, total 420 to 450 d (Fig. 3). No repeat, no synodic-period recurrence, no periodic
  orbit, no cycler. The reference missions are one per opportunity, 1973 to 1990.
- **Venus.** Only the PFB remark (p.338). The Venus-swingby Mars round trips are in Hollister and Prussing (RE-16 and AIAA 65-700), presented at the same conference (AIAA/ION Astrodynamics Specialist Conference, Monterey, 16-17 September 1965; papers 65-700 and 65-695); this chapter adds no Venus numbers.
- **Relation to Ehricke's 1959 AZM-072 report (wanted-list row 40).** AZM-072 is not cited here: the chapter has no reference list and no note names it. Nothing in the text quotes or reproduces it. The row stays open; whether AZM-072 or its
  published part (ASME 60-AV-1) contains related PB or aphelion-brake ideas cannot be told from this chapter.
- **Gate.** There is no gate anchor for this paper in literature_check.py. If one is ever added, the content supports {transfer} (a one-shot Earth-Mars-Earth round trip with heliocentric burns), no cycler, repeated-moon or resonant label.

## 5. Cautions

- The perihelion-brake exchange ratio is the paper's own measure; its definition (delta-V saved at the planetocentric maneuver over the heliocentric delta-V) depends on the fixed entry velocity and the 3-D geometry; the text says it degrades with plane change.
- Eq. 2 and Fig. 2 middle panel are unreliable (sec. 1).
- Fig. 3 has two labelling errors (data file `fig_3.checks`): the 1975 group's Julian date Julian date (printed 2691.5, which is 10-6-75) is one day later than its printed calendar date 10-5, and the Julian date 5740.5 printed on the 1986 group belongs to the 1984 group
  (2-10-84); the 1986 group (3-31-86) should read 6520.5. The calendar dates themselves reproduce Table 1's Mars departure dates when the transfer times are added in all nine cases.
- The "time anomaly" wording of Table 1's footnotes 11 and 13 is a slip for true anomaly (the units are degrees).
- Reading from the text layer alone would have produced several wrong table cells.

## 6. Citation mining

The chapter has no reference list and cites no other work by name (NERVA II is named as an engine). Nothing to mine. No new wanted-list rows. The related wanted-list rows are unchanged: row 40 (Ehricke 1959 AZM-072, still open),
the Deerwester and Ross rows (the Venus-swingby and nonstop round-trip context, from other digests).
