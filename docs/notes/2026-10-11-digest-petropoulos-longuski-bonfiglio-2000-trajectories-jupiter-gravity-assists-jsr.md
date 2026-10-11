# Digest: Petropoulos, Longuski and Bonfiglio 2000, "Trajectories to Jupiter via Gravity Assists from Venus, Earth, and Mars" (JSR 37(6):776-783)

**Digested:** 2026-10-11, owner upload, read in full from the text layer (two extractions). All six tables were parsed by
script from both extractions and read cell by cell against the page renders. Figures were viewed on the page renders.
**Citation:** Anastassios E. Petropoulos, James M. Longuski and Eugene P. Bonfiglio (Purdue University), Journal of
Spacecraft and Rockets 37(6):776-783 (November-December 2000), DOI 10.2514/2.3650. Presented as AIAA Paper 98-4284
(Boston, 10-12 August 1998). Received 19 August 1999, revised 1 March 2000, accepted 13 March 2000.
**Filed as** `cyclers_pdf/papers/petropoulos-longuski-bonfiglio-2000-trajectories-jupiter-gravity-assists-venus-earth-mars-jsr-37-776-doi-10.2514-2.3650.pdf`
(8 pages, text layer, md5 7d3f306548b76524db368e915299d120). Table transcription:
`data/sources/petropoulos-longuski-bonfiglio-2000-tables.yaml`.

## Verdict

A one-way ballistic tour paper. It enumerates Earth-launched gravity-assist paths to Jupiter, 1999 to 2031, and lists
the best case of each path and of each launch year. It holds no cycler, no periodic orbit and no closed sequence of
encounters. It matters to this project for two things: the path-pruning method that precedes a search, and the
statement and the data on resonant legs between repeated flybys of one body. It adds no collision with any catalogue
row, #942 or #943 candidate, or #973 cell.

## What the paper does

- **Problem.** Low launch-energy, low total DeltaV (impulsive model) trajectories to Jupiter, flight time up to about
  seven years, launches 1 January 1999 to 31 December 2031.
- **Tool.** STOUR (Satellite Tour Design Program), an automated patched-conic search over launch date and launch
  V-infinity that finds all ballistic gravity-assist trajectories for a path. It can add DeltaV in three modes:
  powered flyby, broken plane and V-infinity leveraging. One path over three decades takes several days of
  workstation time. MIDAS (JPL) optimises selected cases.
- **Path space.** Venus, Earth, Mars, up to four flybys: 120 sequences.
- **Pruning by gravity assist potential (GAP) plots.** For each path the maximum final aphelion against Earth launch
  V-infinity, assuming ballistic flight, tangential Earth launch, circular coplanar planets, no overturning, flybys
  down to 200 km altitude and phasing ignored. The V-infinity vector is turned positively if the next different planet
  is outward and negatively if inward. Table 2 ranks the 39 paths that reach Jupiter by minimum launch V-infinity
  (VEVE lowest, 2.53 km/s, then VEME 2.60 and VEEM 2.76). The direct-launch minimum is 8.79 km/s.
- **Selected set.** 25 paths (all classical ones and the reasonable nonclassical ones) were run in STOUR.
  - **Classical:** DeltaV-EGA, VEGA, VEEGA, V2 and V3 (Vn is n Venus flybys after Earth launch).
  - **Nonclassical examples:** VMV2, VEME, VMVE, VEVE, VEEM, VMVMV, ME, MM.
  - **No feasible trajectory with flyby altitude of 200 km or more:** VM, MV, VEV, MVM.
- **Findings.**
  - VEEGA, VEGA and V3 dominate. The best Earth-flyby trajectories beat the best no-Earth ones on almost every
    measure.
  - The DeltaV-EGA does not appear in the yearly table.
  - Nonclassical paths sometimes win. VMVE in 2021 reaches Jupiter in 3.6 years at launch V-infinity 4.00 km/s.
    The text compares this with a Hohmann transfer. VEME in 2023 beats the nearest VEEGA. VMV2 in 2012 beats any
    2012 V3. V2E in 2022 fills a year with no VEGA or VEEGA.
  - Mars as a leveraging body seldom helps VEEGAs. It is rare for V3 but very effective when it is available.
  - The Mars geometry recurs every 47 to 49 years. VMVE repeats in similar form about every 45 years.

## Tables (transcribed; three witnesses agree)

All in `data/sources/petropoulos-longuski-bonfiglio-2000-tables.yaml`.

| Table | Content | Rows |
|---|---|---|
| 1 | Legend for the launch-date plots (text) | 7 |
| 2 | GAP rank order of the 39 paths that reach Jupiter, minimum launch V-infinity | 39 |
| 3 | Best candidate for each of the 25 selected paths, unoptimised | 30 (26 cases, 4 with no feasible trajectory) |
| 4 | Promising trajectories with Earth flybys, by launch date | 44 |
| 5 | Promising trajectories without Earth flybys, by launch date | 42 |
| 6 | DeltaV-optimised versions to 2010 plus noteworthy later cases | 35 |

Row format in Tables 3 to 5: path, launch date (YY-MM-DD), launch V-infinity, signed flyby altitude below 200 km or
deep-space DeltaV with the flyby number, arrival V-infinity at Jupiter, time of flight in years, days of each flyby
after launch, with an r where a leg is resonant. Table 6 has launch V-infinity, deep-space DeltaV, arrival
V-infinity, time of flight and flyby days. The text layer prints infinity as "1" and minus as a stray glyph. Both were
read on the page render.

## Resonant and repeated encounters (the part that bears on the project)

**The paper's definition.** Sequential encounters of the same body are resonant if the time between them is an
integer multiple of the body's orbital period. Only same-body consecutive flybys carry the flag. No p:q with p greater
than one is used.

**What the tables show.** Every r in Tables 3 to 6 sits between two flybys of the same body (60 flagged legs, with
repeats across the tables). My own count of the leg lengths, from the printed flyby days against sidereal periods of
224.70 d (Venus), 365.26 d (Earth) and 686.98 d (Mars), is below. It is ours, not the paper's.

| Table | Earth legs, 2 and 3 years | Venus legs, 2 to 5 Venus years | Mars legs |
|---|---|---|---|
| 4 (Earth flybys, unoptimised) | 11 and 4 | none | none |
| 5 (no Earth flybys, unoptimised) | none | 10 (2), 6 (3), 3 (4), 2 (5) | none |

- **Unoptimised legs are exact.** Every flagged leg in Tables 3 to 5 is within 0.02 of an integer number of body
  periods (for example 731 d is 2.00 Earth years and 1096 d is 3.00).
- **Optimised legs are close, not exact.** In Table 6 three Venus legs differ from an integer by 0.02 to 0.08
  periods (VVV 2001-03-04, 2006-01-04 and 2008-12-01), because the optimiser moved the flybys.
- **Paper statements.** Most VEEGAs have about two years between the Earth flybys. Two (launches 1999-06-25 and
  2014-11-23) have a three-year resonance. Roughly half of the V3 cases have about two Venus years between the first
  two Venus flybys and about four between the last two, written [~2; ~4]. About half are [~2; ~3]. Two (2004 and
  2023) reach a five-Venus-year leg. Two (2006 and 2019) are [~3; =3], where = marks an exact resonance. These are
  among the worst listed V3 cases.
- **Alternating sequences.** VEVE is the best GAP path (rank 1, 2.53 km/s). Its listed cases are one-way. For
  example the 2026-08-27 case has Venus flybys at 176 and 1043 days and Earth flybys at 695 and 1590 days, with
  neither leg an integer number of periods. The paper does not claim a repeating Earth-Venus orbit.

## Relevance to this project

- **Two-working-body enumeration and the #973 screens.** This is the Purdue lineage of "which path is energy-feasible
  before phasing". The GAP plot is a Tisserand-style feasibility screen for one-way paths. The project's enumerators
  solve a different problem: periodic returns between two working bodies with a turn gate. The paper's screen is
  background for that gate, not a replacement, and it gives no turn-gate numbers. Strange and Longuski 2002 (held)
  states the Tisserand form of the same idea.
- **Repeated or resonant encounters.** It supports one fact for the Earth-Venus cells: ballistic Earth-launched paths
  use 2-year and 3-year Earth returns and 2 to 5 Venus-year Venus returns, all as same-body legs inside a longer
  path. It lists no closed Earth-Venus cycle, so it does not bear on the 8-year Earth-Venus class of #973.
- **No collision.** All 25 selected paths are one-way and end at Jupiter. No row has a start body equal to its end body
  and no sequence repeats. No catalogue row, no #942 or #943 candidate and no #973 member can match a row here.
- **Tool background.** STOUR with a five-day launch step and 0.25 to 1 km/s V-infinity step is the grid the project's
  own gauntlet resolution should be compared with, if a Jupiter-tour benchmark is ever wanted. Tables 4 to 6 are
  published numbers that a reproduction could target.

## Figures

- **Fig. 1 (page 777).** STOUR launch-date plot for low-energy VEEGAs with arrival V-infinity below 8 km/s. Path
  label 3 2 3 3 5, launch V-infinity list 3.00, 3.25, 3.50, 3.75 km/s (numerals 1 to 4), launch dates searched from
  1999-01-01 to 2032-01-01 in 5-day steps, maximum flight time 2600 days, time-of-flight axis 5.0 to 7.5 years.
  Numerals cluster in separate launch windows.
- **Fig. 2 (page 777).** V3 trajectories with two V-infinity-leveraging maneuvers (combined DeltaV up to 1.25 km/s).
  Path label 3 2 0 2 0 2 0 5, launch V-infinity 4.50 km/s, numerals give total deep-space DeltaV in 1 km/s bands.
- **Fig. 3 a to c (page 778).** GAP plots: final aphelion radius in AU (log axis, Jupiter line marked) against Earth
  launch V-infinity from 2 to 8 km/s. (a) 22 (V2), 23 (VE), 24 (VM) and 222 (V3), where 2 is Venus, 3 Earth and 4
  Mars. (b) 223 (VVE), 233 (VEE), 2233 (V2E2) and 2323 (VEVE). (c) 2334 (VE2M), 2343 (VEME), 2422 (VMV2) and 2423
  (VMVE).
- **Fig. 4 (page 778).** STOUR plot of VEME, path label 3 2 3 4 3 5, launch V-infinity 3.75, 4.00, 4.25, 4.50 km/s.
- **Fig. 5 (page 778).** STOUR plot of VMVE, path label 3 2 4 2 3 5, launch V-infinity 3.50 to 5.50 km/s, maximum
  flight time 2800 days, search minimum altitude -500 km. A single numeral 5 near 3.6 years in 2021
  is the short-flight case.

## Citation mining (references 1 to 25 against CORPUS_INDEX and the wanted list)

Held: Minovitch 1963 (JPL TR 32-464), Sohn 1964, Sturms 1967 (the later JPL report), Gillespie and Ross 1967,
Sims, Longuski and Staugler 1997, Diehl, Kaplan and Penzo 1983 (a different report from the cited Diehl and Myers
1987).
Already on the wanted list: Deerwester 1966 and Hollister 1964 (and, from Sims, Longuski and Staugler 1997,
Hollenbeck 1975 and Sims 1996).
Not held and not wanted before today (one new low-priority row added):
- Stancati, Friedlander and Bender (1976), AIAA 76-797, launch-opportunity classification of VEGA and DeltaV-EGA.
- Diehl, R. E. and Myers, M. R. (1987), JPL Publication D-4677, gravity-assist trajectories to the outer solar system.
- Sauer, C. G. (1989), "MIDAS", J. Astronaut. Sci. 37(3):251-259.
- Patel and Longuski (1993), AAS 93-682, and the Patel (1990) and Williams (1990) Purdue theses on the automated
  design that STOUR implements (the Longuski and Williams 1991 paper is in this same batch).
- Flandro (1966), Astronautica Acta 12(4):329-337, Jupiter flyby reconnaissance missions.
- Niehoff (1966), JSR 3(9):1351-1356, gravity-assisted trajectories to solar system targets.
- Broucke (1988), AIAA 88-4220, celestial mechanics of gravity assist.
Background not listed: Battin 1959 and 1987, Sedov 1960, Ross 1966 (already on the list under the Ross 1963 row),
Rinderle 1986 (JPL D-263 Galileo guide), Staugler 1996 (Purdue report).
