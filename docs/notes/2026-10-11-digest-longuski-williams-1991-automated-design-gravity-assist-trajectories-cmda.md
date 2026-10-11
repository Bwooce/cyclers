# Digest: Longuski and Williams 1991, "Automated design of gravity-assist trajectories to Mars and the outer planets" (CMDA 52:207-220)

**Digested:** 2026-10-11, owner upload, read in full. The PDF is a 300 dpi bilevel scan with a hidden OCR layer from
the source, so every table and figure was read on the page images and the layer was used as a digit check only.
**Citation:** James M. Longuski (Assistant Professor) and Steve N. Williams (graduate student), School of Aeronautics
and Astronautics, Purdue University, Celestial Mechanics and Dynamical Astronomy 52:207-220 (1991), volume and pages
confirmed from the running heads and the first page (issue 3 per Petropoulos et al. 2000, reference 17). Received
7 August 1990, revised 25 March 1991, accepted 5 August 1991. The DOI 10.1007/BF00048484 is not printed on the scan; it
is the one recorded as confirmed on the wanted list.
**Filed as** `cyclers_pdf/papers/longuski-williams-1991-automated-design-gravity-assist-trajectories-mars-outer-planets-cmda-52-207-doi-10.1007-BF00048484.pdf`
(14 pages, scan with OCR layer, md5 4223bc81123b179ffd1c36c703c73337). Table transcription:
`data/sources/longuski-williams-1991-tables.yaml`.

## Verdict

The paper that introduced the Purdue automated conic search (the algorithm of Williams' 1990 thesis, later STOUR) with
two applications: outer planet grand tours and Venus gravity assist to Mars. It holds no cycler and no periodic orbit.
It is lineage for the project's enumerators, because its search is organised by target-planet revolutions, which is a
bookkeeping by integer multiples of the target period. It adds no collision with any catalogue row, #942 or #943
candidate, or #973 cell.

## What the paper does

- **Vis viva matching (section 2).** For a given departure body, date and V-infinity, solve Lambert's problem to the
  next target for a set of discrete flight times and a chosen number of spacecraft revolutions. Each solution implies a
  departure vis viva (C3 = V-infinity dot V-infinity). Plot vis viva against flight time (Fig. 2). Solutions are the
  intersections of the curve with the target vis viva.
- **Target revolutions.** The designer sets the number of target-planet revolutions, which fixes the flight-time range
  to a band between whole periods of the target: 0 revolutions is 0 to 1 period, 1 revolution is 1 to 2 periods. For
  Venus (225 days) the bands are 0 to 225 days, 225 to 450 days, and so on. The band is divided into equal steps.
- **Existence test.** The flyby exists if the required bending angle between the incoming and outgoing V-infinity
  vectors is below delta = 2 asin[1/(1 + Rp V-infinity squared over mu)] for the minimum flyby radius Rp.
- **Search tree (section 3, Fig. 3, Table I).** For the Voyager II launch (20 August 1977, V-infinity 10.2 km/s) all
  Lambert solutions on Earth-Jupiter-Saturn-Uranus-Neptune-Pluto are found. Each node has a four-digit code: target
  planet, spacecraft revolutions, target revolutions, orbit number. Two Earth-Jupiter solutions exist (5001 and
  5002), each with a large tree. The Voyager grand tour reaches Neptune and not Pluto within 20 years.
- **Future grand tours (section 4).** Paths E56789, E5679, E5689, E56879, E5789, E579, E589 and E59 were searched from
  1985 to 2050, at V-infinity 12 km/s and 20-day launch steps.
  - **Classical Voyager tour:** the next opportunity is 9 June 2153, which is 175.8 years after Voyager II. The paper
    ties the repeat to the Uranus-Neptune synodic period, printed as 172 years.
  - **Four-planet tour:** E5789 exists in 1995. The last launch with Pluto is in 1996, with a minimum flight time of
    about 7400 days (20 years) at V-infinity 14 km/s. None to Pluto in 1997.
  - **Fast trips via Jupiter:** Jupiter-Uranus every 14 years (1500 days), Jupiter-Neptune and Jupiter-Pluto every 13
    years (2500 days). Table II lists the three-planet opportunities to 2037.
- **Venus assist to Mars (section 5).** Ballistic and VEGA (Venus-Earth gravity assist, E234) trajectories, 1995 to
  2024, plotted below 12.25 km2/s2 launch vis viva, with one-Venus-flyby swingbys (E24) found only in 2015. Lowest launch
  vis viva: ballistic 7.75 km2/s2 in 2018, VEGA 7.59 km2/s2 in 2010. VEGA is clearly below nearby ballistics in 2002,
  2008, 2010, 2020, 2021 and 2023, and there is no mode between 2004 and 2008. The minimum launch vis viva of the Venus
  swingby (11) is above the VEGA (9) in the same period.

## Tables (transcribed; layer and page image agree)

In `data/sources/longuski-williams-1991-tables.yaml`.
- **Table I.** Four solution blocks for the 1977 launch (5 + 3 + 3 + 5 rows): event, code, planet, flyby altitude, flyby
  angle, V-infinity, orbit period, periapsis and apoapsis, C3, leg flight time and cumulative flight time. The page prints
  it rotated by 90 degrees. Two of the blocks are Jupiter dead ends with flyby altitudes of 1012867 km and 783550 km,
  and a Saturn arrival V-infinity of 5.55 and 5.23 km/s.
- **Table II.** Eleven far-outer-planet opportunity rows, 1993 to 2037 and the periodic fast trips.
- The printed cumulative flight times differ from the sums of the leg times by 0.1 day in several rows, consistent with
  rounding (my own addition check).

## Repeating and resonant content

None of the cycler kind. Three items touch repetition and are all one-way:
- **Target-revolution bands** (above) are the same idea as the project's period classes in the two-working-body
  enumeration: a leg length is an integer number of target periods plus a fraction.
- **The 175.8-year and 172-year repeat** of the grand tour is a planetary-alignment repeat set by a synodic period.
- **The 13-year and 14-year repeat** of the Jupiter fast trips. My own check, not the paper's: the Jupiter synodic
  periods are about 13.8 years with Uranus, 12.8 with Neptune and 12.5 with Pluto, which match the printed 14, 13
  and 13.

## Relevance to this project

- **Lineage of the enumerators.** This is the earliest paper in the corpus that organises a conic search by target
  revolutions and checks every flyby by the bending-angle test. The project's two-working-body generator and turn gate
  solve the periodic version of the same problem. Petropoulos et al. 2000 (digest of today) uses the program that grew
  from it, and Strange and Longuski 2002 (held) gives the graphical form.
- **No collision.** The paper has no closed sequence. No catalogue row, #942 or #943 candidate, or #973 member can match
  Table I or II.
- **Reproduction target.** Table I is a published Lambert and flyby chain for Voyager II with every leg time, period and
  flyby condition. It could serve as a positive control for the project's Lambert and turn-gate code, if one is wanted.
  Block 1 is the actual Voyager II geometry. The values are printed to 0.01 to 0.1 and the paper gives no ephemeris
  model, so a check against them needs a tolerance, not exact agreement.

## Figures

- **Fig. 1 (page 208).** An Earth to Venus search: the Earth orbit with R1, the Venus orbit with a fan of target
  position vectors R2 at equal flight-time steps, and V-infinity at Earth.
- **Fig. 2 (page 210).** Vis viva against flight time (0 to 225 days, log axis 1 to 10000 km2/s2) for the first leg of
  Mariner Venus-Mercury: a Type 1 curve falling through the target vis viva line near 125 to 150 days and a Type 2 curve
  that crosses it twice, near 150 days and again beyond 200 days.
- **Fig. 3 (page 213).** The solution tree for 20 August 1977: Earth, then Jupiter, Saturn, Uranus, Neptune and Pluto
  nodes with four-digit codes and arrival dates in parentheses. Solid lines have cumulative time below 20 years, dashed
  lines above. The labelled arrival dates were not transcribed.
- **Figs. 4 and 5 (pages 214 and 215).** Positions of the outer planets in the ecliptic (axes in AU, vernal equinox of
  1950) for two ten-year spans, 1977 to 1986 and 2153 to 2162, 176 years apart. In the second span Pluto is far from the
  Voyager II geometry.
- **Fig. 6 (page 216).** Flight time (1500 to 7500 days) against launch date, 1984 to 2050, for Jupiter-Uranus,
  Jupiter-Uranus-Neptune and Jupiter-Uranus-Neptune-Pluto. One filled square marks the 1995 four-planet case.
- **Fig. 7 (page 217).** Jupiter-Uranus-Neptune-Pluto flight time against Earth launch date, March to April 1996, at
  V-infinity 14 km/s. A U-shaped curve with its minimum of about 7400 days near 22 to 28 March.
- **Fig. 8 (page 218).** Jupiter-Pluto flight times (2500 to 7500 days) against launch date, 1982 to 2042, with groups of
  short flights (2500 to 4000 days) recurring about every 13 years.
- **Fig. 9 (page 219).** Earth departure vis viva (7.5 to 12.5 km2/s2) for Mars trajectories, launch 2005 to 2015,
  marking VEGA, ballistic and Venus swingby cases.

## Citation mining

Held: Diehl, Kaplan and Penzo 1983 (Galileo tour design).
On the wanted list as of today: Diehl and Myers 1987 and Williams 1990 (row for the STOUR lineage, added with the
Petropoulos digest).
Not held and not wanted before today (one new low-priority row):
- Wallace, R. A., Lane, A. L., Roberts, P. H. and Snyder, G. C. (1981), "Missions to the far outer planets in the 1990's",
  AIAA 19th Aerospace Sciences Meeting, St. Louis. The JPL catalogue of multiple-encounter opportunities that this paper
  extends to 2050.
- Longuski, J. M. and Wolf, A. A. (1986), "The Galileo orbital tour for the 1986 launch opportunity", AIAA/AAS
  Astrodynamics Conference, Williamsburg.
Background not listed: Battin 1987 (book), Broucke 1988 (on the Petropoulos row).
