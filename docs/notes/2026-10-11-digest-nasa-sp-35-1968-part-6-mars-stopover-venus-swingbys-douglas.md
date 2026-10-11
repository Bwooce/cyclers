# Digest: NASA SP-35, Space Flight Handbooks Vol. III (Planetary Flight Handbook), Part 6, "Mars Stopover Missions Using Venus Swingbys" (Douglas Aircraft, 1968; NTRS 19680020626) (#960, owner-supplied)

Prepared by the Missile and Space Systems Division of Douglas Aircraft Company for NASA under contract NAS 2-4175. Study manager J. L. Horsewood; B. A. Glassman,
G. P. Bonham and G. P. Uphoff credited; S. M. Norman (NASA/OART Mission Analysis Division) verified the trajectory data; foreword by Jerry Deerwester, August 1967.
Government Printing Office, 1968. "General Discussion, Contour Charts and Related Graphical Information": 16 front pages, 98 numbered pages, 108 PDF pages.

**This is Part 6. It is not the 1963 Part 1 that the wanted-list row naming "Ross et al. 1963, Space Flight Handbook Vol. 3, NASA SP-35" refers to.** This part cites
Part 1 (1963, prepared for Marshall Space Flight Center) as its Ref. 6, for the orbital elements. Part 1 is not held (sec. 6).

- **Filed as** `cyclers_pdf/papers/nasa-sp-35-1968-space-flight-handbooks-vol-3-part-6-mars-stopover-missions-venus-swingbys-douglas-ntrs-19680020626.pdf`.
  Numbers: `data/sources/nasa-sp35-part6-1968-tables.yaml`.
  - Filed copy: `ocrmypdf --redo-ocr -l eng --rotate-pages -O1` of the received scan, md5 dca8c25680e28c9baf5dd2c9fc4bea37, 5,009,169 bytes, 108 pages. The received scan (a
    read-only copy from the owner's Mac) has md5 e4f5c68d47910f815fe9380695362a17, 4,324,445 bytes, and a poor letter-spaced text layer.
  - A first OCR with `--force-ocr` made a 95 MB file (rasterised at full size) and a Ghostscript recompression ran for more than 8 minutes; both were discarded.
    The `--redo-ocr` copy keeps the original page images. The OCR text is still poor on digits (it reads 4 as h, 7 as T and drops or adds commas), so every number here was read from images.
- **How it was read:** the whole OCR text; every page of the five table pages (printed pp.3, 5, 11, 14, 15) at 130 to 250 dpi; the discussion page of every opportunity (text
  layer, spot-checked on images); the Section 3 introduction; the appendix and the reference list; Figs. 3-1A, 3-1B and 3-1C as images; the sentence that gives the repeat interval on a 200 dpi image.
  The other 45 chart pages (Figs. 3-2 to 3-16 A, B, C) were not read for values.
- **Table check:** three witnesses on Tables 1-1, 1-2, 2-1, 2-2, 2-3 and 2-4 (text layer, images by eye, fresh tesseract on 250 dpi renders of the scan), with the cells where the
  witnesses differ settled on zoomed images. The text layer was wrong in a digit on three cells (below). Table 2-5, an ephemeris of about 910 rows, is not transcribed (sec. 3).

## 0. Verdict

**A systematic set of contour charts and (separately printed) tables for round-trip Mars stopover missions with ONE Venus swingby (outbound or inbound), 16 launch
opportunities from 1975 to 1999. It states, and uses, a 2338-day (about 6.4-year) repeat of the Earth-Mars-Venus configuration. No trajectory in it is periodic or repeating:
each mission is a one-way round trip of 480 to 720 days. The numbers are in Supplements A, B and C, which are NOT in this PDF but are on NTRS.**
- The handbook is the data set behind the Sohn, Hollister, Deerwester, Gillespie-Ross and Deerwester-D'Haem papers that it credits (sec. 4). It is the bridge between the
  held Hollister lineage (R1(b)) and any real-ephemeris check of a Venus-swingby row.
- A two-working-body ev row cannot be compared with it by number from this PDF alone: only contour charts and three minimum-velocity charts per opportunity are printed. The
  tables that carry the V_inf and Venus-passage values are the three supplements.

## 1. Organisation

- **Section 1, General Discussion (printed pp.1-8):** background, handbook content, technical approach, computer programs and models, the magnetic tape.
- **Section 2, Planetary constants and related information (pp.9-24):** Table 2-1 (constants), Fig. 2-1 (velocity increments to leave or brake into a circular orbit, and the
  atmospheric entry speed, against hyperbolic excess speed, for Earth at 485 km and Mars at 500 km), Fig. 2-2 (Venus bend angle, 1 to 5 Venus radii), Tables 2-2 to 2-4 (oppositions of Mars,
  inferior conjunctions of Venus, Mars-Venus alignments), Table 2-5 (positions of Earth, Venus and Mars every 10 days).
- **Section 3, Graphical trajectory data (pp.25-90):** three charts for each of 16 opportunities, with a page of discussion before each set:
  - **A, contour chart:** constant hyperbolic excess speeds at Earth and Mars (outbound leg on the left, return leg on the right), and on the swingby side the constant Venus swingby date and the
    passage distance, on a grid of date at Earth against date at Mars. The direct leg shows the porkchop shape with a ridge separating Type 1 and Type 2.
  - **B, mission velocity:** the minimum total mission velocity against mission duration for 0, 30 and 60 day stopovers, propulsive braking (solid) and aerodynamic braking (dashed) at Mars, with the
    constraint that the Venus passage distance is at least 1.1 Venus radii.
  - **C, Earth entry speed:** for those minimum-velocity missions, mostly for the 30-day stopover only.
  - Only ONE family of solutions per opportunity is on the contour chart. The authors say the charts do not hold all solutions of interest and the tables should be scanned for others.
- **Appendix (pp.91-97):** the magnetic tape format: 53-word swingby records and 26-word direct records, 800 bpi, IBM 7094 Fortran IV, 32 blocks (sec. 3).
- **Supplements A, B, C (not in this PDF):** tabular trajectory data for 0, 30 and 60 day stopovers, 16 opportunities each (contents on printed pp.vii to xi: for example Supplement A has
  sections 4.1 to 4.16, 1975 inbound to 1999 outbound, on pages 1 to 348). "Available upon request" from the Ames Technical Information Division. NTRS carries three records called
  the first, second and third of three supplements: 20070014478, 20070014485 and 20070014479 (A, B, C in that order is my inference from "first, second, third").
- **Method (Sec. 1):** a patched-conic program with fixed, inclined, eccentric ellipses (the Table 2-1 elements), the swingby treated as a rotation of the hyperbolic excess
  velocity vector, and an iteration on the post-encounter transfer time until the departure speed matches the arrival speed to 3 m/s. The Venus swingby date is stepped 2 days
  (Category 5) or 5 days (Category 3), the launch or arrival date 10 days (Table 1-2).

## 2. The mission modes and opportunities

- **Modes charted:** outbound swingby (Earth, Venus, Mars on the way out, direct return) and inbound swingby (direct outbound, Mars, Venus, Earth on the way back). **The double swingby (Venus on both legs) is
  explicitly excluded** (Sec. 1, first paragraph). The authors say swingbys are at best of marginal use for one-way missions, so the handbook is for round trips.
- **Opportunities (Table 1-1, 16 rows):** Earth-Mars opposition year, Mars-Venus alignment, category, profile.

| Opposition year | Opposition JD | Alignment JD | Category | Profile |
|---|---|---|---|---|
| 1975 | 2442762 | 2443016 | 3 | inbound |
| 1978 | 2443530 | 2443367 | 5 | outbound |
| 1978 | 2443530 | 2443685 | 5 | inbound |
| 1980 | 2444295 | 2444034 | 3 | outbound |
| 1982 | 2445060 | 2445366 | 3 | inbound |
| 1984 | 2445832 | 2445693 | 5 | outbound |
| 1984 | 2445832 | 2446032 | 5 | inbound |
| 1986 | 2446622 | 2446362 | 3 | outbound |
| 1988 | 2447433 | 2447701 | 3 | inbound |
| 1990 | 2448223 | 2448028 | 5 | outbound |
| 1990 | 2448223 | 2448370 | 5 | inbound |
| 1993 | 2448995 | 2448694 | 3 | outbound |
| 1995 | 2449761 | 2450027 | 3 | inbound |
| 1997 | 2450525 | 2450378 | 5 | outbound |
| 1997 | 2450525 | 2450694 | 5 | inbound |
| 1999 | 2451293 | 2451047 | 3 | outbound |

- **Venus-passage geometry and categories (Sec. 1, Gillespie-Ross numbering):** the 7 Mars-Venus alignments in one syzygistic period are numbered 1 to 7 from a reference alignment; inbound
  swingbys of Categories 2, 4, 6 have no favourable Venus-to-Earth leg, Categories 1 and 7 are possible but not competitive with direct Mars-to-Earth trajectories, which leaves **Category 3 and
  Category 5**. By the near mirror image between inbound and outbound trajectories the same numbering in reverse order serves outbound. So there are four swingby opportunities per period.
  - Category 3 families: Type 2 between Earth and Venus and Type 1 between Venus and Mars (outbound); the reverse order inbound; contours regular (oval speed at Mars, nearly vertical speed at Earth,
    parabolic passage date).
  - Category 5: irregular contours, smaller departure windows, usually two or more competitive families (Type 1/1, 1/2, 2/1, 2/2 combinations listed in the discussion pages).
  - Per-opportunity family notes are in the tables file. Several charts show the "three-dimensional fold" (VanderVeen): contours of constant passage distance that cross on the plane lie on different
    sheets of a surface with the swingby date as the third axis; two cut-outs (1988 inbound, 1999 outbound) occur where the Venus arrival speed is below the least Venus departure speed.
- **Dates:** departures lie before an alignment and arrivals after (Sec. 1); mission durations in the charts run about 480 to 720 days; opportunity JD ranges are in Table 1-2.

## 3. Numbers (tables file)

- **Transcribed with three witnesses:** Table 1-1 (16 rows), Table 1-2 (16 rows of date ranges), Table 2-1 (16 rows of constants, 3 planets), Table 2-2 (12 oppositions), Table 2-3 (16 inferior
  conjunctions), Table 2-4 (27 alignments with categories). **The text layer was wrong in a digit on at least:** Mars radius (3815 for 3415), 1978 outbound swingby departure range (3240 for 3140), and the
  1992 alignment longitude (294.38 for 294.34). The other differences were glyph noise.
- **Not transcribed:** Table 2-5 (Earth, Venus and Mars heliocentric longitude, latitude, distance, Sun right ascension and declination, and communication distances, every 10 days, 1975-1999). It is a
  derived ephemeris of the Table 2-1 fixed-ellipse model and the text layer cannot be read for it. If it is needed, regenerate it from the Table 2-1 elements and spot-check against a few rows.
- **Chart-read values:** only Figs. 3-1A, 3-1B and 3-1C were looked at for values, one witness each, and at the contour-spacing level (e.g. minimum mission velocity with propulsive braking about
  11 to 13 km/s and with aerodynamic braking about 7.8 to 9.5 km/s for 1975 inbound; Earth entry about 14 km/s at 480 days falling to about 10.4 km/s). Marked in the file.
- **Our checks (all ours):**
  - Table 2-1 is self-consistent: 360 / mean motion gives the printed sidereal periods; the printed synodic periods follow from the mean motions (779.94 and 583.92 days); 1 / sqrt(a) gives the printed EMOS
    speeds. sqrt(GM_Sun / 1 AU) is 29.785 km/s against the printed 1 EMOS = 29.766 km/s (0.06 percent), and the handbook does not define EMOS.
  - Tables 1-1, 2-2 and 2-4 agree: every opposition and every alignment of Table 1-1 appears in Tables 2-2 and 2-4 with the same category in the right column.
  - Repeat interval: 7 x 333.9 (Venus-Mars synodic) = 2337.4 days, 3 x 779.94 = 2339.8, 4 x 583.92 = 2335.7; the printed 2338 days is consistent. The same-category alignment spacing in Table 2-4 is
    2326 to 2350 days, so "almost exactly" is within 1 percent.

## 4. Relevance to this project

- **Which swingby stopover modes it charts:** outbound swingby and inbound swingby, Categories 3 and 5, 16 opportunities (sec. 2). Double swingby is excluded. Stopover times 0, 30, 60 days.
- **Launch opportunities and Venus-passage geometry:** sec. 2 table. The Venus passage constraint is at least 1.1 radii in the charts and at least 0.95 radii on the tape.
- **Does any trajectory repeat or get called periodic?** No trajectory repeats. The text says the Venus-swingby missions "occur periodically and may be predicted far into the future", and gives the reason: the
  Earth-Mars-Venus configuration repeats almost exactly every 2338 days (about 6.4 years), so "similar swingby missions would also repeat every 6.4 years, at least in a qualitative sense".
  - That is a repeating opportunity, not a repeating trajectory: each of the 16 is a one-shot 480 to 720 day round trip. The point matches the Gillespie-Ross digest ("a repeating opportunity, not a repeating trajectory") and
    the Young 1966 digest (6.4 yr near-repeat); the 2338-day figure is the same.
  - **Not a collision with any cycler row.** A repeating E-V-M or Venus-Mars path would sit on the 2338-day (7 Venus-Mars synodic periods) commensurability that this handbook states.
- **How the data relate to the held Hollister lineage (R1(b)):**
  - The Foreword and Sec. 1 credit "the independent efforts of Sohn and Hollister" as the first indications; the references are Sohn 1964 (held), Hollister 1963 Sc.D. thesis (held), Deerwester AIAA 65-89 (wanted),
    Gillespie and Ross 1967 (held), Deerwester and D'Haem 1967 (neither), Ross 1963 Part 1 (neither, sec. 6). Hollister and Prussing 1965 (held) is not cited.
  - **A positive-control link to the held Gillespie-Ross 1967 Table 1.** Four of its entries are the same missions as handbook opportunities, by alignment date: the 31 Jan 1983 Category 3 homebound mission is
    the 1982 inbound opportunity (alignment JD 2445366); 24 Dec 1983 Category 5 outbound is the 1984 outbound (2445693); 27 Nov 1984 Category 5 homebound is the 1984 inbound (2446032); 23 Oct 1985
    Category 3 outbound is the 1986 outbound (2446362). Gillespie-Ross print trip time, stopover and four V_inf values to 0.01 EMOS for those four. The handbook's contour charts and Supplement tables for
    the same four opportunities are a second source for the same missions, so a two-source check on a held digest's numbers is possible once the supplements are fetched.
  - The Category 3 and Category 5 split, the "mirror image" between inbound and outbound, and the E-V-M commensurability are the handbook's restatement of Gillespie and Ross, which Gillespie and Ross
    take from Ross 1963 (Part 1).
  - **Mirror link:** the handbook uses a "near mirror image relationship between inbound and outbound trajectories" to reverse-number the alignments. This is the same time-reversal symmetry as in Miele and Wang 1999
    (digest of 2026-10-11) in a different setting; here it is only used to transfer Category 3 and 5 conclusions from inbound to outbound, not proved.
- **Can a two-working-body ev row be compared with it?**
  - **V_inf values:** the contour charts give speeds at Earth and at Mars in EMOS x 10 (labels about 1.0 to 6.0 in the 1975 inbound chart, i.e. about 3 to 18 km/s with EMOS = 29.766 km/s) on the date grid. A row's
    Earth and Mars V_inf and dates can be placed on the chart. Resolution is the contour spacing, so a row can be shown to lie on the chart's contour field, not matched to a digit.
  - **Venus passage distance:** the charts carry passage-distance contours (1.0 to 6.0 Venus radii) and Venus swingby date contours (for example 2443020 to 2443050 for 1975 inbound). The handbook's passage
    distance is a patched-conic value on the fixed-ellipse model and not a real-ephemeris value.
  - **The model differs from ours:** patched conics, fixed inclined eccentric ellipses (0.1 to 0.2 deg error against time-varying elements in 1999), one swingby per mission, stopover at Mars. An ev row built as a
    cycler (no stopover, repeating) is a different object; the comparison is of leg V_inf at matching dates, not of whole missions.
  - **Numerics are in the Supplements:** NTRS 20070014478, 20070014485 and 20070014479 hold the tabulated V_inf, passage distance, dates and times. Fetch these before any digit-level comparison.
- **Catalogue implication:** none. No row is a cycler.
- **Gate note (project rule):** a comparison against chart-read values is a plausibility check on topology and magnitude, never V1 evidence.

## 5. Other things the handbook says

- Entry speeds at Earth for the minimum-velocity missions vary with stopover time by less than 0.5 km/s in most opportunities, by 1 to 2 km/s in 1975 inbound, 1982 inbound at short durations, 1984 outbound at long durations
  and 1990 inbound (propulsive, 60-day), and the 1997 outbound set has all three stopover times drawn. Cusps in the mission-velocity curves occur at the Type 1 to Type 2 change of the direct leg.
- Aerodynamic braking at Mars: velocity increments for braking into a 500 km Mars orbit are excluded, which is why the dashed curves are about 3 km/s lower on the 1975 inbound chart (plot reading).
- Several minimum-velocity missions use a direct leg at the edge of the date range searched (for example 680-day durations in 1984 outbound and 1990 outbound), so the charts are bounded by the search ranges of Table 1-2.

## 6. The other parts of the handbook

This volume names one other part, Part 1 (its Ref. 6). A NASA Technical Reports Server search on 2026-10-11 (record titles only, nothing downloaded) shows the parts below; none is held and none is
in the wanted list as a row yet (corpus-gobetz-sonnet is adding a row naming Ross et al. 1963, which is Part 1; checked at edit time).

| Part | NTRS record | Title (as NTRS) | Year | Held / wanted / neither |
|---|---|---|---|---|
| 1 | 19640004262 | Speed Contours and Auxiliary Graphs for Manned Missions to Mars and Venus (author Stanley Ross; round-trip missions 1965-1999) | 1963 | neither; a row is being added by corpus-gobetz-sonnet |
| 2 | 19640004263 | Supplementary Trajectory Data: Venus to Earth and Mars to Earth, Part 2 (outbound flights to Venus and Mars, 1960-2000) | 1963 | neither |
| 3 | 20090021916 | Supplementary Trajectory Data, Part 3 (returning flights to Earth, 1960-2000) | 1963 | neither |
| 4 | none found | NTRS did not return a Part 4 for my queries | | neither (not searched under another title) |
| 5 | 19660011756 | Trajectories to Jupiter, Ceres and Vesta (Ross) | 1966 | neither |
| 6 | 19680020626 | Mars Stopover Missions Using Venus Swingbys (this document) | 1968 | held (this digest) |
| 6, Suppl. A, B, C | 20070014478, 20070014485, 20070014479 | the three supplements of tables | 1968 | neither; wanted if the digit-level comparison is to be done |
| 7 | 19690024112 | Direct Trajectories to Jupiter, Saturn, Uranus and Neptune | 1969 | neither |
| 8 | 19690025472 (supplement 19710068714) | Jupiter Swingby Missions to Saturn, Uranus, Neptune and Pluto | 1969 | neither (Jovian swingby lane) |
| 9 | 19700019224 | Direct and Venus Swingby Trajectories to Mercury | 1970 | neither (Mercury swingby lane; compare the held VanderVeen 1969 Bellcomm TM) |

NTRS lists these records as having downloads (flag in the record; not downloaded). The years for Parts 7 to 9 are inferred from the record numbers, not read from the records. The owner's expectation of "chapters of things we want" fits Parts 1 to 3 (Mars and Venus round trips), the Part 6 supplements,
Part 9 (Venus swingby to Mercury) and Part 8 (Jupiter swingby tours). I did not read any of them.

## 7. Citation mining (against CORPUS_INDEX.md and the wanted list, 2026-10-11; proposals only, no row edited)

| Ref. | Item | Status |
|---|---|---|
| 1 | Sohn 1964, J. Spacecraft Rockets 1(5) | held (digest 2026-10-05) |
| 2 | Hollister 1963, MIT Sc.D. thesis | held (digest 2026-10-07) |
| 3 | Deerwester, AIAA 65-89 (Jan 1965) | wanted (row for AIAA 65-89) |
| 4 | Gillespie and Ross 1967, J. Spacecraft Rockets 4(2) | held (digest 2026-10-05) |
| 5 | Deerwester and D'Haem, "Systematic comparison of Venus swingby mode with standard mode of Mars round trips", J. Spacecraft Rockets 4(7):904-911 (July 1967) | neither; **proposed row**: DOI 10.2514/3.28986 (Crossref 2026-10-11), conference form AIAA 67-27 (10.2514/6.1967-27), NTRS 19670039603 and 19670053330; medium, the standard-mode comparison for the R1(b) lineage |
| 6 | Space Flight Handbooks Vol. 3, NASA SP-35 Part 1 (1963) | neither; row being added by corpus-gobetz-sonnet; NTRS 19640004262 |
| 7 | Explanatory Supplement to the Astronomical Ephemeris and the American Ephemeris and Nautical Almanac (1961) | neither; book, low |
| 8 | Clarke, V. C., "Constants and Related Data for Use in Trajectory Calculations", JPL technical report, March 1964 (printed as "32-004" in the OCR; the number is probably 32-604) | neither; low; the source of the planetary elements, useful only to settle the 0.06 percent EMOS difference. NTRS search by title found nothing |
| 9 | Makemson, Baker and Westrom, "Analysis and Standardization of Astrodynamic Constants", UC Technical Report No. 12 (Feb 1961) | neither; low |
| 10 | VanderVeen, A. A., "Venus Swingbys for Manned Mars Missions During the 1978-1986 Time Period", AAS Space Flight Mechanics Specialist Conference, Paper 66-125 (July 1966) | neither; **proposed row** (medium): the three-dimensional fold of the contours and the source for 1978-1986 opportunities; the held VanderVeen 1969 JSR paper is a different, later paper |

- **Proposed wanted rows (digest only):** (a) Part 6 Supplements A, B, C (NTRS 20070014478, 20070014485, 20070014479): medium-high, they hold the numbers; (b) Deerwester and D'Haem 1967;
  (c) VanderVeen AAS 66-125; (d) Parts 2 and 3 (NTRS 19640004263, 20090021916) of the Handbook as the 1960-2000 one-way tables, low to medium.

*Filed as `cyclers_pdf/papers/nasa-sp-35-1968-space-flight-handbooks-vol-3-part-6-mars-stopover-missions-venus-swingbys-douglas-ntrs-19680020626.pdf`. Numbers: `data/sources/nasa-sp35-part6-1968-tables.yaml`.*
