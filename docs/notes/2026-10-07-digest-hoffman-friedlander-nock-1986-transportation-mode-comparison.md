# Digest: Hoffman, Friedlander & Nock 1986, "Transportation Mode Performance Comparison for a Sustained Manned Mars Base" (AIAA 86-2016) (#960 batch 32)

S. J. Hoffman and A. L. Friedlander (Science Applications International Corporation, Schaumburg IL) and K. T. Nock
(Jet Propulsion Laboratory, Pasadena CA), "Transportation Mode Performance Comparison for a Sustained Manned Mars Base",
AIAA/AAS Astrodynamics Conference, Williamsburg VA, 18-20 August 1986, AIAA 86-2016(-CP), doi 10.2514/6.1986-2016.
Proceedings pp. 53-75 (23 pp.). U.S. Government work. Accession stamp "A86-47907" on p.1.
- **Author list confirmed: three authors.** The page-1 image shows Hoffman and Friedlander (SAIC) and Nock (JPL). The
  owner's note "Hoffman & Friedlander" omits Nock. Friedlander et al. 1986 ref. 9 (p.43 image) also gives all three.
- File given: `299af824-hoffman1986.pdf`, md5 `01bcf3789d3f63344e5f80166529c448`, 23 pages, 1,579,814 bytes.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/hoffman-friedlander-nock-1986-transportation-mode-performance-comparison-sustained-manned-mars-base-aiaa-86-2016-doi-10.2514-6.1986-2016.pdf`
- **How I read it:**
  - The scan has an old OCR text layer. The letters are space-separated ("T h i s paper"), so grep is poor, but the
    words are right. It is not image-only, so I did not OCR it. A filer who wants a greppable layer could run
    `ocrmypdf --force-ocr -l eng` on a copy (not `--redo-ocr`).
  - I read every page image at 130 dpi (pp.53-75). I read these on 300 dpi crops: the p.56 orbit sentences (VISIT
    V-infinity, aphelion range, Escalator aphelion), Fig. 3 (p.73) day labels, Table 4 totals (p.68), Table 5 rows (p.69).
  - I read Friedlander, Niehoff, Byrnes & Longuski 1986 (FNBL, held) p.35 and its Table 5 (p.43) on 130 dpi images
    for the comparison.
  - Arithmetic checks: `hoffman_checks.py`, output in `hoffman_checks.out`.
- Wanted list: **row 40**. This file fills the item "Hoffman, Friedlander & Nock (1986), AIAA 86-2016-CP".

## 0. Verdict

**This is the companion systems paper to FNBL 1986 (AIAA 86-2009), from the same conference and the same JPL
contract (957404).** FNBL designs the orbits; this paper costs them. It compares the propellant used over one 15-year
cycle by four ways of keeping 20 people on Mars and 6 at Phobos: conjunction (stop-over) transfers, VISIT orbits,
Up/Down Escalators, and the Down Escalator alone.
- **Cycler content is short and secondary.** Sec. III (p.56) describes the VISIT and Up/Down Escalator orbits in two
  paragraphs, "see Ref. 2" (FNBL) "for a more complete description". Fig. 2 and Fig. 3 are schematics. There is no
  orbit table, no a or e, no flyby altitude and no encounter date list.
- **Cycler numbers it does give** (sec. 2): VISIT period about 1.25 yr, perihelion near 1 AU, aphelion "between 1.38
  and 1.66 AU", 4:5 Earth and 3:2 Mars resonance, **V-infinity 4.2-5.2 km/s at Earth and 3.7-3.9 km/s at Mars**;
  Escalator period 2-1/7 yr, aphelion "about 2.32 AU", 51.4 deg apse precession per orbit, V-infinity 5.4-6.2 km/s at
  Earth and 6.1-11.7 km/s at Mars; Fig. 3 transits Up 147 d, Down 164 d. All check against two-body arithmetic (sec. 3).
- **Main result:** total propellant per 15 years: conjunction 30,205 t; VISIT 21,500 t; Up/Down Escalator 32,830 t;
  Down Escalator only 23,660 t (Table 9). The cyclers save less than expected, because taxis must make large hyperbolic
  rendezvous ΔVs and carry CASTLE refurbishment mass.
- **What it gives the project:** a second, independent-in-text (same team) source for the VISIT V-infinity values that
  both VISIT catalogue rows lack, an aphelion range that bears on the VISIT row's open question, and the early prior art
  for the stop-over architecture. **No new orbit.**
- **Catalogue implication (proposals only):** sec. 4. The main one: fill the null VISIT `vinf_kms` fields from FNBL
  Table 5 (the primary), citing this paper as corroboration.

## 1. Content

- **Abstract and sec. I (pp.53-54).** Steady-state Mars base (about 2035): 20 on Mars, 6 at Phobos, tours of 5-7 yr.
  LOX/LH2 only, made on Earth, the Moon, Phobos and Mars. All vehicles have closed-loop air and water. Propellant for
  Earth surface to LEO is not counted.
- **Sec. II (pp.54-55).** History (Boeing IMIS 1968, North American Rockwell, Case for Mars). CASTLEs (Cycling
  Astronautical Spaceships for Transplanetary Long-duration Excursions, after Hollister, ref. 8) ride the circulating
  orbits; Taxis make the hyperbolic rendezvous.
- **Sec. III, Trajectory options (pp.55-57).** Three options: conjunction round trip, VISIT, Up/Down Escalator.
  - Conjunction (p.56): C3 9-16 km²/s², V-infinity 2.5-4.0 km/s; launches every about 2.14 yr, 7 per 15-yr cycle;
    E-M transit 200-350 d, stopover 330-520 d, M-E 190-360 d; round trip 945-995 d (about 2.7 yr). Fig. 1 is a schematic.
  - VISIT and Escalator: sec. 2 below.
  - Table 1 (p.67): transport ΔV ranges over the 15-yr cycle, km/s (image-read):

    | mode | launch L1 to Mars | Mars orbit capture | launch Mars orbit to Earth | capture at L1 |
    |---|---|---|---|---|
    | Conjunction | 1.12-1.47 | 0.96-1.90 | 1.03-1.65 | 1.12-1.44 |
    | VISIT orbit | 1.91-2.54 | 0.81-1.00 | 2.36-2.55 | 1.91-2.54 |
    | Up/Down Escalator | 2.31-2.74 | 0.80-2.58 | 4.09-10.25 | 2.31-2.74 |

    Footnotes: conjunction capture into a 1.15 x 7.95 R_M parking orbit (0.75-1.05 km/s taxi to Phobos); VISIT and
    Escalator taxis aero-propulsive capture to Phobos; LEO-L1 tanker 4.41 km/s round trip; lunar surface-L1 tanker 5.21;
    Mars-Phobos shuttle 7.56; conjunction OTV return to L1 1.14-1.49, to Phobos 1.54-2.08; **Escalator midcourse
    adjustment 0-1.16 km/s**; midcourse navigation 0.05 km/s per leg (all modes).
  - Table 2 (p.67): ΔV by Earth-Moon staging location (m/s; trans-Mars injection, transport from LEO, transport from the
    lunar surface): LEO 4470/--/2670; GEO 3540/3820/3520; Earth-Moon L1 2050/3670/2510; **Earth-Moon cycler 1408/3058/2550**;
    lunar orbit 2230/3880/1730; lunar surface 3960/5610/--. TMI basis: C3 = 30 km²/s², "typical of up escalator", final
    burn at 6878 km perigee. L1 was chosen over the Earth-Moon cycler for station-keeping and launch-geometry reasons
    (p.58). The Earth-Moon cycler is credited to Aldrin (ref. 10, p.65).
- **Sec. IV, Infrastructure (pp.57-61).** L1 spaceport and a Phobos co-orbiting spaceport. Six vehicle types. Isp 460 s,
  mixture ratio 7 (Table 3). CASTLE dry mass 400 t, crew up to 21; conjunction CASTLE adds a 23 t propulsion system;
  cycling CASTLE adds a 60 t hangar (460 t), needs 35% of its mass delivered over 15 yr, and keeps 6 crew aboard (2
  during flybys).
- **Sec. V, Results (pp.60-64).**
  - Conjunction: 2 CASTLEs, each launched every other opportunity (4.3 yr); 7 sorties (CASTLE 1 four, CASTLE 2 three);
    crew 19 per sortie, tour 4.8 yr (3.2 on Mars). Table 4: 30,205 t propellant, 206 t consumables.
  - VISIT: **3 CASTLEs on separate VISIT orbits**, 21 crew each, 15 for Mars; 8 sorties; tours 5.7-7.9 yr, Mars time
    1.6-5.9 yr. One CASTLE is used only twice, because a flyby sequence can pass Earth twice before Mars. Taxi entry speed
    limit printed "11.1 m/sec" (p.62; surely km/s, a slip). Table 5: 21,500 t, 813 t consumables, "almost 30 percent"
    better than conjunction.
  - Up/Down Escalator: one Up and one Down CASTLE, 7 round trips each; Up carries crews out, Down brings them home.
    Fig. 7: on average 4.1 yr on the base, 0.9 yr in transit, 5.0 yr tour; 20 on Mars 90% of the time. Fig. 8: 17 aboard
    (10 Mars crew, 1 Phobos tender, 4 CASTLE riders, 2 rider/tenders). Table 6: 32,830 t, more than conjunction, because
    of twice as many launches at higher ΔV, plus the periodic orbit-rotation ΔVs.
  - Table 7: sensitivity of Down Escalator sortie 2. The taxi crew module mass dominates (5 t: -20.6%; 15 t: +19.8%).
  - Down Escalator only (Table 8): 23,660 t, 22% below conjunction, about equal to VISIT, one CASTLE.
- **Sec. VI, Summary (pp.64-66).** Tables 9-13 compare the four modes (propellant by source, production rates,
  vehicle and sortie counts, a qualitative chart, CASTLE ΔV). 2,000-3,500 t of LH2 and consumables must be lifted from
  Earth (12-18 HLLV flights of 200 t). Table 13 (CASTLE ΔV per 15 yr, m/s): conjunction 36345, VISIT 1000, Up/Down
  5159, Down 2734; major ΔVs per CASTLE 16/0/3/3 with averages 1298/0/627/627 m/s (min 1011/0/270/270, max 1646/0/1105/1105).
  Short taxi trips (under 7 days) from L1 to a circulating CASTLE "have yet to be found for arbitrary CASTLE flyby time".

## 2. The cyclers in this paper, against FNBL 1986 and the catalogue

### 2.1 VISIT (Niehoff)

Printed (p.56; 300 dpi crop): connecting paths nearly tangent to the planet orbits, "perihelion close to 1 AU and
aphelion placed somewhere between 1.38 and 1.66 AU". "The single orbit shown in Figure 2 has a period of approximately
1.25 years and encounters Mars in the vicinity of its perihelion (1.38 AU)." It revolves four times while Earth makes
five and Mars two; "4:5 resonance with Earth and 3:2 resonance with Mars"; Earth every 5 yr, Mars every 3.75 yr. Swingbys
cause a clockwise drift of several degrees per encounter, reset by "free" gravity assists. **"The relative velocity
characteristics of the VISIT orbit are 4.2-5.2 km/sec at Earth encounters and 3.7-3.9 km/sec at Mars encounters."**

| quantity | this paper (p.56) | FNBL p.35 text (VISIT-1 run, Table 1) | FNBL Table 5 (p.43), VISIT-1 / VISIT-2 | catalogue |
|---|---|---|---|---|
| period | 1.25 yr | 1.25 yr | Earth every 5.0 / 3.0 yr; Mars every 3.75 / 7.5 yr | k=7, 14.95 yr (lines 1638-1639) |
| perihelion / aphelion | near 1 AU / 1.38-1.66 AU | 0.94 / 1.39 AU | not given | VISIT-1 0.94 / 1.40 (lines 1663-1664); VISIT-2 0.95 / 1.67 (lines 1789-1790) |
| V-inf Earth | 4.2-5.2 km/s | 4.2-4.5 km/s | 4.2-4.8 / 3.7-4.0 km/s | null (line 1655; VISIT-2 line 1781) |
| V-inf Mars | 3.7-3.9 km/s | 3.7-3.9 km/s | 3.7-4.1 / 2.6-2.8 km/s | null (line 1658; VISIT-2 line 1784) |

- **The orbit described is VISIT-1.** 1.25 yr, 4:5 Earth, Mars near its perihelion. The paper never says "VISIT-1" and
  uses one VISIT type throughout.
- The 1.38-1.66 AU aphelion range spans VISIT-1 (1.38-1.39) and VISIT-2 (1.66). It brackets the catalogue's Rogers
  2012 values (1.40, 1.67). **It does not support the Wikipedia values 1.89 and 1.45 AU** in the open question of the
  VISIT-1 row note (line 1666).
- **The three V-infinity statements differ.** At Mars, this paper copies FNBL's text range (3.7-3.9). At Earth, it
  gives 4.2-5.2, wider than both FNBL ranges (4.2-4.5 text, 4.2-4.8 Table 5). I checked "5.2" on a 300 dpi crop; it is
  printed. I cannot tell from this paper whether 5.2 comes from a longer run or is a slip. FNBL Table 5 is the fuller
  primary and should be preferred.
- **Loose wording.** "revolves about the Sun four times while Earth completes five revolutions and Mars completes two"
  does not hold for Mars: in 5 yr Mars makes 2.66 revolutions. The 3:2 Mars resonance is 3 VISIT revolutions per 2
  Mars revolutions (3.75 yr), as the next sentence says. FNBL p.35 states it correctly.

### 2.2 Up/Down Escalator (Aldrin)

Printed (p.56): "first proposed by Dr. Buzz Aldrin" (ref. 10); both orbits elliptical with a 2-1/7 yr period and
aphelion beyond Mars; Up oriented for a six-month E-M transfer, Down for a six-month M-E return; on a 15/7 synodic
period orbit the encounter line rotates 1/7 circle per revolution, so a 51.4 deg apse precession per orbit must come from
propulsion or gravity assists (mainly Earth); aphelion "about 2.32 AU" (300 dpi crop); **V-infinity 5.4-6.2 km/s at
Earth and 6.1-11.7 km/s at Mars**. Fig. 3 (p.73, 300 dpi crop): Up: Earth t0, Mars t0+147 d. Down: Mars t0+65 d, Earth
t0+229 d (164 d).

- **FNBL Table 5:** Up 5.7-6.2 (Earth) / 6.1-11.7 (Mars); Down 5.4-6.0 / 6.6-11.6. **This paper's ranges are the union
  of FNBL's Up and Down columns.** FNBL Table 5 gives the E-M / M-E flight times as 0.43 / 1.71 yr (Up) and 1.71 /
  0.43 yr (Down); 0.43 yr = 157 d, between Fig. 3's 147 and 164 d.
- **Catalogue rows:** `aldrin-classic-em-k1-outbound` (line 23, `sense: outbound`, the up escalator) and
  `aldrin-classic-em-k1-inbound` (line 2687, the down escalator). V-infinity 6.5 at Earth and 9.7 at Mars (lines 54, 57;
  2717, 2720; circular-coplanar, Russell 2004). Both inside this paper's ranges.
- **Two framings of the orbit, both valid.** This paper and FNBL take the orbit period equal to the synodic period
  (15/7 yr), so a = 1.662 AU, and with Q = 2.32 AU, q = 1.00 AU and e = 0.396. That orbit returns to the same inertial
  geometry each synodic period; the 51.4 deg apse rotation must then be supplied. The catalogue uses a = 1.60, e = 0.393
  (lines 78, 88), whose orbit period is 2.02 yr (the period note at line 50 says so). Q = 2.23 AU, q = 0.971 AU. Its
  repeat is 2.135 yr through the Earth flyby. These are different idealisations of the same cycler. I do not call either
  wrong.

### 2.3 Conjunction mode (stop-over prior art)

The conjunction mode is the stop-over architecture of Penzo & Nock 2002: two reused CASTLEs, propulsive capture at both
planets, 7 opportunities per 15 years, a 330-520 d wait at Mars. See `digest-penzo-nock-2002.md` sec. 4. Not a catalogue
class (same reasoning).

## 3. Arithmetic checks (`hoffman_checks.py` -> `hoffman_checks.out`)

Constants and models are mine (two-body Sun, Earth circular at 1 AU; Mars a = 1.52368 AU, e = 0.0934 where stated).
- **Synodic period** 2.1353 yr; 7 synodic periods = 14.95 yr. Matches "2.14 years" and the "15 year cycle".
- **Escalator.** P = 15/7 yr gives a = 1.6621 AU. With Q = 2.32: q = 1.004 AU, e = 0.396. Consistent with "aphelion of
  about 2.32 AU" and a perihelion at Earth. 360/7 = 51.43 deg (printed 51.4). Agreement.
- **VISIT.** P = 1.25 yr gives a = 1.1604 AU. With Q = 1.38, q = 0.941; with Q = 1.39, q = 0.931. FNBL prints 0.94 /
  1.39 (sum 2.33 vs 2a = 2.32; within rounding). P = 1.5 yr gives a = 1.3104 AU; with Q = 1.66, q = 0.961. So 1.38-1.66 AU
  spans VISIT-1 and VISIT-2. Agreement.
- **Resonances.** 4 x 1.25 = 5 yr (Earth 4:5). 3.75 yr = 3 VISIT revolutions = 1.994 Mars revolutions (Mars 3:2). See the
  wording note in 2.1.
- **V-infinity, circular-coplanar estimates** (mine, to show the printed ranges are plausible):
  - VISIT-1 (q 0.94, Q 1.39): 4.57 km/s at Earth (printed 4.2-5.2); 3.88 km/s against Mars at its perihelion (printed
    3.7-3.9). Agreement.
  - Escalator (a 1.662, q 1.00): 5.44 km/s at Earth (printed 5.4-6.2). At Mars: 9.87 km/s against circular Mars at
    1.524 AU; 7.6-12.1 km/s against eccentric Mars at 1.524 AU (two signs of Mars radial velocity); 10.1 at Mars
    perihelion, 9.2 at Mars aphelion. The printed 6.1-11.7 is the same spread. Agreement in kind.
  - Catalogue Aldrin (a 1.60, e 0.393): 6.58 km/s at Earth and 9.74 km/s against circular Mars. These reproduce the
    catalogue's 6.5 / 9.7. Agreement.
- **Transit times.** For a = 1.662, q = 1.00: r = 1.0 to 1.524 AU takes 121 d; to 1.666 AU (Mars aphelion) 147 d. Fig. 3's
  147 d Up transit fits an encounter with Mars far out on its orbit. Real-ephemeris geometry, not checked further.
- **Side observation, outside this paper (for the filer to judge).** For the catalogue Aldrin elements (a = 1.60,
  e = 0.393), two-body time from r = 1.0 to r = 1.524 AU outbound is **103 d** (Kepler equation 102.8 d; a leapfrog
  integration in the session gave 102.7 d and an 80.5 deg transfer angle). The catalogue row gives `tof_days: 146`
  (line 74) with the same a/e on the same segment (line 78). The 146 d is the source-attested literature value, and the
  V-infinity values do match a = 1.60, e = 0.393. So the a/e and the 146 d do not belong to one two-body arc. I did not
  trace why. It is not caused by this paper.
- **Table 2 TMI.** C3 = 30 km²/s² from a 6878 km circular orbit: 4.4665 km/s. Printed 4470 m/s. Agreement.
- **Mass tables 4, 5, 6 and 8** (row: LOX + LH2 + Phobos + Mars = total; column sums):
  - **Table 4 (p.68) TOTAL "32,205" is a misprint for 30,205.** The total column sums to 30,205; the abstract, p.62 and
    Table 9 all give 30,205. Read on a 300 dpi crop. Every row sums correctly. (Rows 1C and 1D both total 4,832; both
    sum correctly.)
  - Table 5 (p.69): rows C-1 (2,903 vs printed 2,904) and A-3 (2,356 vs 2,355) are off by 1 t. Re-read on a 300 dpi crop;
    printed as I give them. Rounding. All column totals agree.
  - Tables 6 and 8: every row and column agrees.
  - LOX:LH2 = 7.0 in every row (mixture ratio 7, Table 3).
  - Average daily rates (t/day over 15 x 365.25 d) reproduce the printed values in all four tables (for example Table 5
    LOX 3.104, Table 6 4.011, Table 8 2.960; Table 4 Phobos 0.546 + 0.078 = 0.624, which is 3,419 t / 15 yr).
- **Table 9 (p.71).** Each column sums to its printed total (30,205; 21,500; 32,830; 23,660). Mars and Phobos LOX:LH2
  are 7.0. Its Mars and Phobos totals equal the corresponding Table 4/5/6/8 totals. Table 10 rates in kg/day match.
- **Savings quoted.** VISIT vs conjunction 28.8% ("almost 30 percent"). Down Escalator vs conjunction 21.7% ("22
  percent"). Agreement.
- **Table 13 (p.72).** Per-CASTLE = total / number of CASTLEs: 18172 (36345/2), 333, 2580, 2734. Agreement. Per
  sortie: VISIT 125 (1000/8), Up/Down 369 (5159/14), Down 391 (2734/7). Agreement. Conjunction 4543 = 18172/4 (CASTLE 1's
  four round trips), not 36345/7 = 5192. Footnote b's "OTV total 17592 m/s" plus 18172 is 35,764, not 36,345; I could not
  reconcile footnote b with the column. Not called a misprint.

## 4. Catalogue proposals (proposals only; catalogue is read-only)

1. **VISIT V-infinity gap.** `niehoff-visit1` (lines 1655, 1658) and `niehoff-visit2` (lines 1781, 1784) have
   `vinf_kms` null, with the note "TBD: original Niehoff 1985/1986 conference papers needed". **FNBL Table 5 (p.43,
   held) has the values:** VISIT-1 4.2-4.8 (Earth), 3.7-4.1 (Mars); VISIT-2 3.7-4.0, 2.6-2.8 km/s. Propose: fill from
   FNBL Table 5 (real-ephemeris ranges, 20-yr propagation); cite this paper as corroboration for VISIT-1 Mars (3.7-3.9),
   and record its wider Earth range 4.2-5.2 in the note. My circular-coplanar estimates (4.57 / 3.88 km/s for VISIT-1)
   fall inside the ranges. The existing FNBL digest (`docs/notes/2026-06-22-digest-friedlander-niehoff-byrnes-longuski-1986.md`)
   quotes only the p.35 text range and does not mention Table 5.
2. **VISIT aphelion open question (line 1666).** Add to the note: Hoffman 1986 p.56 gives aphelion 1.38-1.66 AU and
   FNBL p.35 gives 0.94 / 1.39 AU for VISIT-1. Both support Rogers 2012 (1.40, 1.67), not Wikipedia's 1.89 / 1.45.
3. **Source attribution.** FNBL is a `corroborating_sources` entry of the Aldrin outbound row (line 130) and both VISIT
   rows (lines 1696, 1815). **The Aldrin inbound row (line 2687; `first_published` at line 2775) has no FNBL entry.**
   Propose: add FNBL to the inbound row (FNBL Table 5 has a Down Escalator column). This paper could be added as a
   corroborating source to the two Aldrin rows (escalator V-infinity ranges, Fig. 3 transits, 2.32 AU aphelion) and to
   `niehoff-visit1` (VISIT V-infinity, aphelion range). It is secondary to FNBL; FNBL should stay first.
4. **Earth-Moon cycler.** Table 2 lists an "Earth-Moon cycler" staging option credited to Aldrin's 1985 presentation
   (ref. 10). It gives only staging ΔVs, no orbit. The catalogue's Earth-Moon Aldrin-lineage row
   (`genova-aldrin-2015-em-3petal-cycler`, line 9235, `first_published` Genova & Aldrin 2015, line 9405) could note that
   the Earth-Moon cycler idea appears in 1985-86 sources. Low value; optional.
5. **No new rows.** No orbit in this paper is new. The conjunction mode is not a catalogue class (see the Penzo digest).

## 5. Citation mining

| ref | work | status |
|---|---|---|
| [1] | SAIC, "Transportation Mode Performance Comparison ...", presentation to JPL, Contract 957404, 9 Dec 1985 | not held; not on the wanted list. FNBL ref. 8 is the same study as an SAIC Final Report, Dec 1985. **New candidate** (low priority; this paper is its published form) |
| [2] | Friedlander, Niehoff, Byrnes & Longuski 1986, AIAA 86-2009-CP | HELD (`friedlander-niehoff-byrnes-longuski-1986-...`) |
| [3] | Boeing, "Integrated Manned Interplanetary Definition Study (IMIS)", NASA CR 66558-66564, 1968 | not held (the "boeing" corpus hit is Donahue & Duggan 2022, a different work); not on the wanted list. Out of scope |
| [4] | Canetti 1968, MEM tests, NASA CR 65911-65913 | not held; out of scope |
| [5] | Boston (ed.) 1984, The Case for Mars, AAS Sci. Tech. Ser. 57 | not held; out of scope |
| [6] | McKay (ed.) 1985, The Case for Mars II, AAS Sci. Tech. Ser. 62 | not held; out of scope (it may hold early cycler papers; not checked) |
| [7] | Duke & Keaton 1986, Manned Mars Mission working group report | not held; out of scope |
| [8] | Hollister, "Castles in Space", Astronautica Acta, January 1967 | not held. Wanted **row 53**, which gives 1969, 14(2):311-316. This paper and Penzo & Nock 2002 both give 1967. The year needs checking |
| [9] | Niehoff, "Manned Mars Mission Design", Steps to Mars conference, NAS, July 1985 | not held; not on the wanted list as a separate item (wanted row 10 has Niehoff 1986 AAS 86-172 and Niehoff et al. 1991). The catalogue records it as VISIT `first_published` (not online) |
| [10] | Aldrin, "Cyclic Trajectory Concepts", SAIC presentation, JPL, 28 Oct 1985 | not held; not on the wanted list. Catalogue line 115 records it |
| [11] | Tolson et al. 1978, Viking Phobos results, Science 199:61-64 | not held; out of scope |
| [12] | Scott et al. 1985, aerobraking OTV, NASA TM 58264 | not held (the "scott" corpus hit is McAdams et al. 2011); out of scope |
| [13] | Eagle Engineering 1985, Report 85-109B, NAS9-17317 | not held; out of scope |

- **Proposal for row 40:** mark the Hoffman, Friedlander & Nock (1986) item as received.
- **Proposal for row 53:** add the 1967 date note (both papers).
- No new candidates worth fetching for the catalogue. The SAIC Dec 1985 report (ref. 1) would only add detail to the
  same tables.

## 6. Files

- `hoffman_checks.py`, `hoffman_checks.out`: the checks of sec. 3. Table values in the script are from the page images.

*Filed as `cyclers_pdf/papers/hoffman-friedlander-nock-1986-transportation-mode-performance-comparison-sustained-manned-mars-base-aiaa-86-2016-doi-10.2514-6.1986-2016.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
