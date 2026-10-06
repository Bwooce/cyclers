# Digest: VanderVeen 1969, "The 1978 Venus-Swingby-To-Mercury Mission", Bellcomm TM-69-1013-2 (#960 batch 22; #942 R1(b), #951 R5)

A. A. VanderVeen, "The 1978 Venus-Swingby-To-Mercury Mission", Bellcomm Technical Memorandum
TM-69-1013-2, 27 February 1969, Case 103-8. NASA-CR-103637, N69-32484, NTRS 19690023106. This is not
the same paper as VanderVeen 1969 JSR 6(4) "Triple-planet ballistic flybys of Mars and Venus" (held,
digested separately).
- Filed as `cyclers_pdf/papers/vanderveen-1969-1978-venus-swingby-to-mercury-mission-bellcomm-tm-69-1013-2-nasa-cr-103637-ntrs-19690023106.pdf`.
  - This is an `ocrmypdf --force-ocr --oversample 400` copy of the NTRS scan. The source scan is md5
    c2051aeb6788b5435994211f65c85cc0, 17 pp.: a CASI disclaimer, the cover sheet, 6 memo pages, the
    references, and 8 figure pages.
  - Companions: `.txt` holds the reconciled text, read on the page images. `-ocr-disagreement-report.txt`
    holds the arithmetic checks, a per-number witness table and all raw witness text.
- I read the memo text, Tables 1-3, the references and all figure labels on the 400 dpi page images.

## 0. Verdict

**(1) Cyclic considerations (sent to main and twobody-gen2-opus first).** The memo does NOT find the
Earth-Venus 8-yr cycle unreliable.
- It AFFIRMS the cycle for Earth-Venus: "at the end of five such synodic periods, or eight years, the
  slight variation among opportunities will be found to repeat" (memo p.1).
- The "little reliability" conclusion is about the **third planet, Mercury**. Over 8 yr, 25 Mercury
  synodic periods leave "a twenty-five-day Mercury position error". With Mercury's 88-d period,
  inclination and eccentricity, "replicas" of earlier missions cannot be found (pp.2, 5).
- **For ev-A/ev-C (#942 R1(b), 16-yr repeat): no published objection here.** The mean-motion slip,
  which is my computation and not in the memo, is quantified in sec. 2.

**(2) R5 (#951).** One concrete E-V-Me opportunity with tables (sec. 3). It is a one-way mission, not a
cycler. The memo states that any Earth-Venus-third-planet cycle must be a multiple of 8 yr, and that the
exact E-V-Me cycle is 104 yr. That is a useful framing fact for R5.

**(3) References.** Minovitch TR 32-464 is cited with its date (31 Oct 1963) but no NTRS or N-number.
The other five are not held.

## 1. Content (READ, page images)

- **Introduction:**
  - Minovitch (1) first saw that a Venus approach cuts Mercury launch energy.
  - E-V-Me opportunities recur every 1.6 yr (refs 2-6).
  - Ref 3 (Sturms 1966) found no acceptable 1977 or 1978 trajectories below its launch-energy limit.
  - The 1978 opportunity was found from "a predicted 8-year repetitive cycle".
- **Cyclic considerations (pp.1-2):**
  - E-V minimal-energy opportunities recur every 1.6 yr. After 5 synodic periods (8 yr) "the slight
    variation among opportunities will be found to repeat".
  - With a third planet, the absolute period must be a common multiple of both synodic periods and an
    integer number of years, so always a multiple of 8.
  - "In the case of Venus-Mars trajectories it has been demonstrated that the absolute repetitive period
    is 32 years" (no citation given).
  - Manning (ref. 6) found a 13-yr direct-Mercury cycle (41 x 116-d synodic periods), so the exact
    E-V-Me cycle is 104 yr, which is "useless".
  - The approximation: in one Venus synodic period, 5 Mercury synodic periods and 4 V-Me synodic
    periods fit "within about six days". Over 8 yr this gives a 25-d Mercury error.
  - Figure 1 plots longitude, latitude and radius of E, V and Me over JD 244 3000-6400, with alignments
    marked at 1.6-yr steps.
- **Results:** 8 yr "to the day" were added to the 1970 reference dates (ref. 2). That gave a Venus
  energy match. The minimal-energy region spans 50 d of launch and about 20 d of Mercury arrival.
- **Figure 2** (contours vs Earth departure 3700-3750 and Mercury arrival 3860-3880, JD - 244 0000):
  - Earth v_inf .182-.190 EMOS;
  - right ascension 234-242 deg and declination -28 to -29 deg;
  - Venus v_inf .36-.42 EMOS, and Venus passage date 3818-3821;
  - Venus periapsis speed 44,000-50,000 fps, and passage radii .8/1.0/1.2 R;
  - Mercury v_inf .35-.50 EMOS.
  - The surface is bounded below by the 1.0 R Venus impact line and above by a second-leg 180 deg ridge
    near Mercury arrival 3878.
- **Venus flyby (p.4):**
  - It mainly changes inclination: the craft departs Earth at about 2 deg inclination and nearly
    matches Mercury's 7 deg.
  - Periapsis declination is -38 deg, heading due east at periapsis, symmetric about Venus' equator.
  - Almost 4 deg of inclination change.
- **Type I / II alternation (p.4, Table 2, from Minovitch's 1965-1973 list):**
  - First-leg types alternate, which predicts type II in 1978. The 1978 nominal is type I, like 1970.
  - "Both types exist during some opportunities." A type-II family 90 d earlier has only a 10-d window
    with Venus clearance above 50 nm (Table 3).
- **Conclusions:** low launch energy over a wide window, moderate Mercury arrival speeds, short trips.
  The 8-yr E-V-Me cycle is unreliable because of Mercury. A 1977 opportunity should exist but was not
  identified.

## 2. The cycle arithmetic (my computation, mean periods; NOT in the memo)

- Venus synodic period: 583.92 d. 5 x 583.92 = 2919.6 d, against 8 yr = 2922.0 d. The mean E-V geometry
  recurs **2.4 d early**, which rotates the E-V line by **-2.4 deg per 8 yr (-4.8 deg per 16 yr)**.
- Mercury: 5 synodic periods (5 x 115.88 = 579.4 d) fall 4.5 d short of one Venus synodic period.
  4 V-Me synodic periods (4 x 144.57 = 578.3 d) fall 5.7 d short. Both match his "about six days".
  Over 8 yr: 5 x 583.92 - 25 x 115.88 = 22.7 d, against his 25 d.
- 41 x 115.88 = 4751.0 d, against 13 yr = 4748.3 d. Manning's 13-yr cycle checks.
- **For #942 R1(b):** a 16-yr E-V repeat on the real ephemeris has to absorb a mean-phase slip of about
  5 deg plus the eccentricity effects. VanderVeen gives no evidence that the E-V repeat fails. His
  evidence is that the 1978 E-V dates, found by adding 8 yr, worked.

## 3. Tables (READ on the page images; arithmetic checks in the disagreement report)

**Table 1, nominal 1978 E-V-Me mission** (JD - 244 0000; 1 EMOS = 29.8 km/s = 97702 fps)

| Event | Date | v_inf km/s (EMOS) | dV km/s (fps) | Radius km (nm) |
|---|---|---|---|---|
| Earth departure | 10 Aug 78 (3730) | 5.49 (.1843); C3 = 30.2 km^2/s^2 | 4.50 (14,760) | 315 (170) |
| Venus passage | 7 Nov 78 (3820) | 11.94 (.3918) | 14.40 (47,250) | 2950 (1590) |
| Mercury arrival | 28 Dec 78 (3870) | 12.25 (.4020) | 9.84 (32,300) | 463 (250) |

- Footnote c: injection from a 170 nm circular Earth orbit; periapsis speed at Venus; entry into a
  250 nm circular orbit at Mercury.
- The "radius" values are altitudes. Fig. 2d puts the nominal point above 1.2 R, and 6052 + 2950 km =
  1.49 R_Venus.
- **Source inconsistency:** 11.94/.3918 and 12.25/.4020 both equal 30.47 km/s per EMOS, not the 29.8 of
  the footnote. The Earth row matches 29.8. The Fig. 2 contours support the EMOS values. Treat the km/s
  values for Venus and Mercury as uncertain by about 2 %.

**Table 2, E-V-Me sequence 1965-1973** (from Minovitch TR 32-464)

| Launch | 18 Dec 65 | 19 Jun 67 | 23 Jan 69 | 18 Aug 70 | 1 Apr 72 | 4 Nov 73 |
|---|---|---|---|---|---|---|
| Leg 1 angle (deg) | 250 | 108 | 257 | 110 | 265 | 103 |
| Leg 1 days | 170 | 96 | 189 | 101 | 197 | 94 |
| Leg 2 angle (deg) | 228 | 190 | 214 | 146 | 152 | 138 |
| Leg 2 days | 106 | 72 | 108 | 59 | 85 | 58 |
| Total angle (deg) | 478 | 298 | 471 | 256 | 417 | 241 |
| Total days | 276 | 168 | 297 | 160 | 282 | 152 |

Every column sums correctly.

**Table 3, alternate (type II) 1978 family**

| | Earth departure | Venus passage | Mercury arrival |
|---|---|---|---|
| JD - 244 0000 | 3640-3650 | 3809.9-3811.4 | 3878-3882 |
| v_inf (EMOS) | .2349-.2196 | .2785-.2800 | .2916-.3093 |
| v_injection (fps) | 17,220-16,430 | | |
| v_passage (fps) | | 42,900-43,090 | |
| v_entry (fps) | | | 22,210-23,790 |
| Passage distance (nm) | | 105-78 | |

## 4. OCR record (owner's three-witness request)

- Witnesses: the NTRS layer; tesseract 400 dpi (text plus the digit-whitelist pass on table and figure
  pages); macOS Vision (accurate).
- Majority voting failed on the tables. 35 of 93 table numbers had fewer than 3 of 4 readings, and five
  were in no witness at all.
- Every number was therefore read on the image and checked by arithmetic (the Table 2 sums, unit
  conversions, JD to calendar). The details are in the disagreement report.
- Note: the NTRS layer was the best single witness for Table 1's EMOS values (.3918),
  but it missed others.

## 5. Citation mining (references 1-7)

| Ref | Item | Status |
|---|---|---|
| 1 | Minovitch, JPL TR 32-464, 31 Oct 1963 | Not held; wanted rank 1. **No NTRS or N-number in this memo.** |
| 2 | Cutting & Sturms, "Trajectory Analysis of a 1970 Mission to Mercury Via a Close Encounter with Venus", JPL TM 312-505, 7 Dec 1964 | Not held. Added to the wanted list (Tier D, with 3-6) |
| 3 | Sturms, "Earth-Venus-Mercury Mission Opportunities in the 1970's", JPL SPS 37-39 Vol. IV, pp.1-5, 30 Jun 1966 | Not held |
| 4 | Sturms, "Trajectory Analysis of an Earth-Venus-Mercury Mission in 1973", JPL TR 32-1062, 1 Jan 1967 | Not held (NTRS likely) |
| 5 | Wallace, "Trajectory Analysis of a 1975 Mission to Mercury Via an Impulsive Flyby of Venus", AAS 68-113 | Not held |
| 6 | Manning, L. A., "Minimal Energy Ballistic Trajectories for Manned and Unmanned Missions to Mercury", NASA TN D-3900, Apr. 1967 | Not held (NTRS). The 13-yr Mercury cycle source |
| 7 | VanderVeen, "Preliminary Results of an Attractive Earth-Venus-Mercury Mission in 1978", Bellcomm Memorandum for File, 9 Oct 1968 | Not held |

Catalogue implications: none. These are one-way missions, not cyclers.
