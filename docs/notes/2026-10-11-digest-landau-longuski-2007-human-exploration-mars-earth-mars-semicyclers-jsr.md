# Digest: Landau and Longuski 2007, "Human Exploration of Mars via Earth-Mars Semicyclers", JSR 44(1):203-210 (#960)

Damon F. Landau and James M. Longuski, Purdue University. Journal of Spacecraft and Rockets 44(1):203-210,
January-February 2007. DOI 10.2514/1.23189. Received 19 February 2006, accepted 19 July 2006. 8 pdf pages, text layer,
no OCR. Wanted-list row "Earth-Mars semi-cycler sources" (the JSR 2007 form of AAS 05-269; the other items of that row
stay wanted).

- **Filed as** `cyclers_pdf/papers/landau-longuski-2007-human-exploration-mars-earth-mars-semicyclers-jsr-44-1-doi-10.2514-1.23189.pdf`.
  Table transcription: `data/sources/landau-longuski-2007-tables.yaml`.
  - The filed file is the received file, unchanged, md5 ab80be46286976bda97b4962f317159f (775,479 bytes). The pages are
    printed pp.203-210 (PDF page 1 = p.203, page 8 = p.210; p.210 is the last reference page).
  - The DOI, volume, issue and page range were confirmed from the first-page header, the PDF subject field and the
    Crossref lookup of 2026-10-11 in the wanted-list DOI note.
- **How it was read:** all 8 pages on a 130 dpi render. Tables 1-2 and the labels of Figs. 2-5 were read again on 300 dpi
  crops (the Fig. 5 label "10/29/2023" is overlapped by a line at 130 dpi and reads cleanly at 300 dpi).
  - Text layer, witness A (raw order) and witness B (layout mode) were parsed by script for Tables 3-11. The two parses
    are identical (38 rows of Table 3, 8 rows each of Tables 4-7, 22/20/21/19 rows of Tables 8-11).
  - Tables 1-2 were typed from the 300 dpi crops. Every date and value was found in both text extractions, except the
    cell "07 June2029", which the paper prints with no space.
  - No table digit differs between the three witnesses.
  - The text layer garbles the maths: "0.5:0.5" reads "0:5:0:5", "7/5 ~ 1.4" reads "7=5 1:4", and the infinity and delta
    signs drop out. The render was used for these.
  - The figure labels (Figs. 2-5) are drawn graphics with no text layer. They have the render as the only witness, read twice.

## 0. Verdict

**An architecture paper, not a cycler paper. The Earth-Mars semicycler is a reusable vehicle that parks at Earth and flies
by Mars twice. It is not a catalogue cycler, and it is not a ballistic periodic orbit.**
- Sequence per vehicle: Earth departure from a high Earth orbit, Mars flyby, a second Mars flyby about 5.8 years later,
  Earth arrival with propulsive capture into high Earth orbit, about 600 days of refurbishment, then the next departure.
  The crew goes out on the first Mars flyby (taxi detaches and lands) and returns on the second flyby of another vehicle
  (hyperbolic rendezvous at Mars).
- Four trajectory types are found, all Earth-Mars-Mars-Earth (E-M-M-E) with moderate delta-V. Two need four vehicles (about
  7 years, 3.3 synodic periods) and two need three (about 4.8 years, 2.2 synodic periods). The paper recommends the
  four-vehicle types and gives dated itineraries for them (Tables 1-3).
- Headline result: once the four vehicles are in place, the injected mass in low Earth orbit is 10 to 35 percent lower than
  for a semidirect mission (Mars parking orbit for the transfer vehicle). For 40 t of cargo with LH2/LOX and 180 day legs,
  seven missions need 3002 t (repeat) against 3349 t (semidirect).
- **No orbit elements are printed.** There is no semi-major axis, eccentricity, flyby altitude or flyby bend angle anywhere
  in the paper. The only trajectory numbers are dates, V_inf and DSM delta-V (Tables 1-3, Figs. 2-5).

## 1. What the paper does

- **Architecture words (pp.203-204).** The paper sorts Mars architectures by where the transfer vehicle is placed:
  semidirect (parking orbit at Mars arrival, NASA Design Reference Mission), stopover (parking orbits at both ends,
  refs 23-24), cycler (flybys at both, refs 25-29), Mars-Earth semicycler (flyby at Earth, parking orbit at Mars,
  refs 30-31), FLEM (a Mars flyby with limited stay, ref 32, Titus 1966), and the new Earth-Mars semicycler (parking orbit
  at Earth, flyby at Mars). "The inspiration for this architecture is derived from the Mars-Earth semicycler and FLEM
  concepts."
  - **Name collision.** Nock et al. 2003, which our corpus holds, uses "semi-cycler" for the Mars-parking kind (the
    Mars-Earth semicycler of refs 30-31, Bishop et al. 2000 and Aldrin et al. 2001). Landau and Longuski's own
    architecture parks at Earth and is the opposite. The word alone does not say which is meant.
- **Four trajectory types (pp.203-205, Figs. 2-5).** Classified by the ratio of Earth revolutions to spacecraft revolutions.
  All four were found "with moderate delta-V (no other trajectory types were found)".

  | Type | Resonance with Earth | Vehicles | Earth launch to Earth arrival | Body sequence |
  |---|---|---|---|---|
  | Fig. 2 | 7:5 (period about 1.4 yr, about 5 revolutions in 7 yr) | 4 | about 7 yr | E-M-M-E |
  | Fig. 3 | 3:2, 1:1, 3:2 | 4 | about 7 yr | E-M-E-E-M-E |
  | Fig. 4 | 5:4 (about 4 revolutions in 5 yr) | 3 | about 4.8 yr | E-M-M-E with two DSMs |
  | Fig. 5 | 2:1, 0.5:0.5, 2:1 (a half-year inclined Earth-Earth leg) | 3 | about 4.8 yr | E-M-E-E-M-E with two DSMs |

  - Figs. 2-5 are the circular-coplanar optima. The text says they repeat exactly each synodic period in that model.
    In the ephemeris model the 5:4 type's delta-V grows a lot when Mars is near aphelion, and the two four-vehicle types
    need much less delta-V than the three-vehicle ones.
  - Four-vehicle and three-vehicle schemes cannot be mixed: a four-vehicle launch in 2009 and a three-vehicle launch in
    2011 both arrive in 2016 and leave no return in 2018.
- **Method.** Minimise the sum of the Earth and Mars V_inf and the DSM delta-V in a circular-coplanar model, then
  re-optimise with integrated ephemerides for Earth and Mars. The paper cites its AIAA 2004-5095 paper (ref 38) for the
  delta-V data and does not describe the optimiser here.
- **Itineraries (Tables 1-3).** Seven missions, 2009 to 2022, because the trajectories approximately repeat in inertial
  space every seven synodic periods. Transfer time on the Earth-Mars and Mars-Earth legs is at most 180 days in Table 1
  and at most 240 days in Table 2. Table 3 is a timeline of the four vehicles working together.
- **Mission assumptions (24 items, pp.205-206).** Four crew. Taxi capsule 5 t. Transfer vehicle 20 t, 27 percent
  (5.33 t) renewed per mission, each vehicle fully renewed every 15 missions. Consumables 5 kg per person per day. Mars
  landing delta-V 500 m/s. Parking orbit reorientation 300 m/s at Earth and 200 m/s at Mars. Hyperbolic rendezvous
  delta-V at Mars 200 m/s. High-energy parking orbits have a 300 km periapsis and a one-day period. Isp and inert
  fractions: LH2/LOX 450 s and 10 percent, NTR 900 s and 30 percent, methane/LOX 380 s and 10 percent. Boil-off 10 percent.
  ISPP yields 18 t of methane and oxygen per tonne of hydrogen.
- **IMLEO model (Eqs. 1-5).** Launch delta-V with a 5 percent loss factor, upper-stage delta-V to the high orbit, escape
  delta-V from the high orbit, and the rocket equation for n stages. All are standard.
- **Results (Tables 4-11).** The first three semicycler missions cost more than semidirect because an extra Earth-return
  vehicle goes to Mars orbit (assumption 5). From the fourth mission on, the semicycler is cheaper. Savings per mission
  against semidirect, repeat columns: 50 t (Table 4), 48 t (Table 5), 39 t (Table 6), 45 t (Table 7). Relative savings:
  lowest 9.5 percent (Table 6, 40 t cargo), highest 35 percent (Tables 5 and 7, no cargo, ISPP). Moving from 180 to
  240 days saves 58 t (20 percent) per mission with LH2/LOX and 7.0 t (8 percent) with NTR, aerocapture and ISPP.
- **Disadvantages (p.209).** Hyperbolic rendezvous at Mars, with no abort to the surface if docking fails. Continuous
  upkeep of four vehicles. Each vehicle is empty for up to six years between Mars flybys. A variation replaces the four
  reusable vehicles with an expendable one that does not brake into Earth orbit (no number given).

## 2. Tables (all transcribed in `data/sources/landau-longuski-2007-tables.yaml`)

- **Tables 1-2 (p.205):** seven rows each, launch years 2009, 2011, 2014, 2016, 2018, 2020, 2022. Columns are Earth
  launch (V_inf), Mars arrival (V_inf), up to three Earth flyby or DSM cells (Table 2: two), Mars launch (V_inf), Earth
  arrival (V_inf). Ranges of V_inf in Table 1 (ours): Earth launch 3.05 to 4.81, Mars arrival 3.22 to 6.10, Mars
  launch 3.36 to 5.16, Earth arrival 3.07 to 6.07 km/s. The Earth-flyby V_inf values are 3.84 (twice, 2011 row) and
  10.10 and 8.39 (2022 row). The DSMs are 0.31 to 1.10 km/s in Table 1 and 0.19 and 0.23 km/s in Table 2.
- **Table 3 (p.206):** 38 events from 06 Nov. 2009 to 19 March 2025, with V_inf or DSM delta-V in four vehicle columns.
- **Tables 4-7 (p.207):** IMLEO in tonnes per launch year for initial EMSC, repeat EMSC and semidirect, with 40 t of cargo
  and with no cargo.
- **Tables 8-11 (p.208):** IMLEO breakdowns for the 2009 launch year (20 to 22 element rows each).
- **Table 12 (p.209):** a text summary of advantages and disadvantages, transcribed verbatim.

## 3. Independent checks, all ours

- **Column totals.** The seven year rows of each of Tables 4-7 differ from the printed Total row by 0 to 2 t in any
  column (rounding of the printed rows). Table 5 is one tonne below its Total in all six columns.
- **Mass breakdowns.** Sums of the element rows: Table 8, 516.1 against 516; Table 9, 575.9 against 576; Table 10, 196.7
  against 197; Table 11, 253.4 against 253.
- **Prose against tables, all agree:**
  - HEO capture 16.6 + 1.8 = 18.4 t. DSM propulsion 15.6 + 1.7 = 17.3 t. HMO capture 40.0 + 4.4 = 44.4 t. The
    "almost 9 t" eliminated at Mars arrival is 44.4 - 18.4 - 17.3 = 8.7 t.
  - Table 10 taxi feedstock plus inert 3.3 + 5.3 = 8.6 t. Table 11 taxi 1.7 + 2.7 = 4.4 t plus transfer-vehicle
    HMO-to-Earth 22.2 + 2.5 = 24.7 t gives the printed 29.1 t.
  - Savings per mission: (3349 - 3002) / 7 = 49.6 t; (1652 - 1311) / 7 = 48.7 t; (2871 - 2599) / 7 = 38.9 t;
    (1579 - 1262) / 7 = 45.3 t. Relative: 272 / 2871 = 9.5 percent; (971 - 630) / 971 = 35.1 percent; (898 - 581) / 898 = 35.3 percent.
  - 180 to 240 days, no-cargo repeat columns: (2060 - 1657) / 7 = 57.6 t, 19.6 percent of 294 t; (630 - 581) / 7 = 7.0 t, 7.8 percent of 90 t.
    So the printed 20 percent and 8 percent are the no-cargo case, not the 40 t case (13 percent and 5 percent).
  - The renewed fraction 5.33 / 20 = 26.7 percent, printed as 27 percent.
- **One loose figure.** The text says the Mars taxi mass is 65 t in Table 8. The four Mars-taxi lines sum to 60.1 t, and
  65.9 t with the crew, capsule and aerobrake line. The 65 t is the second number, truncated.
- **Transfer times.** Every Earth-Mars and Mars-Earth leg is within its limit (180 d in Table 1, 240 d in Table 2), from
  the printed dates. Shortest and longest: 155 and 180 d (Table 1), 171 and 240 d (Table 2). Earth launch to Earth
  arrival of one vehicle's itinerary is 2463 to 2532 d (6.7 to 6.9 yr) in Table 1 and 2467 to 2563 d in Table 2.
- **Table 3 against Tables 1-2.** Every vehicle entry from Earth launch 2 onward matches the Table 1 row of the same
  Earth launch, by script (TV1: 2009 and 2018 rows; TV2: 2011 and 2020; TV3: 2014 and 2022; TV4: 2016).
- **Crew stay at Mars, ours from Table 3:** 567, 557, 579, 606, 627, 577 and 562 days for missions 1 to 7 (text: 550 days).
  Mars launch to Earth arrival: 180 d except mission 4 (160 d). A vehicle waits 622 to 647 days at Earth between arrival and
  its next launch (text: about 600 days).
- **Seven-synodic-period repeat (the paper's "approximately repeat").** The Earth launches are 760 to 802 days apart
  (0.97 to 1.03 synodic periods). The TV4 events of 2010-2014 against the Table 1 row for 2022 (14 June 2010 and
  05 June 2025; 06 Sept. 2010 and 04 Sept. 2025; 25 May 2011 and 14 May 2026; 26 Dec. 2013 and 09 Dec. 2028; 24 June 2014
  and 07 June 2029) are 5470, 5477, 5468, 5462 and 5462 days apart. Seven synodic periods of 779.94 d are 5459.6 d, so the
  offsets are 10, 17, 8, 2 and 2 days. The V_inf pairs are 10.85 and 10.10, 1.16 and 1.10, 8.86 and 8.39, 4.92 and 5.16, 5.12
  and 5.39 km/s. The repeat is real but loose: the large Earth-flyby V_inf differs by 0.5 to 0.75 km/s per repeat.
- **Earth flyby bend (ours).** At 300 km altitude the largest bend at Earth is 43 degrees for V_inf = 10.10 km/s and
  55 degrees for 8.39 km/s (point-mass, mu = 398600.4 km^3/s^2). The paper prints no flyby altitude, so this tests nothing.
- **Table 3 against assumption 5 (two readings, not adjudicated).**
  - Table 3 shows a steady-state four-vehicle schedule. Its TV3 and TV4 entries before 2014 (TV3 Mars launch 1 and Earth
    arrival 1; TV4 Earth flyby 14 June 2010, DSM 06 Sept. 2010, Earth flyby 25 May 2011, Mars launch 2, Earth arrival 2) are
    tails of earlier-cycle trajectories.
  - Assumption 5 instead returns the crews of missions 1 to 3 on extra vehicles sent to Mars orbit, because the first
    vehicle does not reach Mars again until the end of the fourth mission. The "initial EMSC" columns are built on that.
- **Two rows that print the same V_inf and different dates (recorded as printed).** The 2014 rows of Tables 1 and 2 print the
  same four V_inf values (3.72, 5.81, 3.60, 3.07). The dates differ: Table 1 gives 18 Jan. 2014 and legs of 180 and
  180 d, and Table 2 gives 02 Jan. 2014 and legs of 213 and 177 d. The DSM differs (0.54 and 0.23 km/s). Read on the
  300 dpi crops, so it is not a transcription slip. The pair may reflect one optimum with a time-of-flight constraint that
  is inactive at 240 d; the paper does not say.
- **Typesetting slips.** "07 June2029" (no space) in Table 1; "13 Mar 2014" in Table 3 against "13 March 2014" in Table 1;
  "propellent" in Table 9. None changes a number.

## 4. Figures (described)

- **Fig. 1 (p.204):** a schematic of the mission: Earth parking orbit, taxi, transfer vehicle on the semicycler
  trajectory, entry capsule, hyperbolic rendezvous at Mars, cargo vehicle on a minimum-energy transfer. No numbers.
- **Figs. 2-5 (pp.204-205):** heliocentric plan views in the circular-coplanar model, with Earth and Mars orbits dotted
  and encounter points marked. The labels give dates and V_inf or DSM delta-V. They are in the YAML.
  - **Fig. 2 (7:5):** E1 17 July 2022 (3.7 km/s), M2 12 Feb. 2023 (4.3), M3 8 Jan. 2029 (4.3), E4 6 Aug. 2029 (3.7).
    The V_inf matches at both Earth ends and at both Mars flybys. The intervals are 210, 2157 and 210 days.
  - **Fig. 3 (3:2, 1:1, 3:2):** E1 9 Aug. 2022 (3.4), M2 12 Feb. 2023 (4.6), E3 27 July 2025 and E4 27 July 2026 (3.6),
    M5 8 Jan. 2029 (4.6), E6 13 July 2029 (3.4). E3 to E4 is 365 days (the 1:1 leg).
  - **Fig. 4 (5:4):** E1 1 Aug. 2022 (3.2), M2 27 Feb. 2023 (3.8), DSM 30 Nov. 2023 (0.4 km/s), DSM 3 Feb. 2026 (0.4), M3
    6 Nov. 2026 (3.8), E4 4 June 2027 (3.2). Total 4.84 yr.
  - **Fig. 5 (2:1, 0.5:0.5, 2:1):** E1 6 Aug. 2022 (3.2), M2 21 Feb. 2023 (4.1), DSM 29 Oct. 2023 (2.0 km/s), E3 2 Oct. 2024
    and E4 2 April 2025 (5.4), DSM 6 March 2026 (2.0), M5 12 Nov. 2026 (4.1), E6 30 May 2027 (3.2). Total 4.81 yr.
    The Earth-Earth leg E3 to E4 is 182 days and inclined.
  - The Earth and Mars dates of Figs. 2-5 are the same 2022-2029 window as the 2022 rows of Tables 1-2, but they are
    circular-coplanar values and differ from the ephemeris cells (for example, the 7:5 type has E1 on 17 July 2022 in
    Fig. 2 and the Table 1 2022 row has an Earth launch on 9 Sept. 2022 with the 3:2-1:1-3:2 pattern).

## 5. Relevance to this project

- **Class fit (proposal only, no catalogue edit).** The classes are `cycler`, `quasi_cycler`, `precursor_mga` and
  `mga_tour` (data/README.md lines 620-660; the schema enum has three more, which do not apply).
  - `cycler`: no. A cycler is a periodic orbit that transports between encounters. Here the vehicle stops: it brakes into
    high Earth orbit with real capture propellant (16.6 t propellant and 1.8 t inert in the Table 8 case), waits about
    600 days and is re-injected. The Earth V_inf is 3.0 to 6.1 km/s at both ends. Cycling is a logistics schedule for
    four vehicles, not a closure property of one orbit.
  - `quasi_cycler`: no. The spec asks for ballistic encounters that close up to rotation. Here the encounters at
    Earth include propulsive capture and departure, and the ephemeris itineraries carry DSMs of 0.19 to 1.10 km/s.
  - `precursor_mga`: no. Nothing inserts into a catalogue cycler.
  - `mga_tour`: the closest. A single itinerary (Earth launch, Mars, Mars, Earth arrival, optional Earth flybys and DSMs)
    is an epoch-locked flyby sequence. A Table 1 or Table 2 row could be recorded as one if a record is wanted.
  - **Consistency with precedent.** The Penzo and Nock 2002 digest (stop-over cyclers) and the Nock et al. 2003 digest both
    advise no new row and no new class. The same advice applies here. If an architecture record is wanted, a note in
    docs/notes fits better. This paper supplies the Earth-Mars-semicycler Tables 1-3 for it.
- **Is any semicycler already a catalogue row?** No. `data/catalogue.yaml` has no row whose name, notes or sources mention
  semicyclers, DOI 10.2514/1.23189 or AAS 05-269 (grep 2026-10-11).
  - **Sequence collision check (by script over `sequence_canonical`, 426 rows).** 217 rows have a cyclic sequence with
    two adjacent Mars encounters: 209 are `E-E-M-M` and 8 are `E-E-E-M-M`, all `cycler` class (S1L1, McConaghy and
    Russell-Ocampo lines). No row of another class has one, and no row has a linear `E-M-M-E` form. A semicycler E-M-M-E
    leg is topologically a window of the cyclic E-E-M-M pattern cut at an Earth encounter, but the cycler repeats with Earth
    flybys and the semicycler stops at Earth. The paper prints no semi-major axis or eccentricity, so no orbit-level
    comparison with those rows is possible. No collision is claimed.
- **Should one be proposed?** Proposal for the lead to accept or reject: no catalogue row. Reasons: no orbit elements, the
  closure is propulsive, and the itineraries are epoch-specific with DSMs. One thing the paper does give that the
  catalogue lacks as a sourced control is a published, dated set of seven ephemeris E-M-M-E itineraries with V_inf at every
  encounter (Table 1). That set is usable as a positive control for any E-M-M-E or E-M-E-E-M-E Lambert chain search over
  2009-2025, at the level of plot-free tabulated values. The paper prints no state vectors.
- **Relation to the S1L1 and Aldrin lines.** The semicycler is the Earth-parking mirror of Nock 2003's Mars-parking
  semicycler and an alternative to the Aldrin and S1L1 cyclers that the catalogue holds. Landau and Longuski 2006
  (guidance strategy, held) cites the semicycler as "[4]" in the same lineage and uses a cycler taxi. The 2006 digest
  said the AAS 05-269 source was not held; this paper is the journal form of that source.
- **Literature-gate label the content supports.** No gate anchor for semicyclers exists in literature_check.py (grep for
  "semi", "05-269", "23189", "Landau" on 2026-10-11 finds none). If one is added: topology {mga-tour} for the single
  E-M-M-E itinerary, heliocentric centre, body set {Earth, Mars}, flyby body Mars. "resonant" is only descriptive here
  (the Earth:spacecraft ratios name the types). No "repeated-moon" or "pump-tour" label applies. A candidate
  E-M-M-E sequence would be a published architecture, with new-system status only for a different epoch or leg set.
- **Numbers the project may want to check against its Lambert code** (no more than a V0-level control):
  - Table 1, 2009 row: Earth launch 06 Nov. 2009 V_inf 4.59, Mars 05 May 2010 V_inf 6.10, DSM 06 Oct. 2010 0.59 km/s, DSM
    20 June 2015 0.41 km/s, Mars 16 Feb. 2016 V_inf 3.94, Earth 14 Aug. 2016 V_inf 4.20.
  - Fig. 2: Earth V_inf 3.7 at both ends, Mars V_inf 4.3 at both flybys, legs 210, 2157, 210 days.
- **Risk of misuse.** "Semicycler" has two meanings in the corpus (see sec. 1). The paper's Earth-flyby V_inf of 10.10 and
  8.39 km/s in the 2022 row is its largest value and has no stated flyby altitude. The IMLEO numbers depend on the paper's
  24 assumptions and are not trajectory results.

## 6. Citation mining (checked against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Held:** Byrnes, Longuski and Aldrin 1993 (ref 27); Chen, McConaghy, Landau and Longuski 2005 (ref 28); Penzo and Nock
  2002, AIAA 2002-4424 (ref 24); Hoffman, Friedlander and Nock 1986, AIAA 86-2016 (ref 12); Rall and Hollister 1971,
  AIAA 71-92 (ref 26); Sohn 1964 (ref 9).
- **Already on the wanted list, in the row "Earth-Mars semi-cycler sources":** Bishop et al. 2000, AAS 00-255 (ref 30);
  Aldrin et al. 2001, AIAA 2001-4677 (ref 31); Penzo and Nock 2002, AAS 02-162 (ref 37); Landau and Longuski 2005,
  AAS 05-269 (now attribution only). Also on the list: Niehoff, Friedlander and McAdams 1991, IAF-91-438 (ref 23); Hollister
  1969, "Castles in Space" (ref 25, medium to low); McConaghy et al., AAS 03-509 (ref 29, held as the journal form).
- **Not in the index and not on the wanted list (new candidates):**
  - **Titus, R. R. (1966), "FLEM: Flyby-Landing Excursion Mode", AIAA Paper 66-36 (ref 32).** The paper names it as part of
    the architecture's inspiration. Titus 1965 (AIAA 65-515) was digested the same day. Medium priority: the FLEM is
    the earliest Mars flyby-with-landing scheme in the lineage.
  - **Landau, D. F. and Longuski, J. M. (2004), "Reassessment of Trajectory Options for Human Missions to Mars",
    AIAA 2004-5095 (ref 38), journal form to be published (cited "to be published").** It is the source of the trajectory V_inf
    and delta-V data this paper uses (assumption 24). Medium priority. Not confirmed to be the held "pt1 impulsive" file;
    the 2026-06-04 mining note cites it only as a reference. Do not assume the two are the same.
  - **Landau, D. F. and Longuski, J. M. (2003), "Comparative Assessment of Human Missions to Mars", AAS 03-513 (ref 20).**
    Low priority (an earlier IMLEO comparison; the 2009 comparative assessment is held).
  - **Dixon, F. P. (1963), "EMPIRE Dual Planet Flyby Mission" (ref 8).** Low priority; it is a 1963 Mars-Venus flyby
    study, outside the cycler lineage.
  - **Drake (ed.) 1998 Reference Mission 3.0, Hoffman and Kaplan (eds.) 1997 NASA SP-6107, and the other mission-architecture
    references (refs 13-19, 21-22, 33-36):** not cycler-relevant, not added. Cohen 1989 (the 90 Day Study) is not a row.
- **Wanted-list action:** the row "Earth-Mars semi-cycler sources" stays (Bishop et al., Aldrin et al., Penzo and Nock
  AAS 02-162 remain). The JSR 2007 item is marked received and AAS 05-269 is marked attribution only.

*Filed as `cyclers_pdf/papers/landau-longuski-2007-human-exploration-mars-earth-mars-semicyclers-jsr-44-1-doi-10.2514-1.23189.pdf`. Table transcription: `data/sources/landau-longuski-2007-tables.yaml`.*
