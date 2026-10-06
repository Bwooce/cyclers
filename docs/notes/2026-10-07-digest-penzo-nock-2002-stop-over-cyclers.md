# Digest: Penzo & Nock 2002, "Earth-Mars Transportation Using Stop-Over Cyclers" (AIAA 2002-4424) (#960 batch 32)

P. A. Penzo and K. T. Nock (Global Aerospace Corporation, Altadena CA), "Earth-Mars Transportation Using Stop-Over
Cyclers", AIAA/AAS Astrodynamics Specialist Conference and Exhibit, Monterey CA, 5-8 August 2002, AIAA 2002-4424,
doi 10.2514/6.2002-4424. 7 pp. NIAC (USRA) grants 07600-25 and 07600-58.
- File given: `7e00eada-penzo2002.pdf`, md5 `30baf993d208170a49b06c2601b9b3d0`, 7 pages, 110,944 bytes.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/penzo-nock-2002-earth-mars-transportation-stop-over-cyclers-aiaa-2002-4424-doi-10.2514-6.2002-4424.pdf`
- **How I read it:**
  - The PDF is born-digital (embedded Type 1C fonts). Its text layer is clean and complete. No OCR was needed.
  - I read the whole text layer. I read every page image at 150 dpi (Figs. 1-6 are not in the text layer).
  - I read Table 1 (p.5) on a 300 dpi crop. All 14 rows x 9 cells agree with the text layer.
  - Arithmetic checks: `penzo_checks.py`, output in `penzo_checks.out`.
- Wanted list: **row 40** (Earth-Mars semi-cycler and stopover-cycler sources). This file fills the item
  "Penzo, P. A. (2002), 'Earth-Mars Transportation Using Stop-Over Cyclers', AIAA/AAS Astrodynamics Conf., Monterey".
  - Row 40 gives that item as Penzo alone. The page image has two authors: Penzo and Nock. The Nock et al. 2003 digest
    (ref. "Penzo 2002b") also cites it as Penzo alone.
  - Row 40 lists "Penzo & Nock (2002), AAS 02-162" as a **separate** item. The Nock et al. 2003 digest gives AAS 02-162 the
    title "Hyperbolic Rendezvous for Earth-Mars Cycler Missions". It is a different paper. **It stays open.**
  - The batch also fills the row-40 item Hoffman, Friedlander & Nock 1986 (see `digest-hoffman-1986.md`).

## 0. Verdict

**This paper defines the "stop-over cycler", and it is not a cycler orbit.** A stop-over cycler is a reused vehicle
(the station) that makes direct propulsive transfers. It departs a 4-day high ellipse at one planet, captures into a
4-day high ellipse at the other, waits there for the next opportunity, and goes back. There are no gravity assists, no
flybys and no repeating heliocentric orbit. Each leg is an ordinary minimum-ΔV Lambert transfer at a fixed flight time.
- **The only numbers are Table 1 (p.5):** 14 legs (7 Earth-to-Mars, 7 Mars-to-Earth), 2011-2025, all at 210-day flight
  time, for two vehicles V1 and V2. Each row has the launch and arrival dates, the launch C3, the launch ΔV, the arrival
  V-infinity, the arrival ΔV, the total ΔV and the propellant mass fraction.
- **The paper gives no period, aphelion, perihelion, flyby altitude or orbit elements for any cycler.** The only cycler
  numbers are in one sentence on p.3 about the Aldrin cycler (sec. 2). I did not fill this gap from other sources.
- **Table 1 is internally sound.** I reproduce all 28 burn ΔVs from C3 and V-infinity to 1 m/s or better (sec. 3). One
  misprint: the V1-E 10/28/24 arrival date (sec. 3.1).
- **What it gives the project:** the defining source for the stop-over architecture, a fully checkable table of direct
  E-M and M-E transfers (useful as a conjunction-class Lambert positive control for 2011-2025), and a three-way
  comparison point (cycler, semi-cycler, stop-over) for architecture notes.
- **Catalogue implication (proposal only): no new row and no new class.** See sec. 4.

## 1. Content

- **Abstract and introduction (p.1-2).** Cyclers and semi-cyclers give short (about 6 months or less) flights. They have
  high approach and flyby speeds, long waits between flights, short launch periods for hyperbolic rendezvous, and need
  robotic rendezvous for cargo. The stop-over option claims: lower departure and arrival speeds for a given flight time;
  flexible dates; no hyperbolic rendezvous; the station stays near the planet for refuelling and refurbishment; other
  uses for the station while it waits.
- **Elements (p.2):** Astrotel, Spaceports at each planet, Taxis, Shuttles (surface to LEO or Phobos orbit), tankers,
  and solar-electric cargo vehicles. From the GAC cycler studies (ref. 6).
- **Heliocentric optimisation (p.2).** Point-conic model. Flight time is held fixed and the launch date is varied to
  minimise the sum of the departure and arrival ΔVs. The burns are at 200 km altitude from and into 4-day ellipses:
  "apoapsis of about 200,000 km for Earth, and 95000 km for Mars". Flight times 150-250 days, every opportunity of one
  15-year cycle, 2011-2025.
- **Figs. 1-2 (pp.3-4, read on the image).** Minimum transfer ΔV against flight time, 150-250 d, one curve per launch
  year (2011, 2014, 2016, 2018, 2020, 2022, 2024). Earth-to-Mars ΔV runs from about 1.45 km/s (2018, near 200-210 d) to
  about 5.75 km/s (2011, 150 d). Mars-to-Earth runs from about 1.5 to 6.5 km/s on the same axes. The plot legend repeats
  the orbit sizes: "Apoapsis: 200000km(E) 95000km(M)". I give only these coarse read-offs; I did not digitise the curves.
- **Strategy (p.3).** With a 3.0 km/s vehicle limit (LOX/LH, 50% propellant), flight times are about 170 d (2016, 2022),
  about 200 d (2011, 2014, 2024), and 150 d (2018, 2020, under 2.5 km/s). With 210 d for all legs, the average is
  "2.2 km/s" per trip and propellant is "25% less" than for 3.0 km/s.
- **Table 1 (p.5)** is the 210-day strategy. Two vehicles are needed, because the Mars launch is 1-2 months before the
  Earth launch, so the paths cross. Each vehicle waits "almost 18 months" at a planet.
- **Launch period (p.4).** A 30-day launch period costs 85 to 240 m/s above the Table 1 nominals.
- **Phasing (p.6).** The arrival ellipse will not have the right orientation 18 months later. Options: multi-burn
  correction, a solar-electric stage over months (ref. 9), lunar gravity assist or solar perturbation (ref. 10, the
  Ariane GTO precedent). Not analysed.

## 2. Cycler content (all of it)

The paper names cyclers only to contrast them with the stop-over option.
- p.1-2: "Earth-Mars cyclers and semi-cyclers ... provide fast flight times, of 6 months or less" (refs. 1-5).
- p.3 (image checked): "the Aldrin cyclers, which have flight times near 150 days. There, the flyby velocities are about
  6 km/s for Earth and range from 6 to 12 km/s for Mars".

Comparison with the catalogue (all read-only, `data/catalogue.yaml`):
- `aldrin-classic-em-k1-outbound` (line 23): `tof_days: 146` (line 74), `vinf_kms` 6.5 at Earth (line 54) and 9.7 at Mars
  (line 57). Penzo's "near 150 days", "about 6 km/s" and "6 to 12 km/s" agree in kind. Penzo's 6-12 km/s Mars range
  matches the 6.1-11.7 km/s real-ephemeris range of Friedlander et al. 1986 Table 5 (see the Hoffman digest).
- `aldrin-classic-em-k1-inbound` (line 2687): same values (lines 2717, 2720, 2736).
- `s1l1-2syn-em-cpom` (line 484): not mentioned by Penzo.
- `niehoff-visit1` / `niehoff-visit2` (lines 1616, 1748): not mentioned by Penzo.
- Grep results: "penzo" appears only in the two Voyager rows (Kohlhase-Penzo 1977, lines 50173-50347). "stop-over",
  "stopover", "semi-cycler" and "upescalator": no hits. `2syn`: the S1L1 and Russell rows. `aldrin`: the rows above plus the
  Rogers 2015 establishment rows (lines 2799, 2946) and the Genova-Aldrin Earth-Moon row (line 9235).
- `quasi_cycler`: 20 rows, lines 52926-54967. **All 20 are CR3BP KAM-corridor rows** (Braik-Ross C21/C32 corridors and
  Lyapunov L1 corridors). None is Earth-Mars.

## 3. Arithmetic checks (`penzo_checks.py` -> `penzo_checks.out`)

Constants are standard values chosen by me (μE = 398600.4418, RE = 6378.137 km; μM = 42828.37, RM = 3396.19 km).
None is printed in the paper.
- **4-day orbit size.** A 4-day ellipse with periapsis at 200 km has apoapsis altitude 199,926 km at Earth and 94,213 km
  at Mars (radii 206,304 km and 97,609 km). The paper's "about 200,000 km" and "95000 km" are the apoapsis
  **altitudes**. Agreement.
- **Burn ΔVs.** Each burn is a periapsis burn from or into the 4-day ellipse: sqrt(C3 + 2μ/rp) minus the ellipse
  periapsis speed. Earth μ for Earth launches and Earth arrivals, Mars μ for Mars. **All 28 printed burns are
  reproduced to within 1 m/s.** The largest difference is 1 m/s (V1-E 2011 launch 1.006 vs 1.005; V2-E 2022 arrival 1.008
  vs 1.009; V2-M 2016 launch 0.859 vs 0.860). So the arrival V-infinity column is the hyperbolic excess at the arrival
  planet, and the C3 column is V-infinity squared at departure.
- **Totals.** Total ΔV = launch + arrival in all 14 rows, to the last digit.
- **Propellant fraction.** The column fits 1 - exp(-ΔV/(g0 Isp)) with Isp = 449.6-450.1 s in every row. **The paper
  does not state an Isp. 450 s is my inference.** The text's "3.0 km/s ... implies a propellant loading of 50%" needs
  Isp = 441 s. So the 50% figure is rounded (at 450 s, 3.0 km/s gives 49.3%).
- **Flight time.** 13 of 14 rows are exactly 210 d from launch to arrival. See 3.1.
- **Text against table:**
  - "average velocity per trip will be 2.2 km/s": the 14-row mean of the total ΔV is 2.134 km/s. Rounding up; minor.
  - "propellant usage will be 25% less than for the 3.0 km/s": the mean fraction is 38.0%, which is 24.0% less than 50%.
    Agreement.
  - "it ranges from 30% to47%" (p.4): **the table minimum is 25.51%** (V1-M 03/20/18) and the next is 27.94% (V2-E
    05/10/18). Three rows are below 30%. The text is loose; the table is right (it matches the rocket equation).
  - "Mars launch occurring 1 to 2 months earlier than for Earth launch": 49-65 d in the seven pairs. Agreement.
  - "waiting at alternate planets for almost 18 months": the waits at Mars (arrive, then next Mars launch of the same
    vehicle) are 508-532 d = 16.7-17.5 months. **The waits at Earth are 610-648 d = 20.0-21.3 months.** The "almost 18
    months" holds at Mars only.
  - Earth-launch spacing 767-807 d (25.2-26.5 months): matches "every 26 months, on average".

### 3.1 Misprint: V1-E 10/28/24 arrival date

The cell prints `05/26/26` (300 dpi crop; the text layer agrees). That is 575 d after launch. Launch + 210 d is
**05/26/25**. Every other row is exactly 210 d, the table title says 210 days, and the burn ΔVs of this row reproduce.
So the year digit is a misprint for 25. Keep the printed value in any transcription and add a note.

## 4. Is the stop-over cycler a catalogue class gap? (proposal only)

**My answer: no. It is outside the scope by design, not a gap.**
- The schema enum (`data/catalogue.schema.json` line 480) is: `cycler`, `quasi_cycler`, `precursor_mga`, `mga_tour`,
  `resonant_po`, `torus_homoclinic`, `quasi_periodic_torus`. (README lines 620-660 list only the first four. The schema
  is the real enum.)
- None fits:
  - `cycler`: needs a strictly periodic, non-epoch-locked orbit that transports between encounters. A stop-over leg
    ends in capture.
  - `quasi_cycler`: closes up to rotation over 10-15 system periods with 3-15 returns, by ballistic encounters.
  - `precursor_mga` / `mga_tour`: both are gravity-assist sequences. A stop-over leg has no flyby.
- `docs/notes/2026-06-16-catalogue-scope-taxonomy.md` admits classes that are "mission-actionable, structurally
  searchable (closure equations, not arbitrary flyby chains), and have strong existence priors". A stop-over leg is a
  single Lambert arc. It has no closure condition. The "cycling" is a logistics schedule for one vehicle, not a property
  of the orbit.
- The Nock et al. 2003 digest (`docs/notes/2026-10-06-digest-nock-et-al-2003-interplanetary-rapid-transit-system.md`,
  item 3) reached the same view: do not add semi-cycler or stop-over rows.
- If the project ever wants an architecture record, a note in `docs/notes` (cycler vs semi-cycler vs stop-over, with
  Table 1 of this paper, Hoffman 1986 Tables 4-9, and Nock 2003 Table 1) fits better than a catalogue row.
- **Prior art the paper does not cite:** Hoffman, Friedlander & Nock 1986 (AIAA 86-2016) studied the same architecture
  as their "conjunction" mode: two reused CASTLE vehicles, propulsive capture at both planets, 7 Earth departures per
  15 years, Mars stay 330-520 d. Nock is a co-author of both papers.

## 5. Citation mining

| ref | work | status |
|---|---|---|
| [1] | Hollister, "Castles in Space", Astronautica Acta, 1967 (printed "Holster") | not held. On wanted **row 53**, but row 53 gives "(1969), 14(2):311-316". **Both this paper and Hoffman 1986 give 1967** (Hoffman: "January 1967"). The year needs checking before a fetch |
| [2] | Aldrin, "Cyclic Trajectory Concepts", SAIC presentation, JPL, 1985 | not held; not on the wanted list. Unpublished slides; the catalogue already records it (line 115, date 1985-10-28) |
| [3] | Friedlander, Niehoff, Byrnes & Longuski 1986, AIAA 86-2009 | HELD (`friedlander-niehoff-byrnes-longuski-1986-...`) |
| [4] | Nock & Friedlander, "Elements of Mars Transportation System", printed "Astronautica Acta, Vol. 15, No. 6/7, pp. 505-522, 1987" | HELD (`nock-friedlander-1987-elements-mars-transportation-system-acta-astro-...`). The journal is Acta Astronautica |
| [5] | Bishop, Byrnes, Newman, Carr & Aldrin 2000, AAS 00-255 | not held; on wanted **row 40** |
| [6] | Nock et al., "Cyclical Visits to Mars via Astronaut Hotels", GAC Report 510-04911-007, 2000 (NIAC Phase I) | not held; not on the wanted list. **New candidate** (low priority; the cycler content is likely in Nock 2003, which is held) |
| [7] | Bate, Mueller & White 1971, Fundamentals of Astrodynamics | not held; textbook; not needed |
| [8] | Condon & Pearson 2001, Gateway libration-point missions | not held; not on the wanted list; out of scope |
| [9] | Sweetser, Cheng, Penzo & Finlayson 2001, SEP capture and escape spirals | not held (the one CORPUS_INDEX "sweetser" hit is in a different entry, line 594); out of scope |
| [10] | Penzo 1999, AAS 99-106, Ariane ASAP Mars missions | not held; out of scope |

- **Proposal for row 40:** mark the Penzo (2002) Stop-Over Cyclers item as received. Correct the author to Penzo & Nock
  and add AIAA 2002-4424 / doi 10.2514/6.2002-4424. Keep AAS 02-162 (a different paper) and the other row-40 items open.
- **Proposal for row 53:** add a note that two independent citations give 1967 for "Castles in Space".
- New candidate: ref. [6] (GAC NIAC Phase I report), low priority.

## 6. Files

- `penzo_checks.py`, `penzo_checks.out`: the checks of sec. 3. Table 1 values in the script are from the 300 dpi image.

*Filed as `cyclers_pdf/papers/penzo-nock-2002-earth-mars-transportation-stop-over-cyclers-aiaa-2002-4424-doi-10.2514-6.2002-4424.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
