# Digest: Nock et al. 2003, "An Interplanetary Rapid Transit System Between Earth and Mars" (#960 batch 29)

K. Nock, M. Duke, R. King, M. Jacobs, L. Johnson, A. McRonald, P. Penzo, J. Rauwolf and C. Wyszkowski
(Global Aerospace Corp., Colorado School of Mines, SAIC), in M. S. El-Genk (ed.), *Space Technology and
Applications International Forum (STAIF 2003)*, AIP Conf. Proc. 654:1075-1086 (2003), doi 10.1063/1.1541404.
NIAC-funded (USRA contracts 07600-49 and 07600-59).
- Filed as `cyclers_pdf/papers/nock-et-al-2003-interplanetary-rapid-transit-system-earth-mars-aip-cp-654-1075-doi-10.1063-1.1541404.pdf`.
    - Supplied original: image-only scan, 13 pp. (AIP cover page plus pp. 1075-1086), md5 1835dfe4d057a0af3319318faa2c75f6.
  - The filed copy is the OCR copy (`ocrmypdf --redo-ocr -l eng`), md5 e1ec3bde4ab20cf785b10671c84c8fe3.
- How I read it: the whole OCR text, then the page images.
  - Image-checked items:
    - p.1078 (Aldrin and semi-cycler text), with Fig. 2a enlarged at 400 dpi;
    - p.1079 at 200 dpi (17-month stopover, stopover cyclers 4-7 months, two Astrotels, 7-day taxi limit);
    - p.1080 at 110 dpi (spaceport phasing, Phobos 7.5 h, Astrotel 70 t and 2767 kg);
    - the Fig. 5 mass table (p.1081) at 300 dpi;
    - the Fig. 6 taxi mass tables (p.1082) at 200 dpi;
    - Fig. 8 and Table 1 (p.1083) at 200 dpi.
  - Two items are from the OCR text layer only: L/D 0.63 (p.1081) and the operating costs (p.1085).
    Both are marked below.
- Wanted-list rank 45 ("Cycler architecture background"). Removed in batch 29. (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.)
- Index used: `docs/notes/CORPUS_INDEX.md`. The brief's `cyclers_pdf/CORPUS_INDEX.md` paths do not exist.

## 0. Verdict

**This is an architecture paper, not a trajectory paper. It adds no new cycler orbit and no number that
conflicts with the catalogue.**
- It is the NIAC "Astrotel" architecture. Astrotels are solar-electric cycler ships. Taxis do hyperbolic
  rendezvous. Spaceports sit at lunar orbit radius and near Phobos. It is the systems companion to the held
  Rauwolf, Friedlander & Nock 2002 (AIAA 2002-5046).
- **It chooses the Aldrin cycler (up and down escalators, two Astrotels) as the reference.** It compares
  two other options: semi-cyclers and stopover cyclers. It gives no V-infinity values, no aphelion, no a/e
  and no Aldrin flyby altitudes. The only flyby altitude is the semi-cycler's 15,000 km Earth perigee.
- **The orbit numbers it does give are consistent with the catalogue's Aldrin rows** (see section 2).
  - The mid-course correction on "3 out of 7 orbits" repeats Rauwolf 2002.
  - The 2767 kg of SEP propellant per 15 years is the same about 4% propellant fraction that Rauwolf 2002
    gives.
- **Fig. 2a is reused from Rauwolf 2002 Fig. 2** (the same t0+59d, t0+151d and t0+229d labels).
  - It gives a 151 d up-escalator E->M leg and a 170 d down-escalator M->E leg.
  - The held Rauwolf digest (`2026-06-22-digest-rauwolf-friedlander-nock-2002.md`) does not record these
    two numbers. They are Rauwolf's numbers, not independent evidence.
- **Catalogue proposals (proposals only):**
  1. Add Nock et al. 2003 as a third secondary citation of the Aldrin 1985 SAIC presentation in
     `first_published` of `aldrin-classic-em-k1-outbound` (lines 114-124) and `aldrin-classic-em-k1-inbound`
     (lines 2775-2785).
     - The paper's reference list confirms the date (October 28, 1985), the venue (Interplanetary Rapid
       Transit Study Meeting, JPL) and the title ("Cyclic trajectory concepts").
     - Nock's spelling is "Aldrin, E. E."; the catalogue uses "Aldrin, B.". This is the same person.
  2. **No `fleet_size` proposal.** `data/README.md` (line 225) defines `fleet_size` as the minimum
     number of spacecraft on the cycle needed for the stated cadence. Table 1's "2 Astrotels" counts the
     up and down escalators together, which are two catalogue rows. So it does not give a sourced
     per-row value.
  3. **Do not add semi-cycler or stopover-cycler rows from this paper.** It gives no orbit numbers for them.
     Neither is a strict cycler: the semi-cycler parks at Mars, and the stopover cycler parks at both
     planets.
- **Attribution note:** Nock credits semi-cyclers to "Byrnes, 1993 and Bishop, 2000".
  - Byrnes, Longuski & Aldrin 1993 (held) is the Aldrin-cycler paper. A text search of the held PDF finds
    no "semi" at all.
  - Rauwolf 2002 credits semi-cyclers to its ref. 7 only, which is Bishop, Byrnes, Newman, Carr & Aldrin
    (2000).
  - Landau & Longuski 2006 refs 5-6 credit Mars-Earth semi-cyclers to Bishop et al. 2000 (AAS 00-255)
    and Aldrin, Byrnes, Jones & Davis 2001 (AIAA 2001-4677).
  - So do not copy the "Byrnes 1993" semi-cycler credit.

## 1. Content: the cycler orbits it uses (all read on the page image)

**Aldrin cyclers (p.1078; the reference choice).**
- The period is "approximately equal to the Earth-Mars synodic period (26 months)". The line of apsides
  rotates by "an average of about 51.4°" each orbit.
- There are two types. The "Up Cycler" has the fast leg Earth -> Mars. The "Down Cycler" has the fast
  leg Mars -> Earth.
- With two Astrotels, the transit is "~5 months".
- The paper says "the planetary encounters occur at high relative velocities". It gives no numbers.
- It needs "a modest mid-course correction on 3 out of 7 orbits", flown by SEP.
- Fig. 2a ("Low-thrust Aldrin Cyclers", enlarged on the image):
  - Up: Earth encounter at t0, Mars encounter at t0 + 151 d.
  - Down: Mars encounter at t0 + 59 d, Earth encounter at t0 + 229 d, so the M->E leg is 170 d.
  - The orbit drawing reaches past the 2.0 AU tick. The paper gives no aphelion value, so do not read one
    off the plot.
- Astrotel mass: about 70 t; Fig. 5 total is 68,324 kg. The propellant is "2767-kg ... for all major
  corrections over 15 years" (p.1080).

**Semi-cyclers (pp.1078-1079).**
- Five Earth flybys. "Cycle duration of 78 months (3 synodic periods)".
- A 6-month interplanetary trip. A Mars stopover of about 1.5 yr; elsewhere "about 17-month duration".
- Three Astrotels.
- The Earth flybys are "each 12 months apart and typically at 15,000 km perigee altitude".
- The paper says Mars and Earth flyby velocities are "somewhat lower" than Aldrin's. It gives no numbers.

**Stopover cyclers (p.1079; credited to Penzo 2002b).**
- Direct high-thrust E-M transfers with a stop at each planet. Two Astrotels.
- Flight time 4-7 months. Mars stay about 1.5 yr.

**Table 1, Cycler Orbit Comparison (p.1083; image checked at 200 dpi).**

| Option | Crew trip (months) | Astrotels | Tour of duty (yr) | Gap between crews (months) | Average crew population |
|---|---|---|---|---|---|
| Aldrin Cyclers | 5 | 2 | 5 | 2-3 | 19.0 |
| Semi-cyclers | 6 | 3 | 4.5 | 8-9 | 16.7 |
| Stopover Cyclers | 4-7 | 2 | 4.5 | 8-9 | 16.7 |

- Fig. 8 is a schematic timing chart for 2012-2027. It shows ADC, AUC, SO1, SO2 and SC1-SC3 tracks with
  example crew numbers. It has no table values.

**Other numbers.**
- The taxi must reach the Astrotel within about 7 days. The 3-burn hyperbolic rendezvous is the
  reference (after Penzo & Nock 2002a).
- The Earth spaceport is at lunar orbit radius. Its phasing costs about 1 m/s per degree over a month.
- The Mars spaceport is near Phobos (period about 7.5 h).
- The taxi uses aerocapture. The elliptical raked cone has L/D 0.63 (text layer only).
- Taxi dry mass is 15,581 kg (Fig. 6, image). The crew module total is 4,783 kg.
- Operating cost over 15 years: $23 B plus $14 B, about $2.5 B per year (text layer only).

## 2. Checks (script `cyclers_pdf/papers/<pdf stem>-checks.py`, output reproduced here)

- **Fig. 5 mass table.**
  - Dry column 59,814, consumables 8,510, total 68,324: all sum exactly.
  - The row subtotals 6,618, 9,224 and 1,629 sum exactly.
  - The subtotal column itself sums to 60,324. This is because the Utility Module Base (5,000) and
    Permanent Cargo Bay (3,000) cells are blank in that column. 60,324 + 8,000 = 68,324.
- **Fig. 6 taxi tables.**
  - 7,207 + 1,000 + 4,407 = 12,614, and + 2,967 = 15,581, as printed.
  - The crew-cabin items sum to 4,491, but the printed subtotal is 4,490. The printed total 4,783 equals
    4,491 + 292, so the 4,490 is a one-unit slip.
  - The taxi table's "Crew Module 7,207" is not the crew-module table's 4,783. The text does not explain
    the difference.
- **Aldrin geometry.**
  - 360°/7 = 51.43°, which matches "about 51.4°".
  - 7 × 2.135 yr = 14.95 yr, which matches the paper's "15 years" and "3 out of 7 orbits".
- **Semi-cycler cycle.** 6 (E->M) + 18 (Mars stopover) + 6 (M->E) + 4 × 12 (Earth-flyby loops) =
  78 months = 3 × 26, as printed. With the "about 17-month" stopover, the sum is 77. The paper uses both
  figures.
- **Propellant.** 2767 / 68,324 = 4.05%. Rauwolf 2002 states "about 3 metric tons over a 15-year cycle"
  and a "propellant mass fraction of 4.0%" (for a 75 t Astrotel). These agree.
  - I do not derive a delta-V from 2767 kg. This paper states no Isp, and Rauwolf's mass basis is
    different.

## 3. Comparison with the catalogue (read-only; line numbers in `data/catalogue.yaml`)

- **Search hits:**
  - "Nock", "Rauwolf" and "Interplanetary Rapid Transit" appear only at lines 39 and 2703 (the `dv_band`
    comments cite Rauwolf 2002 for 1.561-1.605 km/s per 15 yr), and at lines 118 and 2779 (the Aldrin 1985
    venue). No row cites Nock 2003.
  - No row exists for semi-cyclers or stopover cyclers. The up escalator is the outbound Aldrin row
    (below).
- **`aldrin-classic-em-k1-outbound` (line 23):**
  - `sense: outbound` with the up-escalator gloss (line 46): matches Nock's "Up Cycler".
  - `period.years: 2.135` (line 50): matches "approximately ... 26 months".
  - The E->M `tof_days: 146` (line 74), with `tof_days_bounds: [161, 172]` from Rogers 2012 ephemeris
    (line 75). Nock and Rauwolf's Fig. 2a gives 151 d, and Table 1 gives "5 months".
    - 151 d lies between the circular-coplanar 146 d and Rogers's 161-172 d.
    - It is a third value, from Rauwolf 2002 Fig. 2, not a correction. Rauwolf does not say which model
      produced it.
    - It could be added as a quoted note, credited to Rauwolf 2002 Fig. 2, with Nock 2003 Fig. 2a as the
      reprint.
  - `dv_band_source: byrnes-longuski-aldrin-1993` and the Rauwolf corroboration (line 39) agree with
    Nock's "3 out of 7 orbits".
- **`aldrin-classic-em-k1-inbound` (line 2687):**
  - `sense: inbound` (line 2710) matches "Down Cycler".
  - The M->E `tof_days: 146` (line 2736). Fig. 2a gives 170 d for the down leg. This is the same kind of
    model difference as above.
- **`s1l1-2syn-em-cpom` (line 484) and the other S1L1 rows:** not used by this paper.
- **Held related files:**
  - `nock-friedlander-1987-...acta-astro...pdf` (index line 131);
  - `rauwolf-friedlander-nock-2002-mars-cycler-low-thrust-AIAA-2002-5046.pdf` (index line 112);
  - `friedlander-niehoff-byrnes-longuski-1986-...pdf` (index line 110).
- **Conclusion:** no catalogue row disagrees with this paper in its source attribution or numbers. The
  proposals are in section 0.

## 4. Citation mining

Held status was checked with `ls cyclers_pdf/papers | grep -i` and with `docs/notes/CORPUS_INDEX.md`.

- **Held:**
  - Byrnes, Longuski & Aldrin 1993 (`byrnes-longuski-aldrin-1993-...`).
  - Nock & Friedlander 1987 (`nock-friedlander-1987-...`).
- **Not held, new wanted candidates (cycler-relevant):**
  - Penzo, P. A. & Nock, K. T. (2002a), "Hyperbolic Rendezvous for Earth-Mars Cycler Missions", AAS 02-162.
    Landau & Longuski 2006 also cites it (their ref. 19). This is the source of Fig. 3's 2-, 3- and 4-burn
    options.
  - Penzo, P. A. (2002b), "Earth-Mars Transportation Using Stop-Over Cyclers", AIAA Astrodynamics
    Conference, Monterey, August 2002. It defines the stopover cycler.
  - Bishop, R. H. et al. (2000), "Earth-Mars Transportation Opportunities: Promising Options for
    Interplanetary Transportation", AAS 00-255. Rauwolf 2002 (ref. 7) and Landau & Longuski 2006
    (ref. 5) credit semi-cyclers to it.
  - Hoffman, S., Friedlander, A. & Nock, K. (1986), "Transportation Performance Comparison for a Sustained
    Manned Mars Base", AIAA 86-2016-CP. This is the same conference as the held Friedlander et al. 1986
    (AIAA 86-2009). It is an early Aldrin-cycler architecture paper. Low to medium priority.
  - None of these four is on the wanted list (checked for Penzo, Bishop, Hoffman, "00-255", "02-162",
    "stopover").
- **Not held, not cycler dynamics, not added:**
  - Aldrin 1985 (the presentation is not online; the catalogue already notes this);
  - Nock 2001 (Princeton High Frontier conference);
  - Johnson, King & Duke 2002 (Colorado School of Mines report);
  - Scott et al. 1985 (aerobrake);
  - the TRW, SCARLET, Brophy, Polk, Colozza and O'Neill power and propulsion references.
