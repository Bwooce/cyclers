# Digest: Landau & Longuski 2006, "Guidance Strategy for Hyperbolic Rendezvous" (#960 batch 29)

D. F. Landau and J. M. Longuski (Purdue), AIAA/AAS Astrodynamics Specialist Conference, Keystone CO,
21-24 August 2006, AIAA 2006-6299, doi 10.2514/6.2006-6299.
- Filed as `cyclers_pdf/papers/landau-longuski-2006-guidance-strategy-hyperbolic-rendezvous-AIAA-2006-6299-doi-10.2514-6.2006-6299.pdf`.
  16 pp., text layer, md5 9b1dc51dc1f5fbfb9da3aa3ad66822c4. Supplied by the owner.
- Do not confuse this paper with the held `landau-longuski-2006-human-mars-trajectories-pt1-impulsive-JSR.pdf`.
  That is a different paper.
- How I read it: the full text.
  - Image-checked:
    - Tables 3-4 (p.4) at 200 dpi;
    - Table 6, Eqs. 4-7 and Table 7 (p.5) and Tables 8-9 (p.6) at 170 dpi, all legible;
    - Table 10 (p.6), Table 11 (p.7), and Tables 12-13 with the docking text (p.8) at 200 dpi.
  - Tables 1, 2, 5 and 14 and the figure captions are from the text layer only. No number in this digest
    depends on them.
- Wanted-list rank 44 ("Cycler taxi rendezvous; operations background"). Removed in batch 29. (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.)
- Index used: `docs/notes/CORPUS_INDEX.md`.

## 0. Verdict

**This is an operations paper for the S1L1 two-synodic cycler of McConaghy, Landau, Yam & Longuski 2006
(held). It gives taxi delta-V and mass, not new cycler orbits.**
- The cycler is "the trajectories presented in Ref. 3" = McConaghy et al. 2006, JSR 43(2):456-465 (the
  "Notable Two-Synodic-Period Earth-Mars Cycler", S1L1). The paper studies crew boardings over seven missions,
  2009-2022: Earth departures onto outbound vehicles 1-2, and Mars departures onto inbound vehicles 3-4. It says Earth-Mars geometry repeats every seven synodic periods.
- **Table 3 restates McConaghy 2006's real-ephemeris flybys as periapsis radius.**
  - The Earth rows come from McConaghy Tables 6/7 (outbound vehicles 1-2). The Mars rows come from
    Tables 8/9 (inbound vehicles 3-4). These are the "total delta-V minimised" itineraries, not Tables 2-5.
  - With R_E = 6400 km and R_M = 3400 km, 11 of 14 rows match McConaghy's altitudes exactly. 2018-E and
    2022-M are each 100 km low. 2011-M is the outlier below. V-infinity
    matches in 13 of 14 rows. Five dates differ by one day.
  - **One row is a real mismatch: the 2011 Mars departure.**
    - Landau: 10/29/2013, V-inf 2.93 km/s, r_p 28,600 km.
    - McConaghy Table 8 Mars-5: 10/29/2013, V-inf 3.00 km/s, altitude 17,500 km (r_p about 20,900 km).
    - I read both on the page images.
    - Landau's own impulsive values for that row (2nd stage 1.034, rendezvous 0.304 km/s) reproduce with
      his 2.93 / 28,600, not with McConaghy's values (which give 1.072 and 0.212).
    - So the paper used a different solution for that flyby. **Unresolved.** It is probably a pre-publication
      itinerary, but this is not stated.
- **Catalogue:** no row cites this paper (no hit for "6299", "Guidance Strategy" or "Hyperbolic
  Rendezvous").
  - The S1L1 rows (`mcconaghy-2006-em-k2`, line 210; `s1l1-2syn-em-cpom`, line 484) are not contradicted.
    Their V-inf 4.7 / 5.0 km/s (lines 275, 278) are the circular-coplanar values. Table 3 is the
    ephemeris itinerary.
  - **Proposal only:** cite it in `mcconaghy-2006-em-k2` as taxi and operations corroboration (crew
    boarding delta-V, Tables 7-8). It is the same research group and the same itinerary, so it is not
    independent evidence. It must not raise the validation level.
- **Use for the project:**
  - It is a sourced cost model for the "taxi" axis that the catalogue keeps out of the maintenance band
    (see the Rauwolf digest). It gives:
    - Eq. 4, a one-day phasing orbit;
    - Eq. 5, injection to the cycler V-inf;
    - Eqs. 6-7, rendezvous delta-V = delta-B / delta-T.
  - I reproduced all of these to 0.001-0.002 km/s (below). They are a cheap closed form for scoring a
    cycler's crew-access cost from its flyby V-inf and r_p.

## 1. Content (READ)

- Architectures (Table 1): the cycler (flyby at Earth and Mars), the M-E semi-cycler (parked at Mars) and
  the E-M semi-cycler (parked at Earth).
- Taxi stages: a 1st upper stage (high orbit, not escape), a 2nd upper stage (to the cycler V-inf), a
  rendezvous engine and docking jets. The engines are dual for redundancy. There is a 30-minute window to
  switch rendezvous engines. The cycler is passive.
- Docking (Fig. 1): the taxi targets 10 km at 40 m/s with the rendezvous engine. It closes along the Sun
  line. Range-rate gates: 10 m/s at 10 km, then 5 m/s, 2 m/s and 0.5 m/s, then 5 cm/s inside 150 m. The
  guidance law ignores gravity (Eqs. 1-3).
- Taxi parameters (Table 4): 1-day rendezvous; 6 t capsule; 200 km circular orbit; Mars descent 0.500 and
  ascent 3.835 km/s; Isp 450 s.
- Engine and navigation errors (Tables 5-6). Inertial 3-sigma errors are 300 m and 1 m/s. Relative errors
  are 3% of range and 10 cm/s per km of range.
- 1000 Monte Carlo runs per mission.

**Table 3, cycler flyby parameters (image-checked).**

| Mission | Earth date | V-inf (km/s) | r_p (km) | Mars date | V-inf (km/s) | r_p (km) |
|---|---|---|---|---|---|---|
| 2009 | 11/27/2009 | 6.19 | 37,600 | 09/10/2011 | 3.05 | 23,900 |
| 2011 | 12/14/2011 | 4.33 | 32,600 | 10/29/2013 | 2.93 | 28,600 |
| 2014 | 01/31/2014 | 4.83 | 40,300 | 12/15/2015 | 2.48 | 23,800 |
| 2016 | 03/29/2016 | 5.29 | 41,600 | 01/15/2018 | 2.90 | 6,900 |
| 2018 | 05/22/2018 | 4.00 | 25,600 | 05/30/2020 | 3.50 | 3,700 |
| 2020 | 07/18/2020 | 3.97 | 11,000 | 07/24/2022 | 3.94 | 3,700 |
| 2022 | 09/19/2022 | 4.65 | 19,900 | 08/25/2024 | 3.74 | 5,600 |

**Results.**
- Table 7 (Earth, image-checked).
  - Integrated, with 3-sigma margin: 1st stage 2.818-3.067 km/s; 2nd stage 1.109-2.008; rendezvous
    0.378-0.851; docking 0.094-0.097.
  - Impulsive: 1st stage 2.787; 2nd stage 1.132-2.058; rendezvous 0.077-0.449.
- Table 8 (Mars, image-checked).
  - Integrated: 1st stage 1.252-1.285; 2nd stage 0.779-1.645; rendezvous 0.312-0.643; docking
    0.093-0.096.
  - Impulsive: 1st stage 1.207; rendezvous 0.001-0.304.
- Taxi IMLEO, averages:
  - 191.6 t with safety features (Table 9), against 148.4 t impulsive;
  - 213.3 t with redundant rendezvous propellant (Table 10);
  - 108.8 t with in-situ propellant (Table 11).
- Docking accuracy (Table 13): sigma 2.9 cm; closing speed -3.9 +- 0.8 cm/s; rendezvous time 24.3 h.
- The paper's conclusion: the taxi docks within 10 cm at 7 cm/s 99% of the time, and the safety features
  add about 10% to IMLEO.

## 2. Checks (`cyclers_pdf/papers/<pdf stem>-checks.py` -> `cyclers_pdf/papers/<pdf stem>-checks.out`)

- **Table 3 against McConaghy 2006 Tables 6-9** (McConaghy values from its text layer; the Table 7 and
  Table 8 cells used were also viewed on its p.462 image):
  - Earth: r_p - 6400 equals McConaghy's altitude in 6 of 7 rows. 2018-E is 100 km low (25,600 against
    25,700).
  - Mars: r_p - 3400 equals McConaghy's altitude in 5 of 7 rows. 2022-M is 100 km low (5,600 against
    5,700). 2011-M does not match (above).
  - Dates differ by one day for 2009-E, 2022-E, 2016-M, 2018-M and 2020-M.
  - With R_E = 6378 km and R_M = 3396 km, every matching row is off by a constant +22 km / +4 km. So
    Landau rounded the radii to 6400 and 3400 km.
- **Eq. 4 (impulsive 1st stage, 200 km circular to a 1-day orbit):** I get 2.787 km/s at Earth (printed
  2.787) and 1.206 at Mars (printed 1.207).
- **Eq. 5 (impulsive 2nd stage):** all 14 rows agree with the printed values to within 0.002 km/s.
- **Eqs. 6-7 (impulsive rendezvous, delta-B / 1 day):** all 14 rows agree to within 0.001 km/s, using
  Table 3's r_p as a radius. This confirms that "Periapsis" in Table 3 is a radius, not an altitude.
- **Mass tables:** Earth taxi + Earth cargo = taxi IMLEO exactly in all rows of Tables 9 and 11.
- **Table 12:** the product of the six mass ratios is 1.604, and 1.604 × 26.7 t = 42.8 t. The text says
  42.8 t; Table 11 gives 43.0 t.

## 3. Citation mining

Held status was checked with `ls cyclers_pdf/papers | grep -i` and with `docs/notes/CORPUS_INDEX.md`.

- **Held:**
  - Byrnes, Longuski & Aldrin 1993 [1];
  - Chen, McConaghy, Landau, Longuski & Aldrin 2005 [2] (`chen-mcconaghy-landau-longuski-aldrin-2005-...`);
  - McConaghy, Landau, Yam & Longuski 2006 [3] (`mcconaghy-landau-yam-2006-...`);
  - Friedlander et al. 1986 [18].
- **Not held, new wanted candidates:**
  - Landau, D. F. & Longuski, J. M., "Mars Exploration via Earth-Mars Semi-Cyclers", AAS 05-269 (2005)
    [4], with a journal version in JSR. I did not verify its DOI. It is the E-M semi-cycler source. Medium
    priority.
  - Bishop et al. 2000, AAS 00-255 [5]. Nock et al. 2003 also cites it.
  - Aldrin, Byrnes, Jones & Davis 2001, "Evolutionary Space Transportation Plan for Mars Cycling Concepts",
    AIAA 2001-4677 [6]. This is the M-E semi-cycler.
  - Penzo & Nock 2002, AAS 02-162 [19]. Nock et al. 2003 also cites it.
  - None of these is on the wanted list (checked for "semi-cycl", "05-269", "00-255", "2001-4677",
    "02-162", Bishop and Penzo).
- **Not held, not added (rendezvous, navigation and engine background):**
  - Young & Alexander 1970 [7];
  - Pearson 1989 [8];
  - Carter 1996 [9];
  - Longuski 1981/1982 [10, 12];
  - Kangas et al. 2004 [11];
  - Lee & Hanover 2005 [13];
  - Larson & Wertz 1999 [14];
  - D'Amario 1997 [15];
  - Kachmar et al. 1992 [16];
  - Zimpfer & Spehar 1996 [17].
