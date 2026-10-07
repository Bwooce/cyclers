# Digest: Takubo, Campagnola, Pellegrini & Anderson 2026, "Preliminary Contingency Trajectory Planning for Europa Clipper's Galilean Moon Tour" (AIAA SciTech 2026-1262) (#960, #943)

Y. Takubo (Stanford), S. Campagnola, E. Pellegrini and B. D. Anderson (JPL), AIAA SciTech 2026 Forum, paper
**AIAA 2026-1262**, **doi 10.2514/6.2026-1262**. Crossref: same title and four authors, "AIAA SCITECH 2026
Forum", issued 2026-01-08, 29 references.
- Journal form, not held: "Preliminary Contingency Trajectory Design for Multi-Moon Tours: Application to Europa
  Clipper", JGCD, doi 10.2514/1.G009868. Crossref: issued 2026-08-20, pp. 1-12, 37 references. It is
  wanted-list row 12.
- Source file: upload `35a00ed4-CL25_4815.pdf`, 16 pp, md5 966251db3633bf3e6c40a09248c98beb. This is the JPL
  clearance copy, CL#25-4815. It was made with pdfTeX on 2025-11-01 and has a digital text layer. It has
  27 references, against Crossref's 29 for the published form. So it is a pre-publication copy and may
  differ a little from the AIAA PDF.
- Proposed filename:
  `cyclers_pdf/papers/takubo-campagnola-pellegrini-anderson-2026-preliminary-contingency-trajectory-planning-europa-clipper-galilean-moon-tour-aiaa-2026-1262-doi-10.2514-6.2026-1262.pdf`.
- How I read it:
  - the full text layer;
  - every page rendered at 110 dpi;
  - Tables 1-4, Figs. 2 and 5-9, and the p.4, p.7 and p.10 statements quoted below, read on the page images.
  - Table 2 and Table 4 have two witnesses each (image and text layer), and Table 4 also has the slides
    (`witness-comparison.tsv`).
  - Table 2 was checked by arithmetic (`check_table2.py`, `check_table2.out`).
- Tables transcribed: `takubo-2026-clipper-contingency-tables.yaml`, for data/sources/.
- Collision check for #943 / #1025: `collision.md`. The companion slide deck has its own digest:
  `digest-takubo-2026-slides.md`.

## 0. Verdict

**This is a contingency-escape method paper for the Clipper tour, not a cycler paper. There is no collision
with gc-1, gc-2, the G-C k = 4-6 class or ge-1..3.**
- **What it is:** a three-step design pipeline for "escape transfers":
  1. a Tisserand graph and V_inf-resonance maps in a circular coplanar model;
  2. a broad search with JPL's Star tool, full-ephemeris patched conic;
  3. COSMIC multiple shooting in a full-ephemeris N-body model.
- The escape goes from a nominal Europa resonance to a low-radiation "safe orbit". The safe orbit is a
  Callisto resonance from Ca 5:3 to Ca 3:1, with perijove near Ganymede's orbit.
- **The three solutions are all one-way.** They are EEECCC, EECCC (EEECCC in N-body) and EEGCCC, each
  170-190 d long. Nothing repeats.
- **The G-C cycler appears only as a named category.** Table 1 (p.3) lists "Ganymede–Callisto Cyclers" as
  an example of the quasi-ballistic different-moon class. The p.4 text says (page image): "subsets of
  free-return trajectories, such as the Ganymede–Callisto cycler, have been explored as candidates for
  inclusion in the nominal tour design [9, 22]. Since no generalized families of such transfers are known,
  extensive numerical searches are typically required". [9] is Campagnola et al. 2019 and [22] is Russell &
  Strange 2009; both are held.
- **The search could not reach the gc range.** Table 3 (p.11) sets the minimum V_inf for flybys 2-6 at
  Eu 3.9, Ga 6.0 and Ca 4.5 km/s. gc-1 (2.40/1.81) and gc-2 (3.62/3.04) lie below these floors.
- **What it gives the project:**
  1. A citable statement from the Clipper mission-design group: no "generalized families" of G-C
     quasi-ballistic transfers are known (p.4). It supports the gc novelty context. It is an assertion,
     not a search.
  2. Table 2: a three-row V_inf conversion, Europa n:1 orbit to Ganymede and Callisto V_inf. I reproduce
     it exactly in a circular coplanar model, so it is a small positive control for any
     Tisserand/V_inf-conversion code.
  3. One period coincidence. The terminal Ca 3:1 orbit (3 T_Ca = 50.07 d, our arithmetic) sits at the
     4 S(G,C) = 50.10 d repeat of the #971 R11 k = 4 class. The paper does not use this. It is lineage,
     not prior art (see `collision.md` item 5).
- **PROPOSALS only:**
  - No catalogue row. The trajectories are one-way mga-tour fragments without published state data.
  - Wanted-list row 12 gives Anderson as "R. L." The paper and Crossref give "Brian D. Anderson" (B. D.).
    Proposal: correct the initials. Keep row 12 open for the JGCD form, which has 37 references and is
    longer.
  - The p.4 sentence could be cited in the gc novelty write-up as context, with attribution.

## 1. Content

### 1.1 Background and nomenclature (pp.2-4)
- V_inf globe with pump angle alpha and crank angle (Fig. 1). The figure labels the crank angle gamma; the
  text uses kappa.
- Resonance rho = n:m: n moon revolutions, m spacecraft revolutions. Non-resonant transfers are short
  n:m- or long n:m+. Exterior means n > m (p.3).
- Table 1 sorts transfers by same or different moon, and by quasi-ballistic or non-ballistic:
  - (i) resonant and pi-transfers: Petal Rotation, COT;
  - (ii) VILT, used in the pump-down;
  - (iii) cyclers and special pi-transfers: Ganymede-Callisto cyclers, switch-flip;
  - (iv) escape transfers: pump-down variants.
- p.4: Europa and Callisto pi-transfers were studied early but are "not included in the final Europa
  Clipper tour design [19]".

### 1.2 Problem (pp.4-5)
- Contingency: an instrument or bus failure while propulsion still works.
- Start state: the nominal closest-approach (C/A) state of the triggering flyby. The next flyby body stays
  as in the nominal tour, and only small targeting changes are allowed. After the second C/A, any moon and
  any burns are allowed.
- Fig. 2 example: E17, then E18c, G01, C01, then the safe orbit (C02, ...).

### 1.3 Energy analysis (pp.5-9)
- Baseline orbits (p.5, from [9]): Eu 4:1 and Eu 5:1 at V_inf,Eu = 4.0 km/s, and Eu 6:1 at 4.5 km/s.
- Table 2: V_inf at Ganymede and Callisto, all three rows checked (section 2).
- Tisserand graph (Fig. 3). The admissible rectangles use minimum flyby altitudes of 25 / 50 / 50 km for
  Europa / Ganymede / Callisto. The maximum-bending equations (1a)-(1b) are on p.6.
- Radiation:
  - The nominal orbits take about 40-45 kRad per revolution.
  - Raising perijove to Ganymede's radius, given as 1.07e6 km, cuts this to below 5 kRad per orbit.
  - The mission TID budget is about 3.0 MRad (p.7, image).
- Safe-orbit set (Fig. 4, purple box): "The resonance ratio ranges from Ca 5:3 to Ca 3:1, and the perijove
  radius spans the vicinity of Ganymede's orbital radius" (p.7).
- V_inf-resonance maps: Callisto (Fig. 5, p.8; V_inf 0-6.5 km/s; Ganymede contours at 6 and 7 km/s; the
  "Ga 3:1" line) and Ganymede (Fig. 6, p.9; V_inf 2.5-8.5 km/s).
- Escape routes inferred from the maps (pp.8-9):
  - from Eu 4:1 (Ca 5.077): two flybys, Ca 6:5 to Ca 5:3;
  - from Eu 5:1 (Ca 5.825): at least three flybys, Ca 3:2, 2:1, 3:1; or via Ga 3:1+ at V_inf,Ga 7.0
    to Ca 2:1;
  - from Eu 6:1 (Ca 6.376): at least three flybys, Ca 5:3, 2:1, 3:1. The final orbit "still intersects
    Ganymede's orbital path".
  - Ganymede is "inefficient for raising the perijove radius" because the Europa orbits have high V_inf,Ga.

### 1.4 Star search (pp.9-12)
- Cases are taken from the 21F31-V5 tour (NAIF SPK named in the footnote, p.9):
  - E19: Petal Rotation, Eu 5:1, near-ecliptic;
  - E12: COT1, Eu 6:1, inclined.
- Table 3 inputs:
  - up to six bodies; B[3-4] is any of E, G, C; B[5,6] is Callisto;
  - minimum altitude 50 / 60 / 60 km;
  - r_min for leg 5 is 9.0e5 km;
  - N_rev at most 5 or 3; prograde only;
  - T[6] within 200 d of t0; TOF per leg 5-100 d, and 30-85 d for the last leg;
  - total Delta-V at most 0.100 km/s (petal) or 0.200 km/s (COT); Delta-V_enc at most 0.010 km/s;
  - V_inf ranges Eu 3.9-4.6, Ga 6.0-8.0, Ca 4.5-7.0 km/s;
  - crank angle in 18-degree steps.
- Results: 16,961 solutions for E19 and 787 for E12. Minimum Delta-V is about 10 m/s for E19 and about
  80 m/s for E12 (p.10; Fig. 7).

### 1.5 COSMIC convergence (pp.10-14)
- MONTE's COSMIC with SNOPT. Point masses: Sun, Earth, Mars, Jupiter and the four Galilean moons. No SRP
  and no J2. Control nodes are kept at least 3 d from any C/A.
- Table 4 (p.13):
  - A: Delta-V 13.89 -> 6.33 m/s, TID 299.13 -> 283.81 kRad, 169.96 -> 169.76 d, ends Ca 2:1.
  - B: 55.30 -> 112.72 m/s because of an untargeted Europa flyby; 144.84 -> 175.48 kRad; ends Ca 2:1+.
  - C: 111.80 -> 43.32 m/s, 241.37 -> 242.04 kRad, ends Ca 4:2-.
- Fig. 9: all three reach the safe-orbit box. Roughly two flybys do the energy change; the others are for
  phasing (p.14).

## 2. Checks (ours)

- **Table 2** (`check_table2.py` / `.out`):
  - Method: a circular coplanar model. Jupiter GM is 126,686,534 km^3/s^2. The moon orbit radii are
    671,100, 1,070,400 and 1,882,700 km. For each Eu n:1 period, I solve for p from the given V_inf,Eu,
    then compute V_inf at Ganymede and Callisto.
  - Result: all six printed values are reproduced to the last printed digit: 6.644/5.077, 7.025/5.825
    and 7.429/6.376.
  - The apojoves are 2.72e6, 3.25e6 and 3.77e6 km. They match the positions of the Fig. 9 diamonds.
  - With these constants the Eu 5:1 perijove is 669,773 km, just inside Europa's orbit.
- **Periods** (same script): T_Ga 7.155 d and T_Ca 16.691 d. The slides print 7.15 and 16.7. S(G,C) is
  12.524 d, 3S = 37.57 d and 4S = 50.10 d. The Ga 3:1 line on the Callisto map sits at rho = 1.286; Fig. 5
  draws it between the 5:4 and 4:3 lines, which agrees.
- **Table 4:** the image, text layer and slides 17-19 agree on all 18 cells.

## 3. Errata and inconsistencies (recorded neutrally)

1. p.10 (image): the final-leg constraint, minimum distance 9.0e5 km, is said to ensure "a perijove radius
   greater than Ganymede's orbital radius". But p.7 gives Ganymede's radius as 1.07e6 km, which is larger
   than 9.0e5 km. p.9 also admits that the Eu 6:1 terminal orbit "still intersects Ganymede's orbital path".
   The 9.0e5 km floor is therefore "near" Ganymede, not above it.
2. p.10: "since Eu 6:5 is the lowest-order resonance transition achievable without a maneuver". From the p.9
   context (Ca 6:5 to Ca 5:3 from Eu 4:1), this probably means Ca 6:5 (INFERRED).
3. p.12: Solution A's "EEECCC" sequence is described as "two Europa and two Callisto flybys". The string
   has three of each.
4. The frame is "Jovicentric CLIPJ2000" in the text (p.12) and "ECLIPJ2000" in the Fig. 8 caption.
5. Fig. 1 labels the crank angle gamma; the text uses kappa.
6. The acknowledgement says "Copyright 2023", while the paper is from 2025-26. The slide deck footer on
   slide 33 says "(392M) Final Presentation", which suggests the work began as a 2023 JPL fellowship project.

## 4. Citation mining

Held status was checked with `ls papers | grep -i` and CORPUS_INDEX.md.

| Ref | Work | Held? |
|---|---|---|
| 1 | Ross, Koon, Lo & Marsden 2003 multi-moon orbiter | held (AAS 03-143) |
| 2 | Strange, Campagnola & Russell 2009, Enceladus orbiter (low-mass moons) | not held; wanted-list row 48 |
| 3 | Campagnola, Strange & Russell 2010 CMDA 108, non-tangent VILT | held as AAS 10-164 (conference form) |
| 4 | Takubo, Landau & Anderson 2024 CMDA 136, "Automated tour design in the Saturnian system" | **not held, not listed** (tour automation; low priority) |
| 5 | Landau et al. 2025 Uranus cruise and tour | held |
| 6 | Buffington, Campagnola & Petropoulos 2012 AIAA 2012-5069 | held |
| 7 | Campagnola, Buffington & Petropoulos 2014 Acta 100 | held |
| 8 | Lam, Buffington & Campagnola 2018 AIAA 2018-0202 | held |
| 9 | Campagnola et al. 2019 JGCD 42(12) | held |
| 10 | Greco, Campagnola & Vasile 2022 JGCD 45(6), belief optimal control | not held, not listed (robust design; not cycler) |
| 11 | Arya, French, Pellegrini & Campagnola 2024 AAS, Delta-V99 optimisation | not held, not listed (not cycler) |
| 12 | Strange, Russell & Buffington 2007 V_inf globe AAS 07-277 | held |
| 13 | Sims, Longuski & Staugler 1997 JGCD 20(3) | held |
| 14 | Campagnola & Russell 2010 JGCD 33(2), Endgame part 1 | held only as AAS 09-224 (conference form) |
| 15 | Buffington, Strange & Campagnola 2012 ISSFD, "Global moon coverage via hyperbolic flybys" | **not held, not listed** (COT/petal origin; low priority) |
| 16 | Anderson, Campagnola & Buffington 2018 JGCD 41(4), petal rotation | held |
| 17 | Campagnola et al. 2023 AAS SFM, "Analysis Delta-V and Eclipse Duration for Crank-over-the-Top and Petal Rotation" | **not held, not listed** (COT/petal as CR3BP periodic orbits, cited p.4 with [18]; medium-low) |
| 18 | Campagnola & Russell 2010 JGCD 33(2), Endgame part 2 | held only as AAS 09-227 (conference form) |
| 19 | Campagnola et al. 2024 ISSFD, 21F31 reference trajectory | held |
| 20 | Lantukh & Russell 2012 AIAA 2012-4749 | held |
| 21 | Scott et al. 2025 JAS 72(6), Clipper pump-down | not held; gc prior-art note row 11 ("unlikely") |
| 22 | Russell & Strange 2009 JGCD 32(1) | held |
| 23 | Landau 2018 JGCD 41(7) | held |
| 24 | Lawden 1963, Optimal Trajectories for Space Navigation | not held (book; background) |
| 25 | Landau, Campagnola & Pellegrini 2022 JAS, "Star Searches for Patched-Conic Trajectories" | **not held, not listed** (the Star tool; method only; low-medium) |
| 26 | Evans et al. 2018 CEAS Space J. 10, MONTE | not held (software; background) |
| 27 | Gill, Murray & Saunders 2005 SIAM Rev., SNOPT | not held (background) |

None of these is likely to hold a G-C or G-E cycler. The cycler sources it cites ([9], [22]) are both held.
If the wanted list is extended, refs 25 and 17 are the only ones with any project use: the Star method, and
COT/petal as periodic orbits.

*Filed as `cyclers_pdf/papers/takubo-campagnola-pellegrini-anderson-2026-preliminary-contingency-trajectory-planning-europa-clipper-galilean-moon-tour-aiaa-2026-1262-doi-10.2514-6.2026-1262.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/takubo-2026-clipper-contingency-tables.yaml`.*
