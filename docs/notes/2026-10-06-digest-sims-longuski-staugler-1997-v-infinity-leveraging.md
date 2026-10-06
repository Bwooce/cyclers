# Digest: Sims, Longuski & Staugler 1997, "V-infinity Leveraging for Interplanetary Missions: Multiple-Revolution Orbit Techniques" (#960 batch 29)

J. A. Sims, J. M. Longuski and A. J. Staugler (Purdue), J. Guid. Control Dyn. 20(3):409-415 (May-June 1997),
doi 10.2514/2.4064. Conference version AAS 95-306 (Halifax, August 1995). Based on Sims's 1996 Purdue PhD thesis.
- Filed as `cyclers_pdf/papers/sims-longuski-staugler-1997-v-infinity-leveraging-interplanetary-multiple-revolution-orbit-techniques-jgcd-20-3-409-doi-10.2514-2.4064.pdf`.
  7 pp., text layer, md5 4f9e0e8fe677d6f202abca9e3049a12a. Supplied by the owner.
- How I read it: the full text.
  - Image-checked at 200 dpi:
    - the K:L(M)± definition and conversion rule (p.409);
    - p.411: the Saturn example, the 3:2/4:3 ranges, the interior-type ± definition, "0.35 year", the
      1:1 cases and the Fig. 4 example;
    - p.413: Tables 1-2, the 6.36 / 3.33 km/s text, "4.92 to 4.23", and the 1.6 AU Mars remark;
    - p.414: Tables 3-4, the Fig. 10 dates and 2.95 km/s, and the Mercury 6.39 km/s;
    - p.415: Table 5.
  - From the text layer only: the algorithm details on p.410 (185 km parking orbit, 200 km minimum flyby,
    1e-5 tolerance). They are marked below.
- Wanted-list rank 57 ("VILT foundations (cited by Campagnola 2010, Lantukh 2015); X1 tooling"). Removed in batch 29. (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.)
- Index used: `docs/notes/CORPUS_INDEX.md`.

## 0. Verdict

**This is the source of the K:L(M)± notation that the catalogue's Rogers-based establishment rows use.
It is also the canonical generalised delta-V-EGA (V-infinity leveraging) paper. It has no cycler content
of its own.**
- **Notation (p.409, image):**
  - K = number of Earth orbit revolutions;
  - L = number of spacecraft orbit revolutions;
  - M = the spacecraft revolution on which the maneuver is made;
  - "+, − = Earth encounter just after (before) the spacecraft passes the line of apsides".
  - For exterior types this means after (before) perihelion. For interior types it means "before or after
    aphelion" (p.411).
  - The nominal orbit period is K/L years. K is also the approximate time in years from launch to the
    Earth flyby.
  - M may be left out when there is one spacecraft revolution. The old "3+ delta-V-EGA" is the new
    3:1(1)+, or 3:1+.
- **Catalogue attribution, PROPOSAL.** In row `mcconaghy-2005-em-case1` (notes, lines 11761-11767) the
  catalogue says the 4:3(2)- prefix "refers
  to the V_inf-leveraging K:L(M)+/- notation introduced in their [Rogers 2012] Nomenclature table".
  - The notation was introduced by Sims, Longuski & Staugler (AAS 95-306, 1995; JGCD 1997). Rogers 2012
    reuses it.
  - Proposed wording: "...notation of Sims, Longuski & Staugler 1997 (JGCD 20(3):409, doi 10.2514/2.4064),
    restated in Rogers 2012's Nomenclature table".
  - The gloss itself is right for exterior types. The gloss in `u0l1-4-3-3-establishment-cc`
    (line 4466, "before line of apsides") matches Sims's wording.
  - The other K:L mentions (lines 75, 80, 3677, 3746 in `visit-1-5-4-3-establishment-cc`, and 4467) make
    no attribution claim and need no change.
  - `grep -n "K:L(M)\|introduced in"` finds no other "introduced in Rogers" claim. `grep -n -i
    nomenclature` hits only SnLm and Russell-Ocampo nomenclature fields.
- **Juno cross-reference.** Under Sims's own conversion rule, Juno's "2+ delta-V-EGA" (Lam 2008 digest,
  `2026-06-17-digest-lam-2008-juno.md` lines 126 and 213) is a **2:1(1)+** exterior delta-V-EGA. That
  digest's action item 6 asked for this paper; it is now supplied.
- **No efficiency table exists.** Leveraging efficiency is shown only as plots (Fig. 4, exterior; Figs. 5-6,
  launch V-inf and deep-space delta-V). The one quoted number is from the text (p.411, image): for a 3:1
  delta-V-EGA, a 1 km/s aphelion delta-V raises the Earth V-inf by about 8 km/s, against 1.7 km/s if the
  same 1 km/s were added at launch. That is a leverage ratio of about 4.7.

## 1. Content (READ)

- **Exterior delta-V-EGA (after Hollenbeck 1975).**
  - Launch tangentially into an orbit with a period slightly longer than the nominal K/L years.
  - Apply a retrograde tangential delta-V at aphelion.
  - Re-encounter Earth non-tangentially with a higher V-inf.
- **Interior delta-V-EGA (new here).**
  - Launch into an orbit with aphelion at 1 AU.
  - Apply a delta-V at perihelion to raise the aphelion.
  - Examples: 1:2, 3:4, 2:3.
  - "Trajectories with periods less than 0.35 year cannot be achieved", so 1:3 and 1:4 are not possible.
  - 1:1− is not possible (interior or exterior). 1:1+ exterior has no advantage over a direct launch.
- **Algorithm (Eqs. 1-17).**
  - Circular Earth orbit at 1 AU; two-body.
  - Iterate the apse delta-V until the true-anomaly mismatch is zero (tolerance 1e-5). Do this separately
    for the + and − branches.
  - Parking orbit 185 km. Minimum flyby altitude 200 km. (Text layer only; my launch-delta-V check below
    reproduces the printed values with 185 km.)
  - Sweetser's Jacobi-integral estimate is used as a cross-check. Aerogravity assist is noted.
- **Performance claims (Fig. 2).**
  - 3:2 types beat single-revolution types for final aphelion of about 1.6-4.3 AU.
  - 3:2(2) is slightly better than 3:2(1) below about 3 AU.
  - 4:3(1) is better still below about 2.4 AU.
  - − types need slightly less total delta-V at the low end.
- **Saturn example (p.411).**
  - Direct launch: V-inf 10.3 km/s, launch delta-V 7.28 km/s.
  - 3:1− : launch V-inf 6.95 km/s (5.23 km/s) plus aphelion delta-V 0.39 km/s = 5.62 km/s total.
  - This saves 1.66 km/s and adds 2.89 yr.

**Table 1, Hygiea, 3:2 delta-V-EGA (analytic / MIDAS; image-checked).**

| | 3:2(1)+ | 3:2(1)− | 3:2(2)+ | 3:2(2)− |
|---|---|---|---|---|
| Launch V-inf (km/s) | 3.60 / 3.56 | 3.51 / 3.46 | 3.47 / 3.43 | 3.36 / 3.33 |
| Aphelion delta-V (km/s) | 0.471 / 0.458 | 0.497 / 0.483 | 0.508 / 0.492 | 0.540 / 0.523 |
| Days from aphelion | 0 / 4.7 | 0 / 5.2 | 0 / 5.1 | 0 / 4.9 |
| E-E time of flight (yr) | 3.13 / 3.13 | 2.86 / 2.87 | 3.14 / 3.13 | 2.86 / 2.86 |

**Table 4, Mercury, 3:4 delta-V-EGA (analytic / MIDAS; image-checked).**

| | 3:4(1)+ | 3:4(1)− | 3:4(4)+ | 3:4(4)− |
|---|---|---|---|---|
| Launch V-inf (km/s) | 4.09 / 4.12 | 3.94 / 3.90 | 3.58 / 3.60 | 3.35 / 3.29 |
| Perihelion delta-V (km/s) | 0.448 / 0.464 | 0.487 / 0.513 | 0.591 / 0.610 | 0.665 / 0.701 |
| Days from perihelion | 0 / 3.8 | 0 / 1.2 | 0 / 2.0 | 0 / 0.7 |
| E-E time of flight (yr) | 3.09 / 3.09 | 2.91 / 2.91 | 3.11 / 3.11 | 2.88 / 2.88 |

- **Table 5** (3:4(1,4)− to Mercury, split maneuver): launch V-inf 3.49 km/s. Perihelion delta-V 0.184 km/s
  on revolution 1 (2.0 d from perihelion) and 0.452 km/s on revolution 4 (1.9 d). E-E time 2.89 yr.
- **Table 2** (Hygiea injected mass):

  | Launch vehicle | Direct (kg) | 3:2(2)− (kg) | Cost |
  |---|---|---|---|
  | Pegasus XL/Star 27 | 40 | 80 | $20 M |
  | Delta II 7925 | 550 | 1000 | $54 M |
  | Atlas IIAS | 1000 | 2200 | $105-145 M (RY 1999) |

- **Mars replaces the deep-space delta-V.** The 3:2 nominal aphelion is 1.6 AU, close to Mars. Fig. 9 is a
  STOUR E-M-E-Hygiea search.
  - Fig. 10 is a 4:3(2)− analogue: launch 3 March 2007, Mars 3 July 2009, Earth 2 January 2011, Hygiea
    2 January 2012. Launch V-inf 2.95 km/s, with no deterministic deep-space delta-V.
  - This Mars-for-delta-V substitution is the mechanism that Rogers 2012/2015 use for cycler
    establishment.

## 2. Checks (`cyclers_pdf/papers/<pdf stem>-checks.py` -> `cyclers_pdf/papers/<pdf stem>-checks.out`)

- **Launch delta-V from a 185 km circular orbit:**
  - V-inf 10.3 gives 7.292 km/s (printed 7.28).
  - V-inf 6.95 gives 5.236 km/s (printed 5.23).
  - 5.23 + 0.39 = 5.62, and 7.28 - 5.62 = 1.66, as printed.
- **Hygiea totals:**
  - Direct (V-inf 6.36): 4.931 km/s (printed 4.92).
  - 3:2(2)− (3.33 + 0.523): 4.243 km/s (printed 4.23).
  - The 0.01 km/s offsets come from my mu and radius choices.
- **Nominal orbits (period K/L yr, circular Earth):**
  - 3:2 aphelion 1.621 AU (paper: "1.6 AU");
  - 3:1 aphelion 3.160 AU;
  - 4:3 aphelion 1.423 AU;
  - 3:4 perihelion 0.651 AU;
  - 1:2 perihelion 0.260 AU.
- **Table 1 / Table 4:** MIDAS (realistic Earth orbit) differs from the analytic launch V-inf by
  -0.06 to +0.03 km/s, and from the analytic apse delta-V by -0.017 to +0.036 km/s. The total delta-V is
  within 0.03 km/s for every type. This supports "agree quite well".
- **Table 5:** 0.184 + 0.452 = 0.636 km/s. This lies between 3:4(1)− (0.513) and 3:4(4)− (0.701), as the
  text says.

## 3. What it gives the project's code

- **The project already has VILT code** (`grep -ril "leverag\|vilt" src`):
  - `src/cyclerfinder/search/vilm.py` (phase-free floor);
  - `leveraging_leg.py` (phase-full single leg);
  - `leveraging_chain.py` (multi-hop endgame);
  - also `endgame_graph.py`, `releg_*`, `tisserand_mga_window.py`, `moon_prune.py`, `asteroid_leveraging.py`.
- **I found no heliocentric K:L(M)± delta-V-EGA solver.**
  - `vilm.py`, `leveraging_leg.py` and `leveraging_chain.py` cite Campagnola & Russell, "The Endgame
    Problem". `leveraging_leg.v_m_kms` reads `core.satellites.SATELLITES`, which holds planetary
    satellites only. So the Sun-Earth cases in Tables 1 and 4 cannot be run through `leveraging_leg` as
    it is now.
  - A grep of the VILT-hit files for `helio`, `"Sun"` and `MU_SUN` finds two uses:
    - `releg_solver.py` (a DSM / low-thrust V-inf retarget, `mu = MU_SUN_KM3_S2` at line 672);
    - `tisserand_mga_window.py` (heliocentric Tisserand windows).
    - I did not assess whether either could solve a K:L(M)± leg.
  - Tests: `tests/data/test_aldrin_establishment.py` pins the Rogers 2012 Table 4 values that the
    catalogue carries. It is a drift guard, not a reproduction. The other leveraging tests
    (`test_golden_multirev_leveraging.py`, `test_leveraging_leg.py` and others) are moon-endgame goldens
    from Campagnola-Russell.
- **There is a notation clash.** `vilm.VilmLeg` (vilm.py lines 54-75) uses Campagnola-Russell's
  n:m_K±:
  - n:m is the same order as Sims's K:L (minor-body revolutions : spacecraft revolutions).
  - Campagnola-Russell's K is the number of full spacecraft revolutions on the H-B arc. Sims's K is the
    number of Earth revolutions.
  - Campagnola-Russell's ± marks long (H−) or short (H+) transfer. Sims's ± marks encounter after or
    before the line of apsides.
  - I did not check whether the two ± conventions map one-to-one.
  - A note in `vilm.py`, or in a docs note, would stop a future reader from taking a Rogers "4:3(2)−" label
    as a Campagnola n:m_K label.
- **Possible use, PROPOSAL only.** Tables 1 and 4 (analytic columns, circular Earth) are sourced goldens
  for a heliocentric delta-V-EGA solver. Such a solver would need the moon-only satellite lookup to be
  generalised to Sun-planet pairs first.
  - The Rogers-based establishment rows (for example `aldrin-4-3-2-establishment`, line 2799, with
    `v_infinity_leveraging_dv_kms: 0.568`) are this same problem with Mars.
  - These goldens would test the algorithm under those rows. They would not test the rows' values, which
    come from Rogers.

## 4. Citation mining

Held status was checked with `ls cyclers_pdf/papers | grep -i` and with `docs/notes/CORPUS_INDEX.md`.
Nothing cited is held, except as noted.

- **Not held, already on the wanted list:** none of the cited works. Row 12 (Strange & Sims 2001,
  AAS 01-437) is a later paper by the same author. It is not cited here.
- **Not held, new wanted candidates (VILT lineage):**
  - Hollenbeck, G. R. (1975), "New Flight Techniques for Outer Planet Missions", AAS 75-087 [1]. This is
    the origin of the delta-V-EGA.
  - Sims, J. A. & Longuski, J. M. (1994), "Analysis of V-infinity Leveraging for Interplanetary Missions",
    AIAA 94-3769 [3]. This is the first general use of the term, with a closed-form approximation.
  - Sims, J. A. (1996), "Delta-V Gravity-Assist Trajectory Design: Theory and Practice", PhD thesis,
    Purdue [7]. This has the full derivations.
  - Sweetser, T. H. (1993), "Jacobi's Integral and Delta-V-Earth-Gravity-Assist Trajectories",
    AAS 93-635 [6]. This is the Jacobi-integral delta-V estimate.
  - Lowest priority: AAS 95-306 [8], the conference version of this paper.
- **Not held, background, not added:**
  - Williams 1990, MS thesis [2] (STOUR automation). This is not the Longuski & Williams 1991 CMDA paper in
    wanted row 69.
  - Beckman & Smith 1973, AAS 73-231 [4] ("orbit pumping").
  - Roberts & Uphoff 1973, JPL EM 393-159 [5] ("orbit cranking"). This is not the Uphoff, Roberts &
    Friedman 1976 paper in wanted row 46.
  - Patel 1993, MS thesis [9]. This is not the held Patel, Longuski & Sims 1998.
  - Sauer 1989, MIDAS, JAS 37(3) [10].
  - The JPL launch-vehicle and cost documents [11, 12].
