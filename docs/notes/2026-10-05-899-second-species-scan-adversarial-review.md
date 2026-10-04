# #899 second-species scan: adversarial review of the results (2026-10-05)

Scope: the full `#899` scan in `data/runlogs/899_scan/` (1060 seeds, resonances 7-3, 2-1, 1-2,
3-2; launched at 2b1924c0, crashed at seed 603 on the divide-by-zero fixed in 86fd5e95, resumed
at 86fd5e95; finished 2026-10-05 05:28 AEDT). Question: which published Earth-Moon orbits did it
reproduce, and is anything new. Prior for novelty: low (Casoliva et al. 2008/2010, Vaquero &
Howell 2013, Ross et al. 2025, Genova & Aldrin 2015, Arenstorf all hold planar Earth-Moon resonant
or cycler orbits).

Model: the planar circular restricted three-body problem at the paper's mass ratio
mu = 0.0121529529 (Casoliva 2010; the catalogue rows use the registry mu 0.0121505844, see the
`#899` STEP 2 entry for why the paper's mu is the reproduction target). Units for this note:
L = 384,400 km (Casoliva's own d_EM), 1 TU = 27.321661 d / 2 pi = 4.34838 d, Earth radius
6378.137 km, Moon radius 1737.4 km. Altitudes below are above the surface. Stability index on the
Casoliva scale, k = lambda + 1/lambda, critical at |k| = 2; k_par is the planar non-trivial pair,
k_perp the vertical pair.

Nothing here edits the catalogue, the runlogs or the source. Scratch scripts were run from the
session scratchpad and deleted afterwards; every number below is from those runs.

## 1. Inventory

Scan outcomes (runlog.jsonl, 1060 records): 811 seeds failed to correct (multiple shooting did not
converge), 1 failed in the propagator, 126 continued but stopped short of the target mu, 122
reached the target mu and were walked in C. Grid seeds (Eq. 14-17) close rarely, as Casoliva
report: of 200 grid seeds per resonance, 16 (1-2), 1 (2-1), 10 (3-2) and 10 (7-3) corrected at
mu = 1e-6 and were continued.

64 resonant hits (T = 2 pi q at the target mu): by resonance 2-1: 22, 3-2: 31, 7-3: 7, 1-2: 4;
by seed type table2: 1, walk (members of the 2008 Table 2 families sampled at mu = 1e-6): 45,
grid: 6, chain (two returning arcs): 12.

Duplicates. The scan's in-memory duplicate list was reset by the resume at seed 603, so its
`duplicate_of` field cannot be trusted across the restart. I re-did it independently: every hit
was re-integrated with `core.cr3bp.cr3bp_eom` (scipy DOP853, rtol = atol = 1e-13), all its y = 0
crossings found by event location, and the crossing farthest from both primaries compared with
every crossing of each group representative and of its mirror image (x, -y, -u, v). All 64 hits
collapse to **7 distinct orbits**; the largest crossing-to-crossing distance inside a group is
9.0e-11, and different groups differ in C by at least 0.11.

| # | orbit | hits | seed types | max in-group distance | mirrored members |
|---|---|---|---|---|---|
| 1 | 2-1 at C 0.48873531 | 22 | table2, walk | 3.4e-11 | 12 |
| 2 | 3-2 at C 0.70893304 | 13 | grid, walk | 1.3e-11 | 6 |
| 3 | 3-2 at C 0.37778607 | 13 | grid, walk | 9.0e-11 | 0 |
| 4 | 1-2 at C 1.47999190 | 2 | grid | 1.0e-12 | 0 |
| 5 | 1-2 at C 2.76298150 | 2 | grid | 6.5e-14 | 0 |
| 6 | 7-3 at C 1.06876239 | 7 | chain | 6.3e-12 | 0 |
| 7 | 3-2 at C 2.61330475 | 5 | chain | 5.2e-12 | 0 |

Orbit 6 is printed twice in Casoliva Table 3 (7-3b and 7-3c): the printed 7-3c crossing is the
mirror image of a 7-3b crossing, so 7-3b and 7-3c are one asymmetric orbit and its mirror, as the
paper itself says ("Flipping cycler 7-3b by 180 deg with respect to x axis yields cycler 7-3c").

## 2. The distinct orbits

Distances are the max-abs state difference between the printed Table 3 crossing (converted to
this project's frame, (-x_i, 0, -u_i, -v_i)) and the nearest y = 0 crossing of the orbit or of its
mirror, all 16 printed rows of every resonance compared, at the paper's mu (where the printed rows
themselves close to about 1e-10). "Independent closure" is my DOP853 re-integration of the stored
state over one period (max-abs end minus start), with the Jacobi drift over the period.

| # | class | C | T (TU / d) | k_par | k_perp | periselene (nd / alt km) | perigee (nd / alt km) | nearest published orbit, distance | independent closure / Jacobi drift |
|---|---|---|---|---|---|---|---|---|---|
| 1 | (a) | 0.4887353098 | 2 pi / 27.32 | 1.2822 | 1.9256 | 0.235136 / 88,649 | 0.187670 / 65,762 | Casoliva 2-1a, 4.2e-11 (catalogued) | 2.1e-9 / 1.0e-12 |
| 2 | (a) | 0.7089330386 | 4 pi / 54.64 | -1.2467 | 1.8752 | 0.219367 / 82,587 | 0.039231 / 8,702 | Casoliva 3-2c, 1.1e-10 (catalogued) | 2.9e-9 / 3.3e-12 |
| 3 | (b) | 0.3777860739 | 4 pi / 54.64 | 414.79 | -118.66 | 0.040613 / 13,874 | 0.073886 / 22,024 | Casoliva 32c family (Fig. 3c), past the plotted end; nearest printed row 3-2c, 0.19 (dC -0.33) | 2.9e-9 / 4.0e-12 |
| 4 | (a) | 1.4799919040 | 4 pi / 54.64 | 0.5333 | 1.9976 | 0.536322 / 204,425 | 0.090672 / 28,476 | Casoliva 1-2b, 4.8e-11 (printed, not catalogued) | 5.0e-13 / 4.0e-13 |
| 5 | (a) | 2.7629814961 | 4 pi / 54.64 | -4.1925 | 1.9998 | 0.698033 / 266,587 | 0.748961 / 281,522 | Casoliva 1-2e, 5.4e-11 (catalogued) | 2.5e-13 / 3.3e-13 |
| 6 | (a) | 1.0687623900 | 6 pi / 81.96 | 57.352 | -2.2709 | 0.034365 / 11,473 | 0.040855 / 9,326 | Casoliva 7-3c 5.2e-11 direct, 7-3b 1.2e-10 (both catalogued) | 1.3e-10 / 1.7e-12 |
| 7 | (c) | 2.6133047523 | 4 pi / 54.64 | -10638.6 | 11.715 | 0.009978 / 2,098 | 0.126868 / 42,390 | none held; nearest printed rows 7-3c 0.11, 3-2d 0.49 at dC +0.96 | 1.6e-10 / 4.5e-13 |

Per-class counts: (a) reproduction of a printed row: **5 orbits, 6 printed rows** (2-1a, 3-2c,
1-2b, 1-2e, 7-3b = 7-3c mirror); (b) member of a published family, not printed: **1** (orbit 3);
(c) not matched to anything held: **1** (orbit 7); (d) artefact: **0**.

Other catalogued Earth-Moon rows (Ross et al. 2025 and Braik & Ross 2026 (k1,k2) cyclers and
corridors, Vaquero 2:1 and 3:1 members, Arenstorf, Genova & Aldrin, Wittal) are not close to any of
the seven: their periods are 2.1 to 31 TU but never 2 pi q with these C values (Ross and Braik-Ross
rows sit at C 3.13 to 3.18, the Vaquero rows at C 1.98 to 3.13 with T 5.7 to 6.5 TU; the nearest
in period is Braik-Ross R52-S, T 12.603 at C 3.1294, a 5:2 orbit), and Arenstorf, Genova-Aldrin
and Wittal hold no state. Different (C, T) at the same mu means a different orbit, so the
catalogue distance is not a state distance there.

## 3. Independent re-verification (orbits 3 and 7, and the five reproductions)

Code: `core.cr3bp.cr3bp_eom` and `cr3bp_stm_eom` with scipy `solve_ivp` DOP853 at rtol = atol =
1e-13; no `second_species_*` module in the verification path. Every lunar pass is above 1e-3 lunar
distances (the closest is 0.00998 on orbit 7), so the KS propagator was not needed.

- Closure over the period: 2.5e-13 to 2.9e-9 for all seven (table). Jacobi drift at most 4e-12.
- Minimal period: the closest return to the start state for t in (0.02 T, 0.98 T) is 0.08 (orbit
  6) to 4.2 (orbit 2) in max-abs state; none is a repeated shorter orbit.
- Monodromy: the trivial pair (the two eigenvalues nearest +1) is 1 +- 1e-4 for all seven; the
  non-trivial planar pair and the vertical pair have products 1 to within 7e-7 (orbit 3, the
  largest k) and 2e-8 elsewhere. k_par and k_perp agree with the scan's values to 6e-4 relative
  on orbit 3 and better than 2e-6 on the others.
- Inertial revolutions about the Earth per period (rotating angle plus t): 2-1a -2, 3-2c -3,
  orbit 3 -3, 1-2b +1, 1-2e +1, 7-3b/c -7, orbit 7 +3. So 2-1a, 3-2c, orbit 3 and 7-3b/c are
  retrograde about the Earth, and orbit 7 is a prograde 3:2 orbit. The labels p-q match the
  winding for all seven.
- Surfaces: every periselene is above 1737 km (lowest 2,098 km altitude, orbit 7) and every
  perigee above 6378 km with margin (lowest 8,702 km altitude, 3-2c).
- Symmetry (time-reversal mirror about the x axis): 2-1a, 3-2c and orbit 3 are symmetric
  (mirror distance 1.5e-12 to 2.7e-11, perpendicular crossings to 1e-11); 1-2b, 1-2e, 7-3b/c and
  orbit 7 are asymmetric (0.10 to 0.27). Casoliva state that 1-2e and 7-3b/c are asymmetric: a
  check passed.

## 4. Orbit 3 (3-2, C 0.3778): class (b)

The walk that finds 3-2c at C 0.709 is one characteristic curve at the target mu. Going down in C
it folds at C -0.372 (T 13.21) and returns on an upper branch, whose T crosses 4 pi again at
C 0.3778: orbit 3. Same seed families, same walks (32a seeds and the grid seed at C -0.267 give
both hits). Casoliva Fig. 3c (page image read 2026-10-05): the dashed 32c segment begins near
C -0.4, at "resonance relation" about 1.42, and rises through 1.5 at C about 0.7. Taking the plotted
ratio as 1.5 times 4 pi / T (it reproduces 3-2a: 1.5 x 12.566 / 12.636 = 1.4917 at the end of the
32a segment, and the 32d segment computed from the printed 3-2d state below), our fold at
C -0.372 has ratio 1.427: the plotted 32c segment stops at this fold, and orbit 3 is past it on
the same curve. The 32b segment (C -0.5 to 0, ratio 1.32 to 1.45) is not this branch; our upper
branch never goes below C -0.372.

So orbit 3 is a member of the published 32c family, not printed and not plotted. It is very
unstable (k_par 414.8, k_perp -118.7: unstable in and out of the plane), symmetric and retrograde,
and so very likely in the Franz & Russell 2022 Earth-Moon symmetric-orbit database (public on
Zenodo, DOI 10.5281/zenodo.6411980, not held as data, mu differs from Casoliva's at about 2e-4
relative). That check was not run. It does not change the (b) label.

## 5. Orbit 7 (3-2, C 2.6133): class (c), not novelty-claimable

Numbers: C 2.6133047523, T 4 pi = 12.566 TU = 54.64 d, prograde, three revolutions about the
Earth and two lunar passes per period (periselene 0.00998 = 2,098 km altitude, and 0.0368 =
12,399 km), perigee 0.1269 (42,390 km altitude), asymmetric, k_par -10638.6, k_perp 11.72. Origin:
the two-arc chain seed (1,1,+1)+(1,2,-1) at mu = 1e-6 (its generating orbit passes the Moon at
5.6e-7, inside the Moon scaled to mu = 1e-6, a mathematical second-species generator), continued
in mu at fixed C 2.611 and walked to T = 4 pi at the target mu; five chain seeds (C 2.611 to 2.833)
reach the same orbit to 5e-12.

Attempts to break it:

1. **Printed rows.** Nearest of all 16 Table 3 rows: 0.11 (7-3c, a different resonance and C).
2. **Casoliva 32d family.** 3-2d (C 1.6506, symmetric, prograde, flies through the Earth) was
   corrected from its printed state and continued in C both ways at the paper's mu with the Earth
   stop disabled. The 32d curve reproduces Fig. 3c's 32d segment (ratio 1.500 at 1.65 rising to
   1.68 near C 2.84, where the paper's curve ends at about 1.67). At C 2.598 to 2.631 it has
   T 12.14 to 12.20 and lies 0.39 from orbit 7. Not the same family.
3. **Branch of a symmetric family.** Asymmetric families are born at pitchforks of symmetric ones,
   so orbit 7 could be a branch of a published symmetric 3:2 family (Ross et al. 2025 (3,2),
   Vaquero 2013 Table 3.5 3:2 member at C 2.85675, T 54.55 d, both symmetric and prograde). Its
   family was continued at the paper's mu in both directions, measuring the mirror distance of the
   crossing set on each member: it stays between 0.089 and 0.61 along the whole curve, which runs
   from Earth collision (perigee 1e-4 at C 1.815) up through a fold at C 2.7434 (k_par passes
   through +2 there: a fold, the curve turns back, it does not split) and down again to Earth
   collision (perigee 3e-4 at C 1.659). No symmetric member, so no pitchfork, on the walked
   curve. The curve never reaches C 3.13 (Ross) or 2.857 (Vaquero), and its members are
   asymmetric. Not checked: the continuation through Earth collision (needs Earth
   regularisation), and the Vaquero 3:2 symmetric family itself (not continued).
4. **2008 seed 32b.** Casoliva 2008 say the 32b seed (C 2.0635 at mu = 1e-6) cannot be continued
   to the Earth-Moon mu. Orbit 7's generating family at mu = 1e-6 was continued in C to below
   2.06: at C 2.056 to 2.068 it is 1.2 from the 32b state. Not 32b.
5. **Literature check.** `search/literature_check.py` with the offline corpus search: for
   sequence ("Moon", "Moon") it returned "published", citing Vasile & Campagnola 2009 (a
   Jupiter-moon tour); for ("E", "Moon") "inconclusive", citing the Aldrin Earth-Mars cycler. Both
   are token matches on body names (the `#647` limitation stated in the module: it is built for
   cycler vocabulary, cannot see states, and does not apply to raw CR3BP resonant periodic
   orbits). It gives no evidence either way. Web search was not run. The Vaquero 2013 and Ross et
   al. 2025 digests and the Ross & Roberts-Tsoukkas 2026 golden table were read for families that
   could contain it: none does, as above (all their tabulated families are symmetric, at C 3.13
   or above for Ross).
6. **Turn gate (the `#888` lesson).** At each lunar pass I compared the Moon-relative inertial
   velocity turn between entry and exit of the lunar sphere of influence (66,000 km) with the
   two-body hyperbolic turn available at the pass radius, 2 asin(1/e). Orbit 7: pass 1 demands
   90.6 degrees against 67.0 available, pass 2 77.8 against 100.0; orbit 3: 23.8 against 12.6.
   **Positive control fails:** the published 7-3b/c orbit, reproduced at 5e-11, demands 32.1
   against 19.1 at its closest pass. So this patched-conic test is not a discriminator for a
   solution of the full restricted problem: the Earth's pull during a multi-day transit of the
   sphere is part of the turn. These orbits are ballistic by construction (closure 1.6e-10,
   Jacobi drift 4.5e-13). The test counts neither for nor against orbits 3 and 7.
7. **Sun sensitivity.** The model is the circular restricted problem only. Orbit 7's largest
   multiplier is 1.06e4 per 54.6 d, an error-doubling time of 4.1 d (orbit 3: 6.3 d; 7-3b/c:
   14.0 d; 1-2e: 27.6 d; 2-1a, 3-2c and 1-2b: no exponential growth). The Sun's tidal
   acceleration in the Earth-Moon system is of order (1/13.4)^2, about 6e-3 of the Earth-Moon
   mutual term, and orbit 7 multiplies a perturbation by about 1e4 per period, so the CR3BP orbit
   says nothing about an operational trajectory. A Sun-forced counterpart (as `#905` does for the
   forced models) is needed before any claim, and with a 4-day doubling time it would probably be
   a different, still more unstable object or none at all.

Verdict on orbit 7: (c), not matched to anything held, **not novelty-claimable**: CR3BP only,
extremely unstable (k about -1e4), asymmetric families of this kind are generic (any p-q
resonance has many second-species families), the literature check is blind to it, and the one
public database that might hold planar Earth-Moon orbits at this level (Franz & Russell 2022)
covers symmetric orbits only, so it can neither confirm nor clear it. Absence from the held
corpus is conditional on the corpus (Broucke 1968, Henon's second-species catalogues and the
Leiva & Briozzo 2006a atlas are not held or not tabulated).

## 6. Published rows: reached and not reached

| row | status after the scan | nearest hit distance | note |
|---|---|---|---|
| 2-1a | reproduced (step 2 and scan) | 4.2e-11 | |
| 3-2c | reproduced (step 2 and scan) | 1.1e-10 | |
| 7-3b / 7-3c | reproduced as one orbit and its mirror | 1.2e-10 / 5.2e-11 | |
| 1-2e | **now reproduced** (grid seeds C 1.402 and 1.590, branch 1) | 5.4e-11 | step 2 had it as not reproduced (inferred impact); the scan overturns that |
| 1-2b | **reproduced** (same walks) | 4.8e-11 | printed, not catalogued |
| 7-3a | not reached | 0.071 (orbit 6) | |
| 1-2c | not reached | 0.62 | the 1-2 walks covered C 1.34 to 3.14 at the target mu on the family through 1-2b and 1-2e; 1-2c and 1-2d are on other segments (Fig. 3a shows five) |
| 1-2d | not reached | 0.56 | as above |
| 2-1b | not reached | 0.69 | the 21a sheet walks span C 0.06 to 1.48 and give only 2-1a |
| 1-2a, 3-2a, 2-1c, 2-1d, 7-3d | not reached (uncatalogued, flagged by Casoliva) | 0.07 to 0.69 | |
| 3-2d | not reached by the scan | 0.49 | corrected only from its printed state in this review (section 5) |

So the scan moved one catalogued row (1-2e) from not reproduced to reproduced and added one
printed, uncatalogued row (1-2b). It did not change the status of 7-3a, 1-2c, 1-2d or 2-1b.
Catalogued Casoliva rows: 6 of 9 now reproduced forward from mu = 1e-6 seeds (2-1a, 3-2c, 7-3b,
7-3c, 1-2e, and 1-2b besides).

Two side findings:

- **1-2a and 1-2b flags.** `TABLE3_ROWS` carries `satisfies_resonance=False` for 1-2a and 1-2b.
  The Table 3 page image does print footnote e ("This cycler trajectory does not satisfy its
  resonance relation") on both, so the transcription is right. But 1-2b closes at T = 4 pi
  exactly (reproduced at 4.8e-11 with one Earth revolution per period), the printed T of both is
  12.5663706144, and the text's own example ("see how the period for 1-2a* is different from
  2 pi q") does not match the printed period. The flag's meaning is unclear from the paper. The
  module was not edited.
- **1-2e is planar-unstable.** k_par -4.19 (planar multipliers -3.94 and -0.254); the printed k
  1.99978 is k_perp, the `#801` block choice. Casoliva list 1-2e among the stable cyclers. Known
  from `#801`, restated here because the scan reproduces it independently.

## 7. Recommendation

- No catalogue action. Orbits 1, 2, 4, 5 and 6 are reproductions; the catalogue rows for 2-1a,
  3-2c, 1-2e, 7-3b and 7-3c could cite this scan as forward reproduction evidence when the rows
  are next touched (owner's call; not done here).
- Orbit 3: record as a member of Casoliva's 32c family past the plotted fold; no row.
- Orbit 7: record as "(c) not matched to anything held, not novelty-claimable". Before anyone
  spends more on it: (i) continue Vaquero's 3:2 symmetric family (Table 3.5, C 2.85675) down in C
  and look for a pitchfork, (ii) query the Franz & Russell database only for symmetric
  neighbours, (iii) run a Sun-forced continuation (`#905` style). Given k about -1e4 the expected
  value of (iii) is low.
- The turn gate needs a positive control before it is used on CR3BP orbits: as built here it
  fails the published 7-3b/c orbit.
