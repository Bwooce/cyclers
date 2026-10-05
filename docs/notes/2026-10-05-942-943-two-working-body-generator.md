# #942 / #943: two-working-body cycler generator (Hollister date corrector), controls and enumeration

Status: BUILD AND CONTROLS DONE; ENUMERATION PENDING the lead's launch (sec. 6). No catalogue writes.
Novelty language: every enumerated cycler is "candidate, pending collision check" until sec. 7 clears it.

## 1. What was built

- `src/cyclerfinder/search/two_working_body.py`: the corrector.
  - Legs:
    - `LambertLeg`: transfers and same-body generic or "symmetric" returns.
    - `ResonantLeg`: an n:m full-revolution return. Its duration is fixed and its V-infinity direction
      is free on a circle.
    - `HalfRevLeg`: an n-pi return. It is analytic and has discrete directions.
  - Unknowns: the start dates of the Lambert legs. Residuals: the V-infinity magnitude difference at
    each massive junction, or the full vector difference at a massless target. This is Hollister &
    Menning 1970, pp.1194-1195. A square system, solved by least squares.
  - Free directions: chosen by minimax of the demanded/available turn ratio over the flybys of a
    block. Ties at inner flybys are broken by spreading the turn evenly. The tie-break cannot change
    the largest ratio, so it cannot change a gate verdict.
  - Every flyby goes through the #888/#937 gate (`verify/turn_gate.py`), on inertial-axis vectors.
    Wrap vectors are rotated in the circular model.
  - Near-180-degree rejection: PROVISIONAL threshold of 175 deg (sec. 8). The largest demanded turn is
    stored, so the number can change without a rerun.
  - Systems:
    - `CircularSystem`: the ideal model.
    - `MeanElementSystem`: Standish-Williams J2000 elements, inclined and elliptic, with optional
      commensurate periods.
  - Independent checks:
    - `encounter_self_consistency`: each leg is re-propagated by its own Kepler step, not the Lambert
      solver, and must land on its body to well inside the sphere of influence (#480 rule).
    - `kepler_step`: an eccentric-anomaly propagator. The core `propagate` fails to converge on some
      ordinary heliocentric arcs; this was reported to the lead with a reproducer.
- `src/cyclerfinder/search/two_working_body_enum.py`:
  - structures (cycle templates).
  - multi-start seeding spread across basins.
  - de-duplication.
  - assessment (gate at the project floor and at H&M's 1.1 radii, extents, sphere of influence).
  - placement of a massless target.
- `src/cyclerfinder/search/hollister_menning_1970.py`: the Table 3 loader, with transcription and print
  fixes (sec. 3.1).
- Scripts:
  - `scripts/run_942_hm_control.py`: the pre-registered H&M control.
  - `scripts/run_942_enumerate.py`: the sharded, resumable driver.
  - `scripts/analyse_942_enumeration.py`: merges physically identical zeros and mirror twins.
- Tests: `tests/search/test_two_working_body.py`.

## 2. Bend capacity ("massive enough to bend", #943 brief)

Available bend (deg) at the project floor. GM and radius are from R-S 2009 Table 2 for the moons and
from the registry for the planets.

| Body | Floor (km) | V_inf 1 | 2 | 3 | 4 | 5 | 6 km/s |
|---|---|---|---|---|---|---|---|
| Europa | 100 | 82.4 | 38.0 | 20.3 | 12.3 | 8.2 | 5.8 |
| Ganymede | 100 | 103.1 | 56.7 | 33.3 | 21.2 | 14.5 | 10.5 |
| Callisto | 200 | 94.4 | 48.1 | 27.1 | 16.9 | 11.4 | 8.1 |

| Body | Floor (km) | V_inf 3 | 4 | 5 | 6 | 7 | 8 km/s |
|---|---|---|---|---|---|---|---|
| Earth | 200 | 121.1 | 104.6 | 90.1 | 77.7 | 67.1 | 58.2 |
| Venus | 300 | 116.5 | 99.2 | 84.4 | 71.9 | 61.4 | 52.7 |
| Mars | 200 | 69.4 | 50.5 | 37.6 | 28.8 | 22.6 | 18.1 |

At the R-S V-infinities (1.7-4.3 km/s), Callisto and Europa bend tens of degrees. Rall's condition holds.

## 3. Positive control 1: Hollister & Menning 1970 Table 3 (pre-registered; FAILED as registered)

### 3.1 Source audit
- Table 3 is the INCLINED-ELLIPTIC 16-yr solution set on real dates (p.1194-1195). It is not
  circular-coplanar. The circular 3.2-yr orbits (Hollister 1969, not held) have no printed numbers.
  The 15 catalogue rows carry `model_assumption: circular-coplanar`; that looks like a mislabel.
  Reported, not edited.
- Transcription errors in `data/sources/hollister-menning-1970-table3.yaml`, checked on 450-dpi page
  images:
  - orbit 1 row 12 is "E 3163" (the YAML has V).
  - orbit 6 row 3 is "995" (the YAML has 993).
- Print errors in the paper (225-d Venus step broken):
  - orbit 2 "5715" (5815).
  - orbit 4 "2573" (2583).
  - orbit 5 "5662" (5652).
  - orbit 6 "2585, 2810" (2360, 2585).
  - orbit 8 "4870" (5870).
  - orbit 5 "977" (inferred 1017: a 265-d first gap gives an 8.9 km/s residual at the printed dates).
- Paper-internal: the text puts the 1.16-radius Earth pass in orbit 5. The table has it in orbit 4.
- 25 printed rows have a (V_r, theta, Rmin) triple that disagrees with r_p = mu/v^2 (1/sin(theta/2) - 1)
  by more than 10 %. In several of them our theta reproduces THEIR Rmin.

### 3.2 Results
Pre-registered criteria are in the header of `scripts/run_942_hm_control.py`:
- dates within 3 d.
- V_r within 0.005 EMOS.
- theta within 3 deg.
- Rmin within 10 %.
- for >= 90 % of encounters.

Stage results:
- Stage F (V_r at the printed dates): failed on all 15 orbits (0.32-0.78 within 0.006 EMOS).
- Stage C, real periods: 0/15. Dates slip about 0.3 d per Venus full revolution (224.70 d against the
  printed 225 d).
- Stage C on H&M's stated "exact periodicity" model, after amendments 1-2 (post hoc, labelled): 4/15
  pass. Earth 5844/16 d and Venus 5844/26 d; anchor JD 2443363.

| Orbit | Earth/Venus returns | Exact zero | Match | Match, source-consistent rows | Note |
|---|---|---|---|---|---|
| 1 | all full-rev | yes | 0.92 | 1.00 | Table 2 1H seed and perturbed seed reach the same zero |
| 2 | Venus FS | yes | 0.92 | 0.92 | |
| 3 | Earth SY x5 | yes | 0.24 | 0.32 | 7 zeros within 35 d; best 0.28 |
| 4 | Earth SY x1 | yes | 0.64 | 0.73 | |
| 5 | Earth SY x2 | yes | 0.44 | 0.48 | |
| 6 | Earth SY x2 | no (fold) | 0.16 | 0.16 | LS minimum 0.0057 EMOS, just above H&M's 0.005 |
| 7 | Earth SY x3 | yes | 0.44 | 0.50 | |
| 8 | Earth SY x4 | yes | 0.76 | 0.76 | |
| 9 | Venus S x1 | no (fold) | 0.92 | 0.96 | LS minimum 0.0015 EMOS, inside H&M's tolerance |
| 10 | Venus S x2 | yes | 0.84 | 0.91 | |
| 11 | Venus S x3 | yes | 0.96 | 0.96 | |
| 12 | Venus S x4 | yes | 0.96 | 0.96 | |
| 13 | Venus S + Earth SY | yes | 0.88 | 0.88 | |
| 14 | | yes | 0.80 | 0.91 | |
| 15 | | yes | 0.88 | 0.96 | |

Reading:
- The family with Earth full-revolution returns (orbits 1, 2, 9-15) reproduces.
- The Earth-symmetric-return family (orbits 3-8, H&M's "sequential modification" set) does not. Our
  model has nearby zeros that differ by about 0.01 EMOS.
- The homotopy from the printed dates folds for orbits 6 and 9.
- The control is FROZEN. A further model change chosen by match would be tuning.
- The open discrepancy goes to the lead and owner.

## 4. Positive control 2: recall in the circular model (production enumerator)

E-V, k = 2 (3.197 yr), production settings (n_phase 36, n_split 12, n_refine 40, basin-spread seeds):
- 1H (full-rev Earth, two full-rev Venus): 3 gate-passing zeros, e.g. V_inf E/V 2.99/3.19 and
  6.14/5.28 km/s.
- 2H (full-rev Earth, full-rev plus symmetric at Venus): 3 gate-passing zeros, e.g. 5.60/6.02 km/s.
- 3H (symmetric Earth, two full-rev Venus): 3 gate-passing zeros, e.g. 4.31/5.51 km/s. The symmetric
  return lasts 1.36 yr. H&M's sequential modification inserted a 1.37-yr one (p.1195); this is a
  consistency check only.

## 5. Positive controls 3-6: Russell & Strange (expected values from the papers' tables)

All four were run BLIND through the production enumerator. In each run the published cycler is the only
gate-passing physical cycler.

| Row | Cell | Ours | Published |
|---|---|---|---|
| GanCal#5 | gc1, k=3 | V_inf 3.238 / 3.34; alt 328.0 km; r 821,913-2,390,844 km | 3.24 / 3.34; 328; 821,915-2,390,844 |
| GanCal#1 | gc1, k=3, f(2:1) leg | 3.180 / 3.255; alt 247.3 km | 3.18 / 3.26; 247 |
| GanEur#43 | ge1, k=2 | 1.872 / 3.893; alt 8862 km; r 564,558-1,072,318 | 1.87 / 3.89; 8861; 564,558-1,072,330 |
| GanEur#5 | ge1, k=5 | 1.658 / 2.573; alt 1819.8 km | 1.66 / 2.57; 1819 |
| VenMar#45 | vm, k=2 | 8.220 / 12.965; alt 19,784 km; r_max 341,571,415 | 8.22 / 12.96; 19,784; 341,571,371 |

GanCal#1 contains a resonant f(2:1) leg, so its 247 km altitude also checks the minimax free-direction
convention against a published value. VenMar#45's published minimum distance (108,067,501 km) is the
conic's perihelion, which the leg does not pass; ours is the along-leg minimum.

## 6. Enumeration (PENDING)

Settings and pruning (first pass):
- ev and em: <= 2 returns per block, 0-rev transfers, generic 1-rev.
- vm: <= 3 Venus returns, transfers 0-1 rev, generic 1-2.
- gc and ge: <= 1 return per block, transfers 0-2 rev, generic 1-2.
- k = 1-3 in every cell.
- Catalogue: resonances (1:1, 2:1, 1:2, 3:2, 2:3); half-revs (1, 0, peri/apo), (3, 1, peri/apo);
  generic n-rev, low and high.

## 7. Literal-collision checks (to be completed per candidate)

Against:
- Hollister & Menning 1970 (catalogue rows hollister-menning-1970-ev-orbit-01..15).
- Jones 2017 VEM.
- R-S 2007 and 2009 rows: EurGan, GanEur, GanCal, GanIo, VenMar.
- Pisarevsky 2008: Table 4 class III, and the Fig. 14 class I.1 candidates for Earth-Mars.
- `data/empty_regions.jsonl`.
- `literature_check.py`.

## 8. Open items for the owner

- The H&M control discrepancy (orbits 3-8).
- The near-180 threshold (provisional 175 deg).
- The catalogue `model_assumption` of the 15 H&M rows.
- The two YAML transcription fixes.
