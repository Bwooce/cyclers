# #1054: the #1000 seed-design fix: pre-registration (2026-10-11)

Owner: earthmoon-opus. Brief (lead, 2026-10-11): perigee at least 1.3 times the floor, omega from the
target lunar-pass distance, a perigee step limit; a test that reproduces the 93-94 % impact rate on the
old seeds and shows the new rate; the expected effect written before the run. No catalogue edits.

## 1. What the old records already say (no solver run; `scripts/run_1054_seed_design.py predict`)

Source: the #1000 grid (5,256 seeds that cross the lunar orbit) and the #1048 gaps (1,728 seeds),
`data/1054_seed_design/predict.json`. "impact@1" means the first arc of the uncorrected skeleton already
hits a body (before any Newton step); "impact@>1" means a later iteration.

| #1000 grid subset | seeds | impact | impact@1 | impact@>1 | converged | candidates per solve | families (of 15) |
|---|---|---|---|---|---|---|---|
| all | 5,256 | 0.590 | 0.286 | 0.304 | 0.286 | 0.0405 | 15 |
| seed perigee 6,578 km (the old "94 %" row) | 720 | 0.936 | 0.849 | 0.087 | 0.031 | 0.0014 | 1 |
| skeleton lunar pass < 5,000 km (the old "93 %" bin) | 393 | 0.934 | 0.789 | 0.145 | 0.053 | 0.0076 | 3 |
| r_p >= 1.3 x floor | 4,536 | 0.535 | 0.197 | 0.338 | 0.327 | 0.0467 | 15 |
| lunar pass 20,000-60,000 km (the #1048 suggestion) | 1,445 | 0.606 | 0.320 | 0.286 | 0.309 | 0.0609 | 9 |
| r_p >= 1.3 x floor and lunar pass >= 5,000 km | 4,201 | 0.504 | 0.153 | 0.352 | 0.348 | 0.0498 | 13 |

Findings that shape the design:
1. **The two 93-94 % bins are iteration-1 impacts.** 85 % (floor row) and 79 % (near-Moon bin) of
   their seeds hit a body on the first arc, before any Newton step. In a sample of 25 floor-row
   iteration-1 impacts all 25 hit the Earth: the lunar perturbation over one arc lowers the skeleton's
   6,578 km perigee below the surface. So the #1048 reading "the Newton steps push such a perigee into
   the Earth" is wrong for these bins; a step limit cannot touch them. (Correction to #1048 note sec. 4.)
2. The second perigee row of the old grid is already 8,578 km = 1.304 x floor. The "1.3 x floor" part
   therefore mostly removes the floor row, and it loses no family.
3. **The #1048 band suggestion (20,000-60,000 km) is a regression**: it keeps 9 of 15 families and
   drops the 7:3 ret wM-1 asymmetric family (the Casoliva 7-3b/7-3c family), 4:3 ret, 5:2 ret wM-1 and
   8:3 ret. Far passes (> 60,000 km) produce candidates too (4:3 ret comes only from them). The omega rule
   below therefore targets a log-spread of pass distances and keeps far passes.
4. After the seed selection, most remaining impacts are later-iteration impacts (0.352 of solves): the
   step limit's target.
5. Denominators: the #1048 "54 %" is impacts over all 5,760 grid seeds (3,101 / 5,760), including
   the 504 that never cross the lunar orbit and are not solved. Over solved seeds it is 59 %. This note
   uses solved seeds.

## 2. The fix (pre-registered)

All new behaviour is behind parameters whose defaults reproduce the old code path exactly
(`run_1000_shooting.shoot`, `solve_seed`; `run_1000_continue` and `run_1000_pitchfork` import `shoot`).

- **(F1) perigee rows:** r_p = geomspace(1.3 x 6,578.137, 42,164, 8) km.
- **(F2) omega from the target lunar-pass distance** (the p:q of `PQ`, whose skeleton crosses the lunar
  orbit): for each (p, q, sense, r_p), scan omega on a 1-degree grid with the vectorised #1048 skeleton
  distance d(omega); for each target d* in {5,000, 10,000, 20,000, 40,000, 80,000} km take every omega
  where d(omega) - d* changes sign (linear interpolation between grid points). A combination with no
  crossing for any target is skipped and counted. `PQ_GAPS` (skeleton never reaches the Moon) keep the
  uniform 10-degree omega grid with the F1 rows.
- **(F3) perigee step limit:** at each Newton iteration, after the existing 0.1 infinity-norm cap, each
  node's osculating Earth perigee (two-body, GM = 1 - mu, Earth-relative inertial velocity) after the
  step must be at least max(floor radius 6,578.137 km, 0.8 x its value before the step). Otherwise the
  step is halved, up to 10 times; after 10 halvings the halved step is taken. Nodes are at apocentres,
  so this is a limit on the osculating, not the actual, perigee.

## 3. Test (fast, deterministic; `tests/scripts/test_1054_seed_design.py`)

- **Reproduction:** a fixed set of old seeds (floor row and lunar pass < 5,000 km, seeded random
  choice) re-solved with the default path give exactly the status recorded in `grid.jsonl`, seed by
  seed.
- **First-arc impact rate** (the 93-94 % mechanism): fraction of seeds whose skeleton's first arcs hit a
  body, over the same (p, q, sense) set: old floor-row seeds vs F1+F2 seeds.
  - Expected: old >= 0.80 (records: 0.849 for impact@1; the test counts any arc of the skeleton, not
    only the first one integrated, so it may be a little higher); new <= 0.25 (records predict 0.153 for
    r_p >= 1.3 x floor with pass >= 5,000 km; F2 oversamples 5,000-20,000 km passes, whose impact@1 is
    0.315, so I allow up to 0.25).
- The 93-94 % bins do not exist under the new design, so the new result is stated as the overall
  per-solve rate plus the yield metrics in sec. 4.

## 4. Pilot ablation (in my session, chunks under 8 min, Pool(2)) and expected effect

A random sample (seed 1054) of 240 seeds per arm, from the #1000 PQ set, same (p, q, sense) mix:

| arm | seeds | step limit | expected impact | expected impact@1 | expected converged | expected candidates per solve |
|---|---|---|---|---|---|---|
| B | old | F3 | 0.45-0.57 (impact@>1 falls by at least a third, from 0.30) | 0.29 (unchanged; F3 cannot act at iteration 1) | 0.29-0.36 (smaller than the impact fall; some become max iterations) | 0.04-0.06 |
| C | F1+F2 | none | 0.50-0.62 | 0.15-0.30 | 0.27-0.38 | 0.05-0.10 |
| D | F1+F2 | F3 | <= 0.45 | as C | >= 0.35 | >= 0.05 |

Arm A (old seeds, old path) is the record itself; the reproduction test ties it to the current code.

- **Positive control (hard stop):** 7:3 seeds at Casoliva 7-3b's C, F1+F2+F3: at least one seed
  converges onto the 7-3b row (old control: 22 of 1,152 seeds); the rate per seed is reported.
- **Families:** the pilot is too small to count families. Losing a family is a regression; the full
  F1+F2+F3 grid (about 5,000-6,000 seeds, hours) with the family count against the 15 of #1000 is
  handed to the lead.
- What would falsify the fix: D's impact rate not below C's and B's; or the control not recovered;
  or (full grid) a #1000 family lost.
