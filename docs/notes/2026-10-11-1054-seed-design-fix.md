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

## 5. Results (2026-10-11, 12:03-13:05 AEDT; load average 35-69 during the runs)

Code: b977513f (F3 in `run_1000_shooting.shoot`/`solve_values`, defaults keep the #1000 path; F1+F2
generator and the arms in `run_1054_seed_design.py`). F2 gives 7,252 seeds for the PQ set (old grid:
5,256 solved), no combination skipped. Data: `data/1054_seed_design/{arm_B,arm_C,arm_D,control}.jsonl`,
`report.json`.

### 5.1 Test (`tests/scripts/test_1054_seed_design.py`)

- Reproduction: 20 of 20 old seeds give the recorded status and iteration count. Pass.
- F3 helper: osculating perigee of a skeleton node equals r_p (rel 1e-12). Pass.
- Skeleton-arc impact rate: old floor row **0.900** (>= 0.80: reproduces the 93-94 % bins); new
  design **0.300**, above the pre-registered bound 0.25. **FAIL.** A larger diagnostic sample (400 new
  seeds, seed 7) gives 0.32: 90 Earth hits and 37 Moon hits. The Earth hits concentrate in the lowest F1
  rows (8,551 km: 26 of 60) and at the 10,000-20,000 km lunar-pass targets. A close lunar pass lowers
  the perigee of a 8,551-13,000 km skeleton below the surface within one arc, so 1.3 x floor is not
  enough when the pass is close. The Moon hits sit at the 5,000 km target (23 of 92).
- The test file is held uncommitted (it fails as pre-registered); ruling requested.

### 5.2 Pilot ablation (240 seeds per arm)

| arm | impact | impact@1 | impact@>1 | converged | other (max iterations) | candidates per solve | families in the sample | mean step halvings |
|---|---|---|---|---|---|---|---|---|
| A: old seeds, old path (records) | 0.496 | 0.225 | 0.271 | 0.333 | 0.171 | 0.0583 (14) | 6 | 0 |
| B: old seeds + F3 | 0.392 | 0.225 | 0.167 | 0.317 | 0.292 | 0.0417 (10) | 4 | 24.8 |
| C: F1+F2, no F3 | 0.683 | 0.312 | 0.371 | 0.237 | 0.079 | 0.0583 (14) | 6 | 0 |
| D: F1+F2+F3 | 0.588 | 0.312 | 0.275 | 0.237 | 0.175 | 0.0583 (14) | 6 | 16.7 |

Against the pre-registered expectations (sec. 4):
- B: impact@1 unchanged (as expected; F3 cannot act at iteration 1); impact@>1 falls by 38 % (expected
  at least a third: met); impact 0.392 (expected 0.45-0.57: lower than expected); converged 0.317 (in
  range, A's 0.333); max iterations rise from 0.171 to 0.292. F3 turns later-iteration impacts
  into max-iteration stalls, not into converged solves. Candidates 10 against A's 14: a difference
  under 1 sigma at these counts, so no detectable yield change.
- C: impact 0.683 (expected 0.50-0.62: **missed**), impact@1 0.312 (expected 0.15-0.30: missed),
  converged 0.237 (expected 0.27-0.38: **missed**), candidates per solve 0.0583 (expected 0.05-0.10:
  met, equal to A).
- D: impact 0.588 (expected <= 0.45: **missed**), converged 0.237 (expected >= 0.35: **missed**),
  candidates per solve 0.0583 (expected >= 0.05: met, equal to A).
- **Falsification criterion met:** D's impact rate (0.588) is not below B's (0.392), a gap of about
  4 sigma at n = 240. The fix as pre-registered does not do what it was meant to do.
- Statistical weight: candidates per solve and family counts show no detectable difference between
  any arm at n = 240 (counts 10-14). The converged gap (C and D 0.237 against A 0.333) is about
  2 sigma: moderate evidence only.

### 5.3 Positive control (7:3 seeds at 7-3b's C, F1+F2+F3)

721 of 1,070 seeds solved (all retrograde seeds of the F1 rows 0-2 and the prograde ones run before the
reorder): 27 converged, **8 converge onto the 7-3b row**. Pass (hard stop cleared). Old control: 22 of
1,152. The remaining 349 seeds are not run; the per-seed rate is not final.

### 5.4 Reading

- The 93-94 % bins are a seed-placement effect at iteration 1 (sec. 1). Removing the floor row (F1)
  helps and loses no family in the old records.
- F2 as designed (log-spread targets from 5,000 km) oversamples close lunar passes. Close passes with a
  low perigee give first-arc Earth hits, and passes at 5,000 km give Moon hits. Under F2 the impact rate rises.
- F3 lowers later-iteration impacts but does not raise convergence. Hypothesis (not measured): B
  averages 24.8 halvings per solve, which fits F3 often using its 10-halving budget within the
  25-iteration cap, so the stalls may come from the limiter and max_it rather than from the Newton
  direction. Untested check: rerun B's stalled seeds with a higher max_it.
- Not tried (would need a new pre-registration): F1 alone at full grid scale; F2 with targets from
  20,000 km plus the uniform grid for far passes; a perigee floor that rises with the closeness of the
  lunar pass.

### 5.5 Hand-off withdrawn

The full F1+F2+F3 grid that sec. 4 handed to the lead is withdrawn: it is not worth the compute
until a ruling on the design.

## 6. Lead ruling 2026-10-11 and pre-registration of F1 alone (written before any F1 run)

Ruling: (a) now, (b) next. (a) The test keeps the reproduction checks and the old floor-row bound
(>= 0.80). The F1+F2 design's 0.300 stays in sec. 5.1 with no bound. F3 stays in the code as an
opt-in, default off, with its own test: three arm-B seeds re-solved with F3 on give their recorded
status, iteration count and halving count. (b) F1 alone is the revised design. The pilot starts
only after the lead reports the catalogue window green.

**Design E (F1 alone):** perigee rows RP_NEW = geomspace(1.3 x 6,578.137, 42,164, 8) km = 8,552,
10,741, 13,490, 16,943, 21,281, 26,728, 33,570, 42,164 km. The uniform 10-degree omega grid. No F2,
no F3.

**Pilot (paired):** arm E is arm A's 240 seeds (the B sample: seed 1054, the #1000 PQ set). Each seed
keeps its (p, q, sense, row index, omega), and its perigee moves from RP[ir] to RP_NEW[ir]. Arm A is
the #1000 record of those seeds (impact 0.496, converged 0.333, candidates per solve 0.0583). A seed
whose F1 skeleton no longer crosses the lunar orbit drops out of the rates and is counted.

**Expected effect, from the old records** (impact rate per row, interpolated in log r_p to the F1
rows; per-row record rates 0.936, 0.801, 0.667, 0.556, 0.480, 0.424, 0.368, 0.368):
- The population model gives A 0.577 and E 0.526: a paired drop of about 0.05. Most of it comes from
  A's 30 floor-row seeds (28 impacts at 6,578 km; at 8,552 km the records give 0.80, about 24).
- A's sample runs below the model (0.496 against 0.577), so the expected E is about 0.44 (range
  0.38-0.50).
- Converged: expected to rise by about the impact drop, to about 0.35-0.40.
- Candidates per solve: no prediction that n = 240 can test (counts 10-14 per arm, sec. 5.2).

**Bound (from the old records):** E impact <= 0.535, the records' rate for r_p >= 1.3 x floor.

**Falsification (any one):**
1. E impact > 0.535.
2. E impact not below A's 0.496 (no paired reduction).
3. E converged < 0.273 (A's 0.333 minus 2 SE, SE 0.030).
4. Positive control: no seed converges onto the 7-3b row.

Power: the expected paired drop of about 0.05 is about 1.5 sigma at n = 240. Criterion 2 can fail by
chance about one time in ten even if F1 works as the records say. If it fails alone, I report it as
weak evidence, not as a refutation.

**Positive control (hard stop):** 7:3 seeds at 7-3b's C, F1 rows x 36 omegas, both senses, no F3:
576 seeds, retrograde first. At least one seed must converge onto the 7-3b row. The per-seed rate is
reported against the old 22 / 1,152.

**After the pilot:** if E passes, the full F1 grid (about 5,000 seeds, hours) goes to the lead,
checkpointed. At full scale the test is that none of the 15 #1000 families is lost.

**Candidate next design (registered, not designed yet):** a perigee floor that rises as the
skeleton's lunar pass gets closer: r_p,min(d) larger for smaller d. This is the principled form of
F2. Sec. 5.1 shows that 8,551-13,000 km perigees with a close lunar pass give first-arc Earth hits.
It is used only if F1 also fails.
