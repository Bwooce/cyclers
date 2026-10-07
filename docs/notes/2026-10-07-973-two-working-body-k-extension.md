# #973: the two-working-body generator at longer periods (gc k = 4-6, ev k = 4-5) and a Europa-Callisto cell (ec k = 1-4)

Status: DONE for gc k = 4 and ec k = 1-4. Pre-registration secs. 1-2 (commit ab24925e); amendment 2.6
(two screens, commit 94dda754); controls and Liang segments sec. 3; gc k = 4 sec. 4, re-screened in sec. 9;
final candidate set sec. 10: 10 clean two-working-body members at k = 4. ec k = 1-4: empty (sec. 5,
stamped). gc k = 5-6 and ev k = 4-5 are launched by the lead (sec. 6). The ev screen question is ruled (amendment
2.7; sec. 9.2). The literature step and rung (d) are #1025.
Novelty language: every gate-passing cycler here is "candidate, literature step deferred". Nothing in
this note is called novel. No catalogue writes. No real-ephemeris runs.

Source of the routes: `docs/notes/2026-10-07-971-fable-corpus-review-2.md` sec. 3, R11 (gc), R12 (ec),
R13 (ev). Generator, settings and gauntlet: `docs/notes/2026-10-05-942-943-two-working-body-generator.md`
secs. 1-6.1 (fixes in 6.3, 6.23, 6.30, 6.34 are in the code used here).

## 1. What runs, and what is deferred

Driver: `scripts/run_973_enumerate.py` (declares `#973` to the preflight gate). Subcommands:
- `enumerate`: the run_942 driver unchanged (`scripts/run_942_enumerate.py::main`, same flags and
  outputs), with this task's preflight and one more cell, `ec`.
- `recall`: the controls of sec. 2.3, solved with the production solver and seeds.
- `liang`: the Liang et al. 2024 open segments of sec. 2.4.
- `gauntlet`: `scripts/gauntlet_942.py`'s DOP853 re-fly, SOI check and literal-collision functions,
  imported unchanged, WITHOUT its literature step, plus a Liang 2024 Tables 3/5/7 per-moon V_inf check
  and a #576 symmetric-closure check.

Where the `ec` cell lives: the brief named `src/cyclerfinder/search/two_working_body_enum.py` for it, but
the cell table is `cell_system` in `scripts/run_942_enumerate.py` (a shared file owned by #942). The
`ec` cell is therefore defined in the #973 driver, which wraps that function. No `src/` file changes in
#973, so the mypy step of the brief does not apply. A control-only cell `ec576` (sec. 2.3) is defined the
same way.

No k cap exists in the generator or the driver (`--k` takes any integer; `structures()` prunes only
templates whose fixed legs plus minimum Lambert times exceed the period). Nothing was raised.

DEFERRED (not run in #973):
- The scripted literature step (`literature_check.py`) and the `#875` policy reading. They run after
  `#972` lands; `literature_check.py` and its callers are being edited for `#972`. The gauntlet module
  is loaded with a stub in place of `literature_check`; the stub raises if called.
- The real-ephemeris rung (d) (`scripts/run_942_realeph_chain.py`, also under edit for `#972`).
- The registry stamp: `data/empty_regions.jsonl` is not edited here; the proposed stamp text goes in
  the results sections for the lead.

## 2. PRE-REGISTRATION (2026-10-07, before any production run)

### 2.1 Cells and settings

Every cell: circular coplanar ideal model, structures `[A-block, A->B, B-block, B->A]`, one visit per
cycle, period k synodic periods. Return catalogue per body as in generator note 6.1 (resonant 1:1, 2:1,
1:2, 3:2, 2:3; half-rev (1,0,p/a), (3,1,p/a), both mirror signs; generic same-body Lambert legs).
Seeds n_phase 36, n_split 12, n_refine 40, minimum separation 0.03 T (the settings at which every
#942/#943 recall control passed).

| Route | Cell | Bodies (both massive) | Constants | k | Returns (A,B) | Transfer revs | Generic revs | Other | Structures |
|---|---|---|---|---|---|---|---|---|---|
| R11 | gc | Ganymede, Callisto | R-S 2009 Table 2 | 4 / 5 / 6 | 1,1 | 0,1,2 | 1,2 | none | 4,068 / 4,589 / 4,844 |
| R13 | ev | Earth, Venus | Venus 0.61520 yr (generator note 6.1) | 4 / 5 | 2,2 | 0 | 1 | `--resonant-only E` | 3,566 / 4,070 |
| R12 | ec | Europa, Callisto | R-S 2009 Table 2 (p.148) | 1 / 2 / 3 / 4 | 1,1 | 0,1,2 | 1,2 | none | 37 / 177 / 470 / 1,330 |

- gc and ev use exactly the #943 cell-3 and #942 cell-5 flags (registry stamps
  `jupiter-ganymede-callisto-two-working-body-rs2009-ideal-k1-3-943` and
  `heliocentric-earth-venus-two-working-body-earth-fullrev-only-ideal-k1-3-942`), only k changes.
- ec uses the gc flags (the gc pattern, as R12 asks). Europa is body A (inner), as Ganymede is in gc.
  Flyby constants: GM and radius from R-S 2009 Table 2 (Europa 3202.74 km^3/s^2, 1561 km; Callisto
  7179.29, 2408 km), floors from the registry (Europa 100 km, Callisto 200 km). Synodic period 4.512 d.
- The 971 note's "variant with 2-3 revolutions on the inter-moon legs (Liang's regime)" is NOT in this
  pre-registration: `--transfer-revs 0,1,2` already includes 2 revolutions; a 3-rev variant would be a
  separate, later run.
- Timing pilots (seeded `--sample 24`/`40 --seed 973`, scratch only, not results): gc k4 0.21 s per
  structure, gc k6 0.78 s, ev k5 0.91 s, ec k4 0.23 s. Each cell/k is split into shards
  (`--shard i/n`, own output directory each) so that no tool call runs longer than 10 min; two shards
  run at a time (`xargs -P 2`, foreground).
- Output: `data/973_<cell>/k<k>/s<i>/` (structures.jsonl, zeros.jsonl, settings.json, log).

### 2.2 Pass definition (unchanged from generator note 6.1)

A gate-passing zero: max |residual| < 1e-8 km/s; #888/#937 gate "pass" at every massive flyby at the
registry floors ("indeterminate" is not a pass); no demanded turn >= 175 deg; independent Kepler-step
re-propagation miss < 1 km at every encounter. Physical cyclers: zeros merged by cyclic flyby sequence
(body, V_inf to 1 m/s, turn to 0.1 deg), mirror twins merged (`scripts/analyse_942_enumeration.py`).

Gauntlet for each gate-passing physical cycler (sec. 1): DOP853 re-fly of every leg (rtol 1e-13);
reported arrival miss, V_inf vector error, junction mismatch and the gate on the integrated vectors;
SOI fraction; literal collisions: R-S 2007/2009 Table 3 rows, Campagnola 2019 GCGC, Hollister 1969
orbits I-III (ev, k = 2 only), catalogue rows on the pair (this covers the 15 H&M 1970 rows and Jones
2017), Liang 2024 Tables 3/5/7 per-moon V_inf (NEAR within 0.3 km/s), the #576 symmetric closures of the
pair (NEAR within 0.3 km/s). Re-fly agreement for the candidate table: arrival miss < 1 km and the gate
still "pass" on the integrated vectors.

ev only: Hollister & Menning 1970 orbits are 16-yr (k = 10) concatenations of 3.2-yr blocks. Before any
k = 4 or k = 5 ev candidate is described as outside H&M, the 15 Table 3 orbits are checked for a 4- or
5-synodic sub-period (`hollister_menning_1970.load_table3` / `block_types`).

Coverage: for each cell/k the number of structures whose zero count reached n_refine (40) is reported.
A structure at that ceiling may have more zeros than were refined. "No gate-passer at k" is conditional
on these seeds and this catalogue.

### 2.3 Controls (recall), tolerances fixed here

Production solver and seeds, each structure solved on its own (`recall` subcommand). A control passes
when a zero of the stated structure has V_inf within the tolerance at both bodies. The gate status of
the recalled zero is reported separately; recall means "reappears as a zero".

| Control | Cell | Structure key | Expected V_inf (km/s) | Source of expected | Tol. (km/s) |
|---|---|---|---|---|---|
| gc-1 | gc | k3\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | G 2.39722, C 1.80669 | `data/943_cell_gc_gauntlet.json` (full precision) | 1e-3 |
| GanCal#5 (as in the #943 run) | gc | k3\|LGanymede>Ganymede/1l\|LGanymede>Callisto/1h\|LCallisto>Ganymede/0s | G 3.23830, C 3.33953 | same | 1e-3 |
| GanCal#5 (published) | gc | same | G 3.24, C 3.34 | R-S 2009 Table 3 | 0.05 |
| ev-A | ev | k2\|LE>V/0s\|RV/1:1\|LV>V/1h\|LV>E/0s | E 4.89283, V 10.36400 | `data/942_cell_ev_gauntlet.json` | 1e-3 |
| ev-C | ev | k2\|LE>V/0s\|LV>V/1h\|LV>E/0s | E 9.07492, V 13.16643 | same | 1e-3 |
| Hollister 1H (orbit I topology) | ev | k2\|RE/1:1\|LE>V/0s\|RV/1:1\|RV/1:1\|LV>E/0s | E 2.99407, V 3.19027 | same | 1e-3 |
| #576 E-C closures, 6 rows | ec576 | k{n}\|LEuropa>Callisto/0s\|LCallisto>Europa/0s and the Callisto-first order | per row, `data/enumerate_576_jupiter_galilean_symmetric_closures.jsonl` | #576 | 1e-3 |
| same 6 rows | ec | same keys | same | #576 | reported, not judged |

Notes on the controls:
- The #942/#943 rows are code-path regressions (expected values are our own earlier output, at 1e-3
  km/s). The only sourced expected value is GanCal#5's published pair. Hollister's published V_inf for
  orbit I are not reproduced from the stated geometry (generator note sec. 4), so 1H is a topology-
  and-regression control.
- The #576 file holds 6 rows that are 3 physical cyclers (n = 3, 4, 5; each listed once from each moon).
  #576 used the registry moon radii and Jupiter GM (`scripts/scan_558_...::residual_at_point`), not
  R-S Table 2: the radii differ by 1.7e-4 (Europa) and 6e-5 (Callisto) relative. So the judged recall is
  in #576's own model (cell `ec576` = `moon_circular`, registry sma and GM), at 1e-3 km/s. In the ec cell
  the closest zero is reported for comparison only. The n = 5 closure lies outside the ec k = 1-4 run
  and is a targeted solve here. The #576 file predates the turn gate (#888); gate status is reported,
  not expected.
- A control that fails to recall stops the route; the numbers go to the lead, nothing is loosened.

### 2.4 Liang et al. 2024 open segments (model check, not a closure)

Liang 2024 Tables 3/5/7 (transcribed in `src/cyclerfinder/search/cge_scaffold.py`, reproduced in their
own model by #222): the C-G-C half (C -> G 31.8973 d, G -> C 18.1697 d for members A/B) and the C-E-C
half (C -> E, E -> C) of each member A, B, C. Each half is flown as an OPEN two-leg segment in our cell
(gc for C-G-C, ec for C-E-C): our moons are rotated so that the start moon is at Liang's angle at the
first flyby and the middle moon at Liang's angle at the middle flyby; the end moon then sits off
Liang's by the period difference of the two models (reported). Lambert legs at Liang's printed times of
flight, revolution counts 0-4, both branches; the solution nearest the printed V_inf is reported
(identification; the residuals are then judged).

Pass bar (fixed before the run): every |V_inf - printed| (start, middle in, middle out, end) within
`cge_scaffold.vinf_print_tolerance_kms` at Liang's epoch of that flyby (the #222 tolerance), and the
middle-flyby |in| - |out| within the same tolerance.

Status of the bar: the bar was coded in the driver before the first and only `liang` run, but that run
(2026-10-07 21:20 AEDT, scratch) came before this note was committed. Its output is kept as is and
reported in sec. 3; the code and the bar are not changed after it.

What the segments show: the C-G-C half (50.067 d) lies at gc k = 4 (50.09 d), so it is evidence that the
k = 4 cell reaches Liang's geometry with the same Lambert machinery. The C-E-C half (about 49.9-50.1 d,
about 11 E-C synodic periods) is far outside ec k = 1-4; it is a model check of the ec constants only,
not evidence of coverage.

### 2.5 Order and stop rules

R11 (gc) first, then R13 (ev), then R12 (ec). Controls of a cell run before its production run is read.
A control that fails to recall, a generator bug, or a cell/k that would take more than 60 min of wall
time (with 2 workers) stops the route and goes to the lead.

### 2.6 AMENDMENT (lead ruling 2026-10-07), written and committed BEFORE the re-screen of secs. 4-5

Two screens join the PASS definition of sec. 2.2 (and of every two-working-body cell from now on). A
gate-passing zero that fails either is not a pass. They are applied in this order; the first that
applies is the recorded status.

1. **Primary screen.** r_min (`leg_extent`, every leg, fixed legs included) must exceed the primary's
   floor. Jovian cells (gc, ge, ec and their control cells): registry Jupiter radius 71,492 km plus the
   registry safe altitude 5,000 km = 76,492 km (`verify/turn_gate.body_constants("Jupiter")`).
   Heliocentric cells: the registry holds no Sun entry, so the IAU 2015 nominal solar radius, 695,700
   km, with no floor. Status on failure: "reject: primary impact".
2. **Unscheduled-pass screen.** Every leg is sampled (3,000 points, Kepler step) and every interior
   local minimum of the distance to each body of the cell is refined (bounded minimisation). The leg
   end points are the scheduled encounters and are not minima of this search.
   - A minimum below the body's radius (the cell's flyby constants: R-S Table 2 for the moons, the
     registry for the planets): "reject: moon impact" ("reject: planet impact" in heliocentric cells).
   - Otherwise a minimum inside the body's Laplace sphere of influence (`sphere_of_influence_km`:
     Ganymede 24,350 km, Callisto 37,681 km): "model-invalid". The leg was propagated as a conic
     while inside that body's SOI. This is distinct from a gate fail.

Scope and limits: only the bodies of the cell are checked (Io and Europa are not in the gc model, for
example). The search is sampled; a pass that does not make a local minimum of the sampled distance
could be missed.

Implementation: `run_973_enumerate.py screen GAUNTLET.json --out F` (also run inside `gauntlet`, whose
records now carry `screen_status`). The re-screen uses the committed shards and gauntlet JSONs; no
re-enumeration. The post-hoc J class and the `passes` flags of sec. 4 are superseded by this
pre-registered definition (they agree with it in substance; the primary floor now includes Jupiter's
5,000 km safe altitude).

Then, as a check on earlier work (not a re-enumeration): the same screens on the #942/#943 candidates
gc-1, gc-2 (`data/943_cell_gc_gauntlet.json`), ev-A, ev-B, ev-C (`data/942_cell_ev_gauntlet.json`).

### 2.7 AMENDMENT (lead ruling 2026-10-07, after 94dda754 and 2fd1c05a; before the 2.7 re-screen)

The ruling was option (c) of the question in sec. 9.2, with (b) registered as the generator task #1027.

Why: the 2.6 screen rejected the Hollister 1H member 2.99/3.19, a member of a published family flown
in the ephemeris (generator note 6.23, D1). It did so on a quantity our own chooser picked freely: the
direction of a full-revolution return. That is a defect of the screen's scope, not evidence against
the orbit. The positive control did its job by exposing it.

The unscheduled-pass screen of 2.6 is split by leg type:
- (i) **Fixed-geometry legs.** These are Lambert legs (direction set by the solve) and half-rev legs
  (discrete directions; HV(3,1,a) is the tilted circle by the 6.34 rule). The screen applies as written:
  impact rejection or model-invalid. A tilted-circle half-rev that meets its body mid-leg is a physical
  rejection in the ideal model.
- (ii) **Free-direction full-rev n:m legs (R legs).** A mid-leg pass inside an SOI is a CONSTRAINT on the
  direction choice, not a rejection. Until the generator re-picks the direction under that constraint
  (#1027: `optimise_block` minimax subject to no unscheduled pass inside any SOI on the leg; the
  expectation is the same or a slightly worse worst ratio), such a candidate is reported as
  "pass, direction-dependent (tilted-circle minimax)". It is listed separately, never dropped silently,
  and never counted as clean.

The ev in-run control is NOT void. Its V_inf recall (sec. 3.1) stands. Under (ii) the screen flags it,
because its Earth 1:1 return at the minimax direction is the tilted circle (e = 3e-16, i = 5.76 deg),
which meets Earth again at half the leg.

Implementation: `screen_candidate` in `scripts/run_973_enumerate.py`; the leg type is read from the leg
key (R = full-rev). The re-screen of secs. 9.1-9.2 is repeated under 2.7. It applies to ev k = 4/5 when
those runs finish.

## 3. Controls and the Liang segments (results, 2026-10-07)

### 3.1 Recall controls (sec. 2.3): ALL RECALLED

Data: `data/973_gc/recall_controls.json`, `data/973_ev/recall_controls.json`,
`data/973_ec/recall_controls.json` (logs beside them). Every hit is an exact zero (max |residual|
< 2e-13 km/s) with Kepler-step encounter miss < 1e-5 km.

| Control | Cell | Zeros of the structure | Recalled zero V_inf (km/s) | Max dV_inf vs expected | Gate at the floor |
|---|---|---|---|---|---|
| gc-1 | gc | 7 | G 2.39722, C 1.80669 | < 1e-5 | pass, worst 0.737 |
| GanCal#5 (run value) | gc | 10 | G 3.23830, C 3.33953 | < 1e-5 | pass, worst 0.940 |
| GanCal#5 (published R-S Table 3, tol 0.05) | gc | 10 | same | 0.0017 | pass |
| ev-A | ev | 13 | E 4.89283, V 10.36400 | < 1e-5 | pass |
| ev-C | ev | 16 | E 9.07492, V 13.16643 | < 1e-5 | pass (2 zeros, one per mirror) |
| Hollister 1H, orbit I topology | ev | 6 | E 2.99407, V 3.19027 | < 1e-5 | pass |
| #576 E-C n = 3, both anchors | ec576 | 6 / 6 | E 6.47945, C 4.05858 | < 1e-5 | FAIL, worst 28.87 |
| #576 E-C n = 4, both anchors | ec576 | 6 / 8 | E 4.02544, C 4.36172 | < 1e-5 | FAIL, worst 8.33 |
| #576 E-C n = 5, both anchors | ec576 | 6 / 6 | E 4.15159, C 4.98224 | < 1e-5 | FAIL, worst 11.43 |
| same six, in the ec cell (R-S constants; reported only) | ec | 6 / 6 / 6 / 8 / 6 / 6 | n=3: E 6.48193, C 4.05888; n=4: E 4.02736, C 4.36188; n=5: E 4.14983, C 4.98188 | 0.0025 / 0.0019 / 0.0018 (the model-constant shift) | FAIL, 28.90 / 8.33 / 11.42 |

- The 6 #576 rows are 3 physical cyclers; the Europa-first and Callisto-first orders give the same
  zero in each case (6 -> 3).
- All three #576 closures FAIL the #888/#937 demanded-turn gate (worst ratios 8-29). #576 judged them
  by bend capacity before the turn gate existed. This is consistent with the #943 note (6.11), which
  found no #576 closure among the gate-passers.
- The n = 3 and n = 4 closures also appear in the ec production run (sec. 5) as zeros of
  `k3|LEuropa>Callisto/0s|LCallisto>Europa/0s` and the k = 4 key, at the ec-cell values above.

### 3.2 Liang 2024 open segments (sec. 2.4): REPRODUCED, 6 of 6

Data: `data/973_liang_open_segments.json` (the single run, 2026-10-07 21:20 AEDT; see 2.4 on its
timing relative to the pre-registration commit).

| Member | Segment | Days | Legs (rev, branch) | Ours: start, middle in/out, end (km/s) | Printed | Max dev. | Tolerance range |
|---|---|---|---|---|---|---|---|
| A | C-G-C (gc) | 50.067 | 1 high, 1 low | 5.6706, 6.9873/6.9878, 5.6676 | 5.6730, 6.9919, 5.6698 | 0.0046 | 0.017-0.058 |
| B | C-G-C (gc) | 50.067 | same | same as A (Liang's first two rows are identical) | same | 0.0046 | same |
| C | C-G-C (gc) | 50.067 | 1 high, 1 low | 7.6394, 10.4857/10.4859, 7.6371 | 7.6433, 10.4922, 7.6409 | 0.0065 | 0.016-0.057 |
| A | C-E-C (ec) | 49.909 | 1 high, 1 low | 5.6723, 4.6822/4.6814, 5.8742 | 5.6698, 4.6685, 5.8721 | 0.0137 | 0.058-0.137 |
| B | C-E-C (ec) | 49.968 | 1 high, 1 low | 5.6720, 4.4978/4.4970, 5.7931 | 5.6698, 4.4853, 5.7914 | 0.0125 | 0.058-0.138 |
| C | C-E-C (ec) | 50.063 | 1 high, 1 low | 7.6275, 11.9807/11.9806, 7.7701 | 7.6409, 12.0213, 7.7838 | 0.0407 | 0.057-0.138 |

- Middle-flyby |in| - |out| is 1e-4 to 8e-4 km/s in all six. The end moon sits 0.039 deg from Liang's
  position (the period difference of the two models over the segment).
- Leg identification is unambiguous: the runner-up Lambert solution is 1.0-3.4 km/s worse.
- Reading: the gc cell's Lambert machinery reaches Liang's C-G-C geometry, at gc k = 4's period
  (50.067 d against 50.09 d), and the ec constants reproduce Liang's C-E-C legs. These are open
  segments of a three-moon cycle, not closures; they are not candidates and not collisions.

## 4. R11 result: gc (Ganymede-Callisto) k = 4 (2026-10-07)

Run: 4,068 structures (all done once, 0 errors; `check` OK), 4,629 exact zeros, 1,439 physical
cyclers, **51 gate-passing** (0 zero-assessment errors). 0 structures reached the n_refine ceiling (no
seed saturation seen). Gauntlet: `data/973_gc/k4_gauntlet.json` (log `k4_gauntlet.log`).

Shard layout as run (the first four 8-way shards ran 8-10 min each under CI load, so the rest was
split 16-way to keep each tool call under 10 min): `s0`-`s3` are shards 0-3 of 8; `s4`, `s5`, `s6`,
`s7`, `s12`, `s13`, `s14`, `s15` are those shards of 16. Together they cover every residue mod 8
exactly once (checked by `run_973_enumerate.py check`). Wall time 41 min with 2 workers under CI load
(load average 17-30); mean 1.07 s per structure, against the 0.21 s pilot at low load.

Gauntlet, all 51:
- DOP853 re-fly: largest arrival miss 0.49 km (gc4-50; all others < 0.03 km), largest V_inf vector
  error 4.8e-6 km/s, gate on the integrated vectors "pass" for all 51. Largest miss / smallest SOI
  2.0e-5.
- Gate at H&M's 1.1-radius floor: 50 of 51 pass.
- Largest demanded turn 54.6 deg; no near-180 demand.
- Literal collisions: none with R-S GanCal#1/#5 (k = 3 rows; no V_inf within 0.3), Campagnola 2019 GCGC
  (none within 0.5 of 3.5/4.5), or Liang 2024 Tables 3/5/7 (no candidate within 0.3 km/s of a member's
  per-moon V_inf). One NEAR: gc4-7 (G 3.823 / C 2.378) is within 0.149 km/s of the #576 n = 1
  Ganymede-Callisto symmetric closure, at another k (1 vs 4) and another structure; recorded, not a
  collision. Catalogue rows on the pair (7: R-S family, GanCal#1/#5, the four Liang rows): no match.
- Literature step: DEFERRED (sec. 1).

What the 51 are. Two checks were added AFTER the run; they are descriptive, they judge nothing, and the
thresholds are mine.

(a) Unscheduled passes (`run_973_enumerate.py passes`, `data/973_gc/k4_unscheduled_passes.json`): every
leg is sampled (3,000 points, Kepler step) and each interior minimum of the distance to either moon is
refined. A minimum inside a moon's Laplace sphere of influence (Ganymede 24,350 km, Callisto 37,681 km)
is an encounter the cycle does not schedule, so the patched-conic cycle is not valid there. Read
against a known state first: on the #943 gc gate-passers (gc-1, gc-2, GanCal#5) it flags none; the
closest is gc-2's Ganymede pass at 33,899 km on its Callisto-Callisto leg (1.4 SOI). On gc k = 4 it
flags **14 of 51**, three of them at 3-32 km from a moon's centre (inside the moon: gc4-12 Ganymede,
gc4-13 and gc4-21 Callisto).

(b) Architecture, by the demanded turns:
- **J, 8: a leg passes inside Jupiter** (r_min < 71,492 km, down to 185 km from the centre). The
  pre-registered pass definition and the #942/#943 gauntlet have no primary-impact screen, so these pass
  as registered; they are not physical trajectories. Papercut filed.
- **P-C, 13: Callisto passes through** (every Callisto turn < 0.01 deg): Ganymede does all the work and
  Callisto is in effect a massless target. This is Russell & Strange 2009's Ganymede -> Callisto
  architecture (their Table 1 lists that system; the GanCal rows). Our R-S digest does not record the
  largest period they searched, so whether k = 4 lies inside their search is not known. 6 of the 13
  have an unscheduled SOI pass.
- **P-G, 5: Ganymede passes through** (Callisto does the work, Ganymede a target): the reverse system,
  Callisto -> Ganymede, which R-S Table 1 does not list. 1 of the 5 has an unscheduled pass (gc4-13,
  inside Callisto).
- **S-C, 5 / S-G, 8: shallow.** Both moons turn, but Callisto (S-C) or Ganymede (S-G) by less than
  1 deg; near the P architecture of the same letter. 2 and 3 of them have an unscheduled SOI pass.
- **W2, 12: both moons turn by at least 1 deg** (two working bodies in the sense of the #943 brief). 2
  have an unscheduled SOI pass: gc4-0 (Ganymede, 9,432 km) and gc4-11 (Ganymede, 21,913 km). **10 W2
  members have neither a Jupiter pass nor an unscheduled SOI pass**: gc4-1, 2, 6, 8, 10, 19, 25, 28, 38,
  41. These ten are the two-working-body candidates of R11 at k = 4.

Three pairs in P (gc4-26/27 P-G, gc4-33/34 and gc4-36/37 P-C) and three in J share r_min and r_max to 1 km (the extremes come from a leg they have
in common) but differ by 0.02 km/s at the pass-through moon; they are distinct zeros and stay as
separate rows. gc4-16 and gc4-17 have the same V_inf and nearly the same extremes; both have an
unscheduled Callisto pass (5,742 and 3,828 km), so neither is a clean representation.

Why k = 4 holds many near-zero-turn members: in this model 4 G-C synodic periods (50.093 d) are within
0.03 d of 3 Callisto periods (50.067 d) and 7 Ganymede periods (50.082 d) (Lynam 2015's 3:4:7 window).
A spacecraft orbit commensurate with that cycle meets both moons again with a small slip, and a small
turn removes it. This is the expected physics of the window, not a solver effect.

Of the ten W2 candidates, eight keep r_min above 0.3 Gm; gc4-38 and gc4-41 dip to 136,000-143,000 km
(1.9-2.0 Jupiter radii, inside Io's orbit). The lowest-V_inf ones (r_min at Ganymede's orbit,
1,070,335 km, i.e. perijove at Ganymede, or 0.94-0.98 Gm):

- gc4-1: G 1.472 / C 1.581 km/s; Ganymede 2 x 36.1 deg (ratio 0.467), Callisto 2 x 42.8 deg (0.677);
  `k4|RGanymede/3:2|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s`.
- gc4-2: G 1.738 / C 1.857; turns 54.6 / 41.1 deg (0.826 / 0.779); two generic returns.
- gc4-8: G 2.419 / C 2.486; turns 13.7 / 10.5 deg (0.306 / 0.293); two generic returns. The largest
  margin of the W2 set.
- gc4-19: G 2.357 / C 4.109; turns 41.8 / 5.0 deg (0.899 / 0.311).

These are "candidate, literature step deferred", NOT novel. Nearest published relatives named by the
#943 prior-art note: the JUICE C-G-C round trip and the Lynam capture windows (one-shot), Liang 2024's
C-G-C half (sec. 3.2: two to four times the V_inf of these four, open, part of a three-moon cycle).

Candidate table, all 51, with both checks (model epoch: both moons at angle 0 at t = 0; period 50.093 d; "Lambert
starts" are the solved start dates of the Lambert legs in key order; flyby count per moon per cycle):

| # | Class | Unscheduled pass inside SOI (moon, km from centre) | Structure (key) | V_inf G / C (km/s) | G: flybys, max turn (deg), max ratio | C: flybys, max turn, max ratio | Worst ratio | r_min / r_max (km) | Lambert starts (d) | Re-fly miss (km) |
|---|---|---|---|---|---|---|---|---|---|---|
| gc4-0 | W2 | G 9,432 | k4\|RGanymede/1:1\|LGanymede>Callisto/1h\|LCallisto>Ganymede/2l | 3.194 / 1.553 | 2, 23.29, 0.767 | 1, 1.25, 0.019 | 0.767 | 759,541 / 1,882,591 | 4.345, 21.408 | 2.5e-06 |
| gc4-1 | W2 | none | k4\|RGanymede/3:2\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 1.472 / 1.581 | 2, 36.13, 0.467 | 2, 42.79, 0.677 | 0.677 | 1,070,335 / 2,243,803 | 10.564, 31.965 | 8.4e-06 |
| gc4-2 | W2 | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 1.738 / 1.857 | 2, 54.58, 0.826 | 2, 41.08, 0.779 | 0.826 | 940,161 / 2,473,504 | 1.044, 11.479, 18.954, 43.662 | 9.4e-06 |
| gc4-3 | P-G | none | k4\|LGanymede>Callisto/2h\|LCallisto>Ganymede/2l | 4.300 / 1.880 | 1, 0.00, 0.000 | 1, 1.28, 0.025 | 0.025 | 795,584 / 1,882,594 | 10.452, 36.615 | 1.0e-04 |
| gc4-4 | S-G | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Ganymede/2h | 4.388 / 1.913 | 2, 0.85, 0.047 | 1, 6.30, 0.124 | 0.124 | 785,058 / 1,915,782 | 2.241, 20.389, 26.406 | 4.6e-05 |
| gc4-5 | P-C | none | k4\|LGanymede>Callisto/1l\|LCallisto>Ganymede/2h | 1.579 / 2.041 | 1, 1.05, 0.014 | 1, 0.00, 0.000 | 0.014 | 1,070,335 / 2,038,405 | 7.474, 24.158 | 7.8e-06 |
| gc4-6 | W2 | none | k4\|RGanymede/3:2\|LGanymede>Callisto/2l\|LCallisto>Ganymede/0s | 5.007 / 2.135 | 2, 3.91, 0.270 | 1, 1.31, 0.030 | 0.270 | 709,202 / 2,056,858 | 1.532, 24.228 | 8.9e-04 |
| gc4-7 | S-G | none | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 3.823 / 2.378 | 2, 0.16, 0.007 | 2, 34.16, 0.894 | 0.894 | 873,051 / 2,422,525 | 0.723, 24.960, 47.860 | 7.7e-06 |
| gc4-8 | W2 | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 2.419 / 2.486 | 2, 13.72, 0.306 | 2, 10.51, 0.293 | 0.306 | 984,807 / 2,264,287 | 2.782, 22.264, 26.635, 48.504 | 1.9e-06 |
| gc4-9 | P-C | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/2h\|LCallisto>Ganymede/0s | 5.004 / 2.511 | 2, 0.14, 0.010 | 1, 0.00, 0.000 | 0.010 | 726,015 / 1,956,299 | 5.078, 22.919, 50.134 | 6.1e-04 |
| gc4-10 | W2 | none | k4\|RGanymede/1:1\|LGanymede>Callisto/2l\|LCallisto>Ganymede/2h | 6.447 / 2.751 | 2, 5.70, 0.620 | 1, 1.40, 0.045 | 0.620 | 464,563 / 1,882,607 | 0.427, 20.778 | 5.3e-04 |
| gc4-11 | W2 | G 21,913 | k4\|RGanymede/3:2\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 1.820 / 2.816 | 2, 51.45, 0.817 | 2, 1.35, 0.045 | 0.817 | 1,070,338 / 2,308,510 | 10.732, 14.218, 35.875 | 9.1e-06 |
| gc4-12 | P-C | G 3 (inside the moon) | k4\|RGanymede/3:2\|LGanymede>Callisto/0s\|LCallisto>Ganymede/1h | 1.855 / 2.922 | 2, 52.46, 0.851 | 1, 0.00, 0.000 | 0.851 | 1,070,334 / 2,327,904 | 10.696, 14.146 | 7.0e-06 |
| gc4-13 | P-G | C 32 (inside the moon) | k4\|LGanymede>Callisto/2l\|RCallisto/1:2\|LCallisto>Ganymede/1h | 6.833 / 2.937 | 1, 0.00, 0.000 | 2, 2.28, 0.082 | 0.082 | 488,606 / 1,882,611 | 0.148, 36.614 | 4.0e-05 |
| gc4-14 | S-G | C 7,912 | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/2l\|LCallisto>Ganymede/1h | 6.963 / 3.006 | 2, 0.51, 0.064 | 1, 6.59, 0.244 | 0.244 | 473,905 / 1,906,050 | 6.120, 25.130, 44.805 | 5.4e-05 |
| gc4-15 | S-C | none | k4\|RGanymede/1:1\|LGanymede>Callisto/2h\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 6.760 / 3.063 | 2, 5.49, 0.652 | 2, 0.61, 0.023 | 0.652 | 438,208 / 1,912,566 | 0.704, 21.711, 37.621 | 4.1e-05 |
| gc4-16 | S-C | C 5,742 | k4\|LGanymede>Callisto/1h\|RCallisto/1:1\|LCallisto>Ganymede/0s | 2.138 / 3.631 | 1, 1.03, 0.020 | 2, 0.02, 0.001 | 0.020 | 1,069,915 / 2,696,073 | 1.755, 48.818 | 2.9e-05 |
| gc4-17 | P-C | C 3,828 | k4\|LGanymede>Callisto/2h\|LCallisto>Ganymede/0s | 2.138 / 3.631 | 1, 1.03, 0.020 | 1, 0.00, 0.000 | 0.020 | 1,070,334 / 2,695,657 | 1.755, 48.818 | 9.6e-05 |
| gc4-18 | P-C | none | k4\|RGanymede/2:1\|LGanymede>Callisto/1h\|LCallisto>Ganymede/0s | 2.253 / 3.889 | 2, 37.84, 0.770 | 1, 0.00, 0.000 | 0.770 | 1,070,334 / 2,872,965 | 3.517, 36.402 | 1.5e-03 |
| gc4-19 | W2 | none | k4\|RGanymede/2:1\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 2.357 / 4.109 | 2, 41.77, 0.899 | 2, 5.01, 0.311 | 0.899 | 969,653 / 3,047,178 | 3.498, 36.486 | 1.5e-04 |
| gc4-20 | P-C | none | k4\|LGanymede>Callisto/0s\|LCallisto>Ganymede/1h | 2.736 / 4.851 | 1, 1.00, 0.026 | 1, 0.00, 0.000 | 0.026 | 1,070,333 / 3,864,538 | 10.866, 13.360 | 1.9e-04 |
| gc4-21 | P-C | C 31 (inside the moon) | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|LCallisto>Ganymede/2l | 9.229 / 5.002 | 2, 0.15, 0.033 | 1, 0.00, 0.000 | 0.033 | 288,525 / 2,087,948 | 7.776, 31.074, 37.671 | 1.2e-04 |
| gc4-22 | P-C | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/1l\|LCallisto>Ganymede/2h | 8.288 / 5.013 | 2, 0.14, 0.024 | 1, 0.00, 0.000 | 0.024 | 414,566 / 2,265,347 | 8.115, 20.043, 31.709 | 1.2e-04 |
| gc4-23 | S-G | C 11,153 | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|RCallisto/1:2\|LCallisto>Ganymede/0s | 9.315 / 5.064 | 2, 0.16, 0.034 | 2, 0.37, 0.033 | 0.034 | 281,595 / 2,092,521 | 7.783, 31.056, 54.343 | 1.3e-04 |
| gc4-24 | S-C | none | k4\|RGanymede/3:2\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 8.664 / 5.109 | 2, 2.28, 0.433 | 2, 0.32, 0.029 | 0.433 | 367,198 / 2,404,135 | 7.115, 15.052, 29.725 | 4.8e-05 |
| gc4-25 | W2 | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 8.152 / 5.228 | 2, 3.46, 0.586 | 2, 7.30, 0.695 | 0.695 | 357,770 / 2,395,646 | 4.216, 20.830, 30.242, 44.898 | 2.1e-05 |
| gc4-26 | P-G | none | k4\|LGanymede>Callisto/0s\|LCallisto>Callisto/2l\|LCallisto>Ganymede/0s | 7.680 / 5.335 | 1, 0.00, 0.000 | 2, 0.21, 0.020 | 0.020 | 527,421 / 2,583,488 | 9.595, 20.607, 50.659 | 4.2e-05 |
| gc4-27 | P-G | none | k4\|LGanymede>Callisto/1l\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1l | 7.704 / 5.335 | 1, 0.00, 0.000 | 2, 0.21, 0.020 | 0.020 | 527,421 / 2,583,488 | 3.318, 19.336, 39.377 | 2.3e-05 |
| gc4-28 | W2 | none | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 8.857 / 5.551 | 2, 4.44, 0.880 | 2, 7.50, 0.797 | 0.880 | 320,944 / 2,601,280 | 7.082, 30.488, 33.880, 53.783 | 2.2e-05 |
| gc4-29 | S-G | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/1h\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 9.128 / 5.598 | 2, 0.49, 0.103 | 2, 1.33, 0.143 | 0.143 | 328,024 / 2,375,598 | 10.214, 22.349, 39.024, 53.648 | 3.1e-05 |
| gc4-30 | P-G | none | k4\|LGanymede>Callisto/1h\|LCallisto>Callisto/1l\|LCallisto>Ganymede/1h | 9.660 / 5.973 | 1, 0.00, 0.000 | 2, 0.29, 0.035 | 0.035 | 297,589 / 2,380,474 | 10.122, 26.965, 41.587 | 1.6e-05 |
| gc4-31 | S-G | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 8.924 / 6.077 | 2, 0.70, 0.141 | 2, 0.92, 0.116 | 0.141 | 397,621 / 2,756,713 | 6.771, 30.798, 34.050, 53.612 | 1.1e-05 |
| gc4-32 | S-C | none | k4\|RGanymede/2:1\|LGanymede>Callisto/1l\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 9.477 / 6.320 | 2, 2.82, 0.636 | 2, 0.23, 0.031 | 0.636 | 354,194 / 3,005,379 | 1.749, 16.781, 36.187 | 9.2e-06 |
| gc4-33 | P-C | C 25,844 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Ganymede/0s | 8.544 / 6.509 | 2, 0.13, 0.024 | 1, 0.00, 0.000 | 0.024 | 508,404 / 3,264,088 | 12.437, 44.079, 47.246 | 6.8e-05 |
| gc4-34 | P-C | C 25,809 | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|LCallisto>Ganymede/1l | 8.544 / 6.527 | 2, 0.13, 0.024 | 1, 0.00, 0.000 | 0.024 | 508,404 / 3,264,088 | 0.117, 18.567, 32.155 | 7.7e-05 |
| gc4-35 | S-G | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 8.615 / 6.570 | 2, 0.13, 0.025 | 2, 0.14, 0.021 | 0.025 | 500,925 / 3,275,187 | 6.416, 24.891, 42.878 | 3.2e-05 |
| gc4-36 | P-C | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/1h\|LCallisto>Ganymede/0s | 10.467 / 7.065 | 2, 0.14, 0.038 | 1, 0.00, 0.000 | 0.038 | 281,779 / 2,835,570 | 7.632, 31.111, 54.803 | 1.8e-04 |
| gc4-37 | P-C | none | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Ganymede/1l | 10.467 / 7.086 | 2, 0.14, 0.038 | 1, 0.00, 0.000 | 0.038 | 281,779 / 2,835,570 | 2.967, 29.582, 39.248 | 8.8e-05 |
| gc4-38 | W2 | none | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1l | 11.460 / 7.217 | 2, 2.75, 0.897 | 2, 4.89, 0.849 | 0.897 | 143,171 / 2,703,105 | 0.961, 16.537, 19.356, 38.175 | 5.9e-05 |
| gc4-39 | P-C | C 21,720 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/2l\|LCallisto>Ganymede/1h | 12.226 / 7.252 | 2, 0.18, 0.068 | 1, 0.00, 0.000 | 0.068 | 96,759 / 2,285,287 | 0.527, 15.873, 35.320 | 9.9e-05 |
| gc4-40 | S-G | C 36,932 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/1h\|RCallisto/1:2\|LCallisto>Ganymede/0s | 12.285 / 7.291 | 2, 0.18, 0.069 | 2, 0.24, 0.042 | 0.069 | 93,608 / 2,286,215 | 9.165, 24.494, 56.461 | 1.0e-04 |
| gc4-41 | W2 | none | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 11.828 / 7.513 | 2, 1.49, 0.516 | 2, 2.32, 0.435 | 0.516 | 135,971 / 2,662,003 | 8.685, 28.884, 36.396, 51.267 | 2.5e-04 |
| gc4-42 | S-C | G 23,690 | k4\|RGanymede/2:1\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 13.421 / 8.882 | 2, 2.12, 0.939 | 2, 0.28, 0.071 | 0.939 | 78,948 / 3,298,961 | 5.747, 15.271, 30.643 | 2.8e-05 |
| gc4-43 | J | none | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/1h\|LCallisto>Ganymede/0s | 15.080 / 9.633 | 2, 0.21, 0.116 | 1, 0.00, 0.000 | 0.116 | 10,157 / 2,679,136 | 10.288, 39.261, 58.015 | 8.1e-03 |
| gc4-44 | J | none | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/2l\|LCallisto>Ganymede/0s | 15.080 / 9.661 | 2, 0.21, 0.116 | 1, 0.00, 0.000 | 0.116 | 10,157 / 2,679,136 | 7.206, 28.326, 49.580 | 2.4e-02 |
| gc4-45 | J | none | k4\|LGanymede>Callisto/1h\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 14.510 / 9.666 | 1, 0.00, 0.000 | 2, 0.29, 0.087 | 0.087 | 36,965 / 3,075,878 | 3.612, 26.416, 42.175 | 5.4e-04 |
| gc4-46 | J | none | k4\|LGanymede>Callisto/1l\|LCallisto>Callisto/2h\|LCallisto>Ganymede/0s | 14.536 / 9.666 | 1, 0.00, 0.000 | 2, 0.29, 0.087 | 0.087 | 36,965 / 3,075,878 | 9.839, 23.290, 57.623 | 1.3e-03 |
| gc4-47 | J | none | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|LCallisto>Ganymede/0s | 14.637 / 9.726 | 2, 0.18, 0.094 | 1, 0.00, 0.000 | 0.094 | 33,310 / 3,082,660 | 3.761, 40.199, 42.500 | 8.2e-05 |
| gc4-48 | J | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/1h\|LCallisto>Ganymede/1l | 14.637 / 9.750 | 2, 0.18, 0.094 | 1, 0.00, 0.000 | 0.094 | 33,310 / 3,082,660 | 2.567, 16.223, 39.067 | 1.3e-04 |
| gc4-49 | J | none | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 14.744 / 9.819 | 2, 0.14, 0.073 | 2, 0.07, 0.022 | 0.073 | 29,981 / 3,085,784 | 5.679, 19.368, 29.648, 45.491 | 1.8e-03 |
| gc4-50 | J | none | k4\|LGanymede>Callisto/2l\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1l | 15.993 / 10.251 | 1, 0.00, 0.000 | 2, 0.46, 0.156 | 0.156 | 185 / 2,692,527 | 3.232, 24.253, 41.025 | 4.9e-01 |


## 5. R12 result: ec (Europa-Callisto, both massive) k = 1-4 (2026-10-07): NO GATE-PASSING CYCLER

Run: 2,014 structures (all done once, 0 errors; `check` OK per k), 462 exact zeros, 215 physical
cyclers (merged within each k: 7 / 31 / 55 / 122), **0 gate-passing**. 0 structures at the n_refine
ceiling. Data: `data/973_ec/k1/s0`, `k2/s0`, `k3/s0`, `k4/s0`, `k4/s1` (k = 4 as 2 shards of 2);
gauntlet summary `data/973_ec/k1-4_gauntlet.json`. Wall time 7 min (mean 0.29-0.34 s per structure).

| k | Period (d) | Structures | Zeros | Gate fail | No minimax directions | Lowest worst ratio (structure) |
|---|---|---|---|---|---|---|
| 1 | 4.51 | 37 | 14 | 7 | 7 | 110.7 (k1\|HEuropa/1,0,a\|LEuropa>Callisto/0s\|LCallisto>Europa/0s) |
| 2 | 9.02 | 177 | 48 | 38 | 10 | 8.10 (k2\|HEuropa/1,0,a\|LEuropa>Callisto/0s\|LCallisto>Europa/0s) |
| 3 | 13.54 | 470 | 112 | 77 | 35 | 4.76 (k3\|REuropa/2:1\|LEuropa>Callisto/0s\|LCallisto>Europa/0s) |
| 4 | 18.05 | 1,330 | 288 | 209 | 79 | 1.36 (k4\|LEuropa>Europa/1l\|LEuropa>Callisto/0s\|LCallisto>Europa/0s; E 4.37, C 3.65 km/s, largest turn 14.3 deg) |

- In-run recall: the #576 n = 3 and n = 4 closures are zeros of the plain E->C->E structures at k = 3
  and k = 4, with the ec-cell values of sec. 3.1. Both fail the gate (28.9, 8.32).
- Reading: in this cell the E-C period at k <= 4 (at most 18 d, about one Callisto period) leaves the
  transfers fast and the V_inf high relative to what Europa and Callisto can turn; the best zero needs
  1.36 times the available turn. The trend with k (110 -> 8.1 -> 4.8 -> 1.36) says the first gate
  passers, if any, lie at k >= 5. That is outside this pre-registration; Liang's C-E-C half sits at
  about 11 synodic periods.
- Proposed registry stamp (for the lead; `data/empty_regions.jsonl` not edited): region
  `jupiter-europa-callisto-two-working-body-rs2009-ideal-k1-4-973`, centre Jupiter, template
  `[A-block, A->B, B-block, B->A]` one visit, k 1-4, returns 1,1, transfer revs 0,1,2, generic revs 1,2,
  catalogue as #943 gc, seeds 36/12/40/0.03, method `two_working_body.correct_dates` + minimax turn
  gate (ballistic, coplanar, patched-conic, circular), points 2,014 structures, zeros 462, physical 215,
  gate-passing 0, errors 0, prune gates as #943 gc. Empty is conditional on this catalogue and these
  seeds.

## 6. Not run by #973: gc k = 5, 6 and ev k = 4, 5 (lead ruling 2026-10-07: the lead launches them)

Each is over 60 min of wall time with 2 workers while CI runs. Controls for both cells are already
recalled (sec. 3.1), so the production runs can start without them.

### 6.1 Launch commands (from the repo root; foreground or as the lead prefers)

One block per cell/k; set the four variables, then run the three commands. Each shard writes its own
directory `data/973_<cell>/k<k>/s<i>/` and log `s<i>.log`.

```sh
# gc k = 5:  CELL=gc K=5 N=16 PILOT=0.5  FLAGS="--max-returns 1,1 --transfer-revs 0,1,2 --generic-revs 1,2"
# gc k = 6:  CELL=gc K=6 N=24 PILOT=0.78 FLAGS="--max-returns 1,1 --transfer-revs 0,1,2 --generic-revs 1,2"
# ev k = 4:  CELL=ev K=4 N=16 PILOT=0.8  FLAGS="--max-returns 2,2 --transfer-revs 0 --generic-revs 1 --resonant-only E"
# ev k = 5:  CELL=ev K=5 N=24 PILOT=0.91 FLAGS="--max-returns 2,2 --transfer-revs 0 --generic-revs 1 --resonant-only E"
export CELL=gc K=5 N=16 PILOT=0.5 FLAGS="--max-returns 1,1 --transfer-revs 0,1,2 --generic-revs 1,2"
mkdir -p data/973_$CELL/k$K
seq 0 $((N-1)) | xargs -n1 -P2 sh -c 'uv run python scripts/run_973_enumerate.py enumerate --cell $CELL --k $K --out data/973_$CELL/k$K/s$0 --shard $0/$N --timing-pilot-s $PILOT --n-phase 36 --n-split 12 --n-refine 40 $FLAGS > data/973_$CELL/k$K/s$0.log 2>&1; echo "shard $0 exit $?: $(tail -1 data/973_$CELL/k$K/s$0.log)"'
uv run python scripts/run_973_enumerate.py check data/973_$CELL/k$K/s*/
uv run python scripts/run_973_enumerate.py gauntlet --cell $CELL data/973_$CELL/k$K/s*/ --out data/973_$CELL/k${K}_gauntlet.json | tee data/973_$CELL/k${K}_gauntlet.log
```

Since amendment 2.6, `gauntlet` also applies both screens; each record carries `screen_status`.
For ev, the tilted-circle question in sec. 9.2 decides how the R-leg rejections are read.

- `xargs -n1` with `sh -c '... $0 ...'` is used because macOS `xargs -I` refuses commands over 255
  bytes. `-P2` = 2 workers; raise it when the machine is idle (shards are independent).
- A subset of shards can be run by replacing `seq 0 $((N-1))` with a list, e.g. `printf "0\n1\n"`.
- The command line was checked with `--count-only` (gc k = 5: 287 structures in shard 0 of 16).

### 6.2 Expected runtime (structure counts from sec. 2.1; per-structure times measured)

| Cell/k | Structures | Pilot s/structure (load 3-4) | Serial, idle | 2 workers, idle | 2 workers, CI load (x5, as gc k4) |
|---|---|---|---|---|---|
| gc k5 | 4,589 | about 0.5 (between the k4 and k6 pilots; not measured) | 38 min | 19 min | 1.6 h |
| gc k6 | 4,844 | 0.78 | 63 min | 32 min | 2.6 h |
| ev k4 | 3,566 | about 0.8 (not measured; ev k5 pilot) | 48 min | 24 min | 2.0 h |
| ev k5 | 4,070 | 0.91 | 62 min | 31 min | 2.6 h |

The gc k4 pilot (0.21 s, 24 structures, low load) underestimated the full run (1.07 s mean, CI load 17-30)
by 5x. Each progress line in a shard log carries the running totals and an ETA; a shard that stops
writing for more than about 5 min while its process is gone has died.

### 6.3 Resume and checkpoint behaviour

- Each shard appends one line per finished structure to `structures.jsonl` (with its zero count and
  time) and one line per zero to `zeros.jsonl`. Re-running the same shard command with the same
  `--out` skips every key already in `structures.jsonl` (errored structures are retried).
- Hazard: the zero lines of a structure are written BEFORE its structure line. A shard killed between
  the two writes solves that structure again on resume and writes its zeros twice. `check` counts
  duplicated zero lines and fails on them; the physical-cycler merge in the gauntlet is not affected,
  only the raw zero count. If `check` reports duplicates after a resume, delete the duplicate zero
  lines of the one structure named in the last `structures.jsonl` line before the kill (or re-run that
  shard into a fresh directory).
- Never mix `--shard i/N` layouts in one cell/k unless the union covers every residue once (sec. 4 did
  8-way plus 16-way; `check` verifies coverage from the structure keys, not from the layout).

### 6.4 Validating a finished cell/k

`uv run python scripts/run_973_enumerate.py check data/973_<cell>/k<k>/s*/` regenerates the full
structure list from the shards' `settings.json` and fails (exit 1) if: settings differ between shards;
any structure is missing, extra, done twice or errored; or a zero line is duplicated. It prints the
number of structures whose zero count reached n_refine (seed saturation; report it). Read a known
state first: on gc k4 it prints `expected 4068 structures, done 4068, zero lines 4629, ... ceiling
(40): 0` and `OK`.

### 6.5 Analysis of those runs (for the fresh agent)

Same as sec. 4: `gauntlet` on the shard dirs; then the classes J (r_min below the primary radius:
Jupiter 71,492 km; for ev the Sun is never reached by these conics, but check r_min against the solar
radius 695,700 km anyway), P (a moon with every turn < 0.01 deg), S (one moon < 1 deg), W2 (both >=
1 deg); the candidate table as in sec. 4 (`flyby_table`, `x_days`, `r_min_km`, `r_max_km`,
`cross_check`); `passes` on the gauntlet JSON (unscheduled passes inside a moon's or planet's SOI);
seed-saturation count from `check`. For ev k = 4/5, the H&M check of sec. 2.2 is done
(below). The literature step and rung (d) stay deferred until #972 lands.

ev, H&M 1970 sub-period check (done, `data/973_ev/hm1970_subperiod_check.json`): every Table 3 orbit is
5 blocks of [E pair, E->V, two Venus returns, V->E], each block 2 synodic periods (3.2 yr); orbits 1-3
repeat block by block (k = 2), orbits 4-15 have no sub-period (k = 10 only). Every H&M orbit has one
E->V transfer per 2 synodic periods; a #973 ev structure at k = 4 or 5 has one per 4 or 5. So no k = 4
or k = 5 ev structure can equal an H&M orbit or a sub-period of one. The catalogue V_inf comparison
(the 15 H&M rows) still runs in the gauntlet.

## 7. Verified and assumed

Verified (by a command whose output is in the data files):
- Every control of sec. 2.3 recalls, at 1e-3 km/s (run values, #576 in its own model) and 0.05 km/s
  (GanCal#5 published).
- Liang 2024 C-G-C and C-E-C halves of members A, B, C reproduce within the #222 tolerance in our cells.
- gc k4 and ec k1-4: every structure of the cell solved exactly once, 0 errors, 0 duplicated zero lines,
  0 seed-saturated structures.
- The 51 gc k4 gate-passers: DOP853 re-fly miss <= 0.49 km, gate pass on integrated vectors, SOI
  fraction <= 2e-5. 8 pass inside Jupiter; 14 have an unscheduled pass inside a moon's SOI (the check
  flags none on the three #943 gc gate-passers).

Assumed or not checked:
- That 36/12/40 seeds are dense enough at k = 4 (free time about 1.3x that of k = 3). No structure hit
  the n_refine ceiling, which is necessary for completeness, not sufficient. A doubled-seed sample was
  not run.
- The ec cell's constants are R-S 2009 Table 2's; R-S treated no Europa-Callisto pair, so no published
  E-C row tests the cell beyond Liang's open segments.
- Novelty: not assessed (literature step deferred).
- The J/P/S/W2 classes are descriptive; their thresholds (R_J, 0.01 deg, 1 deg) were chosen after the
  run and judge nothing. The unscheduled-pass check is a sampled search (3,000 points per leg,
  then bounded refinement); a pass between samples that does not make a local minimum of the sampled
  distance could be missed.

## 8. Papercuts

- `docs/papercuts/2026-10-07-twobody-ext-opus-no-primary-impact-screen.md`: the #942/#943 pass
  definition and gauntlet do not reject conics that pass inside the central body, nor legs that pass
  through a moon's SOI between scheduled encounters.
- `docs/papercuts/2026-10-07-twobody-ext-opus-pilot-under-load.md`: the 24-structure pilot ran at low
  load and underestimated the CI-loaded run by 5x; the 10-min tool-call limit then forced a mixed
  8/16-way shard layout.

## 9. Re-screen under amendment 2.6 (2026-10-07, after commit 94dda754)

### 9.1 gc k = 4 and ec

The re-screen runs on the committed gauntlet JSON (`run_973_enumerate.py screen`); there is no
re-enumeration. Data: `data/973_gc/k4_screen.json` and its log.

| Status | Count | By the sec. 4 class |
|---|---|---|
| pass | 29 | W2 10, P-C 7, P-G 4, S-C 3, S-G 5 |
| model-invalid (unscheduled pass inside a moon's SOI) | 11 | W2 2 (gc4-0, gc4-11), P-C 4, S-C 2, S-G 3 |
| reject: moon impact | 3 | P-C 2 (gc4-12, gc4-21), P-G 1 (gc4-13) |
| reject: primary impact (r_min <= 76,492 km) | 8 | J 8 |

- The primary floor now includes Jupiter's 5,000 km safe altitude, but no candidate moves because of it.
  The three lowest non-J r_min values (78,948, 93,608 and 96,759 km) are all model-invalid on other
  grounds.
- Every gc rejection comes from a Lambert leg, so it does not depend on a free-direction choice (see
  9.2).
- ec k = 1-4 had no gate-passer, and the screens can only remove, so the result stays empty. Registry
  stamp appended: `jupiter-europa-callisto-two-working-body-rs2009-ideal-k1-4-973` in
  `data/empty_regions.jsonl`, with both screens listed among its prune gates.

The 10 clean W2 members all pass the screens. They are the R11 candidate set at k = 4 (sec. 10). The
other 19 screen-passers are pass-through (P) or shallow (S) members, close to the one-working-body
architecture, and are not carried forward as two-working-body candidates.

### 9.2 The #942/#943 candidates, and a finding about tilted-circle fixed legs

Data: `data/973_screen_943_gc.json`, `data/973_screen_942_ev.json`.

- gc-1, gc-2 and GanCal#5 (#943 gc): all pass. gc-2's closest unscheduled pass is Ganymede at
  33,899 km on its Callisto-Callisto leg (SOI 24,350 km).
- ev-A, ev-B and ev-C pass.
  - ev-A: closest Earth 152,190,166 km, Venus 200,339,040 km.
  - ev-B: closest Earth 1,172,392 km, about 1.27 times Earth's Laplace SOI.
  - ev-C: closest Earth 47,110,597 km, Venus 5,723,283 km.
- The full #942 ev gate-passing set has 31 cyclers. The first screen (2.6) gave 22 pass and 9 "reject:
  planet impact". In all 9 a fixed leg's flyby direction is a pure tilted circle (e about 3e-16, a equal
  to the planet's orbit radius, i about 5-6 deg), which meets the planet again mid-leg:
  - on an Earth 1:1 full-rev at its minimax direction, at half the leg, 0.3 km from Earth's centre;
  - on a Venus (3,1,apo) half-rev (the tilted circle by the 6.34 rule), at 1/3 of the leg.
- Under amendment 2.7 (re-screen 2026-10-07; data replaced in place) the result is **22 pass, 4 pass
  flagged "direction-dependent (tilted-circle minimax)", and 5 reject: planet impact**:

| # (gauntlet index) | k | Structure | V_inf E / V (km/s) | 2.7 status | Leg that re-meets the planet |
|---|---|---|---|---|---|
| 0 | 2 | k2\|RE/1:1\|LE>V/0s\|RV/1:1\|RV/1:1\|LV>E/0s | 2.99 / 3.19 | direction-dependent: Hollister 1H, the in-run control | RE/1:1 (Earth, half-leg) |
| 2 | 2 | k2\|RE/1:1\|LE>V/0s\|HV/1,0,a\|LV>V/1h\|LV>E/0s | 3.51 / 4.41 | direction-dependent (the Menning-variation skeleton of 6.12) | RE/1:1 |
| 3 | 2 | k2\|RE/1:1\|LE>V/0s\|LV>V/1l\|RV/1:1\|LV>E/0s | 4.05 / 7.07 | direction-dependent: Hollister 2H member | RE/1:1 |
| 13 | 3 | k3\|RE/1:1\|LE>V/0s\|LV>V/1h\|RV/3:2\|LV>E/0s | 4.30 / 5.82 | direction-dependent | RE/1:1 |
| 1 | 2 | k2\|RE/1:1\|LE>V/0s\|HV/1,0,a\|HV/3,1,a\|LV>E/0s | 2.99 / 3.19 | reject: planet impact (the 1H variant with a half-rev pair at Venus) | HV/3,1,a (Venus, 1/3 of the leg); its RE/1:1 too |
| 12 | 3 | k3\|RE/1:1\|LE>V/0s\|RV/3:2\|HV/3,1,a\|LV>E/0s | 3.34 / 4.04 | reject: planet impact | HV/3,1,a; its RE/1:1 too |
| 16 | 3 | k3\|RE/1:1\|RE/2:3\|LE>V/0s\|HV/1,0,a\|HV/3,1,a\|LV>E/0s | 5.58 / 3.77 | reject: planet impact | HV/3,1,a |
| 17 | 3 | k3\|RE/1:1\|RE/1:1\|LE>V/0s\|RV/1:1\|HV/3,1,a\|LV>E/0s | 5.92 / 3.58 | reject: planet impact | HV/3,1,a |
| 18 | 3 | k3\|RE/2:3\|LE>V/0s\|RV/1:1\|HV/3,1,a\|LV>E/0s | 5.92 / 3.58 | reject: planet impact | HV/3,1,a |

- The ev in-run control is not void: its V_inf recall stands (sec. 3.1). The screen flags it under 2.7
  (ii) because its Earth 1:1 return at the minimax direction is the tilted circle. The other 1H member
  (6.14/5.28) and the other 2H member (5.60/6.02) pass outright. ev-A, ev-B and ev-C pass outright.
- The 4 flagged cyclers wait for #1027 (direction re-picked under the no-unscheduled-pass constraint).
  The 5 rejections are physical in the ideal model.
- idx 1 is a member of a published class: it is the Hollister 1H variant with the Venus half-rev pair.
  Its rejection holds for the IDEAL circular model only. There, the HV(3,1,a) leg is a tilted circle
  and meets Venus again at 1/3 of the leg. The published orbit exists in the real ephemeris, where
  that leg is not a tilted circle. Whenever a published or control member is rejected on fixed-leg
  geometry, the real-ephemeris rung decides, not the ideal-model screen.
- gc k = 4 is unchanged under 2.7 (every gc flag is on a Lambert leg), as are the #943 gc candidates.

### 9.3 Addendum text for the #942/#943 generator note (for the lead to relay; that note is mid-edit)

> Addendum (2026-10-07, from #973): the pass definition of sec. 6.1 had two gaps. It never compared
> r_min with the primary's radius, and it never looked for a leg passing inside a body's SOI between
> the scheduled encounters (`encounter_self_consistency` and the DOP853 re-fly check only the leg end
> points). #973 amendment 2.6 adds both: r_min must exceed the primary radius plus the registry floor
> (Jupiter 71,492 + 5,000 km; Sun 695,700 km), and an interior distance minimum below a body's radius
> is an impact rejection, while one inside its Laplace SOI makes the structure model-invalid. Applied
> to this note's candidates (`data/973_screen_943_gc.json`, `data/973_screen_942_ev.json`), gc-1,
> gc-2, GanCal#5, ev-A, ev-B and ev-C all pass (ev-B's closest unscheduled Earth pass is 1.27 Earth
> SOI). Of the 31 ev gate-passers, 9 have a fixed leg whose direction is the tilted circle (e = 0,
> a = the planet's orbit radius), which meets the planet again mid-leg. Under the lead's ruling
> (#973 amendment 2.7) the 5 with a tilted-circle half-rev (Venus 3,1,apo) are physical rejections.
> The 4 whose only re-encounter is on a free-direction Earth 1:1 full-rev, including the Hollister 1H
> member 2.99/3.19 (the in-run control, whose recall stands), are "pass, direction-dependent" until
> #1027 re-picks those directions. The vm, vm2, em and ge cells have not been re-screened.

## 10. Final candidate table: R11, Ganymede-Callisto, k = 4 (the clean two-working-body members)

Status: "candidate, literature step deferred". The literature step and rung (d) are task #1025, which
waits for #972 and the chain tool. Not novel, and no catalogue writes.

Model: circular coplanar, R-S 2009 Table 2 constants, both moons at angle 0 at t = 0. Period 50.0929 d
(4 Ganymede-Callisto synodic periods). "Lambert starts" are the solved start dates of the Lambert
legs, in key order. Turns are the demanded turns at the minimax directions; ratio = demanded /
available at the registry floor (Ganymede 100 km, Callisto 200 km). Every member is an exact zero
(residual < 1e-8 km/s), passes the gate at every flyby (also on the DOP853-integrated vectors) and
passes both amendment-2.6 screens.

| # | Structure (key) | Lambert starts (d) | V_inf G / C (km/s) | Ganymede flybys: turn deg (ratio, required alt km) | Callisto flybys: turn deg (ratio, required alt km) | Worst ratio | r_min / r_max (km) | Closest unscheduled G / C (km) | Re-fly miss (km) |
|---|---|---|---|---|---|---|---|---|---|
| gc4-1 | k4\|RGanymede/3:2\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 10.5642, 31.9650 | 1.4719 / 1.5814 | 36.13 (0.467, 7,521); 36.13 (0.467, 7,521) | 42.79 (0.677, 2,592); 42.79 (0.677, 2,592) | 0.677 | 1,070,335 / 2,243,803 | 937,659 / 2,861,790 | 8.4e-06 |
| gc4-2 | k4\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 1.0438, 11.4794, 18.9537, 43.6625 | 1.7380 / 1.8572 | 54.58 (0.826, 1,232); 54.58 (0.826, 1,232) | 41.08 (0.779, 1,443); 41.08 (0.779, 1,443) | 0.826 | 940,161 / 2,473,504 | 234,539 / 537,850 | 9.4e-06 |
| gc4-6 | k4\|RGanymede/3:2\|LGanymede>Callisto/2l\|LCallisto>Ganymede/0s | 1.5319, 24.2277 | 5.0073 / 2.1348 | 3.91 (0.270, 8,532); 3.91 (0.270, 8,532) | 1.31 (0.030, 133,459) | 0.270 | 709,202 / 2,056,858 | 780,710 / 1,102,074 | 8.9e-04 |
| gc4-8 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 2.7824, 22.2641, 26.6353, 48.5041 | 2.4191 / 2.4855 | 13.72 (0.306, 9,820); 13.72 (0.306, 9,820) | 10.51 (0.293, 9,118); 10.51 (0.293, 9,118) | 0.306 | 984,807 / 2,264,287 | 106,557 / 705,684 | 1.9e-06 |
| gc4-10 | k4\|RGanymede/1:1\|LGanymede>Callisto/2l\|LCallisto>Ganymede/2h | 0.4268, 20.7783 | 6.4468 / 2.7512 | 5.70 (0.620, 1,917); 5.70 (0.620, 1,917) | 1.40 (0.045, 74,348) | 0.620 | 464,563 / 1,882,607 | 512,601 / 129,111 | 5.3e-04 |
| gc4-19 | k4\|RGanymede/2:1\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | 3.4984, 36.4856 | 2.3565 / 4.1090 | 41.77 (0.899, 580); 41.77 (0.899, 580) | 5.01 (0.311, 6,905); 5.01 (0.311, 6,905) | 0.899 | 969,653 / 3,047,178 | 1,477,091 / 1,465,758 | 1.5e-04 |
| gc4-25 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 4.2162, 20.8302, 30.2416, 44.8978 | 8.1525 / 5.2279 | 3.46 (0.586, 2,139); 3.46 (0.586, 2,139) | 7.30 (0.695, 1,458); 7.30 (0.695, 1,458) | 0.695 | 357,770 / 2,395,646 | 693,108 / 1,386,805 | 2.1e-05 |
| gc4-28 | k4\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s | 7.0815, 30.4882, 33.8800, 53.7826 | 8.8570 / 5.5509 | 4.44 (0.880, 491); 4.44 (0.880, 491) | 7.50 (0.797, 924); 7.50 (0.797, 924) | 0.880 | 320,944 / 2,601,280 | 1,009,509 / 185,038 | 2.2e-05 |
| gc4-38 | k4\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1l | 0.9608, 16.5373, 19.3562, 38.1754 | 11.4600 / 7.2168 | 2.75 (0.897, 424); 2.75 (0.897, 424) | 4.89 (0.849, 689); 4.89 (0.849, 689) | 0.897 | 143,171 / 2,703,105 | 196,345 / 285,872 | 5.9e-05 |
| gc4-41 | k4\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 8.6852, 28.8845, 36.3961, 51.2666 | 11.8283 / 7.5126 | 1.49 (0.516, 2,732); 1.49 (0.516, 2,732) | 2.32 (0.435, 3,744); 2.32 (0.435, 3,744) | 0.516 | 135,971 / 2,662,003 | 151,840 / 1,401,694 | 2.5e-04 |

Notes on the set:
- gc4-6 and gc4-10 meet Callisto once per cycle with a turn of 1.3-1.4 deg (ratio 0.03-0.05). They are
  W2 by the 1-deg line, but Callisto barely works.
- gc4-38 and gc4-41 reach 136,000-143,000 km from Jupiter (1.9-2.0 Jupiter radii, inside Io's orbit;
  Io is not in this model).
- The lowest V_inf is gc4-1 (G 1.47 / C 1.58 km/s, perijove at Ganymede's orbit). The largest margin
  is gc4-8 (ratios 0.31 / 0.29).

## 11. R11 result: gc (Ganymede-Callisto) k = 5 (lead launch; analysed 2026-10-08)

Run: launched by the lead with the sec. 6.1 commands (16 shards of 16). Results:
- 4,589 structures, all done once, 0 errors (`check` OK).
- 9,699 exact zeros, 3,106 physical cyclers, **8 gate-passing**, 0 zero-assessment errors.
- 0 structures at the n_refine ceiling.

Gauntlet and screens: `data/973_gc/k5_gauntlet.json` and its log. The gauntlet applies the 2.6/2.7
screens; the sec. 4 classes were added by the analysis.
- DOP853 re-fly: largest miss 1.1e-4 km, gate "pass" on the integrated vectors for all 8.
- Screens: 7 pass, 1 model-invalid. gc5-2 makes an unscheduled Ganymede pass at 3,562 km on its 2-rev
  Ganymede-Ganymede leg, a Lambert leg, so the rejection is fixed-geometry. No primary impact, no
  moon impact, no direction-dependent flag.
- Classes: W2 7 (6 pass, 1 model-invalid) and P-C 1 (gc5-1: Callisto passes through, the R-S Ganymede
  -> Callisto architecture). No J, no S.
- Literal collisions: none. Checked against R-S GanCal#1/#5, Campagnola 2019 GCGC, Liang 2024 Tables
  3/5/7 per-moon V_inf, the #576 G-C closures and the 7 catalogue rows on the pair. None is within
  the NEAR bands.
- Literature step: deferred (#1025).

**Clean two-working-body members at k = 5: 6** (gc5-0, 3, 4, 5, 6, 7).

Why k = 5 differs from k = 4. In this model (R-S 2009 Table 2 periods) k G-C synodic periods are:

| k | Period (d) | in Callisto periods | in Ganymede periods |
|---|---|---|---|
| 3 | 37.570 | 2.251 | 5.251 |
| 4 | 50.093 | 3.002 | 7.002 |
| 5 | 62.616 | 3.752 | 8.752 |
| 6 | 75.139 | 4.502 | 10.502 |

At k = 4 the cycle is within 0.2 % of 3 Callisto and 7 Ganymede periods (the 3:4:7 window). A
spacecraft orbit commensurate with the cycle meets both moons again with a small slip, so many
near-zero-turn members pass the gate (sec. 4: 18 P and 13 S of 51).

At k = 5 the cycle is about three quarters of a period off a whole revolution of each moon
(3.752 and 8.752; the moons' phases repeat with respect to each other, but not with respect to
inertial space). The nearest whole-revolution counts, 4 Callisto (66.76 d) and 9 Ganymede (64.39 d),
miss the cycle by 4.1 d and 1.8 d. No spacecraft orbit commensurate with the cycle returns to both
moons with a small slip, so the moons must turn the orbit by large angles:
- All 7 W2 members demand 10-68 deg at each moon. Their worst ratios are 0.76-0.97, against 0.01-0.94
  at k = 4.
- Only 1 P member survives.
- The gate-passer count drops from 51 to 8, although the zero count roughly doubles (4,629 -> 9,699).

At k = 6 the cycle is within 0.05 % of 4.5 Callisto and 10.5 Ganymede periods: half-period
commensurate. This is a structure the half-rev returns can use; it is noted here as a prediction
for the k = 6 run, not as a result.

The k = 5 set is the lowest-V_inf of the route so far:
- gc5-0: G 1.490 / C 1.255 km/s, with turns of 68 / 64 deg (ratios 0.89 / 0.82).
- gc5-3: G 1.54 / C 1.83, the largest margin (worst ratio 0.761).
All 6 clean members keep r_min at or above 678,771 km (outside Europa's orbit) and r_max at or below
2.93 Gm.

Candidate table, all 8 (period 62.6162 d; model epoch as in sec. 10; the screen status is the
2.6/2.7 status):

| # | Class | Screen | Structure (key) | V_inf G / C (km/s) | G: flybys, max turn (deg), max ratio | C: flybys, max turn, max ratio | Worst ratio | r_min / r_max (km) | Lambert starts (d) | Closest unscheduled G / C (km) | Re-fly miss (km) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gc5-0 | W2 | pass | k5\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1h | 1.490 / 1.255 | 2, 68.21, 0.891 | 2, 64.45, 0.816 | 0.891 | 945,432 / 2,085,249 | 5.5332, 23.3428, 28.6227, 51.2839 | 124,906 / 370,431 | 1.2e-05 |
| gc5-1 | P-C | pass | k5\|LGanymede>Ganymede/2h\|LGanymede>Callisto/1h\|LCallisto>Ganymede/1l | 1.499 / 1.578 | 2, 68.78, 0.903 | 1, 0.00, 0.000 | 0.903 | 1,032,361 / 1,938,679 | 5.3970, 32.5034, 51.5014 | 868,341 / 350,185 | 8.6e-06 |
| gc5-2 | W2 | model-invalid | k5\|LGanymede>Ganymede/2h\|LGanymede>Callisto/1h\|LCallisto>Ganymede/2h | 3.292 / 1.645 | 2, 27.77, 0.958 | 1, 35.19, 0.581 | 0.958 | 704,892 / 1,892,585 | 0.0926, 16.8994, 34.4372 | 3,562 / 766,149 | 1.4e-05 |
| gc5-3 | W2 | pass | k5\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 1.536 / 1.831 | 2, 31.90, 0.428 | 2, 40.81, 0.761 | 0.761 | 1,053,110 / 2,464,863 | 1.3327, 23.7138, 31.4858, 56.1769 | 964,228 / 532,023 | 5.4e-06 |
| gc5-4 | W2 | pass | k5\|RGanymede/3:2\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 1.745 / 1.902 | 2, 55.77, 0.848 | 2, 39.56, 0.773 | 0.848 | 1,058,116 / 2,488,277 | 10.7318, 18.9386, 43.6776 | 983,320 / 557,262 | 6.6e-06 |
| gc5-5 | W2 | pass | k5\|LGanymede>Callisto/0s\|LCallisto>Callisto/2l\|LCallisto>Ganymede/0s | 1.936 / 2.760 | 1, 49.86, 0.848 | 2, 16.82, 0.546 | 0.848 | 1,070,338 / 2,654,042 | 12.5232, 22.5166, 65.1460 | 362,630 / 547,077 | 1.9e-05 |
| gc5-6 | W2 | pass | k5\|RGanymede/2:1\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 2.191 / 3.206 | 2, 45.27, 0.890 | 2, 13.98, 0.573 | 0.890 | 1,070,338 / 2,925,405 | 0.8929, 12.2320, 37.8609 | 660,651 / 445,194 | 1.8e-05 |
| gc5-7 | W2 | pass | k5\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/1l | 4.701 / 4.597 | 2, 15.63, 0.967 | 2, 10.43, 0.788 | 0.967 | 678,771 / 2,926,921 | 3.4184, 21.4057, 24.9696, 45.4855 | 516,757 / 201,809 | 1.1e-04 |

Clean members, per-flyby detail (turn deg, ratio, required altitude km; registry floors Ganymede 100 km,
Callisto 200 km):

| # | Ganymede flybys | Callisto flybys |
|---|---|---|
| gc5-0 | 2 x 68.21 (0.891, 854) | 2 x 64.45 (0.816, 1,583) |
| gc5-3 | 2 x 31.90 (0.428, 8,419) | 2 x 40.81 (0.761, 1,594) |
| gc5-4 | 2 x 55.77 (0.848, 1,062) | 2 x 39.56 (0.773, 1,472) |
| gc5-5 | 1 x 49.86 (0.848, 985) | 2 x 16.82 (0.546, 3,093) |
| gc5-6 | 2 x 45.27 (0.890, 657) | 2 x 13.98 (0.573, 2,636) |
| gc5-7 | 2 x 15.63 (0.967, 208) | 2 x 10.43 (0.788, 988) |

Status: "candidate, literature step deferred" (#1025). Not novel. Raw data (lead rule 2026-10-08, from k = 5 on):
- Committed: each shard's `settings.json` and `structures.jsonl` (`data/973_gc/k5/s*/`), plus the
  gauntlet JSON and log.
- Not committed: `zeros.jsonl`.
- The full shard directories (zeros and logs) are archived outside the repo at
  `~/dev/references/cyclers-runs/973/gc/k5/` (14 MB; `diff -rq` against the working copy was
  clean at archive time).
- The gauntlet can be reproduced either from the archive, or by re-running the committed driver with
  the committed settings (`--shard i/16`, the sec. 6.1 gc k = 5 line).

md5 of each `zeros.jsonl` (lines total 9,699):

| Shard (of 16) | Zero lines | md5 |
|---|---|---|
| s0 | 575 | 7c1d0a26706142378a11418d67ca36f4 |
| s1 | 593 | 03ee13931ff81b827268d8752da9a9fe |
| s2 | 628 | ae0aeaf13fa24d6fbc44e1cafa1cd377 |
| s3 | 584 | 3a13849930190dde754c3e3eaaea0bbc |
| s4 | 625 | 12c0d7ec9ed703325f4a920b9be81867 |
| s5 | 565 | 8f354d4f57d1258f9811f1cfab905880 |
| s6 | 637 | 009d4fe9640dba54a99ec72dc0cafe19 |
| s7 | 608 | cbf438647c292a0f0052012329a06dda |
| s8 | 680 | 563af646966ebb8a4c879d72a5e2ff9a |
| s9 | 644 | 58599ce1dec683c91105218025a02249 |
| s10 | 585 | 63add7b154667be5cbbadcff1592d715 |
| s11 | 573 | ccf5bc395d69e3ec949cdc1074915a30 |
| s12 | 609 | acd80ded40e88de31a88206cba1bede8 |
| s13 | 553 | 8c08456006dfbfcd0c5fb2f8dd3b7c96 |
| s14 | 650 | 7351077a5bc3eeff47b2e08a922f09a1 |
| s15 | 590 | 01fb38c41896f58f1b2279a5d5f572d8 |

## 12. R11 result: gc (Ganymede-Callisto) k = 6 (lead launch; analysed 2026-10-08)

Run: 24 shards of 24, launched by the lead.
- 4,844 structures, all done once, 0 errors (`check` OK).
- 13,897 exact zeros, 4,806 physical cyclers, **4 gate-passing**, 0 zero-assessment errors.
- 0 structures at the n_refine ceiling.

Gauntlet and screens (`data/973_gc/k6_gauntlet.json`):
- DOP853 re-fly: largest miss 1.1e-4 km; gate "pass" on the integrated vectors for all 4.
- Screens (2.6/2.7): all 4 pass. No impact, no model-invalid, no direction-dependent.
- Classes: all 4 are W2.
- Literal collisions: one NEAR. gc6-0 (G 1.673 / C 1.257) is within 0.253 km/s of the #576 n = 3 G-C
  symmetric closure, at k = 3 against 6, with a different structure (#576 has no returns; gc6-0 has a
  2-rev Ganymede return). It is recorded, not a collision. No match with R-S GanCal, GCGC, Liang 2024
  or the catalogue rows.
- Literature step: deferred (#1025).

**Clean two-working-body members at k = 6: 4** (gc6-0, 1, 2, 3).

Resonance structure: 6 synodic periods (75.139 d) are 4.502 Callisto and 10.502 Ganymede periods.
That is half-period commensurate to 0.05 %. Sec. 11 predicted that this would favour the half-rev
(n-pi) returns. **The prediction is not borne out.** None of the 4 passers uses a half-rev leg: all
are generic Lambert returns (one has a 2-rev Ganymede return and 2-rev transfers). All need large
turns (worst ratios 0.77-0.92), as at k = 5.

The gate-passer count keeps falling with k: 51 (k = 4), 8 (k = 5), 4 (k = 6). The zero count keeps
rising: 4,629, 9,699, 13,897. Only k = 4's near-integer window (3 Callisto and 7 Ganymede periods)
produces the shallow-turn population.

Notable members:
- gc6-0 has the lowest Callisto V_inf of the route, 1.257 km/s. Callisto barely works (one flyby,
  9.7 deg, ratio 0.12), while Ganymede turns 63 deg twice (ratio 0.92).
- gc6-1 (G 1.607 / C 1.837) has the most margin: worst ratio 0.768, both moons turning 41-46 deg.
- All 4 keep r_min at or above 879,637 km and r_max at or below 2.89 Gm.

Candidate table, all 4 (period 75.1394 d; columns as in sec. 11; A = Ganymede, B = Callisto):

| # | Class | Screen | Structure (key) | V_inf A / B (km/s) | A: flybys, max turn (deg), max ratio | B: flybys, max turn, max ratio | Worst ratio | r_min / r_max (km) | Lambert starts (d) | Closest unscheduled A / B (km) | Re-fly miss (km) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gc6-0 | W2 | pass | k6\|LGanymede>Ganymede/2h\|LGanymede>Callisto/2h\|LCallisto>Ganymede/2h | 1.673 / 1.257 | 2, 63.08, 0.919 | 1, 9.68, 0.123 | 0.919 | 879,637 / 1,883,137 | 10.2125, 27.3572, 56.3546 | 134,590 / 779,588 | 1.0e-05 |
| gc6-1 | W2 | pass | k6\|LGanymede>Ganymede/2h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 1.607 / 1.837 | 2, 46.03, 0.645 | 2, 40.99, 0.768 | 0.768 | 1,062,677 / 2,466,968 | 1.1855, 36.3842, 44.0069, 68.7023 | 872,129 / 532,898 | 6.3e-06 |
| gc6-2 | W2 | pass | k6\|LGanymede>Ganymede/2l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/0s | 2.468 / 2.012 | 2, 32.46, 0.742 | 2, 41.06, 0.860 | 0.860 | 1,058,639 / 2,892,497 | 0.5986, 36.9711, 43.9480, 68.7611 | 1,014,134 / 579,791 | 1.1e-04 |
| gc6-3 | W2 | pass | k6\|LGanymede>Ganymede/1h\|LGanymede>Callisto/0s\|LCallisto>Callisto/1l\|LCallisto>Ganymede/1h | 1.664 / 2.192 | 2, 61.37, 0.889 | 2, 30.75, 0.720 | 0.889 | 1,043,144 / 2,584,336 | 3.8525, 24.1611, 32.7416, 57.6765 | 835,935 / 669,337 | 5.1e-06 |

Per-flyby detail (turn deg, ratio, required altitude km):

| # | Ganymede flybys | Callisto flybys |
|---|---|---|
| gc6-0 | 2 x 63.08 (0.919, 587) | 1 x 9.68 (0.123, 46,920) |
| gc6-1 | 2 x 46.03 (0.645, 3,330) | 2 x 40.99 (0.768, 1,540) |
| gc6-2 | 2 x 32.46 (0.742, 1,552) | 2 x 41.06 (0.860, 876) |
| gc6-3 | 2 x 61.37 (0.889, 794) | 2 x 30.75 (0.720, 1,734) |

Status: "candidate, literature step deferred" (#1025). Not novel.

Raw data (lead rule of 2026-10-08):
- Committed: the shards' `settings.json` and `structures.jsonl` (`data/973_gc/k6/s*/`), plus the
  gauntlet JSON and log.
- Not committed: `zeros.jsonl`.
- Full shard directories archived at `~/dev/references/cyclers-runs/973/gc/k6/` (21M; `diff -rq`
  clean). md5 of each `zeros.jsonl` (lines total 13,897):

| Shard (of 24) | Zero lines | md5 |
|---|---|---|
| s0 | 556 | 4ebdbf79cf95f8e995b7f90f8ebc7cb2 |
| s1 | 553 | 76eda3b037f1392a1298b895b27c7093 |
| s2 | 535 | 7e05daf37c7b4321b5188135b0a754f0 |
| s3 | 570 | bb02fb6f0ffaedacaf529743eeeb566d |
| s4 | 560 | 83d50bb29dcc40beb3f7f1e5879c81d2 |
| s5 | 577 | 8f6f6235e2011dcde6bc59c5f2ab7c32 |
| s6 | 554 | 7ef3c6d11e1ced1140a5a498bf1bb798 |
| s7 | 615 | 6c0b33af564129e48563797f20f9e243 |
| s8 | 579 | 15ef566aed16ebc25c1d9ace267a9048 |
| s9 | 615 | e023db8909d2e29d4c45a8fb62381a23 |
| s10 | 569 | e1e878abe302ca6059e0f0ce2627ae93 |
| s11 | 612 | e69a2fafad2634dc4276ed1fcfe552bc |
| s12 | 584 | 07f7086e35ca4f684ae71376378413f8 |
| s13 | 624 | 7a820d193280138c985e3ed3b3de3c1a |
| s14 | 579 | cf6f550670ea31f7a8a85c2981548451 |
| s15 | 596 | 7c08d8d5b7e8cff0297474af2510f3fd |
| s16 | 571 | 278ba3949f212ce530a9b7ea47e1ec73 |
| s17 | 648 | fd7a52300afa68bc05bdb297608da080 |
| s18 | 575 | 26259eecb8040a51fe3671c78e54a370 |
| s19 | 595 | 922a3e55336074774f6d8714b528afb9 |
| s20 | 515 | cdc105cf1c02819ddd51ff36d402bfed |
| s21 | 574 | ae3ab3dfb813c3957827fdb032f4a94c |
| s22 | 570 | e8f460b25931352ac1fd585268ffece6 |
| s23 | 571 | 05e4a4448372efbc7361ffde0052efbe |

## 13. R13 result: ev (Earth-Venus) k = 4 (lead launch; analysed 2026-10-08)

Run: 16 shards of 16, with the cell-5 flags (Earth returns full-rev only).
- 3,566 structures, all done once, 0 errors (`check` OK).
- 22,586 exact zeros, 6,206 physical cyclers, **12 gate-passing**, 0 zero-assessment errors.
- 0 structures at the n_refine ceiling.

Gauntlet and screens (`data/973_ev/k4_gauntlet.json`):
- DOP853 re-fly: largest miss 1.8e-3 km; gate "pass" on the integrated vectors for all 12.
- Screens (2.6/2.7): 11 pass, 1 "reject: planet impact", 0 direction-dependent. The rejected ev4-2 has
  a Venus HV(3,1,a) leg, the tilted circle, which meets Venus at 1/3 of the leg (2.7 (i), physical in
  the ideal model). It is not a published member.
- Classes: all 12 are W2. In 4 of them (ev4-4, 5, 8, 9) one of the three Earth flybys has a 0.0-deg turn (an Earth
  full-rev pair joined without a turn), but every member has at least one Earth and one Venus turn of
  1 deg or more.
- Literal collisions: none.
  - Hollister orbits I-III are k = 2 topologies only.
  - No catalogue row is within the bands. That includes the 15 H&M 1970 rows, which are k = 10, so
    only the 0.1 km/s "CATALOGUE V_INF" band applies to them.
  - Nearest H&M row by per-body V_inf, per transfer skeleton: orbit 5 (E 5.15 / V 5.6) at 0.69-1.74
    km/s.
  - Structurally no k = 4 one-visit structure can equal an H&M orbit (sec. 6.5 check).
- Literature step: deferred (#1025). Sources to attribute as in R13: Hollister 1969, H&M 1970,
  VanderVeen 1969.

**Clean members: 11.** They share only 5 distinct transfer skeletons (V_inf E / V; members differing
only in the order or type of the Earth and Venus return blocks):

| Skeleton (E / V km/s) | Members | Worst ratio (best member) | r_min / r_max (AU) of the best member |
|---|---|---|---|
| 4.462 / 5.197 | ev4-0 | 0.858 | 0.616 / 1.641 |
| 4.595 / 4.529 | ev4-1 | 0.835 | 0.672 / 1.641 |
| 6.806 / 4.017 | ev4-3 | 0.990 | 0.722 / 2.196 |
| 6.830 / 4.038 | ev4-4, 5, 6, 7 | 0.779 (ev4-4, 6, 7) | 0.496 / 1.168 (ev4-4) |
| 6.893 / 3.968 | ev4-8, 9, 10, 11 | 0.939 (all four; the Venus turn binds) | 0.495 / 1.173 (ev4-8) |

Resonance structure (model periods: Earth 1 yr, Venus 0.61520 yr; synodic 1.5988 yr):

| k | Period (yr) | in Earth periods | in Venus periods |
|---|---|---|---|
| 2 | 3.1975 | 3.198 | 5.198 |
| 3 | 4.7963 | 4.796 | 7.796 |
| 4 | 6.3950 | 6.395 | 10.395 |
| 5 | 7.9938 | 7.994 | 12.994 |

At k = 4 the cycle is 0.4 of a revolution off for both planets, with no near-commensurability. As at
gc k = 5, every passer turns hard: Earth 37-84 deg, Venus 31-94 deg, worst ratios 0.78-0.99.

Candidate table, all 12 (period 2,335.78 d = 6.395 yr; A = Earth, B = Venus; Lambert starts in days
from the model epoch, both planets at angle 0 at t = 0):

| # | Class | Screen | Structure (key) | V_inf A / B (km/s) | A: flybys, max turn (deg), max ratio | B: flybys, max turn, max ratio | Worst ratio | r_min / r_max (km) | Lambert starts (d) | Closest unscheduled A / B (km) | Re-fly miss (km) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ev4-0 | W2 | pass | k4\|RE/3:2\|RE/1:1\|LE>V/0s\|LV>V/1l\|RV/1:1\|LV>E/0s | 4.462 / 5.197 | 3, 83.80, 0.858 | 3, 32.81, 0.401 | 0.858 | 92,198,207 / 245,490,544 | 438.5278, 600.1102, 1151.7230 | 7,060,609 / 21,981,230 | 3.4e-04 |
| ev4-1 | W2 | pass | k4\|RE/1:1\|RE/3:2\|LE>V/0s\|RV/1:1\|RV/1:1\|LV>E/0s | 4.595 / 4.529 | 3, 79.94, 0.835 | 3, 66.70, 0.732 | 0.835 | 100,516,501 / 245,455,614 | 438.5278, 1100.6184 | 14,518,057 / 70,245,874 | 2.3e-04 |
| ev4-2 | W2 | reject: planet impact | k4\|RE/1:1\|RE/2:3\|LE>V/0s\|RV/3:2\|HV/3,1,a\|LV>E/0s | 5.148 / 4.214 | 3, 74.87, 0.849 | 3, 92.51, 0.965 | 0.965 | 78,700,757 / 175,861,400 | 547.8750, 1673.4678 | 118,677,025 / 1 | 7.9e-04 |
| ev4-3 | W2 | pass | k4\|RE/2:1\|RE/3:2\|LE>V/0s\|LV>E/0s | 6.806 / 4.017 | 3, 68.35, 0.990 | 1, 31.40, 0.317 | 0.990 | 107,988,489 / 328,536,056 | 329.1806, 583.9444 | 295,145,901 / 57,341,351 | 2.8e-04 |
| ev4-4 | W2 | pass | k4\|RE/2:3\|RE/2:3\|LE>V/0s\|LV>V/1l\|RV/1:1\|LV>E/0s | 6.830 / 4.038 | 3, 38.32, 0.557 | 3, 76.86, 0.779 | 0.779 | 74,275,175 / 174,747,357 | 325.3446, 570.9186, 1120.1041 | 44,329,531 / 16,983,253 | 3.5e-04 |
| ev4-5 | W2 | pass | k4\|RE/2:1\|RE/2:1\|LE>V/0s\|RV/1:1\|LV>V/1l\|LV>E/0s | 6.830 / 4.038 | 3, 57.79, 0.840 | 3, 76.86, 0.779 | 0.840 | 95,753,330 / 330,718,685 | 551.7110, 856.4308, 1180.9146 | 44,329,531 / 16,983,253 | 3.8e-04 |
| ev4-6 | W2 | pass | k4\|RE/1:1\|RE/3:2\|LE>V/0s\|LV>V/1l\|RV/1:1\|LV>E/0s | 6.830 / 4.038 | 3, 37.68, 0.548 | 3, 76.86, 0.779 | 0.779 | 95,753,330 / 254,757,029 | 325.3446, 570.9186, 1120.1041 | 44,329,531 / 16,983,253 | 2.5e-04 |
| ev4-7 | W2 | pass | k4\|RE/3:2\|RE/1:1\|LE>V/0s\|LV>V/1l\|RV/1:1\|LV>E/0s | 6.830 / 4.038 | 3, 37.68, 0.548 | 3, 76.86, 0.779 | 0.779 | 95,753,330 / 254,757,029 | 325.3446, 570.9186, 1120.1041 | 44,329,531 / 16,983,253 | 2.5e-04 |
| ev4-8 | W2 | pass | k4\|RE/2:3\|RE/2:3\|LE>V/0s\|RV/1:1\|LV>V/1h\|LV>E/0s | 6.893 / 3.968 | 3, 38.30, 0.562 | 3, 93.67, 0.939 | 0.939 | 74,105,790 / 175,511,189 | 322.9502, 802.9221, 1110.8676 | 37,916,347 / 15,920,721 | 4.2e-04 |
| ev4-9 | W2 | pass | k4\|RE/2:1\|RE/2:1\|LE>V/0s\|RV/1:1\|LV>V/1h\|LV>E/0s | 6.893 / 3.968 | 3, 56.67, 0.831 | 3, 93.67, 0.939 | 0.939 | 92,289,835 / 330,932,279 | 322.9502, 802.9221, 1110.8676 | 37,916,347 / 15,920,721 | 1.8e-03 |
| ev4-10 | W2 | pass | k4\|RE/1:1\|RE/3:2\|LE>V/0s\|RV/1:1\|LV>V/1h\|LV>E/0s | 6.893 / 3.968 | 3, 36.87, 0.541 | 3, 93.67, 0.939 | 0.939 | 92,289,835 / 255,020,599 | 322.9502, 802.9221, 1110.8676 | 37,916,347 / 15,920,721 | 4.2e-04 |
| ev4-11 | W2 | pass | k4\|RE/3:2\|RE/1:1\|LE>V/0s\|RV/1:1\|LV>V/1h\|LV>E/0s | 6.893 / 3.968 | 3, 36.87, 0.541 | 3, 93.67, 0.939 | 0.939 | 92,289,835 / 255,020,599 | 322.9502, 802.9221, 1110.8676 | 37,916,347 / 15,920,721 | 4.2e-04 |

Status: "candidate, literature step deferred" (#1025). Not novel.

Raw data (lead rule of 2026-10-08):
- Committed: the shards' `settings.json` and `structures.jsonl` (`data/973_ev/k4/s*/`), plus the
  gauntlet JSON and log.
- Not committed: `zeros.jsonl`.
- Full shard directories archived at `~/dev/references/cyclers-runs/973/ev/k4/` (36M; `diff -rq`
  clean). md5 of each `zeros.jsonl` (lines total 22,586):

| Shard (of 16) | Zero lines | md5 |
|---|---|---|
| s0 | 1407 | 58a128671615661e42ba451ee3f089f6 |
| s1 | 1386 | 57b11781f7b9848da90d95046de0a945 |
| s2 | 1458 | d003a7a3decd5519b190a9f83153b27c |
| s3 | 1416 | 8bed2b19eead766c3298e42380b838f1 |
| s4 | 1374 | 73f7e714297a66f3fdd70c85bf9c2a98 |
| s5 | 1431 | 33ea59d7f53714e422457f588c424ad3 |
| s6 | 1442 | 67455c8d7f0ddddba89c969f7bf5545e |
| s7 | 1404 | b0b6bd58113f90958dd7fd6ec58fa07a |
| s8 | 1436 | 7871bc0a4874c36216975f458ac11966 |
| s9 | 1373 | 31427b1e8b415c1c2abea3854326d760 |
| s10 | 1431 | 0af440a43ad98c3f79fffceb65b52533 |
| s11 | 1420 | 4ece325c485b04a19d740aa4644c9ceb |
| s12 | 1356 | b84fb4904d61621479d0c18b7c7b44bb |
| s13 | 1441 | e6c33a84b0796cc7732dffd406e76ad5 |
| s14 | 1426 | dde010805d48cb4d22f01976fb8972ab |
| s15 | 1385 | 3585354a11f1344fc803092dee83b8f0 |
