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
