# #973: the two-working-body generator at longer periods (gc k = 4-6, ev k = 4-5) and a Europa-Callisto cell (ec k = 1-4)

Status: PRE-REGISTRATION (secs. 1-2), written and committed before any production run. Results follow
in secs. 3-6 as each route finishes.
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
