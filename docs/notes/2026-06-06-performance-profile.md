# Performance profile — evidence-based optimisation plan (Task #128 Stage 1)

**Date:** 2026-06-06
**Author:** profiling pass (read-only; no source/test/data changed)
**Scope:** measure, don't change. Targets: one `ballistic_correct` solve, the
ephemeris layer, scan scaling, fast-suite time, the maintenance solve.

## Load-contamination caveats (READ FIRST)

This box (16 cores) was **shared with two other compute agents** for the entire
run. `uptime` load average ranged **2.8 → 20.6** across measurements (16-core
box; >16 = oversubscribed). Mitigations applied:

- Single-solve / per-call timings: `nice -n 19 taskset -c 0-3` + warm-up +
  repeated medians (min reported alongside).
- Scan-scaling: `taskset -c 0-7`, medians of 3, run when load had dropped to
  ~6–9. **Even so the absolute scan numbers are noisy** (see Target 3) — the
  *shape* of the scaling curve is the trustworthy signal, not the absolutes.
- cProfile numbers carry profiler overhead (≈4–2× wall) and are used only for
  **relative** attribution, never as wall-clock truth. Wall-clock medians are
  reported separately and are the authority for "how long".
- A fast-suite **FAILED** test (`test_schema_v45_fields.py::test_live_aldrin_outbound_is_v1_everything_else_v0`)
  appeared on the second full run. The working tree changed *during* this
  profiling pass — another agent began editing `data/catalogue.yaml`,
  `src/cyclerfinder/data/validate.py`, and that exact test file. **The failure
  is in-flight WIP from another agent, NOT a baseline regression and NOT caused
  by this read-only profiling.** First full suite run (before those edits
  landed) had no such failure.

---

## Target 1 — one `ballistic_correct` solve (S1L1 E-M-E-E cell)

Invocation: the `test_correct_s1l1.py` cell — `sequence=("E","M","E","E")`,
`per_leg_revs=(0,1,2)`, `branch=("single","single","low")`, `slack_leg=2`,
`tof_seed=(154,379)`. Both backends converge (`vinf_mars` ≈ 10.5–10.8).

### Wall-clock medians (warm, `taskset -c 0-3`, nice 19)

| backend  | reps | median | min | max |
|----------|------|--------|-----|-----|
| circular | 9    | **31.6 ms** | 31.3 | 31.9 |
| astropy  | 5    | **360.9 ms** | 343.4 | 370.0 |

`state()` calls per solve = **124** (circular) / **120** (astropy) — same
optimiser path; the 11x wall difference is entirely the per-call ephemeris cost.

### cProfile cumulative top (verbatim, profiler overhead included)

**circular** (4.189 s under profiler; state is negligible — Lambert + scipy jacobian dominate):
```
   ncalls  tottime  cumtime  filename:lineno(function)
        1    0.000    4.189  search/correct.py:302(ballistic_correct)
        1    0.000    4.182  scipy/optimize/_lsq/least_squares.py:266(least_squares)
       31    0.072    3.563  search/correct.py:86(_vinf_nodes)
       93    0.213    3.489  core/lambert.py:511(lambert)
        7    0.201    2.913  scipy/optimize/_numdiff.py:278(approx_derivative)   <- finite-diff jacobian
        1    0.000    2.878  {scipy.optimize._minpack._lmder}
       93    0.381    1.375  core/lambert.py:170(_solve_single_rev_newton)
        1    0.206    1.287  scipy/optimize/_differentiable_functions.py:542(__init__)
```
Circular state total ≈ **0.5 ms / solve** (124 calls × 3.8 µs). The circular
solve is **Lambert-bound + scipy-finite-difference-jacobian-bound**, NOT
ephemeris-bound.

**astropy** (1.569 s under profiler; ephemeris dominates):
```
   ncalls  tottime  cumtime  filename:lineno(function)
        1    0.000    1.569  search/correct.py:302(ballistic_correct)
       30    0.002    1.559  search/correct.py:86(_vinf_nodes)
       30    0.000    1.402  search/correct.py:117(<listcomp>)        <- the state() list-comp
      120    0.007    1.400  core/ephemeris.py:252(state)             <- 89% of the solve
        7    0.000    1.356  scipy/optimize/_numdiff.py:278(approx_derivative)
      240    0.001    0.604  astropy .../get_body_barycentric_posvel  <- 2 posvel calls per state()
1800/1560    0.008    0.512  astropy/units/core.py:2018(__call__)     <- units conversion churn
```
**~89% of the astropy solve is `Ephemeris.state()`.** Lambert and the scipy
jacobian structure are identical to circular but are now in the noise.

---

## Target 2 — the ephemeris layer (the real lever)

### Per-call cost (warm process, kernel mmapped)

| measurement | median |
|---|---|
| COLD first `state()` in a warm process (kernel already on disk) | 11.2 ms |
| WARM `state(E, **SAME** epoch)` repeated | 2.73 ms |
| WARM `state(E, **DISTINCT** epoch)` | 2.67 ms |
| WARM `state(M, distinct)` | 2.32 ms |

**astropy does NOT cache per-(body,epoch) state.** Same-epoch repeat costs the
same as a fresh epoch (2.73 vs 2.67 ms) — there is no kernel-segment / result
memoisation we benefit from at this call granularity.

The "18.5 s cold vs 12 s warm" in the brief is **process-level** (fresh
interpreter: import astropy + first DE440 mmap), not per-call. Within a warm
process the marginal `state()` is ~2.7 ms.

### Per-call breakdown (each `state()` does **two** `get_body_barycentric_posvel`)

| component | median |
|---|---|
| `Time(unix, scale="tdb")` construction | 0.73 ms |
| `get_body_barycentric_posvel("earth")` (Time fixed) | 0.73 ms |
| `get_body_barycentric_posvel("sun")` | 0.42 ms |

Plus the `.xyz.to("km")` units conversions and the ICRS→ecliptic rotation. The
`state()` body in `_AstropyBackend.state` (ephemeris.py:252) calls posvel
**twice** (body + sun) and constructs a fresh `Time` each call.

From the maintenance `tottime` view (Target 5), the irreducible Chebyshev math
(`jplephem/spk.py:202 generate`) is **3.89 s of 31 s** — astropy's `Time` +
units wrapper machinery around it (`day_frac`, `_get_time_fmt`, ply/yacc
`parseopt_notrack`, units `_new_view`/`__array_finalize__`/`_expand_and_gather`)
costs **roughly 2× the actual ephemeris math**. Most of the per-call cost is
astropy framework overhead, not DE440 evaluation.

### Memoisation ceiling

One S1L1 solve makes **120 `state()` calls but only 70 distinct (body,epoch)
pairs** → **42% are exact duplicates** within a single solve (the scipy
finite-difference jacobian holds `t0` fixed and perturbs only the ToFs, so the
`t0` encounter and the unperturbed columns recompute the same epochs
repeatedly; the slack-leg reconstruction also re-hits shared epochs).

Estimated ceilings on the astropy solve (state = 89% of 360 ms ≈ 320 ms):
- **Exact-(body,epoch) memo within a solve:** 120→70 calls ⇒ ≈ −42% of state ⇒
  solve ≈ **360 → ~225 ms (~−37% wall)**.
- **Hoist the Sun state per epoch** (Sun is fetched once per `state()` but is the
  same for body-E and body-M at the same epoch; with an epoch grid the Sun
  posvel is shared) ⇒ removes the 0.42 ms Sun posvel from ~half the calls.
- **Batch/vectorise `Time` + posvel over the whole epoch grid** (astropy
  accepts array `Time`): collapses 70 `Time()` constructions + per-call units
  wrapping into a handful of vectorised calls ⇒ the framework overhead (≈2× the
  Chebyshev math) largely disappears; only the ~0.5–1 ms/epoch Chebyshev math
  remains. Plausible solve **360 → ~80–120 ms (~−70%)**, ceiling set by the
  irreducible SPK math.

---

## Target 3 — scan scaling (`scan_parallel` vs `scan_serial`)

`build_epoch_branch_grid`, astropy backend, `taskset -c 0-7`, medians of 3,
measured when load had fallen to ~6–9. The grid API exposes `max_workers` but
**no chunksize knob** (`pool.map` default chunking) — chunked task sizing is not
tunable without a source edit, so not attempted.

| workers | median | speedup vs serial | efficiency | ideal |
|---|---|---|---|---|
| serial  | 3.08 s (193 ms/pt) | 1.00x | — | — |
| 2 | 1.64 s | **1.88x** | 94% | 2x |
| 4 | 1.01 s | **3.05x** | 76% | 4x |
| 8 | 0.93 s | **3.30x** | 41% | 8x |

This **reproduces the brief's sub-linear ~4x/8-worker result** (I see 3.30x; the
shortfall vs 3.98x is the two other agents contending for my 8 pinned cores).

**Same curve on the circular backend** (32-pt grid, pure CPU, no astropy I/O):
2→1.82x, 4→2.73x, 8→**3.16x (39% eff)**. The sub-linearity is **backend-
independent** ⇒ it is NOT astropy-I/O or kernel-contention; it is general CPU /
memory-bandwidth contention on an oversubscribed shared box.

**Worker-startup cost is cheap and NOT the bottleneck:** a 1-point parallel job
at w=4/8/16 = **0.14–0.16 s flat** (pool spin-up + per-worker
`Ephemeris('astropy')` construction). `Ephemeris('astropy')` ctor in a warm
process ≈ **0 ms** — `ProcessPoolExecutor` forks, so the DE440 mmap is inherited,
not reloaded per worker. The "per-worker Ephemeris construction" feared in
scan.py's docstring is essentially free under fork.

**Bottleneck verdict:** not startup amortisation, not task granularity (each
task ≈ 360 ms, far above the ~150 ms pool overhead). It is **core/memory-
bandwidth contention** — here aggravated by the two co-tenant agents. On an idle
box this should approach the brief's ~4x; the real ceiling is hardware
parallelism, not the scan code.

---

## Target 4 — fast suite (`-m 'not slow' -n auto`)

Default `addopts`: `-ra -q --strict-markers --timeout=600 -n auto -m 'not slow'`.
No WIP test files needed scoping at the time of the first run (only source files
were dirty). Total wall ≈ **3.4 min** (matches the ~3 min brief; xdist
`-n auto` already on).

`--durations=25` (verbatim head):
```
108.63s  test_optimise_ephemeris_mode.py::test_maintenance_mode_result_unchanged_for_aldrin_cell
 82.98s  test_catalogue_rediscovery_tagging.py::test_rediscovered_2syn_em_tagged_known
 72.54s  test_discover.py::test_discover_em_k2_yields_known_for_2syn
 64.00s  test_sequence.py::test_deepening_frontier_no_repeats
 34.12s  test_discover.py::test_discover_accepts_multirev_params
 28.84s  test_vem_rediscovery.py::test_all_multibody_rows_period_round_trips_via_loader
 28.80s  (setup) test_crosscheck.py::test_crosscheck_cycler_represents_single_rev_legs
 28.13s  test_real_closure.py::test_2syn_em_cpom_periodic_over_2_cycles_astropy
 27.16s  (setup) test_crosscheck.py::test_v1_lambert_crosscheck_aldrin
 27.05s  (setup) test_crosscheck.py::test_crosscheck_leg_missing_endpoint_raises
 26.47s  test_real_closure.py::test_real_closure_uses_m6a_machinery
 25.73s  test_real_closure.py::test_aldrin_powered_cycler_solver_and_drift_floor_on_de440
 24.85s  test_real_closure.py::test_resolve_real_t_start_picks_low_mismatch_window
 21.80s  test_real_closure.py::test_aldrin_ballistic_closure_fails_because_powered
 ... (then VEM-rediscovery folds 17–19 s, optimize 13–17 s, loader 12–13 s)
```
The top ~15 are **all astropy-backed** maintenance/discover/real-closure/
crosscheck tests. The same ephemeris per-call cost (Target 2) drives suite time.
A `state()` epoch-grid memo / batch (and/or a session-scoped astropy
`Ephemeris` fixture so the kernel mmaps once per xdist worker) would compound
across these dozen tests. The three crosscheck **setup** entries (27–29 s each)
suggest an expensive per-test fixture rebuilding astropy state — a session/
module-scoped fixture is a likely cheap win there.

---

## Target 5 — the maintenance solve (Aldrin E-M-E, astropy)

`optimise_aldrin_maintenance_dv(Ephemeris("astropy"), real_window_priority_date=…)`.

| backend | wall |
|---|---|
| circular | **0.54 s** (median of 3) |
| astropy | **16.75 s** (median, warm) — matches the brief's 12–18 s |

### cProfile cumulative top (verbatim, 31.2 s under profiler):
```
   ncalls   cumtime  filename:lineno(function)
        1    31.247  maintain.py:717(optimise_aldrin_maintenance_dv)
    18351    27.062  core/ephemeris.py:332(state)              <- 87% ephemeris-bound
        1    16.289  maintain.py:478(optimise_maintenance_dv)  (DE + SLSQP phase)
     3192    16.086  maintain.py:387(_objective)
     3193    15.831  search/construct.py:40(construct_cycler)
        1    14.958  maintain.py:441(_resolve_aldrin_real_t0_guess)  (phase-match window scan)
        1    14.958  search/phase_match.py:567(find_candidate_windows)
     4386    14.939  search/phase_match.py:273(_mismatch_at_date)
        1    11.238  scipy .../_differentialevolution.py(differential_evolution)
       60    11.057  _differentialevolution.py:1601(__next__)   <- 60 DE generations
    36702     9.702  astropy get_body_barycentric_posvel        <- 2x per state()
```
### tottime top (where the CPU actually burns):
```
   ncalls   tottime  filename:lineno(function)
   142422    3.891   jplephem/spk.py:202(generate)          <- irreducible Chebyshev SPK math
    91755    0.718   astropy/time/utils.py:19(day_frac)       \
  1101060    0.658   astropy/units/quantity.py __array_finalize__ |
   403722    0.611   astropy/units/quantity.py _new_view          | astropy framework
   110106    0.599   astropy/units/core.py _expand_and_gather     | overhead ≈ 2x the
    36702    0.578   astropy _get_body_barycentric_posvel         | real math above
    18351    0.524   astropy/extern/ply/yacc.py parseopt_notrack /
   410787    0.574   core/lambert.py:144(_dt_dz)            <- Lambert: small slice
  1071588    0.452   core/_stumpff.py stumpff_s
```
**Verdict: decisively EPHEMERIS-bound, not optimiser-bound.** 18,351 `state()`
calls = 87% of cumulative. Two roughly-equal sub-phases, *both* ephemeris-bound:
1. `_resolve_aldrin_real_t0_guess` → phase-match window scan (`_mismatch_at_date`
   × 4386) ≈ **15 s** — a launch-window grid that hammers `state()`.
2. DE (60 generations × popsize, 3192 `_objective` → `construct_cycler` → Lambert
   + `state()`) + SLSQP polish ≈ **16 s**.

The Lambert/Stumpff inner math is a small slice (~1.5 s). Of the ephemeris time,
the actual SPK Chebyshev (`spk.generate`, 3.9 s) is **half** of the astropy
per-call cost; the other half is `Time`/units wrapper churn.

---

## Ranked optimisation plan (gain × risk × algorithm-survivability)

> "Survives #122" = survives the residual-mode redesign (magnitude→vector). All
> ephemeris-layer wins are in `core/ephemeris.py`, **below** the
> residual/optimiser layer, so they are orthogonal to #122 and to any optimiser
> retune.

### 1. Ephemeris epoch-grid memoisation / batching (BIGGEST, SAFEST)
- **Evidence:** state() = 89% of the astropy solve (Target 1), 87% of the
  maintenance solve (Target 5). 42% of per-solve calls are exact (body,epoch)
  duplicates (Target 2). astropy does no caching of its own.
- **Two tiers:**
  - (a) **memoise exact (body,epoch)** results in `_AstropyBackend` (small LRU
    keyed on `(body, round(t_sec, ~6))`): est. **−37% on a single solve**, and
    larger on the maintenance/phase-match scans which revisit the same launch
    epochs across DE generations and the window grid.
  - (b) **vectorise `Time` + `get_body_barycentric_posvel` over the epoch grid**
    and **share the Sun state per epoch**: collapses the ~2× astropy `Time`/units
    framework overhead; ceiling ≈ **−70% (360→~100 ms/solve; 16.75→~5–7 s
    maintenance)**, limited by irreducible SPK Chebyshev math (3.9 s tottime).
- **Risk:** LOW. Pure provider-layer change; result values bit-identical (memo)
  or within float tolerance (vectorised path returns the same DE440 states).
  Must respect the existing ICRS→ecliptic rotation and the `(r,v)` tuple
  contract. Memo cache must be per-`Ephemeris` instance (per worker process) so
  it never crosses the pickling boundary — fits scan.py's existing
  construct-in-worker design.
- **Survivability:** HIGH — below the residual layer; survives #122 entirely and
  every optimiser retune. Compounds across Target 4's dozen astropy tests
  (consider also a session/worker-scoped `Ephemeris` fixture).

### 2. Single-`state()` Sun-state reuse (CHEAP SUBSET OF #1)
- **Evidence:** every `state()` calls `get_body_barycentric_posvel` **twice**
  (body 0.73 ms + Sun 0.42 ms). The Sun term is identical for all bodies at a
  given epoch but is recomputed per body. With the E/M encounters at shared
  epochs the Sun posvel is recomputed redundantly.
- **Gain:** removes ~0.42 ms from a large fraction of the 120/18,351 calls —
  order ~10–15% of state time on its own; a stepping-stone if full batching (1b)
  is deferred.
- **Risk:** VERY LOW. **Survivability:** HIGH (provider layer).

### 3. Phase-match window-scan epoch memo (maintenance-specific)
- **Evidence:** `_resolve_aldrin_real_t0_guess` → `_mismatch_at_date` × 4386 ≈
  15 s = ~half the maintenance solve, all `state()`. The DE pass then re-queries
  overlapping epochs.
- **Gain:** a shared epoch cache across the window scan **and** the subsequent DE
  pass (i.e. #1a with a cache that lives for the whole `optimise_*` call) could
  remove a large share of the 14.9 s scan and the DE re-queries — plausibly
  **−30–40% of the maintenance wall** beyond what a per-solve memo gives.
- **Risk:** LOW–MED (need to confirm the window scan and DE query the *same*
  epochs often enough; cache hit-rate not directly measured here — measure
  before committing). **Survivability:** HIGH (provider-layer cache; orthogonal
  to optimiser tuning).

## DO-NOT-BOTHER (measured, found cheap or not the bottleneck)

- **Lambert micro-optimisation** (`_solve_single_rev_newton`, `_stumpff`,
  `_dt_dz`): on astropy it is **in the noise** (state = 89%); on the maintenance
  solve Lambert+Stumpff ≈ 1.5 s of 16.75 s. Churn-risk against #122's residual
  redesign, tiny payoff. Skip. (Only the *circular* path is Lambert-bound, and
  circular solves are already 31 ms — not a target.)
- **scipy `least_squares` / `differential_evolution` overhead:** the optimiser
  scaffolding (minpack `_lmder`, `approx_derivative` structure) is small vs the
  user objective. The maintenance solve is ephemeris-bound *inside* the
  objective, not optimiser-bound. Don't swap solvers for speed.
- **Worker / pool startup, per-worker `Ephemeris` construction:** measured flat
  ~0.15 s for a 1-point job at any worker count; ctor ≈ 0 ms warm (fork shares
  the DE440 mmap). The scan.py docstring's worry about per-worker Ephemeris cost
  is unfounded under fork. Not a lever.
- **scan task chunking / granularity:** tasks are ~360 ms each, well above the
  ~150 ms pool overhead; the sub-linear scaling is core/bandwidth contention
  (identical on the I/O-free circular backend), not short tasks. (Also: no
  chunksize knob without a source edit.) Tuning chunksize won't fix
  contention-bound scaling.
- **Ephemeris construction caching of the ICRS→ecliptic rotation matrix:**
  already precomputed once in `_AstropyBackend.__init__`. Nothing to do.
- **`Ephemeris("circular")` anything:** 31 ms/solve, 0.54 s maintenance. Fast
  enough; it's the test-surface fast path and not in any hot loop that matters.

---

## Reproduction (scripts were in /tmp, removed after; recreate to re-measure)

- T1/T2: warm + median `ballistic_correct` on the `test_correct_s1l1` cell for
  both backends; wrap `ephem.state` in-harness to count calls; per-call timing of
  `Time()` / `get_body_barycentric_posvel` in isolation.
- T3: `build_epoch_branch_grid` → `scan_serial` vs `scan_parallel(max_workers=…)`,
  medians of 3, both astropy and circular backends; 1-point parallel job for
  startup cost.
- T5: `optimise_aldrin_maintenance_dv(astropy, real_window_priority_date=…)` under
  cProfile (cumulative + tottime) and wall-clock median.
- T4: `pytest --durations=25 -p no:cacheprovider`.
All under `nice -n 19` + `taskset` to limit cross-talk with co-tenant agents.
