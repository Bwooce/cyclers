# Multi-rev Lambert, Full-3D, VEM Enumeration & Real-Ephemeris Discovery — Master Roadmap

> **For agentic workers:** This is a MASTER ROADMAP spanning four subsystems, not a single executable plan. Each milestone below gets its own detailed `docs/superpowers/plans/YYYY-MM-DD-<milestone>.md` written with the `superpowers:writing-plans` skill before execution. M-L is ready for task-level planning now; M-3D and M-ED need a design pass (brainstorming) first.

**Goal:** Make the ~217 currently-unreachable catalogue rows (multi-encounter E-E-M-M families, VEM triples, and families the circular-coplanar model cannot host) constructible, optimisable, verifiable, and discoverable on real ephemeris.

**Architecture:** Four dependent milestones. Multi-rev Lambert (M-L) is the only genuinely-unbuilt *core* capability and is foundational. Full-3D geometry (M-3D) lifts the coplanar assumption. N-encounter + VEM enumeration (M-N) opens the discrete search/loader layer. Real-ephemeris blind discovery (M-ED) implements the `optimise_cell_ephemeris` stub + M7 TCM-budget machinery and wires the discovery loop. Critical path: `M-L → (M-3D ∥ M-N) → M-ED`.

**Tech stack:** Python 3.11, numpy, scipy (`differential_evolution` + SLSQP), astropy DE440 ephemeris, pytest + xdist, uv-managed venv, ruff + mypy.

---

## Scope decisions (locked 2026-06-02)

| Decision | Choice | Consequence |
|---|---|---|
| Real-ephemeris scope | **Full blind discovery** | M-ED implements `optimise_cell_ephemeris` + M7 TCM budgeting + `enable_v3` discovery wiring, not just reproduction of known rows. |
| Geometry | **Full 3D inclination** | M-3D lifts the coplanar assumption across frames / construct / propagate / closure. Idealized circular model may stay 2D as a fast pre-filter. |
| VEM / moon rows | **Include VEM now** | M-N adds multi-body beat period + multi-synodic enumeration. (Non-heliocentric moon tours stay under separate task #76 — planet-centric, orthogonal.) |
| Ephemeris source | **Stay on DE440** | No change needed. DE440 already supplies full 3D states; pin + document the kernel and its 1550–2650 validity window. See "Ephemeris note" below. |

### Ephemeris note (resolves "are we getting fresh ephemeris from JPL?")

We already use the authoritative JPL ephemeris: `ephemeris.py:117` calls `solar_system_ephemeris.set("de440")`; states come from `get_body_barycentric_posvel` (line 144), heliocentric, rotated ICRS→ecliptic (line 151). astropy fetches the DE440 SPK kernel from the JPL/NAIF mirror on first use and caches it under `~/.astropy/cache` (the in-code comment "bundled with astropy" is inaccurate — only the smaller `de432s` ships in-wheel; `de440` is downloaded). A planetary DE is a versioned release, not a live feed, so there is no staleness to chase. DE441 (longer span, same accuracy) and live Horizons (identical DE440 numbers for planets, just slower) add nothing for our epochs/tolerance. **The full-3D work needs no fresher data** — DE440 already returns 3D states and we currently keep the z-component; the coplanar limitation is downstream, not in the data. Action item (M-3D): explicitly pin the expected kernel and assert it at backend init for reproducibility.

---

## Verified current state ("is it already built?")

Investigated read-only on 2026-06-02. Summary: **only multi-rev Lambert is genuinely unbuilt at the core**; the N-encounter compute stack and the real-ephemeris maintenance optimiser already exist and are general — gated, stubbed at the edges, or Aldrin-locked, not missing.

### Multi-rev Lambert
- **EXISTS:** `lambert(r1, r2, tof, *, mu, prograde, max_revs)` signature returning `list[LambertSolution]` with `n_revs`/`branch`/`v1`/`v2` (`core/lambert.py:166`). Universal-variable + Stumpff (`core/_stumpff.py`). `Leg.n_revs`/`branch` fields (`model/cycler.py:139-140`). `construct_cycler` already accepts + threads `max_revs_per_leg`/`branch_per_leg` (`search/construct.py:45,141,198`). Catalogue segments carry `n_revs`.
- **STUBBED/GATED:** `lambert()` validates but **ignores** `max_revs` and always returns the single-rev solution (`core/lambert.py:217,383-390`). Downstream raises `MultiRevLambertRequiredError` when `n_revs>0`: `verify/real_closure.py:476-481` (→ status `"v3-skipped-multirev"` at :749), `verify/crosscheck.py:128-131`. `s1l1-2syn-em-cpom` is in EXPECTED_SKIPS (`real_closure.py:241`).
- **MISSING:** the solver math (z-domain to `(2π(n+1))²`, per-n bracketing, low/high branch) + multi-solution selection in `construct.py`.

### N-encounter construction / optimisation
- **ALREADY ARBITRARY-N (verified):** `construct_cycler` (`search/construct.py:40`, validates only `n>=2` + monotonic times); `Cycler`/`Encounter`/`Leg` (`model/cycler.py:144`, `len(legs)==len(bodies)-1`, closure uses `encounters[0]`/`[-1]`); `enumerate_cells` (`search/sequence.py:145`, `for length in range(2, l_max+1)`); `optimise_cell_idealized` (`search/optimize.py:960`, free-vars `n_interior = len(cell.sequence) - 2` at :1053); `_build_chain`/`optimise_maintenance_dv` (`search/maintain.py:303,389`).
- **HARDCODED-TO-2-BODY:** `_target_period_sec` uses only `cell.bodies[0:2]` (`search/optimize.py:229`) — fine for 2-body E-E-M-M, **breaks for VEM** (needs multi-body beat period; `search/resonance.py` has `multi_body_beat_days`). Loader gate `_is_two_body_alternation` enforces exactly 2 encounters (`tests/_catalogue_loader.py:134-152`), classifying multi-encounter rows as `MULTI_ENCOUNTER_SEQUENCE` and excluding them.
- **MISSING:** structural inference mapping a catalogue `E-E-M-M @ k-synodic` row → a `Cell` with correct per-leg `n_revs`/`branch`.

### Real-ephemeris mode
- **EXISTS & GENERAL:** `optimise_maintenance_dv` is sequence/body-agnostic and runs on `Ephemeris("astropy")` (`search/maintain.py:389`); `construct_real_ephemeris_cycler`/`verify_real_closure` general for single-rev chains (`verify/real_closure.py:445,642`); DE440 returns full 3D states (`core/ephemeris.py:129`).
- **EXISTS-BUT-ALDRIN-LOCKED:** public BVP `solve_powered_periodic_cycler` pins E-M-E via `optimise_aldrin_maintenance_dv` (`search/bvp.py:92,149`).
- **MISSING:** `optimise_cell_ephemeris` is a `NotImplementedError` stub (`search/optimize.py:1144-1178`); `enable_v3` discovery path guarded off (`data/discover.py:94`); no real-ephemeris driven discovery; coplanar construction despite 3D states (no inclination handling).

---

## Milestones

### M-L · Multi-rev Lambert (foundational, no deps)

**Objective:** `lambert(..., max_revs=N)` returns the full solution set (`n=0` single + `n∈[1,N]` low/high branches); construction selects the catalogue-specified `(n_revs, branch)`; the `n_revs>0` rejection gates are removed.

**Touch points:**
- `core/lambert.py` — extend the z-bracket scan beyond `z_high_single_rev = 4π²` to per-revolution domains `z ∈ [(2πn)², (2π(n+1))²]`; for each `n` the time-of-flight curve has a minimum splitting low/high branches — bracket and Newton-solve each; return all converged solutions.
- `search/construct.py:141-199` — replace the "≤1 solution" assumption with selection of the requested `(n_revs, branch)` from the returned list.
- `verify/real_closure.py:476-481,528` — delete the `MultiRevLambertRequiredError` raise; thread `max_revs` from the leg into `lambert()`.
- `verify/crosscheck.py:128-131,177` — allow `n_revs>0` legs through the crosscheck.

**Design questions to settle in the detailed plan:** branch-selection contract (how the catalogue's `branch` string maps to low/high); minimum-tof feasibility check per `n` (some `n` have no solution for a given tof — return empty, don't error); whether to keep the existing single-rev fast path untouched for `max_revs=0` (preserve current bit-for-bit results).

**Test gates that flip green:** `tests/verify/test_real_closure.py::test_construct_raises_on_multi_rev_leg` becomes a *construction* test; `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (currently XFAIL) and `test_2syn_em_cycler_periodic_over_3_laps_astropy` (XFAIL) target passing; new unit tests for `lambert(max_revs=1..2)` against Vallado worked examples (sourced golden values only — never values our own solver computed).

**Ready for task-level planning:** YES.

---

### M-3D · Full 3D / inclination geometry (deps: none, but pairs with M-L)

**Objective:** Lift the coplanar assumption so construction, frame transforms, propagation, and drift all respect real planetary inclination from DE440's 3D states.

**Touch points (to be confirmed during design):**
- `core/ephemeris.py` — pin/assert the DE440 kernel; keep full 3D (already does).
- `core/frames.py` — `synodic_omega` + rotating-frame transforms currently assume the ecliptic plane is the orbit plane; generalise the frame definition to a 3D angular-velocity vector / inclined reference.
- `search/construct.py` — Lambert is already 3D-native (takes 3D `r1,r2`); verify encounter/leg geometry carries z.
- `verify/propagate.py` + `verify/real_closure.py` — 3D multi-lap drift in the dynamic rotating frame; revisit `REAL_DRIFT_TOLERANCE_KM`.
- Decision: keep the idealized circular model 2D as a fast pre-filter, or add inclination there too.

**Design pass required first** (brainstorming): the frame-definition change is the crux and has multiple valid formulations. Write a design doc before the TDD plan.

**Test gates:** new 3D-state assertions (sourced inclinations from DE440); existing coplanar tests must either still pass or be consciously re-baselined with documented reasons.

**Ready for task-level planning:** NO — needs design first.

---

### M-N · N-encounter + VEM enumeration & rediscovery (deps: M-L)

**Objective:** Surface multi-encounter (E-E-M-M) and VEM (3-body) catalogue rows to the optimiser and rediscovery gauntlet.

**Touch points:**
- `tests/_catalogue_loader.py:134-152` — generalise `_is_two_body_alternation` to admit N-encounter sequences; reclassify `MULTI_ENCOUNTER_SEQUENCE`/`NOT_TWO_BODY` rows that are now constructible. Keep the coverage-audit census invariant intact.
- New: structural inference — given a catalogue `sequence_canonical` + `period.k`, infer per-leg `n_revs`/`branch` (the intermediate same-body loop legs are multi-rev) to build the matching `Cell`.
- `search/optimize.py:229` — `_target_period_sec`: add multi-body beat period (use `search/resonance.py::multi_body_beat_days`) for VEM; keep the 2-body fast path.
- `search/sequence.py` — confirm/extend `enumerate_cells` produces VEM cells with correct per-leg revs/branch defaults.

**Test gates:** the `MULTI_ENCOUNTER_SEQUENCE` exclusion count drops; new rediscovery cases for at least one E-E-M-M family and one VEM triple (skipped/xfailed until M-ED can actually converge them on ephemeris — wire as XFAIL with documented reasons, flip on M-ED).

**Ready for task-level planning:** YES once M-L lands (the loop legs need multi-rev to construct).

---

### M-ED · Real-ephemeris blind discovery (deps: M-L, M-3D, M-N)

**Objective:** Implement `optimise_cell_ephemeris` and the M7 TCM-budget machinery; wire the `enable_v3` discovery loop so `find_cyclers` can discover families directly on DE440.

**Touch points:**
- `search/optimize.py:1144-1178` — implement the stub: real-ephemeris periodic BVP over a `Cell` (the general `optimise_maintenance_dv` already does the heavy lifting on astropy; this wraps it with the cell constraints + TCM budget objective).
- `search/maintain.py` / `search/bvp.py` — expose the general `optimise_maintenance_dv` / generalise `solve_powered_periodic_cycler` off the Aldrin lock.
- New: M7 TCM-budget machinery (per-family maintenance-ΔV bounds) — referenced by the stub's own error message as not-yet-shipped.
- `data/discover.py:94-100` — flip `enable_v3` on and route discovery through the ephemeris optimiser.

**Design pass required first** (brainstorming): TCM-budget definition + the discovery objective/constraints. Largest and least-specified milestone.

**Test gates:** a previously-unreachable published family (e.g. S1L1 5.65/3.05) rediscovered on real ephemeris to its sourced anchors; the M5 `test_2syn_em_rediscovers_5_65_kms_earth` XFAIL flips to a passing ephemeris-mode test; `test_ephemeris_mode_stubbed_until_m6` retired.

**Ready for task-level planning:** NO — needs design first (TCM budget).

---

## Sequencing & handoff

1. **M-L** — detail + execute now (foundational, well-understood, ~3 files).
2. **M-3D** — design pass (brainstorming) → detailed plan → execute. Can start design in parallel with M-L execution.
3. **M-N** — detail + execute after M-L.
4. **M-ED** — design pass → detailed plan → execute last, after M-L/M-3D/M-N.

**Discipline (carried from project memory):** golden/validation tests assert only source-attested expected values (never values our own code computed); run `uv run ruff check .` + `ruff format --check .` + `mypy` before each commit; never commit without explicit user request; no AI attribution in commit messages.
