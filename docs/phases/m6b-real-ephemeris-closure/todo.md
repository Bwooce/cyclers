# M6b — Real-ephemeris closure verification (todo)

Working checklist for the M6b milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: predecessor recap → catalogue loader + tests → phase-match plumbing → `verify/real_closure.py` skeleton → helpers (each with paired unit tests) → `verify_real_closure` + gate tests → regression set → `optimise_cell_ephemeris` error-message edit → local green → CI → closeout.

## Predecessor recap

- [ ] Re-read spec §12(a) (idealised vs ephemeris-mode optimisation) — M6a built the verification half; **M6b uses it via the ephemeris backend** but does NOT wire `optimise_cell_ephemeris`.
- [ ] Re-read spec §12(c) (dynamic ephemeris frame + tolerant verification) — M6b reuses M6a's dynamic frame without changes.
- [ ] Re-read spec §12.1 (phase-match bridge) — M6b's construction path *composes* `phase_match.find_real_windows`; the M6-slice (commit `9b2611d`) already implemented it.
- [ ] Re-read spec §12.2 (the three-representation framework) — **M6b is the V1-idealised → V2-real-ephemeris-instance promotion gate**; do NOT flatten representations.
- [ ] Re-read spec §14 (V0–V5 gauntlet) — focus on V2 (real-ephemeris periodicity, M6b's gate) and V3 (TCM budget, M7's gate, deferred). Note the V2-real vs V2-idealised distinction.
- [ ] Re-read spec §16.1 (catalogue schema v2 fields) — M6b reads `model_assumption`, `trajectory_regime`, `primary`, `priority_date`, `bodies`, `legs[].tof_days`, `legs[].n_revs`, `vinf_kms_at_encounters[].vinf_kms`.
- [ ] Re-read spec §16.2 (canonical signature) — M6b's real-ephemeris instances **inherit their idealised parent's signature**; do NOT compute a new signature.
- [ ] Re-read M6a's `phases/m6a-idealized-closure-verification/plan.md` §3.3 (`StabilityReport` shape — M6b composes it), §4.3 (50,000 km tolerance derivation — M6b's 200,000 km tolerance follows the same logic but for real-eccentricity case), §5 risks (#1, #2, #3 carry over).
- [ ] Re-read M6a's hand-off note in `phases/m6a-idealized-closure-verification/todo.md` — note the frame-bodies decision, the actual `max_drift_km` reproduced, any escalations.
- [ ] Re-read M6a's `verify/propagate.py` end-to-end: `verify_long_term_stability` is the M6b binding composition target; `_resolve_frame_bodies` is the policy helper M6b inherits.
- [ ] Re-read M6 slice's `search/phase_match.py` end-to-end: `PhaseSignature`, `phase_signature_from_catalogue_entry`, `find_real_windows`, `LaunchWindow`. These are the M6b-consumed surfaces.
- [ ] Re-read Pascarella 2024 entry in `docs/v2-future-references.md` — the patched-conic → medium-fidelity pipeline; M6b is the medium-fidelity output.
- [ ] Sanity-check: `verify_long_term_stability(cycler, n_laps=2, ephem=Ephemeris("circular"))` works (M6a tested at `n_laps=3`; M6b defaults to `n_cycles=2`).
- [ ] Sanity-check: `Ephemeris("astropy").state("E", 0.0)` returns a real heliocentric position (≈ 1.496e8 km).
- [ ] Confirm M5 binding gate brokenness (task #54) — `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)` returns either an empty list or a non-closing cycler. M6b's construction path is **Lambert-chain from catalogue**, NOT via `find_cyclers`. Documented in plan §3.1.
- [ ] Confirm multi-rev Lambert blocker — `core/lambert.py::lambert(r1, r2, tof, max_revs=1)` raises `LambertGeometryError` for non-zero `max_revs` (per the M1 stub). M6b documents this as the binding limitation; `EXPECTED_SKIPS` carries the affected entries.
- [ ] Confirm parallel Schema-v2 backfill / other agents are NOT touching: `verify/propagate.py`, `core/frames.py`, `core/ephemeris.py`, `search/phase_match.py`, `data/seed_cyclers.yaml` (mid-write reads should fail gracefully; sequential coordination is the parent agent's job).

## Catalogue loader (test infrastructure)

- [ ] Create `tests/data/__init__.py` (empty).
- [ ] Create `tests/data/_catalogue_loader_m6b.py` with module docstring noting: "test infrastructure; M7's full catalogue loader (`data/catalog.py`) supersedes this."
- [ ] Implement `load_m6b_entries() -> list[dict]`:
  - [ ] Read `data/seed_cyclers.yaml` via `pyyaml`.
  - [ ] Filter to `model_assumption in (None, "circular-coplanar")`.
  - [ ] Filter to `trajectory_regime in (None, "ballistic")`.
  - [ ] Filter to `primary in (None, "Sun")`.
  - [ ] Apply V1-pass filter if `validation.gates.V1.pass` field exists; else accept all (M6a may not have written this field yet).
  - [ ] Use `entry.get(field, default)` per schema-v2 default rules (missing field = circular-coplanar / ballistic / Sun).
- [ ] Define `M6B_REGRESSION_IDS: Final[tuple[str, ...]]` per plan §3.2 (5 entries: aldrin-classic-em-k1-outbound, aldrin-classic-em-k1-inbound, mcconaghy-2006-em-k2, russell-ocampo-2.1.1+2-case2, russell-ocampo-2.5.1+0).
- [ ] Create `tests/data/test_catalogue_loader_m6b.py`:
  - [ ] `test_loader_filters_v1_pass_circular_coplanar_ballistic_sun_only`: returned count is between 180 and 215; every entry's `model_assumption in (None, "circular-coplanar")`.
  - [ ] `test_loader_excludes_cr3bp_entries`: no entry has `model_assumption == "cr3bp"`.
  - [ ] `test_loader_excludes_analytic_ephemeris_entries`: no entry has `model_assumption == "analytic-ephemeris"`.
  - [ ] `test_loader_excludes_non_sun_primaries`: no entry has `primary in ("Earth", "Jupiter", "Saturn")`.
  - [ ] `test_m6b_regression_ids_all_in_loader`: every id in `M6B_REGRESSION_IDS` is present in `load_m6b_entries()` output.
- [ ] Confirm `uv run pytest tests/data/` green.
- [ ] Confirm `uv run mypy tests/data/` clean.

## Phase-match plumbing test

- [ ] Add `test_find_real_windows_for_aldrin_signature_within_priority_window` to `tests/search/test_phase_match.py`:
  - [ ] Build `PhaseSignature` for Aldrin (`bodies=("E","M")`, `leg_durations_s=(146 * SECONDS_PER_DAY,)`, `vinf_target_kms=(5.5, 5.5)` — Aldrin's published V∞).
  - [ ] Call `find_real_windows(sig, Ephemeris("astropy"), (datetime(1980,1,1,tzinfo=UTC), datetime(1995,1,1,tzinfo=UTC)), n=3, mismatch_cap_kms=10.0)`.
  - [ ] Assert returned list is non-empty.
  - [ ] Assert at least one window's `departure_date` is within ±5 years of `datetime(1985, 10, 28, tzinfo=UTC)`.
- [ ] Confirm `uv run pytest tests/search/test_phase_match.py` green (existing tests remain unchanged).

## `verify/real_closure.py` — skeleton

- [ ] Create `src/cyclerfinder/verify/real_closure.py`:
  - [ ] Module docstring referencing spec §12.1, §12.2, §14 V2-real, §16.1; plan path; the "Pascarella 2024 medium-fidelity stage" architectural framing.
  - [ ] Imports: `Ephemeris`, `Cycler`, `verify_long_term_stability`, `multi_lap_propagation`, `lap_to_lap_drift`, `lambert`, `LambertConvergenceError`, `LambertGeometryError`, `phase_signature`, `phase_signature_from_catalogue_entry`, `find_real_windows`, `LaunchWindow`.
  - [ ] Define `REAL_DRIFT_TOLERANCE_KM: Final[float] = 200_000.0` with the §4.3 rationale in the docstring.
  - [ ] Define `N_CYCLES_DEFAULT: Final[int] = 2` with the trade-off rationale.
  - [ ] Define exception types: `MultiRevLambertRequiredError(Exception)`, `RealClosureConstructionError(Exception)`.
  - [ ] Define `RealClosureResult` frozen dataclass per plan §3.3 with all 11 fields:
    - `cycler_id: str | None`
    - `n_cycles_propagated: int`
    - `max_drift_km: float`
    - `per_cycle_drift_km: tuple[float, ...]`
    - `per_encounter_vinf_mismatch_kms: tuple[float, ...]`
    - `closes: bool`
    - `v3_status: str`
    - `horizon_tcm_mps: float`  *(M6b: 0.0; M7 populates)*
    - `per_cycle_tcm_mps: tuple[float, ...]`  *(M6b: zeros; M7 populates)*
    - `frame_used: str`  *(M6b: always "dynamic")*
    - `t_start_sec: float | None`
  - [ ] Define `EXPECTED_SKIPS: Final[dict[str, str]]` registry per plan §4.4: keys = catalogue ids, values = human-readable skip reason.
  - [ ] Stub `verify_real_closure`, `construct_real_ephemeris_cycler`, `_resolve_real_t_start`, `_construct_cycler_from_lambert_chain`, `_check_vinf_continuity` with full docstrings + `NotImplementedError` bodies.
- [ ] Update `src/cyclerfinder/verify/__init__.py` to re-export M6b's public surface per plan §1.2.
- [ ] Confirm `uv run mypy src/cyclerfinder/verify/real_closure.py` clean on the skeleton.
- [ ] Confirm `uv run mypy src` overall clean — the new module's exports don't break any downstream type checks.

## Helpers (one at a time, paired tests)

### `_resolve_real_t_start`

- [ ] Implement per plan §3.1.2: signature → `find_real_windows(signature, ephem, date_range_around_priority, n=5, mismatch_cap_kms=20.0)`; return the lowest-mismatch window's `_dt_to_t_sec(departure_date)`. Return `None` if no window beats the cap.
- [ ] Helper for date-range computation: ±10 years around `priority_date` (broad enough to catch any window; M6a's xfail test uses 8-year range).
- [ ] Test `test_resolve_real_t_start_prefers_priority_window` — Aldrin signature + 1985-10-28 priority → returned `t_sec` within ±5 years of priority.
- [ ] Test `test_resolve_real_t_start_returns_none_when_no_window` — synthesise a signature with absurd `vinf_target_kms = (50.0, 50.0)` → returns `None`.

### `construct_real_ephemeris_cycler`

- [ ] Implement per plan §3.1.1:
  - [ ] Extract `bodies`, `legs`, `period_years` from catalogue entry dict.
  - [ ] Compute encounter epochs from `t_start_sec + cumsum(legs[].tof_days * SECONDS_PER_DAY)`.
  - [ ] For each leg: read real planet positions at the two endpoint epochs via `ephem.state(body, t)`.
  - [ ] Check `leg.get("n_revs", 0)`; raise `MultiRevLambertRequiredError(catalogue_id, leg_index=j)` if > 0.
  - [ ] Call `lambert(r1, r2, tof, max_revs=0)`; catch `LambertConvergenceError` / `LambertGeometryError` and raise `RealClosureConstructionError(catalogue_id, leg_index=j, cause=e)`.
  - [ ] Build `Encounter` / `Leg` objects per the M3 dataclass shape.
  - [ ] Construct and return `Cycler(...)` with `period = period_years * SECONDS_PER_YEAR`.
- [ ] Test `test_construct_real_ephemeris_cycler_aldrin`:
  - [ ] Load Aldrin entry; call `_resolve_real_t_start` for a `launch_window`.
  - [ ] Call `construct_real_ephemeris_cycler(aldrin_entry, ephem, launch_window)`.
  - [ ] Assert returned `Cycler.bodies == ("E", "M", "E")` (Aldrin's encoded sequence — verify against catalogue).
  - [ ] Assert `cycler.encounters[0].t == t_start_sec` (the launch epoch).
  - [ ] Assert `_check_vinf_continuity(cycler, ephem)` returns small mismatches (< 2 km/s; Aldrin has no interior encounter so returns `()`).
- [ ] Test `test_construct_raises_on_multi_rev_leg`:
  - [ ] Build a synthetic catalogue entry with `legs[0].n_revs = 1`.
  - [ ] Assert `construct_real_ephemeris_cycler(...)` raises `MultiRevLambertRequiredError`.
- [ ] Test `test_construct_raises_on_lambert_geometry_error`:
  - [ ] Build an entry whose leg geometry is degenerate (zero TOF, or 180° transfer).
  - [ ] Assert `construct_real_ephemeris_cycler(...)` raises `RealClosureConstructionError`.

### `_check_vinf_continuity`

- [ ] Implement per plan §3.1: for each interior encounter (index 1..n-2), compute `|cycler.encounters[i].vinf_in| - |cycler.encounters[i].vinf_out|` and return the tuple of absolute differences.
- [ ] Test `test_check_vinf_continuity_diagnostic`:
  - [ ] For Aldrin (2-encounter chain), returns `()` (no interior encounter).
  - [ ] For a 3-encounter cycler (e.g. McConaghy 2-syn with E-M-E sequence), returns length-1 tuple.

## `verify_real_closure` — compose helpers

- [ ] Implement per plan §3.1.2 pipeline:
  - [ ] Resolve `t_start` if `None`: signature → `_resolve_real_t_start(signature, ephem, signature_priority_date)`. If still `None`, return `RealClosureResult(..., v3_status="v3-no-real-window", closes=False, max_drift_km=inf, ...)`.
  - [ ] If `cycler` is a dict, call `construct_real_ephemeris_cycler(...)`. Catch `MultiRevLambertRequiredError` → return result with `v3_status="v3-skipped-multirev"`. Catch `RealClosureConstructionError` → return result with `v3_status="v3-construction-error"`.
  - [ ] Delegate to `verify_long_term_stability(cycler, n_laps=n_cycles, ephem, t_start, frame_bodies, cycler_id, use_uniform_frame=False)`.
  - [ ] Compute `closes = stability.max_drift_km < REAL_DRIFT_TOLERANCE_KM`.
  - [ ] Compute `v3_status = "v3-real-closure-pass" if closes else "v3-real-closure-fail"`.
  - [ ] Compute `_check_vinf_continuity` diagnostics.
  - [ ] Build and return `RealClosureResult` with V3 placeholders (`horizon_tcm_mps=0.0`, `per_cycle_tcm_mps=(0.0,)*n_cycles`).
- [ ] Test `test_real_closure_uses_m6a_machinery` (**M6b composition gate** — fail-fast for reimplementation):
  - [ ] Use `unittest.mock.patch` on `cyclerfinder.verify.real_closure.verify_long_term_stability`.
  - [ ] Configure return value to a hand-built `StabilityReport`.
  - [ ] Call `verify_real_closure(...)`; assert the patched function was called exactly once with `n_laps=2`.
- [ ] Test `test_aldrin_cycler_periodic_over_2_cycles_astropy` (**M6b BINDING GATE — spec §8**):
  - [ ] Load `aldrin-classic-em-k1-outbound` from catalogue.
  - [ ] Call `verify_real_closure(entry, n_cycles=2, ephem=Ephemeris("astropy"), signature_priority_date=datetime(1985,10,28,tzinfo=UTC), cycler_id="aldrin-classic-em-k1-outbound")`.
  - [ ] Assert `result.closes == True`, `result.max_drift_km < REAL_DRIFT_TOLERANCE_KM`, `result.n_cycles_propagated == 2`, `result.v3_status == "v3-real-closure-pass"`, `result.horizon_tcm_mps == 0.0`, `result.per_cycle_tcm_mps == (0.0, 0.0)`, `result.frame_used == "dynamic"`, `result.cycler_id == "aldrin-classic-em-k1-outbound"`.
  - [ ] **If this fails, escalate per plan §5 risk #1.** DO NOT widen `REAL_DRIFT_TOLERANCE_KM`. Try: (a) loosening `find_real_windows` mismatch_cap to 15.0; (b) `frame_bodies=("E",)` per M6a hand-off ambiguity; (c) sanity check on circular ephemeris (where drift should be near-zero). Document findings.
- [ ] Test `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (**xfail**):
  - [ ] Decorate with `@pytest.mark.xfail(strict=False, reason="multi-rev Lambert blocker — see §5 risk #2")`.
  - [ ] Load `s1l1-2syn-em-cpom`.
  - [ ] Same gate-assertion shape as the Aldrin test.
  - [ ] Test must run; xfail will mark it as `XFAIL` if it fails (expected) or `XPASS` if multi-rev Lambert was secretly implemented (informational, non-failing).
- [ ] Test `test_real_drift_rejects_open_trajectory`:
  - [ ] Construct an Aldrin-derived cycler.
  - [ ] Manually rotate one `vinf_out` by 5 degrees (build a new `Cycler` with the perturbed `Encounter`).
  - [ ] Call `verify_real_closure(perturbed_cycler, n_cycles=2, ephem)`.
  - [ ] Assert `result.closes == False`, `result.max_drift_km > 5 * REAL_DRIFT_TOLERANCE_KM`, `result.v3_status == "v3-real-closure-fail"`.
- [ ] Test `test_real_closure_result_frozen_and_v3_fields_locked`:
  - [ ] Build a result via the Aldrin gate path.
  - [ ] Assert `result.closes = ...` raises `FrozenInstanceError`.
  - [ ] Assert `result.horizon_tcm_mps == 0.0`.
  - [ ] Assert `result.per_cycle_tcm_mps == (0.0,) * result.n_cycles_propagated`.
  - [ ] Assert `result.v3_status in ("v3-real-closure-pass", "v3-real-closure-fail", "v3-skipped-multirev", "v3-no-real-window", "v3-construction-error")`.
  - [ ] Assert `result.frame_used == "dynamic"`.

## Regression set

- [ ] Implement `test_real_closure_regression_set` parametrised over `M6B_REGRESSION_IDS`:
  - [ ] For each id: load entry; if id in `EXPECTED_SKIPS`, `pytest.skip(EXPECTED_SKIPS[id])`; else call `verify_real_closure(entry, n_cycles=N_CYCLES_DEFAULT, ephem=Ephemeris("astropy"), signature_priority_date=_parse_date(entry["priority_date"]))`.
  - [ ] Assert `result.closes` is `True`.
  - [ ] Assert `result.max_drift_km < REAL_DRIFT_TOLERANCE_KM`.
- [ ] Run the regression set; for each failure, decide:
  - [ ] If `v3_status == "v3-skipped-multirev"`: move id to `EXPECTED_SKIPS` with reason `"multi-rev Lambert blocker"`. Document follow-up.
  - [ ] If `v3_status == "v3-no-real-window"`: investigate per risk #3; if a per-entry `mismatch_cap_kms` relax works, encode it; else move to `EXPECTED_SKIPS`.
  - [ ] If `v3_status == "v3-real-closure-fail"`: investigate per risk #1; if drift 200k–500k, document as borderline (may indicate REAL_DRIFT_TOLERANCE_KM needs raising for this milestone; consult plan §4.3); if > 500k, escalate as gate failure.
- [ ] Finalise `EXPECTED_SKIPS` registry: the M6b-shipping set with the discovered blockers.
- [ ] Document any per-entry `mismatch_cap_kms` overrides in the regression test (parametrise per id).

## `optimise_cell_ephemeris` error message update

- [ ] Edit `src/cyclerfinder/search/optimize.py::optimise_cell_ephemeris` body's `NotImplementedError` message per plan §3.4:
  - [ ] Old text: `"requires M6 ephemeris backend (multi-lap propagator + TCM budget). Stub locked in M5; body lands in M6b."`
  - [ ] New text: `"requires M6b real-ephemeris closure (shipped: verify.real_closure.verify_real_closure) AND M7 TCM budget machinery (not yet shipped). M6b's verify_real_closure is the right drift-feasibility check; full wiring lands in M7."`
- [ ] Confirm `tests/search/test_optimize.py::test_optimise_cell_ephemeris_raises_not_implemented` still passes (the test should only assert the exception class, not the message text; if it asserts message text containing "M6", the assertion is still satisfied).

## Sanity: M0–M6a regressions

- [ ] `uv run pytest tests/core/` green — M0–M3 core unchanged.
- [ ] `uv run pytest tests/model/` green — M3 Aldrin gate unchanged.
- [ ] `uv run pytest tests/search/` green — M4 / M5 / M6-slice unchanged.
- [ ] `uv run pytest tests/verify/test_propagate.py` green — M6a unchanged; the xfailed `test_2syn_em_cycler_periodic_over_3_laps_astropy` still xfailed (M6b's xfail for a similar reason; M6a's astropy test is on a slightly different fixture).

## Local green

- [ ] `uv run pytest` → green (M0–M6b all passing; xfail tests still xfail; expected skips still skipped).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean.
- [ ] `uv run mypy src tests` → clean under strict mode.

## CI

- [ ] Commit: `m6b: real-ephemeris closure verification (verify_real_closure + Aldrin gate over 2 cycles)`.
- [ ] Push; confirm GitHub Actions runs and all checks pass.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M6b status → `completed`; M7 row → `planned`. Add a §2 / §4 note: "M6b adds real-ephemeris closure verification (`verify.real_closure.verify_real_closure`, V2-real gate per spec §14). TCM budget computation (V3 gate) deferred to M7. Aldrin closes over 2 cycles on DE440 within 200,000 km drift (the M6b binding gate); 5 regression entries close per spec; multi-rev-Lambert and VEM entries documented in `EXPECTED_SKIPS`."
- [ ] (Optional) Update `docs/overview.md` §5 anchors table: add row "Aldrin real-ephemeris closure over 2 cycles, drift < 200,000 km, M6b gate."
- [ ] Append a `## Hand-off to M7` section to this `todo.md` (below).

## Hand-off to M7

*(Filled in at M6b closeout. Placeholder structure below; replace bullets with measured values from the M6b test pass.)*

- **Actual `max_drift_km` reproduced for the Aldrin gate:** `<TBD>` km (vs the 200,000 km bound). Cycle pair with worst drift: `<TBD>`. The headroom `(200,000 - max_drift_km)` is M7's pre-TCM budget; M7's TCM optimiser should drive drift below `POST_TCM_DRIFT_TOLERANCE_KM` (M7 defines).
- **Regression-set per-entry outcomes:**
  - `aldrin-classic-em-k1-outbound`: `<closes? max_drift_km>`
  - `aldrin-classic-em-k1-inbound`: `<closes? max_drift_km>`
  - `mcconaghy-2006-em-k2`: `<closes? max_drift_km>`
  - `russell-ocampo-2.1.1+2-case2`: `<closes? max_drift_km>`
  - `russell-ocampo-2.5.1+0`: `<closes? max_drift_km>`
- **Final `EXPECTED_SKIPS` list:** `<TBD>` (catalogue ids + reasons). This is the input M7's "outstanding entries" follow-up consumes. Multi-rev Lambert is the dominant cause — implementing it is the highest-leverage M-future task.
- **`find_real_windows` mismatch_cap sufficiency:** `<TBD: default cap (5.0 km/s) was sufficient for N of 5 entries; M of 5 needed per-entry overrides>`. Document per-entry overrides (id → cap).
- **Frame-bodies decision for the Aldrin gate:** the M6a `_resolve_frame_bodies(cycler, None)` policy returned `<TBD>`; this was `<TBD: validated as correct | overridden because …>`. M7's batch runner should default to the same policy.
- **Per-test wall-clock runtime:**
  - `test_aldrin_cycler_periodic_over_2_cycles_astropy`: `<TBD>` s.
  - `test_real_closure_regression_set` (full): `<TBD>` s.
  - This informs M7's batch-runner compute budget. Full catalogue (200+ entries) over n_cycles=5 would take ~`<estimated total>` minutes.
- **`n_cycles=2` vs `n_cycles=3` anomaly check:** running `verify_real_closure(aldrin, n_cycles=3, ...)` produced `<TBD>` km drift vs `<TBD>` at n_cycles=2. If the ratio is roughly linear, drift accumulates as expected; if it jumps, escalate per risk #13.
- **`v3_status` distribution across regression set:**
  - `v3-real-closure-pass`: `<TBD>` entries
  - `v3-real-closure-fail`: `<TBD>`
  - `v3-skipped-multirev`: `<TBD>`
  - `v3-no-real-window`: `<TBD>`
  - `v3-construction-error`: `<TBD>`
- **Locked `RealClosureResult` shape** (M7 inherits exactly this):
  - `cycler_id: str | None`
  - `n_cycles_propagated: int`
  - `max_drift_km: float`
  - `per_cycle_drift_km: tuple[float, ...]`
  - `per_encounter_vinf_mismatch_kms: tuple[float, ...]`
  - `closes: bool`
  - `v3_status: str`
  - `horizon_tcm_mps: float`  *(M7 populates)*
  - `per_cycle_tcm_mps: tuple[float, ...]`  *(M7 populates)*
  - `frame_used: str`
  - `t_start_sec: float | None`
- **Contract for M7's catalogue writer:**
  - MUST check `result.v3_status == "v3-real-closure-pass"` before writing `validation.gates.V2.pass = true` and `validation.gates.V2.max_drift_km = result.max_drift_km`.
  - MUST NOT write any `validation.gates.V3.*` field until M7's TCM-budget runner produces a new `RealClosureResult` with populated TCM fields.
  - Per spec §16.2: do NOT compute a new canonical signature for the real-ephemeris instance — inherit the idealised parent's signature.
- **API surfaces M7 will consume immediately:**
  - `from cyclerfinder.verify.real_closure import verify_real_closure, RealClosureResult, REAL_DRIFT_TOLERANCE_KM, N_CYCLES_DEFAULT, EXPECTED_SKIPS, construct_real_ephemeris_cycler`
  - `from cyclerfinder.verify.propagate import verify_long_term_stability, StabilityReport` — M7's TCM optimiser uses `StabilityReport.max_drift_km` as drift-feasibility check inside its objective.
  - `from cyclerfinder.search.optimize import optimise_cell_ephemeris` — M6b left this raising; M7 fills the body using both `verify_real_closure` (drift feasibility) and the new M7 TCM-cost objective.
- **Spec ambiguities encountered during M6b implementation:**
  - **TCM application semantics.** Spec §12(a) is silent on *when* (cycle-start? mid-cycle? at flybys?) TCMs are applied. M7's plan should decide explicitly.
  - **Per-cycle TCM vs per-encounter TCM.** `RealClosureResult.per_cycle_tcm_mps` is the locked field shape; M7 may need finer granularity. Options: (a) extend with a `per_encounter_tcm_mps` field with default `()`, backward-compatible; (b) add a parallel `TCMReport` dataclass.
  - **5-cycle vs 2-cycle horizon.** Spec §12(a) says 3–5 lap horizon. M6b binds at 2; M7 must extend. Verify that drift remains < tolerance at n_cycles=5 on the Aldrin gate before committing M7 to the 5-cycle horizon.
- **Recommended M7 first steps:**
  1. Read this hand-off + spec §12(a) (ephemeris-mode optimisation) + spec §14 V3 (TCM budget gate).
  2. Write the M7 plan, deciding: TCM application policy, per-cycle vs per-encounter, n_cycles horizon.
  3. Extend `RealClosureResult` with whatever shape M7 needs (additive only; do not reshape the locked fields).
  4. Wire `optimise_cell_ephemeris` against `verify_real_closure` as drift-feasibility check + new TCM-cost objective.
  5. Build catalogue YAML writer that consumes `RealClosureResult` and populates `validation.gates.V2.*` (always) and `validation.gates.V3.*` (only after TCM-budget runs).
  6. Run the full-catalogue batch verification (200+ entries) offline; identify the closes / open / multi-rev-blocked / no-window subsets.
  7. Implement multi-rev Lambert per `core/lambert.py` (the M-future task that unlocks `s1l1-2syn-em-cpom` and similar entries). This may be a parallel M7 sub-task or a dedicated M-future.

M7's first task is the M7 plan doc itself; this hand-off is the input it consumes.
