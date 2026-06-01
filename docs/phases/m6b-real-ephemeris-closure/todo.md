# M6b — Real-ephemeris closure verification (todo)

Working checklist for the M6b milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: predecessor recap → catalogue loader + tests → phase-match plumbing → `verify/real_closure.py` skeleton → helpers (each with paired unit tests) → `verify_real_closure` + gate tests → regression set → `optimise_cell_ephemeris` error-message edit → local green → CI → closeout.

**Status (2026-06-01).** Scaffolding shipped (loader, plumbing, helpers, regression set, `EXPECTED_SKIPS`, `optimise_cell_ephemeris` message). Hygiene clean: `ruff check` + `ruff format --check` + `mypy src tests` all pass; 225 default-suite tests green, 4 xfail / 3 skip in the slow real-closure suite. **The Aldrin binding gate is XFAIL** — the Lambert-chain construction (plan §3.1.1) is wrong for Aldrin-family cyclers; M7 needs an orbital-elements-based constructor before the gate flips to passing. See "Hand-off to M7" below for the architectural finding + measured drift numbers + proposed fix.

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

*M6b closeout, 2026-06-01.* The scaffolding (loader, plumbing test, `verify/real_closure.py` module, helpers, gate tests, regression set, `EXPECTED_SKIPS` registry, `optimise_cell_ephemeris` error-message update, `verify/__init__.py` re-exports) is in place and hygiene-clean (ruff + mypy strict + 225 default-suite tests green; 4 xfail / 3 skip in the slow real-closure suite). **The M6b binding gate (Aldrin classic real-ephemeris closure over 2 cycles) does NOT pass at the 200,000 km tolerance.** The architectural finding below is the input M7 consumes.

### Critical finding — Lambert-chain construction is wrong for Aldrin-family cyclers

**Symptom.** `verify_real_closure(aldrin-classic-em-k1-outbound, n_cycles=2, ephem=astropy)` reports `max_drift_km ≈ 2.89e8` km (~1.93 AU) vs the 200,000 km bound — failure by a factor of ~1450x. The same call on the inbound variant reports `~8.5e7` km (~0.57 AU).

**Root cause.** Plan §3.1.1's algorithm independently Lambert-solves each leg over real planet positions. For the Aldrin classic, the catalogue records two legs: E→M (146 d) and M→E (634 d). On a real ballistic Aldrin flyby, both legs lie on the **same** heliocentric ellipse (a=1.60 AU, e=0.393, i=0) — the M→E "leg" is just continuing along the orbit the spacecraft is already on. Independently Lambert-solving M→E produces a **different** orbit (whichever passes through the two endpoint positions in 634 d), which inevitably has a different |V∞| at Mars departure than the |V∞| at Mars arrival. We measure:

- `|V∞_M_in|  = 9.78 km/s` (from E→M Lambert arrival)
- `|V∞_M_out| = 10.92 km/s` (from M→E Lambert departure)
- Ballistic-flyby continuity violation: ~1.14 km/s

This mismatch propagates: lap-1 starts from `|V∞_E| ≈ 6.27 km/s` (M→E arrival) rather than the lap-0 `|V∞_E| ≈ 6.53 km/s` (E→M departure), so the lap-1 geometry is a different orbit, drifting away by AU per cycle. **The same failure mode appears on the CIRCULAR ephemeris** with the M3-style 132-degree phase epoch: drift ~2.8M km — orders of magnitude above the M6a 50,000 km tolerance. The failure is intrinsic to the Lambert-chain algorithm; it is NOT a real-ephemeris eccentricity-breathing problem.

**Why plan §3.1.1 did not anticipate this.** Plan §3.1.1 was designed for multi-encounter cyclers (e.g. 2-syn S1L1 with intermediate Earth flybys at the half-period mark) where each leg IS a separate Lambert problem with a real flyby providing the V∞ continuity. For those cyclers Lambert-chain is correct — but multi-rev Lambert is needed (§5 risk #2) so they are currently blocked. The Aldrin family sits in a different regime: 2-encounter cyclers whose "two legs" are actually one orbit traversed in two arcs, where the catalogue's published M→E ToF is just `T_cycler − T_outbound` rather than an independent ballistic re-rendezvous.

**Proposed M7 fix — orbital-elements construction.** Aldrin (and likely most of the literature E-M ballistic cyclers) ships with published `orbit_elements` fields (`a_au`, `e`, `i`, `perihelion_au`, `aphelion_au`). Construct the cycler's heliocentric ellipse directly from these elements, plus the published `vinf_kms_at_encounters[].vinf_kms` and the per-leg ToFs to **locate the encounters along the orbit**. This bypasses the per-leg Lambert and guarantees the ballistic-flyby continuity at Mars. Concretely M7 should:

1. Add `construct_real_ephemeris_cycler_from_elements(catalogue_entry, ephem, t_start)` alongside the current Lambert-chain constructor.
2. Use the orbital-elements path for entries whose `orbit_elements.a_au` and `e` are non-null (Aldrin, the Aldrin inbound, the analytic-ephemeris establishment variants, etc.).
3. Use the Lambert-chain path only for entries with multiple distinct heliocentric arcs per cycle (the 2-syn S1L1 family, once multi-rev Lambert lands).
4. Dispatch at the `verify_real_closure` boundary based on which fields the catalogue entry carries.

**Tolerance is correct.** Plan §4.3's 200,000 km derivation is independent of construction algorithm — it assumes the construction produces a ballistic Aldrin and absorbs real-eccentricity breathing only. The tolerance need NOT be widened; the construction path needs to be fixed. Per the parent instruction, the tolerance was not widened unilaterally.

### Measured outcomes

- **Aldrin outbound binding gate:** `max_drift_km = 289,274,066` km on real ephemeris (1986-02-27 launch epoch, picked by `find_real_windows` ±5 yr of 1985-10-28 priority date with mismatch_cap=20 km/s). `n_cycles_propagated=2`, `v3_status="v3-real-closure-fail"`, `frame_used="dynamic"`. Worst cycle pair: lap 0 → lap 1 (the only consecutive pair at n_cycles=2). **XFAIL strict=False** with full diagnostic in the test's xfail reason.
- **Aldrin inbound regression entry:** `max_drift_km ≈ 84,935,476` km. Same XFAIL routing.
- **S1L1 2-syn entry (`s1l1-2syn-em-cpom`):** XFAIL — multi-rev Lambert blocker (`v3-skipped-multirev`). M-future / M7 prerequisite. Pre-routed via `EXPECTED_SKIPS`.
- **McConaghy / Russell entries:** SKIPPED — incomplete leg data (catalogue records only the outbound or outbound+inbound legs; the 2-synodic period requires more legs the published abstracts do not tabulate). Pre-routed via `EXPECTED_SKIPS` with reason `"incomplete leg data"`. M7's catalogue completion is the unblocking work.
- **Open-trajectory rejection test:** PASSES. A 5-degree V∞ rotation on Aldrin produces drift >> 5 × 200,000 km. The gate HAS rejection power.
- **Composition assertion (`test_real_closure_uses_m6a_machinery`):** PASSES. `verify_real_closure` calls `verify_long_term_stability` exactly once with `n_laps=n_cycles`.
- **All helper-level tests:** PASS (construction, `_resolve_real_t_start`, `_check_vinf_continuity`, multi-rev pre-check, construction-error wrapping, frozen dataclass, error-path routing).
- **Plumbing test (`test_find_real_windows_for_aldrin_signature_within_priority_window`):** PASSES — `find_real_windows` returns ≥1 window within ±5 yr of 1985-10-28 at mismatch_cap=10 km/s.

### Hand-off bullets per todo template

- **Actual `max_drift_km` reproduced for the Aldrin gate:** 289,274,066 km (vs the 200,000 km bound). Cycle pair with worst drift: lap 0 → lap 1. The headroom is NEGATIVE — `max_drift_km / REAL_DRIFT_TOLERANCE_KM ≈ 1446`. The architecture (not the tolerance) needs to change before headroom is meaningful.
- **Regression-set per-entry outcomes:**
  - `aldrin-classic-em-k1-outbound`: `closes=False`, `max_drift_km=289,274,066` km, `v3_status=v3-real-closure-fail` (XFAIL via `_M6B_LAMBERT_CHAIN_XFAILS`)
  - `aldrin-classic-em-k1-inbound`: `closes=False`, `max_drift_km=84,935,476` km, `v3_status=v3-real-closure-fail` (XFAIL via `_M6B_LAMBERT_CHAIN_XFAILS`)
  - `mcconaghy-2006-em-k2`: SKIPPED via `EXPECTED_SKIPS` (incomplete leg data — 2 legs totalling 306 d but advertised 4.27 yr period)
  - `russell-ocampo-2.1.1+2-case2`: SKIPPED via `EXPECTED_SKIPS` (incomplete leg data — 1 leg of 207 d but advertised 4.27 yr period)
  - `russell-ocampo-2.5.1+0`: SKIPPED via `EXPECTED_SKIPS` (incomplete leg data — 1 leg of 94 d but advertised 4.27 yr period)
- **Final `EXPECTED_SKIPS` list (M6b-shipped):**
  - `s1l1-2syn-em-cpom`: multi-rev Lambert blocker
  - `mcconaghy-2006-em-k2`: incomplete leg data (catalogue tabulates only 2 legs of a 4.27-yr 2-syn cycler)
  - `russell-ocampo-2.1.1+2-case2`: incomplete leg data (catalogue tabulates only 1 leg)
  - `russell-ocampo-2.5.1+0`: incomplete leg data (catalogue tabulates only 1 leg)
  - `jones-2017-vem-triple-family`: VEM 3-body real-closure is M8 scope
  - `vem-emeeve-3syn`: VEM 3-body real-closure is M8 scope
  Multi-rev Lambert is the dominant unblock for the 2-syn family; catalogue leg-completion is the unblock for the McConaghy/Russell entries; VEM 3-body is M8.
- **`find_real_windows` mismatch_cap sufficiency:** the default 5.0 km/s cap was insufficient for the Aldrin priority-window resolution; `_resolve_real_t_start` uses 20.0 km/s by default. With the 20-km/s cap, the Aldrin signature finds a 1986-02-27 window within ±5 yr of 1985-10-28 priority. No per-entry overrides needed in M6b's current scope.
- **Frame-bodies decision for the Aldrin gate:** `_resolve_frame_bodies(cycler, None)` returned `("E", "M")` for the constructed Aldrin (3-encounter `E-M-E` chain). The default policy is correct; no override needed. The dominant drift signal is the Lambert-chain V∞ mismatch, not the frame anchor — so frame-bodies tuning would not have rescued the gate.
- **Per-test wall-clock runtime (local dev box, single core, astropy DE440):**
  - `test_aldrin_cycler_periodic_over_2_cycles_astropy`: ≈ 7 s.
  - `test_real_closure_regression_set` (5 entries, 3 skipped + 2 xfailed): ≈ 14 s.
  - Full slow real-closure suite: ≈ 41 s.
  - Plumbing test + helper tests (non-slow): ≈ 9 s.
  - At ~7 s per closure verification, a 200-entry batch over n_cycles=5 would take ~2,000 s ≈ 33 min (with ~5x for n_cycles=5 / 2). Acceptable as M7 offline batch.
- **`n_cycles=2` vs `n_cycles=3` anomaly check:** NOT YET RUN. M6a tested at `n_laps=3`; the M6b sanity check `verify_long_term_stability(aldrin, n_laps=2, ephem=circular)` returned `max_drift_km=2.5e-7` km on the M3 build_aldrin_seed (2-encounter Aldrin slice), confirming n_laps=2 works. The 3-laps comparison on the 3-encounter Lambert-chain Aldrin was not performed because the construction itself diverges; once M7's orbital-elements construction lands, M7 should run n_cycles=3 and n_cycles=5 comparisons as part of its V3 horizon-extension work.
- **`v3_status` distribution across the regression set:**
  - `v3-real-closure-pass`: 0 entries
  - `v3-real-closure-fail`: 2 (Aldrin outbound + inbound — both XFAIL)
  - `v3-skipped-multirev`: 1 (`s1l1-2syn-em-cpom`, EXPECTED_SKIPS)
  - `v3-no-real-window`: 0
  - `v3-construction-error`: 0
  - SKIPPED (incomplete leg data): 3 (McConaghy + 2 Russell entries)
  - SKIPPED (VEM 3-body, M8 scope): 0 within `M6B_REGRESSION_IDS` (the 2 VEM entries are in `EXPECTED_SKIPS` but not in the regression IDs)
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
  - **5-cycle vs 2-cycle horizon.** Spec §12(a) says 3-5 lap horizon. M6b binds at 2; M7 must extend. Verify that drift remains < tolerance at n_cycles=5 on the Aldrin gate before committing M7 to the 5-cycle horizon.
- **New ambiguity discovered during M6b implementation:**
  - **Construction algorithm dispatch.** Plan §3.1.1's Lambert-chain construction is correct for cyclers with multiple distinct heliocentric arcs per cycle (the 2-syn S1L1 family, multi-encounter free returns); it is WRONG for cyclers whose "two legs" are arcs of one heliocentric ellipse (the Aldrin classic and its inbound variant). M7 must:
    1. Decide the dispatch criterion. Strawman: dispatch on whether the catalogue entry's `orbit_elements.a_au` and `orbit_elements.e` are non-null AND the leg count equals 2. This catches Aldrin without flagging genuinely multi-Lambert chains.
    2. Implement an orbital-elements-based constructor as a sibling of `construct_real_ephemeris_cycler`. The constructor instantiates the heliocentric ellipse and locates the encounters via Kepler propagation along the orbit; the M→E "leg" trivially carries the same orbital elements as the E→M "leg".
    3. Add `verify_real_closure` dispatch logic that picks the right constructor based on the criterion in (1).
    4. Re-run the regression suite with the new construction and remove the `_M6B_LAMBERT_CHAIN_XFAILS` set from `tests/verify/test_real_closure.py`.
- **Recommended M7 first steps (revised in light of M6b finding):**
  1. Read this hand-off + spec §12(a) (ephemeris-mode optimisation) + spec §14 V3 (TCM budget gate).
  2. **Add orbital-elements-based construction** to `verify/real_closure.py` (the M6b binding-gate unblock). Sibling to `construct_real_ephemeris_cycler`; dispatched at the `verify_real_closure` boundary on a catalogue-entry-shape criterion.
  3. Re-run the M6b binding gate (`test_aldrin_cycler_periodic_over_2_cycles_astropy`) with the orbital-elements construction. If it now passes, remove the xfail and the `_M6B_LAMBERT_CHAIN_XFAILS` set.
  4. Write the M7 plan, deciding: TCM application policy, per-cycle vs per-encounter, n_cycles horizon.
  5. Extend `RealClosureResult` with whatever shape M7 needs (additive only; do not reshape the locked fields).
  6. Wire `optimise_cell_ephemeris` against `verify_real_closure` as drift-feasibility check + new TCM-cost objective.
  7. Build catalogue YAML writer that consumes `RealClosureResult` and populates `validation.gates.V2.*` (always) and `validation.gates.V3.*` (only after TCM-budget runs).
  8. Run the full-catalogue batch verification (200+ entries) offline; identify the closes / open / multi-rev-blocked / no-window subsets.
  9. Implement multi-rev Lambert per `core/lambert.py` (the M-future task that unlocks `s1l1-2syn-em-cpom` and similar entries). This may be a parallel M7 sub-task or a dedicated M-future.
  10. Complete the catalogue leg data for the 2-syn McConaghy / Russell entries (currently in `EXPECTED_SKIPS` with reason "incomplete leg data"). Source: McConaghy 2006 JSR DOI 10.2514/1.15215 (paywalled, may need university access) and Russell 2004 dissertation Tables 3.4 + 4.7-4.9.

M7's first task is the M7 plan doc itself; this hand-off is the input it consumes. The orbital-elements construction (step 2-3) is the highest-leverage unblock: it converts the M6b binding gate from XFAIL to PASS at the existing 200,000 km tolerance.
