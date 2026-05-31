# M5 — Optimisation (todo)

Working checklist for the M5 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: predecessor recap → skeleton → helpers (each with paired unit tests) → idealized mode + gate tests → ephemeris stub → top-level `find_cyclers` → local green → CI → closeout.

## Predecessor recap

- [ ] Re-read M4's `## Hand-off to M5` (in `docs/phases/m4-enumeration-scoring/todo.md`): noted Aldrin regression anchor (`composite_score = 4.239371`, `max_vinf_kms = 9.743359`, `taxi_cost_kms = 6.530070`); noted M5's primary consumed surfaces (`Cell`, `feasible_cells`, `Score`, `score()`, `composite_score`, `rank`); noted M3's open-sequence boundary convention (`vinf_in == vinf_out` at first/last encounters by construction).
- [ ] Re-read spec §12(a) (two optimisation modes) and §12(d) (hard inequality constraints) once more — these are the binding architecture decisions for M5.
- [ ] Re-read spec §13.4 (structured inner search — free-return seed + multi-start grid + local polish) — the single most important paragraph in M5's design.
- [ ] Confirm `construct_cycler(sequence, encounter_times_sec, ephem, ...)` signature in `search/construct.py` — encounter epoch is `seconds from t=0`.
- [ ] Confirm `synodic_period_days(body_a, body_b)` in `search/resonance.py` and `SECONDS_PER_DAY`, `DAYS_PER_JULIAN_YEAR` in `core/constants.py` for the `target_period_sec` derivation.
- [ ] Confirm `Cycler.closure_residual(omega_rad_per_s)` API in `model/cycler.py` — M5 calls this with `omega = 2*pi/T` to evaluate the idealized-mode residual.
- [ ] Confirm `is_ballistic_feasible` is a bool and **cannot** be used directly as a smooth SLSQP constraint (this is why plan §3.3 derives the smooth `_r_p_required` helper).
- [ ] Confirm `scipy>=1.13` is in `pyproject.toml` — no dependency additions needed.

## `search/optimize.py` — skeleton

- [ ] Create `src/cyclerfinder/search/optimize.py` with module docstring referencing spec §12(a), §12(d), §13.4, §9 (degenerate-solution guard).
- [ ] Define `_StartRecord` frozen dataclass (module-private) per plan §3.4: `start_index`, `x0`, `x_final`, `objective_value`, `constraints_satisfied`, `nit`, `success`.
- [ ] Define `OptimisationResult` frozen dataclass per plan §3.4: `cell`, `best_cycler`, `best_score`, `closure_residual_kms`, `optimiser_history`, `converged`, `constraints_satisfied`.
- [ ] Stub `optimise_cell_idealized`, `optimise_cell_ephemeris`, `find_cyclers` with full docstrings + `raise NotImplementedError(...)` bodies.
- [ ] `uv run mypy src` clean on the skeleton (catches `tuple[_StartRecord, ...]` typing, scipy `OptimizeResult` interop, `Iterator[Cell]` consumption).

## Helpers (one at a time, paired tests)

### `_target_period_sec`

- [ ] Implement: `synodic_period_days(cell.bodies[0], cell.bodies[1]) * cell.period_k * SECONDS_PER_DAY`. Document M8 multi-body follow-up (VEM beat from `resonance.find_beats(...)`).
- [ ] Test `test_target_period_sec_em_2syn` — for the 2-syn E-M cell, returns `2 * 2.135 yr * SECONDS_PER_DAY * DAYS_PER_JULIAN_YEAR` within 0.1%.

### `_free_return_seed`

- [ ] Implement: equispaced interior epochs `t_i = i * T / (N-1)` for `i = 0 … N-1`. Returns `tuple[float, ...]` of length `N` (including the pinned endpoints).
- [ ] Test `test_free_return_seed_em_2syn` — returns exactly `(0.0, T/2, T)` for the 2-syn E-M-E cell (`N=3`).
- [ ] Test `test_free_return_seed_vem_3syn` — returns `(0, T/5, 2T/5, 3T/5, 4T/5, T)` for a hypothetical 6-encounter VEM cell (`N=6`). (Unit test on the helper; not exercised by M5's gate.)

### `_multi_start_grid`

- [ ] Implement: deterministic perturbation table; `numpy.random.default_rng(seed).permutation(table)` to pick `n_starts` entries. Start 0 is always the free-return seed exactly.
- [ ] Test `test_multi_start_grid_deterministic` — two calls with `seed=42` produce bitwise-identical lists.
- [ ] Test `test_multi_start_grid_distinct` — `n_starts=5` returns 5 vectors with `len(set(map(tuple, starts))) == 5`.
- [ ] Test `test_multi_start_grid_includes_free_return` — `starts[0]` equals `_free_return_seed(cell, ephem, T)`.

### `_build_cycler_from_x`

- [ ] Implement: unpack `x` into the full encounter-time list (`[0, x_0, x_1, …, T]`), call `construct_cycler(cell.sequence, times, ephem, max_revs_per_leg=list(cell.per_leg_revs), branch_per_leg=list(cell.per_leg_branch))`.
- [ ] Wrap call in `try/except ValueError`; on failure return `None` (the `_objective` / `_constraints` callers handle `None`).
- [ ] Test `test_build_cycler_from_x_aldrin` — hand-built `x` for the Aldrin cell reproduces `build_aldrin_seed(ephem)`'s `max_vinf()` and `radial_span()` within float tolerance.
- [ ] Test `test_build_cycler_from_x_returns_none_on_pathology` — `x` with non-monotonic times → returns `None`, does not raise.

### `_r_p_required`

- [ ] Implement: given `vin, vout, mu`, compute `delta = bend_angle(vin, vout)` and `V_inf = mean(||vin||, ||vout||)`, then `r_p = (1/sin(delta/2) - 1) * mu / V_inf**2`. Document derivation from the `max_bend` formula.
- [ ] Test `test_r_p_required_known_aldrin_flyby` — Aldrin Earth flyby's (vin, vout) → `_r_p_required` exceeds `SAFE_PERIHELION_KM['E']`.
- [ ] Test `test_r_p_required_overbent_returns_small` — synthetic (vin, vout) with 178° angle → `_r_p_required` near 0 (effectively infeasible).

### `_objective`

- [ ] Implement: build cycler via `_build_cycler_from_x`; if `None`, return `1e6`. Else return `closure_residual(omega=2*pi/T) + sum(flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out) for enc in cycler.encounters)`.
- [ ] Test `test_objective_zero_at_closed_aldrin` — `x` for Aldrin's known geometry → objective is a small finite number (< 20 km/s composite — Aldrin is naively closed, not optimised).
- [ ] Test `test_objective_never_raises` — parametrised over pathological `x` (NaN, non-monotonic, all-zero) → returns a finite float, never raises.

### `_constraints`

- [ ] Implement: return a `list[dict]` of `{type: 'ineq', fun: callable}` per encounter — two constraints (V∞ cap, r_p floor) per encounter. Each `fun(x)` builds the cycler internally and returns the slack value. On `None` from `_build_cycler_from_x`, return `-1.0` (constraint violated).
- [ ] Test `test_constraints_ineq_form` — `_constraints(x_aldrin, …)` returns a list of dicts; each `dict["fun"](x_aldrin)` is positive (Aldrin satisfies V∞ ≤ 12 cap and r_p ≥ rp_min).
- [ ] Test `test_constraints_reject_at_cap` — hand-built `x` producing V∞ ≈ 12 km/s with `vinf_cap=7.0` → V∞-cap constraint's `fun` returns negative.
- [ ] Test `test_constraints_per_encounter` — `len(constraints) == 2 * N` for an N-encounter cell.

### `_polish`

- [ ] Implement: call `scipy.optimize.minimize(_objective, x0, args=(...), method='SLSQP', constraints=_constraints(...), options={'maxiter': 200, 'ftol': 1e-6})`. Wrap return in a `_StartRecord`.
- [ ] Evaluate `constraints_satisfied` at `x_final` by re-checking each constraint's `fun(x_final) >= -1e-6` (SLSQP allows small numerical violations at the boundary).
- [ ] Test `test_polish_returns_start_record` — `_polish(_free_return_seed(...), ...)` returns a `_StartRecord` with all six fields populated; `nit > 0`.
- [ ] Test `test_polish_respects_constraints` — starting from a feasible point, `_polish` never returns an `x_final` whose V∞ > `vinf_cap + 0.01` (small numerical slack allowed).

### `_de_pass`

- [ ] Implement: `scipy.optimize.differential_evolution(_objective, bounds, args=(...), constraints=scipy.optimize.NonlinearConstraint(...), seed=seed, polish=False, maxiter=50, popsize=8, tol=1e-4)`. Bounds: each free `t_i` ∈ `(0.01*T, 0.99*T)` strictly inside the period. Return `_StartRecord(start_index=-1, x0=x_initial_de_population_mean, x_final=result.x, ...)`.
- [ ] Test `test_de_pass_obeys_seed` — two calls with `seed=7` → bitwise-identical `x_final`.
- [ ] Test `test_de_pass_returns_in_bounds` — `x_final` strictly inside the bounds box.

## `optimise_cell_idealized` — compose helpers

- [ ] Implement per plan §3.2 pseudocode:
  - [ ] Resolve `target_period_sec` via `_target_period_sec(cell)` if not supplied.
  - [ ] Build `n_starts` initial points via `_multi_start_grid`.
  - [ ] `_polish` each start → list of `_StartRecord`.
  - [ ] If `use_de=True`, also `_polish(_de_pass(...))` and append to records.
  - [ ] Pick best via `_composite_with_constraints` helper (infeasibles → `+inf`).
  - [ ] Build final cycler via `_build_cycler_from_x(best.x_final, ...)`; if `None`, the optimisation has failed — return an `OptimisationResult` with `converged=False, constraints_satisfied=False` and a sentinel `best_cycler` (the seed cycler from start 0).
  - [ ] Score the final cycler with `score(...)`.
  - [ ] Return `OptimisationResult(...)`.
- [ ] Gate test `test_2syn_em_rediscovers_5_65_kms_earth` (**M5 BINDING GATE**) — cell `("E","M","E"), period_k=2, per_leg_revs=(0,0), per_leg_branch=("single","single")`, `vinf_cap=7.0`, `seed=0` → result has Earth V∞ ≈ 5.65 km/s (±0.2), Mars V∞ ≈ 3.05 km/s (±0.2), `constraints_satisfied=True`, `closure_residual_kms < 0.5`.
- [ ] Gate test `test_2syn_em_rejects_high_vinf_degenerate` (**M5 BINDING GATE**) — feed deliberately bad initial guess (`t_1 = 0.05 * T`) via a custom override (test-only `_seed_override` kwarg, or use `n_starts=1, seed=…` chosen to land at that point); assert the per-start constraint check rejects high-V∞ "closures" AND the overall result picks a feasible alternative start.
- [ ] Test `test_aldrin_regression_anchor` — Aldrin cell at `vinf_cap=12.0`; `composite_score(result.best_score) ≤ 4.239371 * (1 + 1e-3)`; `max_vinf_kms ≈ 9.74 ±0.1`; `taxi_cost_kms ≈ 6.53 ±0.1`.
- [ ] Test `test_optimisation_result_frozen_and_seeded` — `result.best_cycler = …` raises `FrozenInstanceError`; two runs with `seed=42` produce bitwise-identical `result.best_score.max_vinf_kms`.
- [ ] Test `test_optimisation_result_carries_history` — `len(result.optimiser_history) == n_starts` (or `n_starts + 1` if `use_de=True`); every entry is a `_StartRecord`.
- [ ] `uv run pytest tests/search/test_optimize.py::test_2syn_em_rediscovers_5_65_kms_earth` → green. **If this fails, escalate per plan §5 risk #1** (raise n_starts, tighten grid, add trust-constr fallback, re-examine the free-return anchor — do NOT soften the V∞ tolerance).

## `optimise_cell_ephemeris` — stub

- [ ] Implement: body is `raise NotImplementedError("requires M6 ephemeris backend")`. Signature locked per plan §3.1 so M6 fills only the body.
- [ ] Test `test_ephemeris_mode_stubbed_until_m6` — calling the function raises `NotImplementedError`; the message contains the substring `"M6 ephemeris"`.

## `find_cyclers` — compose M4 + M5 + M4 rank

- [ ] Implement per plan §3.5 pseudocode:
  - [ ] Default `ephem = Ephemeris("circular")` if `None`.
  - [ ] `cells = list(feasible_cells(bodies, L_max, k_synodic, N_max, vinf_cap, ephem, branch_set))` — materialised (we iterate twice).
  - [ ] `results = [optimise_cell_idealized(c, ephem, vinf_cap=vinf_cap, n_starts=n_starts, seed=seed, use_de=use_de, rp_factors=rp_factors) for c in cells]`.
  - [ ] Filter: `[r for r in results if r.constraints_satisfied and r.best_score.hard_constraints_pass]`.
  - [ ] Sort ascending by `composite_score(r.best_score)`.
  - [ ] Return first `n_keep`.
- [ ] Gate test `test_find_cyclers_em_top_level` (**M5 BINDING GATE**) — `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, n_keep=5, seed=0)`; `len(results) >= 1`; the top result's Earth/Mars V∞ within ±0.2 km/s of the 5.65/3.05 published anchors.
- [ ] Test `test_find_cyclers_empty_when_caps_too_low` — `vinf_cap=1.0` → empty list.
- [ ] Test `test_find_cyclers_n_keep_truncation` — `n_keep=2` → at most 2 results.
- [ ] Test `test_find_cyclers_results_sorted` — output sorted ascending by `composite_score(r.best_score)`.
- [ ] Test `test_find_cyclers_all_results_feasible` — every result has `r.constraints_satisfied=True` and `r.best_score.hard_constraints_pass=True`.

## Optional re-exports

- [ ] Edit `src/cyclerfinder/search/__init__.py` to re-export the M5 public names if helpful. Decision: skip unless test verbosity warrants it (tests use the fully-qualified import per M4 convention).

## Local green

- [ ] `uv run pytest` → green (M0–M5 all passing; M4 Aldrin anchor still produces `composite_score = 4.239371`).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean.
- [ ] `uv run mypy src tests` → clean under strict mode (scipy `OptimizeResult` is `Any` at the boundary — explicit cast or `# type: ignore[no-any-return]` documented where used).

## CI

- [ ] Commit: `m5: per-cell inner-timing optimiser + find_cyclers pipeline (rediscovers 2-syn E-M)`.
- [ ] Push; confirm GitHub Actions runs and all four checks pass.

## Closeout

- [ ] Update `docs/overview.md` §2 deferred-decisions table: move "Global optimiser — scipy DE vs pygmo (decided in M5)" out of the deferred list into the main decisions table with the chosen value (`scipy.optimize.differential_evolution`; pygmo deferred to stretch) and rationale (already a dep; the §13.4 structured grid carries most coverage; pygmo's archipelago belongs at the M7 work-queue layer, not the inner optimiser).
- [ ] Update `docs/overview.md` §4 milestone table: M5 status → `completed`; M6a row → `planned`.
- [ ] Append a `## Hand-off to M6a` section to this `todo.md` (below).

## Hand-off to M6a

*(Filled in at M5 closeout. Placeholder structure below; replace bullets with measured values from the M5 test pass.)*

- **Actual V∞ values reproduced for the 2-syn E-M-E gate:** Earth V∞ = `<TBD>` km/s (vs published 5.65); Mars V∞ = `<TBD>` km/s (vs published 3.05); the gap is the circular-coplanar → real-ephemeris bridge that M6a/M6b absorbs.
- **Closure residual at convergence for the gate cell:** `<TBD>` km/s. This is the budget M6a inherits for its "verify periodic over ≥3 laps" gate — anything M6a's `verify_long_term_stability` reports as drift must be consistent with this idealised residual.
- **Empirical SLSQP iteration count (`nit`) for the 2-syn case:** `<TBD>` per start; total convergence cost `<TBD>` seconds. **DE generation count (if `use_de=True`):** `<TBD>`. Informs M8's compute-budget default for the VEM campaign.
- **Cells in `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, L_max=4, N_max=0)` that failed to converge:** `<TBD>`. Each failed cell is logged with `(cell.id, reason)` for M6a's verification ordering — known-bad cells deprioritised at verify time.
- **`_build_cycler_from_x` try/except trap fires:** `<count TBD>` during the test pass. Suggests robustness of the parameter parameterisation for M6a's ephemeris mode (a higher count means tighter bounds are needed in M6).
- **`use_de=True` vs `use_de=False` quality comparison:** did the multi-start grid + SLSQP alone hit the V∞ anchors, or was the DE wrapper needed? `<TBD>`. Informs whether M8 can run VEM cells with `use_de=False` to keep per-cell budget low.
- **Spec ambiguities resolved during M5 implementation:**
  - `target_period_sec` for VEM cells (`len(bodies) >= 3`) — M5 punts to `synodic_period_days(bodies[0], bodies[1]) * cell.period_k`, which is single-pair-correct. M8's VEM campaign will need a multi-body beat dispatch (`resonance.find_beats`). M5 documents this in `_target_period_sec`'s docstring; no M6 action needed.
  - Per-cell `_StartRecord` history is module-private (`_` prefix). M7's ledger may want a public diagnostic accessor; document the intent at M7's catalogue schema design.
- **API surfaces M6a will consume immediately:**
  - `from cyclerfinder.search.optimize import OptimisationResult, optimise_cell_idealized, optimise_cell_ephemeris`
  - `from cyclerfinder.search.optimize import find_cyclers` (top-level)
  - `optimise_cell_ephemeris` is M6a/M6b's primary deliverable — the signature is locked; M6 fills the body.

M6a's first task is the `verify/propagate.py` multi-lap propagator, against which M5's gate cyclers are checked for bounded drift in the dynamic rotating frame (spec §12(c)). The M5 contract is: any `OptimisationResult` with `constraints_satisfied=True ∧ converged=True ∧ closure_residual_kms < 0.5` is a candidate for M6a verification.
