# M4 — Enumeration + scoring (todo)

Working checklist for the M4 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: Cell + unpruned enumerator → Tisserand gate → feasible/frontier → score → rank → re-exports → CI.

## Predecessor recap

- [ ] Re-read M3's `## Hand-off to M4` (in `docs/phases/m3-model-construct/todo.md`): noted Aldrin solver outputs (`a=1.6017 AU, e=0.3929, V∞_E=6.53, V∞_M=9.74`), 2-synodic E-M-E `closure_residual=8.21 km/s` (naive — M5 fixes), `synodic_omega("V")` raises NotImplementedError.
- [ ] Re-read M2 plan §3.2 + the M2 `linkable` signature: `linkable(body_a, body_b, vinf_kms, tol_au=0.01, tol_e=0.01, a_range_au=(0.3, 5.0), n_points=200) -> bool`. Confirmed coplanar-only.
- [ ] Re-read M3 `Cycler` API: `max_vinf()`, `radial_span()`, `encounters[i].vinf_in/out/r/v_planet`, `Cycler.maintenance_dv()` (which M4's `score()` deliberately does **not** call — `score()` uses `flyby_dv_for` instead).
- [ ] Re-read spec §13.1–§13.8 once more; confirm `Cell.id` format and the iterative-deepening contract.

## `search/sequence.py` — skeleton

- [ ] Create `src/cyclerfinder/search/sequence.py` with module docstring referencing spec §13.1–§13.3, §13.8.
- [ ] Define `Cell` frozen dataclass with the five fields per plan §3.1.1: `bodies`, `sequence`, `period_k`, `per_leg_revs`, `per_leg_branch`.
- [ ] Implement `Cell.id` as a `@property` returning the spec §13.8 deterministic string. Branch alphabet: `{"single":"s","low":"l","high":"h"}`.
- [ ] Stub `enumerate_cells`, `tisserand_feasible`, `feasible_cells`, `deepening_frontier` with `raise NotImplementedError` bodies and full docstrings per plan §3.1.2–§3.1.5.
- [ ] `uv run mypy src` clean on the skeleton (catches `Iterator[Cell]` typing issues before any logic).

## Cell tests

- [ ] Create `tests/search/__init__.py` (empty — first occupant under `tests/search/`).
- [ ] Create `tests/search/test_sequence.py` with the Cell dataclass tests:
  - [ ] `test_cell_is_frozen` — `FrozenInstanceError` on assignment.
  - [ ] `test_cell_id_format` — exact-string match against the spec §13.8 worked example `VEM|E-V-M-E-M-E|k3|r00101|bsslsl` (with the branch-alphabet mapping documented in plan §3.1.1).
  - [ ] `test_cell_hashable` — `{Cell(...): "x"}` works.
- [ ] `uv run pytest tests/search/test_sequence.py::test_cell_id_format` → green.

## Enumerate cells

- [ ] Implement `enumerate_cells(body_set, L_max, k_max, N_max, branch_set=("single",))`:
  - [ ] Loop `L in range(2, L_max+1)`, `k in range(1, k_max+1)`.
  - [ ] Sequence generation: `itertools.product(body_set, repeat=L)` then filter to adjacency-distinct (`s[i] != s[i+1]`).
  - [ ] Per-leg revs: `itertools.product(range(N_max+1), repeat=L-1)`.
  - [ ] Per-leg branches: for each `revs` tuple, choose `"single"` for 0-rev legs and cycle through `branch_set ∩ {"low","high"}` for ≥1-rev legs. If no valid branch exists, skip the rev tuple silently (documented).
  - [ ] Yield `Cell(bodies=body_set, sequence=tuple(seq), period_k=k, per_leg_revs=tuple(revs), per_leg_branch=tuple(branches))`.
- [ ] Enumeration tests:
  - [ ] `test_enumeration_count_em_l2_k1` — exact count **2** per plan §4.4.
  - [ ] `test_enumeration_count_em_l4_k2` (**gate**) — exact count **12** per plan §4.4.
  - [ ] `test_enumeration_excludes_consecutive_same_body`.
  - [ ] `test_enumeration_iterator_not_list` — `isinstance(it, Iterator) and not isinstance(it, list)`.
  - [ ] `test_cell_id_uniqueness` — all 12 ids in the L=4, k=2 enumeration are distinct.
- [ ] `uv run pytest tests/search/test_sequence.py -k enumeration` → green.

## Tisserand pruning gate

- [ ] Implement `tisserand_feasible(cell, vinf_cap, ephem=None)`:
  - [ ] If `vinf_cap <= 0`, return False.
  - [ ] If `len(cell.sequence) < 2`, return False.
  - [ ] Build 24 V∞ samples in `np.linspace(0.5, vinf_cap, 24)`.
  - [ ] For each consecutive pair `(cell.sequence[i], cell.sequence[i+1])`: return False if no V∞ in the grid has `tisserand.linkable(pair[0], pair[1], vinf) is True`.
  - [ ] Return True if every pair found a linkable V∞.
  - [ ] Wrap the whole body in a `try/except Exception` that returns False — never raise (mirrors `tisserand.linkable`'s contract).
- [ ] Tests:
  - [ ] `test_tisserand_pruning_rejects_low_vinf_em` (**gate**) — `(E,M)` cell at `vinf_cap=2.0` returns False.
  - [ ] `test_tisserand_pruning_accepts_aldrin_vinf_em` (**gate**) — `(E,M)` cell at `vinf_cap=8.0` returns True.
  - [ ] `test_tisserand_pruning_propagates_through_sequence` — 3-encounter cell with one infeasible pair returns False.
  - [ ] `test_tisserand_feasible_never_raises` — parametrised over `vinf_cap ∈ {0.0, -1.0, 1e6}`; all return a `bool`.
- [ ] `uv run pytest tests/search/test_sequence.py::test_tisserand_pruning_rejects_low_vinf_em` → green.
- [ ] `uv run pytest tests/search/test_sequence.py::test_tisserand_pruning_accepts_aldrin_vinf_em` → green.

## feasible_cells and deepening_frontier

- [ ] Implement `feasible_cells` as the documented one-liner generator.
- [ ] Implement `deepening_frontier`:
  - [ ] Track yielded cell ids in an in-memory `set`.
  - [ ] Outer loop over `(L, k, N)` raised by `(L_step, k_step, N_step)` each tier.
  - [ ] At each tier, call `feasible_cells` with the current caps and yield cells whose id is not in the set; add to the set.
  - [ ] Document in the docstring that this is **single-process, in-memory** dedup — replaced by M7's persistent ledger.
- [ ] Tests:
  - [ ] `test_feasible_cells_subset_of_enumerate` at `vinf_cap=8.0`.
  - [ ] `test_feasible_cells_strict_subset_at_low_cap` at `vinf_cap=1.0`.
  - [ ] `test_deepening_frontier_yields_in_increasing_complexity` — first 20 cells, sequence lengths non-decreasing within tier.
  - [ ] `test_deepening_frontier_no_repeats` — first 50 ids distinct.
- [ ] `uv run pytest tests/search/test_sequence.py` → green (all sequence tests pass).

## `model/score.py` — skeleton

- [ ] Create `src/cyclerfinder/model/score.py` with module docstring referencing spec §5 step 4 and §12(d).
- [ ] Define `Score` frozen dataclass with the six fields per plan §3.2.1.
- [ ] Module-level constant `DEFAULT_WEIGHTS` per plan §3.2.4.
- [ ] Stub `score`, `taxi_cost_kms`, `composite_score`, `rank` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src` clean on the skeleton.

## taxi_cost_kms

- [ ] Implement `taxi_cost_kms(cycler)`:
  - [ ] Collect `||enc.vinf_in||` for encounters with `enc.body == "E"`.
  - [ ] Return `max(...)` if non-empty, else `0.0`.
- [ ] Tests:
  - [ ] `test_taxi_cost_earth_only` — Earth-only cycler with two encounters of V∞ 5 and 7 → returns 7.
  - [ ] `test_taxi_cost_zero_when_no_earth` — V-M-only cycler → returns 0.
  - [ ] `test_taxi_cost_aldrin` — `build_aldrin_seed(ephem)` → ≈ 6.5 km/s (±0.1).

## score()

- [ ] Implement `score(cycler, ephem, vinf_cap, target_period_sec=None, rp_factors=None)`:
  - [ ] `total_maintenance_dv_kms = sum(flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out) for enc in cycler.encounters)`.
  - [ ] `max_vinf_kms = cycler.max_vinf()`.
  - [ ] `radial_span_au = cycler.radial_span()`.
  - [ ] `period_error_yr = abs(cycler.period - target_period_sec) / (DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY)` if `target_period_sec` given, else `0.0`.
  - [ ] `taxi_cost_kms = taxi_cost_kms(cycler)`.
  - [ ] `hard_constraints_pass`: True iff for every encounter, `max(||vinf_in||, ||vinf_out||) <= vinf_cap` AND `is_ballistic_feasible(vinf_in, vinf_out, mu_planet=PLANETS[enc.body].mu_km3_s2, rp_min=rp_factor*SAFE_PERIHELION_KM[enc.body])` returns True (using `rp_factors.get(enc.body, 1.0)` as the multiplier).
  - [ ] Return `Score(...)`.
- [ ] Tests:
  - [ ] `test_hard_constraints_pass_ballistic` — hand-built feasible cycler.
  - [ ] `test_hard_constraints_fail_on_vinf_cap` — V∞ exceeds cap.
  - [ ] `test_hard_constraints_fail_on_overbent_pair` — bend > max at Mars.
  - [ ] `test_score_aldrin_passes_hard_constraints` (**gate**) — Aldrin at `vinf_cap=12.0` passes, `max_vinf_kms ≈ 9.7 ± 0.1`.

## composite_score and rank

- [ ] Implement `composite_score(s, weights=None)`:
  - [ ] If `s.hard_constraints_pass is False`, return `math.inf`.
  - [ ] Else return `sum(weights.get(field, 0.0) * getattr(s, field) for field in DEFAULT_WEIGHTS)`.
- [ ] Implement `rank(cyclers, ephem, vinf_cap, n_keep=20, target_period_sec=None, weights=None)`:
  - [ ] Score each cycler.
  - [ ] Filter on `hard_constraints_pass`.
  - [ ] Sort ascending by `composite_score`.
  - [ ] Return first `n_keep` as `[(score, cycler), ...]`.
- [ ] Tests:
  - [ ] `test_composite_finite_and_reproducible` — bitwise-identical floats on repeat call.
  - [ ] `test_composite_infinite_on_hard_constraint_fail` — `math.inf` returned.
  - [ ] `test_rank_orders_low_dv_first`.
  - [ ] `test_rank_filters_infeasible`.
  - [ ] `test_rank_empty_input` — `rank([], ...) == []`.
  - [ ] `test_rank_n_keep_truncation`.

## Re-exports

- [ ] Edit `src/cyclerfinder/model/__init__.py` to re-export `Score` alongside the existing `Cycler`, `Leg`, `Encounter`. Confirm `from cyclerfinder.model import Score` resolves.

## Local green

- [ ] `uv run pytest` → green (no M0–M3 regressions; all new tests pass).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean (run `uv run ruff format .` first if needed).
- [ ] `uv run mypy src tests` → clean under strict mode. The `Iterator[Cell]` annotation and the `dict[str, float]` weight-bag are the most likely friction points — resolve in module, not via `# type: ignore`.

## CI

- [ ] Commit: `m4: cell enumeration with Tisserand pruning + per-cycler scoring`.
- [ ] Push; confirm GitHub Actions runs and all four checks pass.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M4 status `planned` → `completed`; M5 row `not yet planned` → `planned`.
- [ ] Append a `## Hand-off to M5` section to this `todo.md` covering:
  - The exact enumeration count produced for the `L=4, k=2, single-branch, (E,M)` case (12 expected per plan §4.4) and any deviation with rationale.
  - The observed Aldrin `Score` under default weights — `total_maintenance_dv_kms`, `max_vinf_kms`, `taxi_cost_kms`, `composite_score`. M5 will need this as a regression anchor when its optimiser starts wiring `construct → score → rank`.
  - Whether `deepening_frontier`'s in-memory dedup was tested at a meaningful cap and how many cells it covered before slowing down — informs the M7 ledger urgency.
  - Any tisserand.linkable surprises noticed during the 24-sample-grid sweep (false-negatives at low V∞? grid coarseness at high V∞?) — feedback for M2 if needed.
  - Confirmation that nothing about M3's `Cycler.maintenance_dv()` was changed; `score()`'s use of `flyby_dv_for` is purely additive.

## Hand-off to M5

_(To be filled in at the close of M4.)_
