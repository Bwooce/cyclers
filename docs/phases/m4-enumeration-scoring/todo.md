# M4 — Enumeration + scoring (todo)

Working checklist for the M4 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: Cell + unpruned enumerator → Tisserand gate → feasible/frontier → score → rank → re-exports → CI.

## Predecessor recap

- [x] Re-read M3's `## Hand-off to M4` (in `docs/phases/m3-model-construct/todo.md`): noted Aldrin solver outputs (`a=1.6017 AU, e=0.3929, V∞_E=6.53, V∞_M=9.74`), 2-synodic E-M-E `closure_residual=8.21 km/s` (naive — M5 fixes), `synodic_omega("V")` raises NotImplementedError.
- [x] Re-read M2 plan §3.2 + the M2 `linkable` signature: `linkable(body_a, body_b, vinf_kms, tol_au=0.01, tol_e=0.01, a_range_au=(0.3, 5.0), n_points=200) -> bool`. Confirmed coplanar-only.
- [x] Re-read M3 `Cycler` API: `max_vinf()`, `radial_span()`, `encounters[i].vinf_in/out/r/v_planet`, `Cycler.maintenance_dv()` (which M4's `score()` deliberately does **not** call — `score()` uses `flyby_dv_for` instead).
- [x] Re-read spec §13.1–§13.8 once more; confirm `Cell.id` format and the iterative-deepening contract.

## `search/sequence.py` — skeleton

- [x] Create `src/cyclerfinder/search/sequence.py` with module docstring referencing spec §13.1–§13.3, §13.8.
- [x] Define `Cell` frozen dataclass with the five fields per plan §3.1.1: `bodies`, `sequence`, `period_k`, `per_leg_revs`, `per_leg_branch`.
- [x] Implement `Cell.id` as a `@property` returning the spec §13.8 deterministic string. Branch alphabet: `{"single":"s","low":"l","high":"h"}`.
- [x] Stub `enumerate_cells`, `tisserand_feasible`, `feasible_cells`, `deepening_frontier` with `raise NotImplementedError` bodies and full docstrings per plan §3.1.2–§3.1.5. (Skipped the stub step in favour of implementing directly; mypy strict was verified on the full file.)
- [x] `uv run mypy src` clean on the skeleton (catches `Iterator[Cell]` typing issues before any logic).

## Cell tests

- [x] Create `tests/search/__init__.py` (already existed under M2's test layout).
- [x] Create `tests/search/test_sequence.py` with the Cell dataclass tests:
  - [x] `test_cell_is_frozen` — `FrozenInstanceError` on assignment.
  - [x] `test_cell_id_format` — exact-string match against the spec §13.8 worked example `VEM|E-V-M-E-M-E|k3|r00101|bsslsl` (the branch-alphabet mapping documented in plan §3.1.1). Also added `test_cell_id_format_all_low_branches` covering the spec's literal `blllll` example and `test_cell_id_format_em_2syn_direct` for the M4 native E-M-E case.
  - [x] `test_cell_hashable` — `{Cell(...): "x"}` works.
- [x] `uv run pytest tests/search/test_sequence.py::test_cell_id_format_spec_worked_example` → green.

## Enumerate cells

- [x] Implement `enumerate_cells(body_set, L_max, k_max, N_max, branch_set=("single",))`:
  - [x] Loop `L in range(2, L_max+1)`, `k in range(1, k_max+1)`.
  - [x] Sequence generation: `itertools.product(body_set, repeat=L)` then filter to adjacency-distinct (`s[i] != s[i+1]`).
  - [x] Per-leg revs: `itertools.product(range(N_max+1), repeat=L-1)`.
  - [x] Per-leg branches: for each `revs` tuple, choose `"single"` for 0-rev legs and cycle through `branch_set ∩ {"low","high"}` for ≥1-rev legs. If no valid branch exists, skip the rev tuple silently (documented).
  - [x] Yield `Cell(bodies=body_set, sequence=tuple(seq), period_k=k, per_leg_revs=tuple(revs), per_leg_branch=tuple(branches))`.
- [x] Enumeration tests:
  - [x] `test_enumeration_count_em_l2_k1` — exact count **2** per plan §4.4.
  - [x] `test_enumeration_count_em_l4_k2` (**gate**) — exact count **12** per plan §4.4.
  - [x] `test_enumeration_excludes_consecutive_same_body`.
  - [x] `test_enumeration_iterator_not_list` — `isinstance(it, Iterator) and not isinstance(it, list)`.
  - [x] `test_cell_id_uniqueness_in_em_l4_k2` — all 12 ids in the L=4, k=2 enumeration are distinct.
  - [x] `test_enumeration_multirev_requires_branch` — bonus: multi-rev rev tuples need a non-single branch in `branch_set` to yield.
- [x] `uv run pytest tests/search/test_sequence.py -k enumeration` → green.

## Tisserand pruning gate

- [x] Implement `tisserand_feasible(cell, vinf_cap, ephem=None)`:
  - [x] If `vinf_cap <= 0.5` (the sampling floor), return False.
  - [x] If `len(cell.sequence) < 2`, return False.
  - [x] Build 24 V∞ samples in `np.linspace(0.5, vinf_cap, 24)`.
  - [x] For each consecutive pair `(cell.sequence[i], cell.sequence[i+1])`: return False if no V∞ in the grid has `tisserand.linkable(pair[0], pair[1], vinf) is True`.
  - [x] Return True if every pair found a linkable V∞.
  - [x] Wrap the whole body in a `try/except Exception` that returns False — never raise (mirrors `tisserand.linkable`'s contract).
- [x] Tests:
  - [x] `test_tisserand_pruning_rejects_low_vinf_em_gate` (**gate**) — `(E,M)` cell at `vinf_cap=2.0` returns False.
  - [x] `test_tisserand_pruning_accepts_aldrin_vinf_em_gate` (**gate**) — `(E,M)` cell at `vinf_cap=8.0` returns True.
  - [x] `test_tisserand_pruning_propagates_through_sequence` — 3-encounter cell with infeasible pair at `vinf_cap=1.0` returns False.
  - [x] `test_tisserand_feasible_never_raises` — parametrised over `vinf_cap ∈ {0.0, -1.0, 1.0e6}`; all return a `bool`.
- [x] `uv run pytest tests/search/test_sequence.py::test_tisserand_pruning_rejects_low_vinf_em_gate` → green.
- [x] `uv run pytest tests/search/test_sequence.py::test_tisserand_pruning_accepts_aldrin_vinf_em_gate` → green.

## feasible_cells and deepening_frontier

- [x] Implement `feasible_cells` as the documented one-liner generator.
- [x] Implement `deepening_frontier`:
  - [x] Track yielded cell ids in an in-memory `set`.
  - [x] Outer loop over `(L, k, N)` raised by `(L_step, k_step, N_step)` each tier.
  - [x] At each tier, call `feasible_cells` with the current caps and yield cells whose id is not in the set; add to the set.
  - [x] Document in the docstring that this is **single-process, in-memory** dedup — replaced by M7's persistent ledger.
  - [x] Add `max_tiers` kwarg so callers (and tests) have a clean stop without `itertools.islice` indirection.
- [x] Tests:
  - [x] `test_feasible_cells_subset_of_enumerate` at `vinf_cap=8.0`.
  - [x] `test_feasible_cells_strict_subset_at_low_cap` at `vinf_cap=1.0`.
  - [x] `test_deepening_frontier_yields_in_increasing_complexity` — first 20 cells, sequence lengths non-decreasing (with `k`/`N` pinned via large steps so only `L` rises across tiers; see test docstring for the rationale on why strict length-monotonicity holds only under that restriction).
  - [x] `test_deepening_frontier_no_repeats` — first 50 ids distinct.
  - [x] `test_deepening_frontier_step_validation` — all-zero steps raise `ValueError`.
- [x] `uv run pytest tests/search/test_sequence.py` → green (all 22 sequence tests pass).

## `model/score.py` — skeleton

- [x] Create `src/cyclerfinder/model/score.py` with module docstring referencing spec §5 step 4 and §12(d).
- [x] Define `Score` frozen dataclass with the six fields per plan §3.2.1.
- [x] Module-level constant `DEFAULT_WEIGHTS` per plan §3.2.4.
- [x] Stub `score`, `taxi_cost_kms`, `composite_score`, `rank` with full docstrings + `NotImplementedError` bodies. (Implemented directly; mypy strict was verified on the full file.)
- [x] `uv run mypy src` clean on the skeleton.

## taxi_cost_kms

- [x] Implement `taxi_cost_kms(cycler)`:
  - [x] Collect `||enc.vinf_in||` for encounters with `enc.body == "E"`.
  - [x] Return `max(...)` if non-empty, else `0.0`.
- [x] Tests:
  - [x] `test_taxi_cost_earth_only` — Earth-only cycler with two encounters of V∞ 5 and 7 → returns 7.
  - [x] `test_taxi_cost_zero_when_no_earth` — V-M-only cycler → returns 0.
  - [x] `test_taxi_cost_aldrin` — `build_aldrin_seed(ephem)` → ≈ 6.5 km/s (±0.1).

## score()

- [x] Implement `score(cycler, ephem, vinf_cap, target_period_sec=None, rp_factors=None)`:
  - [x] `total_maintenance_dv_kms = sum(flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out) for enc in cycler.encounters)`.
  - [x] `max_vinf_kms = cycler.max_vinf()`.
  - [x] `radial_span_au = cycler.radial_span()`.
  - [x] `period_error_yr = abs(cycler.period - target_period_sec) / (DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY)` if `target_period_sec` given, else `0.0`.
  - [x] `taxi_cost_kms = taxi_cost_kms(cycler)`.
  - [x] `hard_constraints_pass`: True iff for every encounter, `max(||vinf_in||, ||vinf_out||) <= vinf_cap` AND `is_ballistic_feasible(vinf_in, vinf_out, mu_planet=PLANETS[enc.body].mu_km3_s2, rp_min=rp_factor*SAFE_PERIHELION_KM[enc.body])` returns True (using `rp_factors.get(enc.body, 1.0)` as the multiplier).
  - [x] Return `Score(...)`.
- [x] Tests:
  - [x] `test_hard_constraints_pass_ballistic` — hand-built feasible cycler.
  - [x] `test_hard_constraints_fail_on_vinf_cap` — V∞ exceeds cap.
  - [x] `test_hard_constraints_fail_on_overbent_pair` — bend > max at Mars (180° at 8 km/s).
  - [x] `test_score_aldrin_passes_hard_constraints_gate` (**gate**) — Aldrin at `vinf_cap=12.0` passes, `max_vinf_kms ≈ 9.7 ± 0.1`.

## composite_score and rank

- [x] Implement `composite_score(s, weights=None)`:
  - [x] If `s.hard_constraints_pass is False`, return `math.inf`.
  - [x] Else return `sum(weights.get(field, 0.0) * getattr(s, field) for field in DEFAULT_WEIGHTS)`.
- [x] Implement `rank(cyclers, ephem, vinf_cap, n_keep=20, target_period_sec=None, weights=None)`:
  - [x] Score each cycler.
  - [x] Filter on `hard_constraints_pass`.
  - [x] Sort ascending by `composite_score`.
  - [x] Return first `n_keep` as `[(score, cycler), ...]`.
- [x] Tests:
  - [x] `test_composite_finite_and_reproducible` — bitwise-identical floats on repeat call.
  - [x] `test_composite_infinite_on_hard_constraint_fail` — `math.inf` returned.
  - [x] `test_rank_orders_low_dv_first` — implemented via `target_period_sec` axis (period_error is monotone in the input ordering).
  - [x] `test_rank_filters_infeasible`.
  - [x] `test_rank_empty_input` — `rank([], ...) == []`.
  - [x] `test_rank_n_keep_truncation`.

## Re-exports

- [x] Edit `src/cyclerfinder/model/__init__.py` to re-export `Score` alongside the existing `Cycler`, `Leg`, `Encounter`. Confirm `from cyclerfinder.model import Score` resolves.

## Local green

- [x] `uv run pytest` → green (153 tests pass; no M0–M3 regressions).
- [x] `uv run ruff check .` → clean.
- [x] `uv run ruff format --check .` → clean.
- [x] `uv run mypy src tests` → clean under strict mode.

## CI

- [x] Commit: `m4: cell enumeration + Tisserand pruning + scoring + ranking`.
- [x] Push; confirm GitHub Actions runs and all four checks pass.

## Closeout

- [x] Update `docs/overview.md` §4 milestone table: M4 status → `completed`; M5 row → `planned`.
- [x] Append a `## Hand-off to M5` section to this `todo.md` (below).

## Hand-off to M5

**Enumeration count gate (plan §4.4):** `enumerate_cells(("E","M"), l_max=4, k_max=2, n_max=0, branch_set=("single",))` yields **exactly 12** cells, matching the §4.4 derivation. No deviation.

**Aldrin Score (regression anchor for M5):**

Built via `build_aldrin_seed(Ephemeris(model="circular"))`, scored with `vinf_cap=12.0`, default weights:

| Field | Value |
|---|---|
| `total_maintenance_dv_kms` | `0.000000` km/s — M3 boundary convention forces `vinf_in == vinf_out` at the open endpoints, so `flyby_dv_for` is zero at both encounters by construction |
| `max_vinf_kms`             | `9.743359` km/s — matches the M3 hand-off (V∞_M ≈ 9.7 km/s) |
| `radial_span_au`           | `(0.9724, 2.2310)` — matches the Aldrin literature anchors (perihelion 0.97 AU, aphelion 2.23 AU) |
| `period_error_yr`          | `0.000000` — no target supplied |
| `taxi_cost_kms`            | `6.530070` km/s — matches the M3 hand-off (V∞_E ≈ 6.5 km/s) |
| `hard_constraints_pass`    | `True` |
| `composite_score(DEFAULT_WEIGHTS)` | `4.239371` — purely ordinal; use as a regression anchor when M5's optimiser starts running cells through `construct → score → rank` |

M5 should be able to reproduce these to ≤ 1e-3 relative tolerance if it leaves the scoring layer untouched.

**`deepening_frontier` in-memory dedup behaviour at the test cap (vinf_cap=8.0, body_set=("E","M"), l_initial=2, k_initial=1, n_initial=0, max_tiers=10):** the test materialises 50 cells with distinct ids in well under a second; the dedup `set` never exceeded ~50 entries in practice during M4's test pass. The persistent-ledger urgency for M7 is **not** elevated — the in-memory set is adequate for single-process exploration at M4/M5 scales (E-M only). M7 will need to revisit once VEM enumeration at higher caps lands in M8.

**`tisserand.linkable` surprises during the 24-sample-grid sweep:** none — the grid resolution is comfortably above the linkable bands' widths in the Aldrin neighbourhood. At `vinf_cap=2.0` the E-M pair correctly fails at every sample; at `vinf_cap=8.0` the pair succeeds on the first sample inside the band. No false-negative regressions, no need for a denser grid in M4. M2's coplanar-only restriction is inherited and re-stated in the module docstring; M6 will add the ephemeris-aware (3-D) variant when astropy lands.

**M3 `Cycler.maintenance_dv()` untouched.** `score()` uses `flyby_dv_for` for the bend-and-magnitude decomposition; M3's naive velocity-discontinuity sum is unchanged and its existing tests continue to consume it. Both methods exist side by side; documented in the M4 plan §3.2 risk table (item 7).

**API surfaces M5 will consume immediately:**
- `from cyclerfinder.search.sequence import Cell, enumerate_cells, feasible_cells, tisserand_feasible, deepening_frontier`
- `from cyclerfinder.model import Score`
- `from cyclerfinder.model.score import score, composite_score, rank, taxi_cost_kms, DEFAULT_WEIGHTS`

M5's optimiser is the first thing that wires the enumerator to the constructor and through the score layer: take a feasible `Cell`, search the timing DOF inside it, build a `Cycler` via `construct_cycler`, `score()` it, and feed many of them to `rank()`. The M4 contract is that those calls compose without further glue.
