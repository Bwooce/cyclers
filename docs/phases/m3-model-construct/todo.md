# M3 — Model + construction (todo)

Working checklist for the M3 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Frame transforms (land first — the spec §10 risk gate)

- [ ] Create `src/cyclerfinder/core/frames.py` with module docstring that states the uniform-frame scope and **explicitly notes the dynamic frame (spec §12(c)) is deferred to M6**.
- [ ] Implement `to_rotating(r_inertial, v_inertial, t, omega_rad_per_s) -> tuple[NDArray, NDArray]`.
- [ ] Implement `from_rotating(r_rot, v_rot, t, omega_rad_per_s) -> tuple[NDArray, NDArray]` as the algebraic inverse, not by computing a matrix inverse.
- [ ] Implement `synodic_omega(body: str) -> float`:
  - [ ] Returns Earth's mean motion (rad/s) for `"E"` and `"M"`, derived from `constants.PLANETS["E"].mean_motion_deg_day`.
  - [ ] Raises `NotImplementedError` for `"V"` (Venus-anchored frames are M8 work).
- [ ] Create `tests/core/test_frames.py` with named tolerances `TOL_ROUND_TRIP_KM = 1e-10` and `TOL_ROUND_TRIP_KMS = 1e-13` at module top.
- [ ] `test_round_trip_identity`: random (r, v) at 1 AU scale, various t and omega, assert `(from_rotating ∘ to_rotating) = identity` within tolerance. **This is the risk-mitigation test for spec §10.**
- [ ] `test_omega_zero_is_identity`: at ω = 0, `to_rotating` is the identity for any t.
- [ ] `test_circular_orbit_is_stationary_in_its_own_frame`: a circular orbit at radius a, transformed at its mean motion, has `‖v_rot‖ < 1e-6 km/s`. (Catches the sign-error bug the round-trip test would miss.)
- [ ] `test_synodic_omega_earth_matches_constants`: agrees with `PLANETS["E"].mean_motion_deg_day` conversion to within 1e-12 rel.
- [ ] `test_synodic_omega_venus_raises`: `pytest.raises(NotImplementedError)`.
- [ ] `uv run pytest tests/core/test_frames.py` green.

## Cycler dataclasses (model layer)

- [ ] Create `src/cyclerfinder/model/__init__.py` re-exporting `Cycler`, `Leg`, `Encounter`.
- [ ] Create `src/cyclerfinder/model/cycler.py` with module docstring describing the closure convention (residual = velocity mismatch after rotating the endpoint state by ω·T).
- [ ] Define `Vec3 = NDArray[np.float64]` alias at module top.
- [ ] Implement `@dataclass(frozen=True) class Encounter` with the 6 fields from plan §3.1.
- [ ] Implement `@dataclass(frozen=True) class Leg` with the 8 fields, defaults `n_revs=0`, `branch="single"`.
- [ ] Implement `@dataclass(frozen=True) class Cycler` with 4 fields and 4 methods.
- [ ] `Cycler.maintenance_dv()`: `sum(‖vinf_out - vinf_in‖)` over encounters; **no `flyby_dv` call in M3** (the bend-feasibility version arrives in M4).
- [ ] `Cycler.closure_residual(omega=None)`: defaults to `synodic_omega("E")` from `frames.py`; returns the velocity residual `‖v_first - rot(v_last, ω·T)‖` in km/s.
- [ ] `Cycler.radial_span()`: returns `(min_perihelion_AU, max_aphelion_AU)` computed analytically from each leg's `(r, v_depart)`.
- [ ] `Cycler.max_vinf()`: `max(‖enc.vinf_in‖)`; uses `vinf_in` only, not both, to avoid double-counting in a steady cycler.
- [ ] Create `tests/model/__init__.py` (empty).
- [ ] Create `tests/model/test_cycler.py` with the six tests from plan §4.2.
- [ ] `test_encounter_is_frozen`: assigning to a field raises `FrozenInstanceError`.
- [ ] `test_leg_defaults`: omitted `n_revs` and `branch` use the documented defaults.
- [ ] `test_cycler_max_vinf_returns_largest_magnitude`.
- [ ] `test_cycler_maintenance_dv_zero_when_vinf_in_equals_vinf_out` (within `1e-12 km/s`).
- [ ] `test_cycler_radial_span_two_body_circular`: hand-crafted Hohmann transfer; assert departure/arrival radii within `1e-6 AU`.
- [ ] `test_closure_residual_uses_default_earth_omega`: `omega=None` matches `omega=synodic_omega("E")`.
- [ ] `uv run pytest tests/model/test_cycler.py` green.

## Constructor + Aldrin gate

- [ ] Create `src/cyclerfinder/search/construct.py` with module docstring stating "no optimisation in M3 — encounter times are inputs."
- [ ] Implement `construct_cycler(sequence, encounter_times_sec, ephem, mu_sun=MU_SUN_KM3_S2, max_revs_per_leg=None, branch_per_leg=None) -> Cycler` per plan §3.3.
  - [ ] Validate: `len(sequence) == len(encounter_times_sec) >= 2`; strictly increasing times; every body in `constants.PLANETS`.
  - [ ] Default `max_revs_per_leg = [0] * (n-1)`, `branch_per_leg = ["single"] * (n-1)`.
  - [ ] Look up each planet's `(r, v_planet)` from `ephem.state(...)`.
  - [ ] Lambert-solve each leg; select the solution matching the requested `(n_revs, branch)`; raise `ValueError` if no match.
  - [ ] Build `Encounter` objects with the boundary-encounter convention `vinf_in = vinf_out` at the first/last node of an open sequence.
  - [ ] Build `Leg` objects with `v_depart`/`v_arrive` from the Lambert solution.
  - [ ] `period = encounter_times_sec[-1] - encounter_times_sec[0]`.
- [ ] Implement `build_aldrin_seed(ephem, t_start_sec=0.0) -> Cycler`:
  - [ ] Calls `construct_cycler(["E", "M"], [t_start_sec, t_start_sec + 146*SECONDS_PER_DAY], ephem)`.
  - [ ] Returns the result directly.
- [ ] Create `tests/model/test_aldrin.py` with named tolerances at module top:
  - [ ] `TOL_A_AU = 0.01`, `TOL_E = 0.02`, `TOL_PERI_AU = 0.02`, `TOL_APO_AU = 0.02`, `TOL_TOF_DAYS = 1.0`.
- [ ] Write `orbital_elements_from_rv(r, v, mu) -> (a_AU, e)` as a test-local helper (vis-viva + eccentricity vector; ~10 lines; not promoted to module).
- [ ] `test_aldrin_E_to_M_leg_elements`: build seed, extract `(a, e, perihelion, aphelion, tof)`, assert all five within tolerance of spec §9 anchors (1.659 AU, 0.41, 0.98 AU, 2.34 AU, 146 d).
- [ ] `test_aldrin_closure_residual_small`: call `cyc.closure_residual()`, assert finite and `>= 0` (the closed-loop residual is exercised in `test_construct.py`, not here).
- [ ] `uv run pytest tests/model/test_aldrin.py -v` green. **This is the M3 gate.**

## Construct + 2-synodic E–M–E test

- [ ] Create `tests/model/test_construct.py` with named tolerance `TOL_CLOSURE_KMS = 0.05`.
- [ ] `test_construct_two_synodic_em_cycler`:
  - [ ] `T_syn = synodic_period_seconds("E", "M")` from M2.
  - [ ] `times = [0.0, 146*SECONDS_PER_DAY, T_syn]`.
  - [ ] `cyc = construct_cycler(["E", "M", "E"], times, ephem)`.
  - [ ] Assert `len(encounters) == 3`, `len(legs) == 2`, `period == pytest.approx(T_syn)`.
  - [ ] Assert `cyc.closure_residual() < TOL_CLOSURE_KMS`.
- [ ] `test_construct_validates_input_lengths`: length mismatch raises `ValueError`; non-monotonic times raise `ValueError`.
- [ ] `test_construct_unknown_body_raises`: `["E", "X"]` raises `ValueError`.
- [ ] `test_construct_multi_rev_branch_selection`: a long E→E transit where `n_revs=0` has no solution; passing `n_revs=1, branch="low"` returns a valid `Leg` with the right labels. (If the M1 Lambert solver does not yet support multi-rev, this test is marked `xfail` and noted in hand-off to M4.)
- [ ] `uv run pytest tests/model/test_construct.py` green.

## Local green (full sweep)

- [ ] `uv run pytest` green across the whole suite (M0 + M1 + M2 + M3 tests).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean (run `uv run ruff format .` if it isn't).
- [ ] `uv run mypy src tests` clean. Resolve any `NDArray[np.float64]` strict-mode complaints in the modules themselves; do **not** sprinkle `# type: ignore`.

## CI

- [ ] Push branch; confirm GitHub Actions runs all four checks (ruff lint, ruff format, mypy, pytest) and is green.
- [ ] Commit message: `m3: cycler model, rotating frame, patched-conic construct; reproduce Aldrin`.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M3 status `planned` → `completed`; M4 row `not yet planned` → `planned`.
- [ ] Append a `## Hand-off to M4` section to this todo.md noting:
  - [ ] Which Lambert branch the 2-synodic E–M–E test actually used (n_revs=0? n_revs=1, low? high?).
  - [ ] Observed `closure_residual` value for the Aldrin and 2-synodic constructions (purely informational; informs M5's residual-minimisation target).
  - [ ] Whether any test tolerance had to be loosened from the plan, and why.
  - [ ] Any API surface that turned out wrong (e.g. `branch="single"` vs other vocabulary choices) and what M4 should know before extending it.
  - [ ] Confirmed: ready to write the M4 plan doc.

## Hand-off to M4

_To be filled in when M3 completes._
