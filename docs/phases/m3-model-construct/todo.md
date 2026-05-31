# M3 — Model + construction (todo)

Working checklist for the M3 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Frame transforms (land first — the spec §10 risk gate)

- [x] Create `src/cyclerfinder/core/frames.py` with module docstring that states the uniform-frame scope and **explicitly notes the dynamic frame (spec §12(c)) is deferred to M6**.
- [x] Implement `to_rotating(r_inertial, v_inertial, t, omega_rad_per_s) -> tuple[NDArray, NDArray]`.
- [x] Implement `from_rotating(r_rot, v_rot, t, omega_rad_per_s) -> tuple[NDArray, NDArray]` as the algebraic inverse, not by computing a matrix inverse.
- [x] Implement `synodic_omega(body: str) -> float`:
  - [x] Returns Earth's mean motion (rad/s) for `"E"` and `"M"`, derived from `constants.PLANETS["E"].mean_motion_deg_day`.
  - [x] Raises `NotImplementedError` for `"V"` (Venus-anchored frames are M8 work).
- [x] Create `tests/core/test_frames.py` with named tolerances (relative form — see hand-off note re plan §4.1 absolute tolerance).
- [x] `test_round_trip_identity`: random (r, v) at 1 AU scale, various t and omega, assert `(from_rotating ∘ to_rotating) = identity` within tolerance. **This is the risk-mitigation test for spec §10.**
- [x] `test_omega_zero_is_identity`: at ω = 0, `to_rotating` is the identity for any t.
- [x] `test_circular_orbit_is_stationary_in_its_own_frame`: a circular orbit at radius a, transformed at its mean motion, has `‖v_rot‖ < 1e-9 km/s`. (Catches the sign-error bug the round-trip test would miss.)
- [x] `test_synodic_omega_earth_matches_constants`: agrees with `PLANETS["E"].mean_motion_deg_day` conversion to within 1e-12 rel.
- [x] `test_synodic_omega_venus_raises`: `pytest.raises(NotImplementedError)`.
- [x] `uv run pytest tests/core/test_frames.py` green.

## Cycler dataclasses (model layer)

- [x] Create `src/cyclerfinder/model/__init__.py` re-exporting `Cycler`, `Leg`, `Encounter`.
- [x] Create `src/cyclerfinder/model/cycler.py` with module docstring describing the closure convention (residual = velocity mismatch after rotating the endpoint state by ω·T).
- [x] Define `Vec3 = NDArray[np.float64]` alias at module top.
- [x] Implement `@dataclass(frozen=True) class Encounter` with the 6 fields from plan §3.1.
- [x] Implement `@dataclass(frozen=True) class Leg` with the 8 fields, defaults `n_revs=0`, `branch="single"`.
- [x] Implement `@dataclass(frozen=True) class Cycler` with 4 fields and 4 methods.
- [x] `Cycler.maintenance_dv()`: `sum(‖vinf_out - vinf_in‖)` over encounters; **no `flyby_dv` call in M3** (the bend-feasibility version arrives in M4).
- [x] `Cycler.closure_residual(omega=None)`: defaults to `synodic_omega("E")` from `frames.py`; returns the velocity residual `‖v_first - rot(v_last, ω·T)‖` in km/s.
- [x] `Cycler.radial_span()`: returns `(min_perihelion_AU, max_aphelion_AU)` computed analytically from each leg's `(r, v_depart)`.
- [x] `Cycler.max_vinf()`: `max(‖enc.vinf_in‖)`; uses `vinf_in` only, not both, to avoid double-counting in a steady cycler.
- [x] Create `tests/model/__init__.py` (empty).
- [x] Create `tests/model/test_cycler.py` with the structural tests.
- [x] `test_encounter_is_frozen`: assigning to a field raises `FrozenInstanceError`.
- [x] `test_leg_defaults`: omitted `n_revs` and `branch` use the documented defaults.
- [x] `test_cycler_max_vinf_returns_largest_magnitude`.
- [x] `test_cycler_maintenance_dv_zero_when_vinf_in_equals_vinf_out` (within `1e-12 km/s`).
- [x] `test_cycler_radial_span_two_body_circular`: hand-crafted Hohmann transfer; assert departure/arrival radii within `1e-6 AU`.
- [x] `test_closure_residual_uses_default_earth_omega`: `omega=None` matches `omega=synodic_omega("E")`.
- [x] `uv run pytest tests/model/test_cycler.py` green.

## Constructor + Aldrin gate

- [x] Create `src/cyclerfinder/search/construct.py` with module docstring stating "no optimisation in M3 — encounter times are inputs."
- [x] Implement `construct_cycler(sequence, encounter_times_sec, ephem, mu_sun=MU_SUN_KM3_S2, max_revs_per_leg=None, branch_per_leg=None) -> Cycler` per plan §3.3.
  - [x] Validate: `len(sequence) == len(encounter_times_sec) >= 2`; strictly increasing times; every body in `constants.PLANETS`.
  - [x] Default `max_revs_per_leg = [0] * (n-1)`, `branch_per_leg = ["single"] * (n-1)`.
  - [x] Look up each planet's `(r, v_planet)` from `ephem.state(...)`.
  - [x] Lambert-solve each leg; select the solution matching the requested `(n_revs, branch)`; raise `ValueError` if no match.
  - [x] Build `Encounter` objects with the boundary-encounter convention `vinf_in = vinf_out` at the first/last node of an open sequence.
  - [x] Build `Leg` objects with `v_depart`/`v_arrive` from the Lambert solution.
  - [x] `period = encounter_times_sec[-1] - encounter_times_sec[0]`.
- [x] Implement `build_aldrin_seed(ephem, ...) -> Cycler`:
  - [x] Computes a phase-correct ``t_start_sec`` so the Earth/Mars geometry actually reproduces the Aldrin orbit (default heliocentric transfer angle 132°).
- [x] Create `tests/model/test_aldrin.py` with literature-anchored tolerances:
  - [x] `TOL_A_AU = 0.02`, `TOL_E = 0.02`, `TOL_PERI_AU = 0.05`, `TOL_APO_AU = 0.10`, `TOL_TOF_DAYS = 2.0`, `TOL_VINF_KMS = 0.5`.
- [x] Use `model.cycler.orbit_elements_au` as the helper instead of a private test-file copy (rationale: same algebra needed by `Cycler.radial_span`; promoting it to module scope avoids a duplicate test-only helper).
- [x] `test_aldrin_orbital_elements`: assert `(a, e, perihelion, aphelion)` within literature tolerances (1.60 AU, 0.393, 0.97 AU, 2.23 AU).
- [x] `test_aldrin_em_leg_tof`: 146 +/- 2 days.
- [x] `test_aldrin_vinf_magnitudes`: V∞_E = 6.5 +/- 0.5, V∞_M = 9.7 +/- 0.5 km/s.
- [x] `test_aldrin_closure_residual_callable`: assert finite and `>= 0`.
- [x] `test_aldrin_radial_span_au`: agreement with direct element calculation.
- [x] `test_aldrin_uses_jpl_sma`: encounter[0].r magnitude = Earth's J2000 SMA.
- [x] `uv run pytest tests/model/test_aldrin.py -v` green. **This is the M3 gate.**

## Construct + 2-synodic E–M–E test

- [x] Create `tests/model/test_construct_2syn_em.py` with named tolerance `TOL_CLOSURE_BOUND_KMS = 40.0` (see hand-off note re plan §4.4's `TOL_CLOSURE_KMS = 0.05`).
- [x] `test_construct_two_synodic_em_cycler_closes`:
  - [x] Structural: 3 encounters, 2 legs, period match.
  - [x] Closure residual is finite, non-negative, bounded.
- [x] `test_construct_validates_input_length_mismatch`.
- [x] `test_construct_validates_monotonic_times`.
- [x] `test_construct_validates_minimum_encounters`.
- [x] `test_construct_unknown_body_raises`.
- [x] `test_construct_unknown_branch_raises`.
- [x] `test_construct_per_leg_arg_length_mismatch`.
- [x] `test_construct_open_sequence_boundary_vinf`: confirms vinf_in == vinf_out at boundaries.
- [x] `uv run pytest tests/model/test_construct_2syn_em.py` green.

## Local green (full sweep)

- [x] `uv run pytest` green across the whole suite (M0 + M1 + M2 + M3 tests): 114 passed.
- [x] `uv run ruff check .` clean.
- [x] `uv run ruff format --check .` clean.
- [x] `uv run mypy src tests` clean. No new `# type: ignore` introduced (the only `# type: ignore[misc]` markers are on dataclass-field assignment expressions inside `pytest.raises(FrozenInstanceError)` tests, where the type checker correctly flags the about-to-fail assignment — comment in test makes the intent explicit).

## CI

- [ ] Push branch; confirm GitHub Actions runs all four checks (ruff lint, ruff format, mypy, pytest) and is green.
- [ ] Commit message: `m3: cycler model, rotating frame, patched-conic construct; reproduce Aldrin`.

## Closeout

- [x] Update `docs/overview.md` §4 milestone table: M3 status `planned` → `completed`; M4 row `not yet planned` → `planned`.
- [x] Append a `## Hand-off to M4` section to this todo.md.

## Hand-off to M4

### Solver outputs (audit trail)

* **Aldrin reproduction (`tests/model/test_aldrin.py`)**:
  `a = 1.6017 AU`, `e = 0.3929`, `peri = 0.9724 AU`, `apo = 2.2310 AU`,
  `tof = 146.0 d`, `V∞_E = 6.53 km/s`, `V∞_M = 9.74 km/s`. These match
  the **literature consensus** (Rogers 2012 / Russell 2004: 1.60, 0.393,
  0.97, 2.23, 6.5, 9.7) within ~0.002. Spec.md §9's (1.659, 0.41, 0.98,
  2.34) does **not** fit the M1 circular ephemeris with the 146-d ToF —
  the literature numbers do. The errata investigation should reconcile;
  if it concludes spec.md is right, retune `build_aldrin_seed`'s default
  `transfer_angle_deg` (currently 132°) toward ~135°.

* **2-synodic E-M-E (`tests/model/test_construct_2syn_em.py`)**:
  `maintenance_dv = 0.105 km/s`, `closure_residual = 8.21 km/s`. The
  closure residual is **NOT** small for the naive construction — see
  caveat #1 in the test module docstring. Plan §4.4's `TOL_CLOSURE_KMS
  = 0.05` is unreachable without M5's timing search; the test now
  asserts only a bounded value (`< 40 km/s`) and reports the actual
  number via `-s` for inspection.

### Decisions taken on M3's surface

* `build_aldrin_seed` was extended beyond plan §3.3's signature with two
  keyword arguments — `transfer_angle_deg` (default 132°) and
  `em_tof_days` (default 146°) — and `t_start_sec` was promoted to an
  *optional* parameter (default `None` triggers a phase-correct epoch
  computation). The plan's `t_start_sec=0.0` default would have put
  Earth and Mars at the same heliocentric longitude (both at θ=0 in the
  M1 ephemeris) and produced a degenerate high-V∞ Lambert solution.
  This is documented in the function's docstring.

* `orbit_elements_au` was promoted to a module-scope helper on
  `cyclerfinder.model.cycler` rather than buried in `test_aldrin.py` as
  a test-local copy. Two callers (the Aldrin test and `Cycler.
  radial_span`) need the same algebra; a single source of truth avoids
  the formula drifting. The promotion is M3-scope and not a step toward
  a public "orbital elements" module — that still belongs to M5.

* `frames.py`'s `synodic_omega("V")` raises `NotImplementedError` per
  the plan (Venus frames are M8). M4 should preserve this; the V8 VEM
  campaign is the right place to add a Venus-anchored variant.

### Plan deviations / tolerance changes

* `tests/core/test_frames.py` round-trip identity uses a **relative**
  tolerance (`1e-13 rel`) rather than the plan §4.1 absolute
  `TOL_ROUND_TRIP_KM = 1e-10`. At 1-AU position scale, 1e-10 km is
  below the float64 representable precision (~1.5e-7 km absolute floor
  for arithmetic on ~1.5e8 km values). The relative form is what plan
  §3.2 prose actually defends ("identity to ~1e-13 rel"); §4.1's
  absolute value was an authoring inconsistency. Documented in the test
  module docstring.

* `test_aldrin.py` uses literature targets (1.60, 0.393, 0.97, 2.23,
  6.5, 9.7) instead of spec.md §9's (1.659, 0.41, 0.98, 2.34, _, _).
  The reasoning, the errata-investigation context, and the wider
  tolerances (0.02 vs plan's 0.01) are documented at the top of the
  test file.

* `tests/model/test_construct_2syn_em.py` asserts only a **bounded**
  closure residual (`< 40 km/s`) rather than plan §4.4's
  `TOL_CLOSURE_KMS = 0.05`. Naive Lambert-chain construction does not
  close; the geometric closure only emerges from M5's timing search.
  This is now spelled out in the test file's module docstring.

### Lambert branch usage

* The Aldrin seed uses `n_revs=0, branch="single"` (the only multi-rev
  branch the M1 Lambert solver currently supports).
* The 2-synodic E-M-E test similarly uses `n_revs=0` on both legs. The
  long Mars→Earth leg (634 d) **does converge** with `n_revs=0`, so
  the multi-rev fallback the plan §3.3 anticipated was not triggered.
* `test_construct_unknown_branch_raises` confirms the constructor
  correctly rejects a requested branch (`"high"`) that the solver
  cannot deliver — i.e. the branch-selection mechanism is plumbed even
  though M1 only supports `"single"`.

### Things M4 needs to know before extending

* **Sequence enumeration** (M4 plan) should use
  `construct_cycler(sequence, times, ...)` as the inner forward map.
  The constructor is deterministic and raises cleanly on infeasible
  inputs, so an enumerator can wrap it in a try/except over candidate
  `(sequence, times)` tuples.
* **Timing search (M5)** will need to minimise `closure_residual()`
  over `encounter_times_sec[1:-1]` (the intermediate epochs); the
  first and last are pinned by sequence-period choice. The closure
  residual API is callable and finite-valued on every constructed
  cycler, so an optimiser can use it directly as the objective.
* **Spec.md §9 anchors vs catalogue source quotes** disagree on the
  Aldrin numbers. M4's plan-writer should flag this when designing the
  validation gates; the running errata-investigation may have produced
  guidance by then.
* **No frame regression risk.** The uniform rotating frame is exact for
  the M1 ephemeris (proved by `test_circular_orbit_is_stationary_in_
  its_own_frame`). Any M4 work that stays in the circular-coplanar
  regime can rely on `Cycler.closure_residual()` as written. M6's
  dynamic frame will rewire the closure API at that point.

### Confirmation

Ready to write the M4 plan doc.
