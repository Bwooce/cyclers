# M6a — Idealized closure verification (todo)

Working checklist for the M6a milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: predecessor recap → dynamic-frame extension to `core/frames.py` → `verify/` subpackage scaffold → `verify/propagate.py` skeleton → helpers (each with paired unit tests) → `verify_long_term_stability` + gate tests → local green → CI → closeout.

## Predecessor recap

- [ ] Re-read spec §12(c) (dynamic ephemeris frame + tolerant verification) — **the binding architecture decision for M6a**.
- [ ] Re-read spec §10 (closure-frame correctness risk) — the binding unit-test requirement on the rotating-frame transform.
- [ ] Re-read spec §8 (M6 milestone: "best E–M cycler verified periodic over ≥3 laps") and §14 V2 (multi-lap periodicity gate) — the gate semantics.
- [ ] Re-read M5's `## Hand-off to M6a` in `docs/phases/m5-optimisation/todo.md`: note actual V∞ reproduced, closure residual at convergence, whether `use_de` was needed, any failed cells. Use the reported `closure_residual_kms` as a sanity floor — M6a's drift cannot reasonably be smaller than this idealised residual integrated over a lap.
- [ ] Re-read M3's `phases/m3-model-construct/plan.md` §3.2 (uniform-frame design) — M6a extends, does NOT replace.
- [ ] Confirm `core/ephemeris.py::Ephemeris("astropy")` works on the dev environment (M6 slice `9b2611d` introduced it). Spot-check `state("E", 0)` returns a real heliocentric position ≈ 1 AU.
- [ ] Confirm `core/kepler.py::propagate(r0, v0, dt)` signature — used by `propagate_lap`.
- [ ] Confirm `model/cycler.py::Cycler.legs` carries `(t_depart, t_arrive, v_depart, v_arrive)` per leg — these are the per-leg states the multi-lap propagator consumes.
- [ ] Confirm `search/phase_match.py::find_real_windows` is independent of M6a (it produces inputs M6a consumes — `t_start` epoch — not the reverse).
- [ ] Confirm `search/optimize.py::optimise_cell_ephemeris` raises `NotImplementedError("requires M6 ephemeris backend")` — M6a leaves this in place; M6b fills it.
- [ ] Confirm parallel Schema-v2 backfill agent is NOT touching `core/frames.py`, `tests/core/test_frames.py`, `tests/_catalogue_loader.py`, or anything under `src/cyclerfinder/verify/` (the M6a workspace).

## `core/frames.py` — dynamic frame extension (strictly additive)

### Skeleton + module docstring update

- [ ] Update `core/frames.py` module docstring to note coexistence: the M3 uniform frame stays for circular-coplanar work (M3 Aldrin, M5 optimiser inner loop), the M6a dynamic frame is for real-ephemeris verification. Cite spec §12(c).
- [ ] Add `from cyclerfinder.core.ephemeris import Ephemeris` import to `core/frames.py` — this is the only new module dependency (frames was previously constants-only).
- [ ] Stub `synodic_omega_dynamic`, `to_rotating_dynamic`, `from_rotating_dynamic` with full docstrings + `NotImplementedError` bodies. Cite spec §12(c) on each.
- [ ] `uv run mypy src/cyclerfinder/core/frames.py` clean on the skeleton (catches `tuple[str, ...]` typing, NDArray shape annotations).

### `synodic_omega_dynamic`

- [ ] Implement per plan §3.1: compute `(r_b0, v_b0) = ephem.state(bodies[0], t_sec)` then `omega = |r_b0 × v_b0| / ||r_b0||²` with sign from `(r_b0 × v_b0)[2]` to match prograde convention. For `len(bodies) >= 2`, optionally average with `bodies[1]`'s rate if both are needed (M6a's 2-body case uses `bodies[0]` only; document the extension point for M8 VEM).
- [ ] Test `test_synodic_omega_dynamic_matches_uniform_for_circular_earth` — at `t ∈ {0, 1 yr, 2 yr, 5 yr}` on the circular backend, the returned rate equals `synodic_omega("E")` to 1e-12 rad/s.
- [ ] Test `test_synodic_omega_dynamic_varies_with_eccentricity_for_astropy` — at Earth's perihelion vs aphelion (separated by ~6 months on the astropy backend), rates differ by ≥ 2%.

### `to_rotating_dynamic` and `from_rotating_dynamic`

- [ ] Implement `to_rotating_dynamic` per plan §3.1.1: read `(r_b0, _v_b0) = ephem.state(bodies[0], t_sec)`; compute frame angle `theta = atan2(r_b0[1], r_b0[0])` (read, not integrated); compute instantaneous `omega = synodic_omega_dynamic(t_sec, bodies, ephem)`; apply Coriolis correction `v - omega × r`; rotate `r` and corrected `v` by `R(-theta)`. Return `(r_rot, v_rot)`.
- [ ] Implement `from_rotating_dynamic` as the algebraic inverse: read the **same** `theta` from `ephem.state(bodies[0], t_sec)` (this is what makes the round-trip exact); rotate by `R(+theta)`; add back `omega × r_inertial`. Return `(r_inertial, v_inertial)`.
- [ ] Test `test_dynamic_frame_roundtrip_identity` (**spec §10 BINDING GATE**) — across a 5-year grid of `t_sec` at 1-month spacing, for a synthetic state at 1.5 AU position + 30 km/s velocity, `from_rotating_dynamic(*to_rotating_dynamic(r, v, t, ("E","M"), ephem), t, ("E","M"), ephem)` recovers `(r, v)` to ≤ 1e-10 relative per-component. **If this fails, halt M6a — every drift measurement downstream depends on it.**
- [ ] Test `test_dynamic_frame_reduces_to_uniform_for_circular` — for `Ephemeris("circular")`, `bodies=("E","M")`, agreement with `to_rotating(r, v, t, synodic_omega("E"))` to 1e-8 km position / 1e-12 km/s velocity across a representative grid.
- [ ] Test `test_to_rotating_dynamic_anchor_x_axis` — `to_rotating_dynamic` of body[0]'s own state has y-component ≤ 1e-6 km (body[0] sits on +x by construction).
- [ ] Test `test_to_rotating_dynamic_inertial_at_t_zero_matches_uniform` — at `t=0` on circular backend, dynamic = uniform exactly.

### Module-level checks

- [ ] Confirm `uv run pytest tests/core/test_frames_dynamic.py` green.
- [ ] Confirm `uv run pytest tests/core/test_frames.py` (M3's existing uniform-frame tests) still green — M3 must not regress.
- [ ] Confirm `uv run mypy src/cyclerfinder/core/frames.py` clean.

## `verify/` subpackage scaffold

- [ ] Create empty `src/cyclerfinder/verify/__init__.py`.
- [ ] Create empty `tests/verify/__init__.py`.
- [ ] Confirm `uv run mypy src` clean — empty `__init__.py` is fine.
- [ ] Confirm `uv run pytest tests/verify` runs (no tests yet, exits cleanly with collected=0).

## `verify/propagate.py` — skeleton

- [ ] Create `src/cyclerfinder/verify/propagate.py` with module docstring referencing spec §5 step 7, §8 M6, §12(c), §14 V2.
- [ ] Define `DRIFT_TOLERANCE_KM: Final[float] = 50_000.0` with the §4.3 rationale in the docstring.
- [ ] Define `StabilityReport` frozen dataclass per plan §3.3 with all eight fields (`cycler_id`, `n_laps_propagated`, `max_drift_km`, `max_drift_lap_index`, `per_lap_drift_km`, `stable`, `per_lap_dv`, `total_tcm_dv`, `frame_used`).
- [ ] Stub `propagate_lap`, `multi_lap_propagation`, `lap_to_lap_drift`, `verify_long_term_stability` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src/cyclerfinder/verify/propagate.py` clean on the skeleton.

## Helpers (one at a time, paired tests)

### `_resolve_frame_bodies`

- [ ] Implement per plan §3.4: if `frame_bodies` is supplied, return it; else look at `cycler.bodies`, dedupe in order, take first two. For Aldrin `["E","M","E"]` returns `("E","M")`.
- [ ] Test `test_resolve_frame_bodies_aldrin` — `cycler.bodies = ["E","M","E"]` → `("E","M")`.
- [ ] Test `test_resolve_frame_bodies_2syn_em` — `cycler.bodies = ["E","E","M","M"]` → `("E","M")`.
- [ ] Test `test_resolve_frame_bodies_explicit_override` — caller passes `("E",)`; helper returns `("E",)` unchanged.

### `propagate_lap`

- [ ] Implement per plan §3.2.1: build `t_sample = np.linspace(t_start, t_end, n_samples)`; for each `t`, find the leg containing it, propagate from leg start state via `core.kepler.propagate(r_depart, v_depart, t - t_depart)`. Return `(n_samples, 7)` matrix `[t, x, y, z, vx, vy, vz]`.
- [ ] Test `test_propagate_lap_matches_construct_at_encounters` — for `build_aldrin_seed(Ephemeris("circular"))`, every `Encounter.r` appears in the sampled output to ≤ 10 km.
- [ ] Test `test_propagate_lap_n_samples_respected` — `n_samples=100` → returned shape `(100, 7)`.
- [ ] Test `test_propagate_lap_continuous_across_legs` — at a leg boundary, consecutive sample positions are within `||v|| * dt` where `dt = (t_end - t_start) / n_samples`.

### `multi_lap_propagation`

- [ ] Implement per plan §3.2.2: continuous propagation through `n_laps` laps from `t_start`. For lap `i`, propagate from the lap-`i-1` end state through to `t_start + (i+1) * cycler.period`, sampling `n_samples_per_lap` points. Return dict with `"samples"`, `"lap_indices"`, `"lap_start_times"` keys.
- [ ] **Critical:** lap `i+1` starts from the propagated end state of lap `i`, NOT from a fresh re-instantiation of the cycler. This is the "continuous" word in spec §12.
- [ ] Test `test_multi_lap_propagation_lap_count` — `n_laps=3, n_samples_per_lap=100` → `samples.shape == (300, 7)`.
- [ ] Test `test_multi_lap_propagation_lap_indices_monotone` — `lap_indices` strictly increasing; `lap_indices[0] == 0`; `lap_indices[-1] == 300`.
- [ ] Test `test_multi_lap_propagation_continuity_at_lap_boundaries` — sample at the end of lap `i` and the start of lap `i+1` have identical `(x, y, z, vx, vy, vz)` to numerical precision (they're the same propagated state).

### `lap_to_lap_drift`

- [ ] Implement per plan §3.2.3: for each row index `i`, transform `samples_0[i]` and `samples_n[i]` into the dynamic rotating frame via `to_rotating_dynamic` at their respective `t`s, compute `||r_0_rot - r_n_rot||`. Return max across all `i`.
- [ ] Test `test_lap_to_lap_drift_same_lap_is_zero` — `lap_to_lap_drift(samples_0, samples_0, ("E","M"), ephem) == 0.0` exact.
- [ ] Test `test_lap_to_lap_drift_zero_for_circular_aldrin` — Aldrin circular-coplanar, lap-0 vs lap-2: drift < 100 km (numerical-noise floor).
- [ ] Test `test_lap_to_lap_drift_translation_in_inertial_frame_does_not_lie` — a pure inertial translation of `samples_lap_n` produces nonzero rotating-frame drift; the frame does not absorb it.

## `verify_long_term_stability` — compose helpers

- [ ] Implement per plan §3.3 pseudocode:
  - [ ] Resolve `bodies = _resolve_frame_bodies(cycler, frame_bodies)`.
  - [ ] Propagate continuously via `multi_lap_propagation(cycler, ephem, n_laps, t_start, n_samples_per_lap)`.
  - [ ] For each consecutive pair `(i, i+1)`, compute `lap_to_lap_drift(samples_lap_i, samples_lap_{i+1}, bodies, ephem)`. Track `max_drift_km` and `max_drift_lap_index`.
  - [ ] Build `per_lap_drift_km[i] = lap_to_lap_drift(samples_lap_0, samples_lap_{i+1}, bodies, ephem)` for `i ∈ range(n_laps)`.
  - [ ] Set `stable = max_drift_km < DRIFT_TOLERANCE_KM`.
  - [ ] Populate `per_lap_dv = (0.0,) * n_laps`, `total_tcm_dv = 0.0` — M6b's territory.
  - [ ] Set `frame_used = "uniform" if use_uniform_frame else "dynamic"`.
  - [ ] Return `StabilityReport(...)`.
- [ ] Test `test_aldrin_cycler_periodic_over_3_laps_circular` (warm-up) — `build_aldrin_seed(Ephemeris("circular"))`, `n_laps=3`, `use_uniform_frame=True` → `stable=True`, `max_drift_km < 1.0` km, `frame_used="uniform"`.
- [ ] **Gate test** `test_2syn_em_cycler_periodic_over_3_laps` (**M6a BINDING GATE — spec §8 M6**):
  - [ ] Build the 2-syn E-M cycler. Preferred: `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)[0].best_cycler`. Fallback: construct from `s1l1-2syn-em-cpom` catalogue entry via `tests/_catalogue_loader.py`.
  - [ ] `report = verify_long_term_stability(cycler, n_laps=3, ephem=Ephemeris("astropy"))`.
  - [ ] Assert `report.stable == True`, `report.max_drift_km < 50_000`, `report.n_laps_propagated == 3`, `report.frame_used == "dynamic"`, `report.per_lap_dv == (0.0, 0.0, 0.0)`, `report.total_tcm_dv == 0.0`.
  - [ ] **If this fails, escalate per plan §5 risk #2.** DO NOT widen `DRIFT_TOLERANCE_KM` — investigate (a) M5's `closure_residual_kms`, (b) alternative `frame_bodies` like `("E",)`, (c) 2-lap vs 3-lap behaviour. Document findings in the hand-off note.
- [ ] Test `test_2syn_em_cycler_drifts_predictably_per_lap` — same cycler; `per_lap_drift_km` is monotone non-decreasing or grows linearly (`per_lap_drift_km[i+1] >= per_lap_drift_km[i] * 0.9`).
- [ ] Test `test_5_lap_drift_stays_under_tolerance` — same cycler, `n_laps=5` → `max_drift_km < 50_000`. M6b needs 5-lap horizons; M6a confirms the geometry holds.
- [ ] Test `test_stability_report_frozen_and_fields_locked` — `report.max_drift_km = …` raises `FrozenInstanceError`; `per_lap_dv == (0.0,) * n_laps`; `total_tcm_dv == 0.0`; `frame_used in ("dynamic", "uniform")`.
- [ ] Test `test_verify_long_term_stability_deterministic` — two calls with same inputs produce bitwise-identical `max_drift_km`.
- [ ] Test `test_verify_long_term_stability_independent_of_t_start_on_circular` — circular backend is time-translation invariant; `t_start=0` and `t_start=1e6` produce the same `max_drift_km` to numerical precision.

## Sanity: M3 frame and M5 optimiser still pass

- [ ] `uv run pytest tests/core/test_frames.py` green — M3's uniform-frame tests are not affected by the additive M6a extension.
- [ ] `uv run pytest tests/model/test_aldrin.py` green — Aldrin gate from M3 unchanged.
- [ ] `uv run pytest tests/search/test_optimize.py` green — M5's optimiser inner loop still uses `to_rotating` (uniform), not `to_rotating_dynamic`.
- [ ] `uv run pytest tests/search/test_phase_match.py` green — the M6 slice's phase-match tests are independent of M6a.

## Local green

- [ ] `uv run pytest` → green (M0–M6a all passing).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean.
- [ ] `uv run mypy src tests` → clean under strict mode.

## CI

- [ ] Commit: `m6a: dynamic rotating frame + verify/propagate.py — 2-syn E-M periodic over 3 laps`.
- [ ] Push; confirm GitHub Actions runs and all checks pass.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M6a status → `completed`; M6b row → `planned`.
- [ ] (Optional) Update `docs/overview.md` §2 with a note on the dynamic/uniform frame coexistence: "M6a adds a non-uniform dynamic rotating frame (spec §12(c)) to `core/frames.py` for real-ephemeris verification. M3's uniform frame stays for circular-coplanar work (M3 Aldrin, M5 optimiser); the dynamic frame is the verification frame for M6a/M6b/M7 V2 batch validation."
- [ ] Append a `## Hand-off to M6b` section to this `todo.md` (below).

## Hand-off to M6b

*(Filled in at M6a closeout. Placeholder structure below; replace bullets with measured values from the M6a test pass.)*

- **Actual `max_drift_km` reproduced for the 2-syn E-M gate:** `<TBD>` km (vs the 50,000 km bound). The headroom `(50,000 - max_drift_km)` is the budget M6b's ephemeris-mode TCM optimiser inherits — TCMs reduce drift further, so M6b's `POST_TCM_DRIFT_TOLERANCE_KM` should be set well below `max_drift_km`.
- **Worst lap pair (`max_drift_lap_index`):** `<TBD>`. If `max_drift_lap_index == 0` the cycler shows immediate drift; if `>= 1`, drift accumulates as expected for an idealised cycler on a real ephemeris.
- **Frame-bodies decision:** the gate's `_resolve_frame_bodies` policy returned `<TBD>`; this was `<TBD: validated as correct | overridden because …>`. M6b's `optimise_cell_ephemeris` should default to the same policy.
- **Per-test wall-clock runtime for the M6a binding gate:** `<TBD>` seconds. M6b's ephemeris-mode optimiser calls `verify_long_term_stability` inside its objective; the per-call cost is the unit of M6b's compute budget. Roughly `<TBD>` calls per ephemeris-mode optimisation × `<TBD>` seconds = `<TBD>` seconds per cell.
- **5-lap drift result:** `test_5_lap_drift_stays_under_tolerance` reported `max_drift_km = <TBD>` at `n_laps=5`. If `< 50,000`, M6b can safely set its horizon to 5 laps. If borderline, M6b may need `n_laps=3` or to tighten the M5-fed initial guess.
- **Risks fired:**
  - Risk #2 (drift > 50k): `<TBD: no | yes — resolution was …>`
  - Risk #3 (propagator noise): `<TBD: no | yes — Kepler tolerance tightened to …>`
  - Risk #5 (frames module import of astropy): `<TBD: no | yes — guarded by passing Ephemeris as parameter>`
- **Whether M5's `find_cyclers` was usable as the gate fixture, or whether the catalogue fallback was needed:** `<TBD>`. If the catalogue fallback was needed, document why (M5 reproduction outside tolerance? M5 result not seeded reproducibly?) — this is a signal for M5 hand-off accuracy that future milestones can use.
- **Locked `StabilityReport` shape** (M6b inherits exactly this):
  - `cycler_id: str | None`
  - `n_laps_propagated: int`
  - `max_drift_km: float`
  - `max_drift_lap_index: int`
  - `per_lap_drift_km: tuple[float, ...]`
  - `stable: bool`
  - `per_lap_dv: tuple[float, ...]`  ← M6b populates with per-lap TCM ΔV (km/s).
  - `total_tcm_dv: float`  ← M6b populates with `sum(per_lap_dv)`.
  - `frame_used: str`  ← M6b uses `"dynamic"` (the default) throughout.
- **API surfaces M6b will consume immediately:**
  - `from cyclerfinder.verify.propagate import verify_long_term_stability, StabilityReport, DRIFT_TOLERANCE_KM, multi_lap_propagation, lap_to_lap_drift`
  - `from cyclerfinder.core.frames import to_rotating_dynamic, from_rotating_dynamic, synodic_omega_dynamic` — M6b's ephemeris-mode optimiser objective transforms candidate trajectories through this frame to compute the drift cost.
  - `from cyclerfinder.search.optimize import optimise_cell_ephemeris` — M6a leaves this raising `NotImplementedError`; M6b fills the body.
- **Spec ambiguities encountered during M6a implementation:**
  - **Frame anchor: Sun-Earth vs Earth-Mars pair.** Spec §12(c) wording "Sun–Earth (or Earth–Mars)" leaves both viable. M6a's `_resolve_frame_bodies` defaults to `("E", "M")` (Earth-Mars pair); the helper is overrideable. M6b plan should revisit if its TCM-budget converges differently in the two frames.
  - **Per-lap-DV semantics.** M6a locks `per_lap_dv` as a `tuple[float, ...]` of length `n_laps_propagated`. M6b should clarify: is this the **TCM at lap boundary** or the **integrated TCM through the lap**? M6a doesn't pin this; M6b decides at implementation.
- **Recommended M6b first steps:**
  1. Read this hand-off + spec §12(a) (ephemeris-mode optimisation) + spec §12.1 (phase-match bridge — already implemented in the M6 slice via `find_real_windows`).
  2. Wire `optimise_cell_ephemeris` against `verify_long_term_stability` as the drift-feasibility check inside its objective.
  3. Use `phase_match.find_real_windows(signature, ephem, date_range)` to pick the launch epoch; pass it through to `verify_long_term_stability` as `t_start`.
  4. M6b gate: best 2-syn E-M cycler's `total_tcm_dv` over 5 laps stays within a documented budget (e.g. 200 m/s, citing McConaghy 2006).

M6b's first task is the M6b plan doc itself; this hand-off is the input it consumes.
