# Modelling Upgrades — Merged Build Plan (4 Stages)

Date: 2026-06-03
Status: planning (some stages large / likely-partial — flagged inline)

This plan merges four stage blueprints into a single dependency-ordered build
plan. Two independent tracks run in parallel:

- **Track A (real-ephemeris closure):** STAGE 3 epoch resolver → STAGE 1 multi-rev
  optimiser. (epoch → optimiser)
- **Track B (3D geometry):** STAGE 2 inclined ephemeris → STAGE 4 3D Tisserand.
  (ephem3d → tiss3d)

The two tracks are independent and may be built concurrently. The only shared
file is `core/constants.py` (PlanetData), touched by STAGE 2 and STAGE 4 — see
the file/conflict map.

## GOLDEN-TEST DISCIPLINE (binding for every test below)

Any test's EXPECTED side must trace to a PUBLISHED source or be explicitly
labelled `# COMPUTED` / `# INVARIANT` / `# BOUND-ONLY`. Never fabricate a number
and never re-assert a value our own code produced as if it were sourced.

Sourced anchors used across this plan (inputs only, unless noted as EXPECTED):

- S1L1 V_inf: **5.65 km/s (Earth), 3.05 km/s (Mars)** — spec.md §9; catalogue
  `s1l1-2syn-em-cpom` source_quotes. (EXPECTED only in the STAGE 1 idealised
  multi-rev test; everywhere else an input.)
- S1L1 E→M ToF **154 d** — Spreen et al. 2020 / Rogers et al. 2012 (input only).
- Aldrin V_inf: **6.5 km/s (E), 9.7 km/s (M)** — Rogers et al. 2012 Table 1 /
  McConaghy 2002 (input only).
- Aldrin E→M ToF **146 d** — McConaghy 2002 Table 4 (input only).
- Aldrin orbit **a=1.60 AU, e=0.393** — Rogers et al. 2012 Table 1; Russell 2004
  Table 3.4 cycler 1.0.1.-1 (EXPECTED — existing golden anchors, must stay
  bit-identical).
- Aldrin Earth flyby **84° required / 72° achievable** — McConaghy/Longuski/
  Byrnes 2002 AIAA 2002-4420 Table 4 (EXPECTED — existing anchors).
- Planet J2000 elements (inc_deg, lan_deg, ecc) — Standish & Williams,
  "Approximate Positions of the Planets", JPL SSD, Table 1 (J2000):
  Venus inc=3.39458°, lan=76.67984255°, ecc=0.00677727;
  Mars inc=1.84969°, lan=49.55953891°, ecc=0.09340065;
  Earth inc≈0.00005°, lan=0.0, ecc=0.01671123.
  (EXPECTED in STAGE 2/4 tests; cite inline in constants.py.)

---

## Build order (recommended)

1. **STAGE 3 — Robust real-ephemeris epoch/family search** (Track A, foundation)
2. **STAGE 1 — Multi-rev + multi-encounter real-ephemeris optimiser** (Track A,
   depends on STAGE 3 epoch resolver being complete)
3. **STAGE 2 — Inclined analytic ephemeris** (Track B, foundation — independent
   of Track A, may run in parallel)
4. **STAGE 4 — 3D Tisserand feasibility** (Track B, depends on STAGE 2 PlanetData
   fields)

Rationale: STAGE 1 consumes the epoch chosen by `_resolve_real_t_start`, so
STAGE 3 must land first to avoid re-touching the resolver. STAGE 4 reads
`PlanetData.inc_deg` / `ecc`, so STAGE 2 (which introduces those fields) must
land first within Track B. The tracks share only `constants.py`; build STAGE 2
before STAGE 4 and the PlanetData field block is added once (STAGE 2 adds
`inc_deg`/`lan_deg`; STAGE 4 adds `ecc` — see conflict map for the single merge
point).

---

## STAGE 3 — Robust real-ephemeris epoch/family search (Track A, FIRST)

### What it unlocks
Replaces the single-signature, calendar-proximity-biased `_resolve_real_t_start`
with a ranked multi-seed candidate pool ranked by V_inf mismatch. This is the
**bug fix that lets the real-ephemeris path land in the correct basin** for
asymmetric families — a precondition for S1L1 closure and for the SnLm /
Hollister rediscovery to find real windows at all.

### Size / partial-ness
Medium. Self-contained: two new pure functions + one resolver-body swap + one
optimiser helper. Low risk. **Does NOT by itself close S1L1** — see risk note.

### Files
- Modify: `src/cyclerfinder/search/phase_match.py`
- Modify: `src/cyclerfinder/verify/real_closure.py` (`_resolve_real_t_start` @321)
- Modify: `src/cyclerfinder/search/optimize.py` (`optimise_cell_ephemeris` @1169)
- Modify: `tests/search/test_s1l1_real_rediscovery.py` (xfail-reason update only)
- New: `tests/search/test_phase_match_robust.py`

### TDD task list
Phase 1 — new pure helpers in `phase_match.py` (no callers yet):
1. RED: write `tests/search/test_phase_match_robust.py`:
   - `test_leg_duration_seeds_returns_symmetric_plus_perturbations` — 2-leg
     (146 d, 634 d), 5 perturb_fracs → 5 PhaseSignatures, every leg > 0,
     `sum(leg_durations_s) == period_s` for all. `# INVARIANT` (period
     conservation is geometry, not physics).
   - `test_leg_duration_seeds_clips_to_min_leg` — −50% on a 60-d leg clips to
     `min_leg_days=30` rather than going negative. `# COMPUTED` (clip boundary).
   - `test_find_candidate_windows_deduplicates_close_windows` — two seeds yield
     windows within 20 d → 1 survivor (lower mismatch_kms). `# COMPUTED`.
   - `test_find_candidate_windows_ranks_by_mismatch_not_calendar` — **primary
     regression test for the proximity bug**: seed_A window near priority with
     mismatch 3.5, seed_B window 4 yr away with mismatch 0.8 → `windows[0]` is
     seed_B (0.8). `# COMPUTED` (the 3.5/0.8 are synthetic inputs; assertion is
     the ranking contract, not a sourced value).
   - `test_find_candidate_windows_aldrin_merges_seeds` — Aldrin sig (E-M,
     tof 146 d, vinf 6.5/9.7) + 2 perturbed variants over 2026–2036 →
     `len(windows) >= 3` and mismatch_kms ascending. `# COMPUTED` (count/order,
     not magnitudes).
2. GREEN: implement `leg_duration_seeds(...)` then `find_candidate_windows(...)`
   in `phase_match.py`; add both to `__all__` (create `__all__` at module
   bottom — none currently). Signatures:
   ```python
   def leg_duration_seeds(bodies, primary_leg_durations_s, vinf_target_kms,
       period_s, *, perturb_fracs=(0.0, 0.10, -0.10, 0.20, -0.20),
       min_leg_days=30.0, max_leg_days_frac=0.95) -> list[PhaseSignature]
   def find_candidate_windows(signatures, ephem, date_range, n=10,
       step_days=5.0, mismatch_cap_kms=20.0, dedup_window_days=30.0)
       -> list[LaunchWindow]
   ```
   Period conservation: 2-leg `[d1,d2]`, fraction f → `[d1*(1+f), d2 - d1*f]`.
   N-leg: apply f to leg 0, `-f*d0/(N-1)` to each remaining leg; clip all to
   `(min_leg_days*SPD, max_leg_days_frac*period_s)` after.
   Dedup: O(N²) sweep, `|dep_A - dep_B| < dedup_window_days` keeps lower mismatch.

Phase 2 — swap `_resolve_real_t_start` body (`real_closure.py:321`):
3. RED: rename existing `test_resolve_real_t_start_prefers_priority_window` →
   `test_resolve_real_t_start_picks_low_mismatch_window`; keep the ±5-yr band
   check (inherited bound) but change the assertion to "resolved window mismatch
   below 20 km/s". `# COMPUTED` bound — NOT a sourced V_inf. Keep
   `test_resolve_real_t_start_returns_none_when_no_window` unchanged (absurd
   vinf=50 km/s → None). `# INVARIANT` (no fabricated window).
4. GREEN: replace body (lines ~356–374): build primary signature (unchanged),
   call `leg_duration_seeds(...)` with `period_s=sum(signature.leg_durations_s)`,
   call `find_candidate_windows(seeds, ephem, (priority-Δ, priority+Δ), ...)`,
   return `windows[0]` epoch (lowest mismatch) or None. Keep the 4-param public
   signature byte-identical. Add imports of the two new helpers.
   ```python
   def _resolve_real_t_start(signature, ephem, priority_date, *,
       window_years=10.0, n_candidates=5, mismatch_cap_kms=20.0) -> float | None
   ```

Phase 3 — `optimise_cell_ephemeris` helper (`optimize.py:1169`, call @1309):
5. GREEN: add `_resolve_t0_multi_seed(cell, seed_days, priority_date, ephem,
   vinf_targets_kms, target_period_sec, n_candidates=5) -> float | None` above
   `optimise_cell_ephemeris`; it builds the primary PhaseSignature from
   `seed_days`, derives period_s from `target_period_sec`, calls
   `leg_duration_seeds` then `find_candidate_windows(..., mismatch_cap_kms=20.0)`,
   returns best-window t0. Replace the single `_resolve_real_t_start` call at
   ~1309 with this helper. `tof_seed_days` is retained but now seeds the
   perturbation grid rather than being passed raw.

Phase 4 — regression + lint:
6. `uv run pytest tests/search/test_phase_match_robust.py
   tests/verify/test_real_closure.py tests/search/test_s1l1_real_rediscovery.py
   tests/search/test_optimize_ephemeris.py -x`
7. `uv run ruff check src/cyclerfinder/search/phase_match.py
   src/cyclerfinder/verify/real_closure.py src/cyclerfinder/search/optimize.py`
8. Update `test_s1l1_real_rediscovery.py` xfail reason from "seeding" to
   "topology confirmed" (the resolver still finds no E→M→E window matching
   V_inf_M=3.05 — that is a topology, not seeding, problem). STAGE 1 then
   re-points it.

### Risks (honest)
- **S1L1 stays xfail after STAGE 3.** The E→M→E cell cannot match
  V_inf_M=3.05 km/s on real ephemeris regardless of ToF seeding; the robust
  resolver confirms the diagnosis but cannot fix topology. Update xfail reason
  to "topology confirmed", do not re-investigate seeding.
- Performance: 5 seeds × ~1460 steps ≈ 7 s/call vs ~1.5 s now. Mitigate with
  `step_days=10` for perturbed seeds, `5` for primary. Mark slow tests.
- Period conservation only clean for 2-leg; N>2 spreads across N−1 legs and
  clipping breaks it — document and test with an N=3 case.
- `_vinf_at_lambert` is single-rev (max_revs=0): asymmetric multi-rev return
  legs still score high mismatch here. NOT fixed in STAGE 3 — that is STAGE 1.

---

## STAGE 1 — Multi-rev + multi-encounter real-ephemeris optimiser (Track A, SECOND)

### What it unlocks
Threads per-leg revolution/branch metadata (`Cell.per_leg_revs` /
`per_leg_branch`) through the whole maintenance solve chain
(`optimise_cell_ephemeris → optimise_maintenance_dv → _build_chain`) and adds
multi-rev-aware ToF floors. This makes the **~200 SnLm rows** and the
**Hollister-15 family** reachable on real ephemeris, and provides the multi-rev
plumbing needed for **S1L1 closure** once the winning topology is pinned.

### Size / partial-ness
Large and load-bearing. The scipy `args=` positional coupling makes the
`_objective` change subtle (a missed call-site silently passes wrong revs).
S1L1 idealised test starts xfail pending topology characterisation. **Likely
partial** until the topology script names the winner.

### Files
- Modify: `src/cyclerfinder/search/maintain.py`
  (`_build_chain` @303, `_objective` @349, `optimise_maintenance_dv` @389,
  `optimise_aldrin_maintenance_dv` @608 — audit only)
- Modify: `src/cyclerfinder/search/optimize.py`
  (`_ephemeris_tof_seed_and_bounds` @1144, `optimise_cell_ephemeris` @1169)
- Modify: `tests/search/test_maintain.py`
- Modify: `tests/search/test_optimize_ephemeris.py`
- Modify: `tests/search/test_s1l1_real_rediscovery.py`
- New: `tests/search/test_s1l1_idealised_multirev.py`

### TDD task list
1. RED `test_maintain.py::test_build_chain_passes_revs_and_branch` — `_build_chain`
   with non-None per_leg_revs/per_leg_branch returns a Cycler, no raise; circular
   ephem; no velocity golden. `# COMPUTED`.
2. GREEN: `_build_chain(x, sequence, ephem, *, per_leg_revs=None,
   per_leg_branch=None)` — None → all-single-rev (must leave the construct_cycler
   call byte-identical to today for the Aldrin path).
3. RED `test_maintain.py::test_optimise_maintenance_dv_multirev_param_threaded` —
   3-leg sequence, per_leg_revs=(0,1,0), per_leg_branch=('single','low','single')
   on circular ephem; `converged=True`, finite maintenance_dv_kms; never asserts a
   sourced V_inf. `# COMPUTED`.
4. GREEN: extend `_objective(x, sequence, ephem, per_leg_revs, per_leg_branch)`
   and `optimise_maintenance_dv(..., *, per_leg_revs=None, per_leg_branch=None,
   ...)`. **AUDIT EVERY scipy `args=` site** (`differential_evolution`,
   `minimize`) in `optimise_maintenance_dv` so the args tuple exactly matches the
   non-x params: `(sequence, ephem, per_leg_revs, per_leg_branch)`. A missed site
   is a silent wrong-revs bug.
5. RED `test_optimize_ephemeris.py::test_ephemeris_seed_bounds_multirev_floor` —
   Cell with per_leg_revs=(0,1); leg-1 lower bound ≥ 1-rev t_min floor;
   `bounds[1][0] > bounds[0][0]`. `# COMPUTED`.
6. GREEN: `_ephemeris_tof_seed_and_bounds` — per-leg lower bound
   `max(0.1*share, estimated_t_min(n_revs))` using Hohmann semi-major axis proxy
   (`t_min ≈ π√(a³/μ)` scaled by revs). Too high → no convergence; too low →
   LambertConvergenceError caught by objective penalty (acceptable).
7. RED `test_optimize_ephemeris.py::test_optimise_cell_ephemeris_threads_revs_to_build_chain`
   — Cell ('E','M','E') per_leg_revs=(0,0) matches current behaviour
   (`result.best_cycler.legs[0].n_revs == 0`); a spied `_build_chain` confirms
   max_revs_per_leg passed through. Circular ephem. `# COMPUTED`.
8. GREEN: `optimise_cell_ephemeris` extracts `cell.per_leg_revs` /
   `cell.per_leg_branch` and passes to `optimise_maintenance_dv`. Verify the N>3
   path: for a 4-encounter cell the vinf_targets_kms dict lookup over
   `cell.sequence` yields 4 entries; PhaseSignature requires
   `len(leg_durations_s) == len(bodies) - 1` (3 legs / 4 bodies — consistent).
9. RED `test_s1l1_idealised_multirev.py::test_s1l1_4enc_multirev_closes_to_published_vinf_anchors`
   — **SOURCED EXPECTED**: V_inf_E=5.65, V_inf_M=3.05 km/s (spec §9 / catalogue),
   tol ±0.4 km/s. Cell bodies=('E','M'), sequence=('E','M','E','E'), period_k=2,
   per_leg_revs=(0,0,1) or (0,0,2) — whichever the characterisation script names.
   `Ephemeris('circular')`. Marked `@pytest.mark.slow` and initially
   `@pytest.mark.xfail(strict=False)` pending topology. Flips to strict pass once
   topology pinned.
10. AUDIT: `optimise_aldrin_maintenance_dv` (@608) must pass None defaults so the
    Aldrin golden anchors (a=1.60/e=0.393, 84°/72°) stay bit-identical. Run
    `test_maintain.py` golden assertions after the change.
11. Re-point `test_s1l1_real_rediscovery.py` xfail to reference the new topology:
    "topology mismatch resolved by STAGE 1 multi-rev plumbing; flip to pass once
    topology confirmed". SOURCED anchors 5.65/3.05 remain the only targets.
12. Lint: `uv run ruff check src/cyclerfinder/search/maintain.py
    src/cyclerfinder/search/optimize.py`; `uv run mypy src`.

### Risks (honest)
- scipy `args=` positional coupling — the single most error-prone change; a
  missed site silently corrupts revs and is invisible without the multi-rev
  golden test. Audit all DE + SLSQP sites in lockstep.
- Aldrin bit-reproducibility — None defaults must leave construct_cycler call
  identical. Verify against golden anchors.
- Multi-rev ToF floor calibration — too high blocks convergence.
- Census ratchet: `test_catalogue_rediscovery.py` EXPECTED_COVERAGE
  MULTI_ENCOUNTER_SEQUENCE:202 must remain 202 (gauntlet only builds 3-encounter
  cells; S1L1 covered by the dedicated test). Do not touch the gauntlet loader.

---

## STAGE 2 — Inclined analytic ephemeris (Track B, FIRST)

### What it unlocks
Adds orbital inclination + longitude of ascending node to the analytic
`_CircularBackend` so Venus (~3.39°) and Mars (~1.85°) have real 3D states. This
is the **VEM 3D foundation**: it lifts the coplanar/idealised limitation that
cannot host real cyclers (the 154-d E→M leg is near-hyperbolic in-plane). Also
the prerequisite for STAGE 4's 3D Tisserand (PlanetData fields).

### Size / partial-ness
Small, purely additive, low risk. The i=0 path stays byte-identical via an exact
`inc_deg == 0.0` short-circuit. No public interface change. Astropy backend
untouched.

### Files
- Modify: `src/cyclerfinder/core/constants.py` (`PlanetData` @78 — add
  `inc_deg: float = 0.0`, `lan_deg: float = 0.0`; set values in PLANETS dict)
- Modify: `src/cyclerfinder/core/ephemeris.py` (`_CircularBackend` @72 — delegate
  to new private `_InclinedCircularBackend` when `inc_deg != 0.0`)
- New: `tests/core/test_ephemeris_inclined.py`

### TDD task list
1. RED `test_inclined_circular_i0_byte_identical_to_flat` — for inc_deg=0 the
   inclined path equals the flat path by exact numpy array equality at
   t=0/73d/365d. Guards the short-circuit. `# INVARIANT`.
2. RED `test_existing_planar_tests_unchanged` — re-run `test_planar` semantics:
   r[2]==0, v[2]==0 for i=0 bodies. Regression guard. `# INVARIANT`.
3. GREEN: add `PlanetData.inc_deg`/`lan_deg` (zero-default, **after** existing
   fields per 3.11 dataclass ordering). Cite Standish & Williams inline. Confirm
   PLANETS is the only construction site and uses keyword args (verified).
4. GREEN: `_InclinedCircularBackend.state(body, t_sec) -> (Vec3, Vec3)`; rotate
   the in-plane state by `R_z(-lan) @ R_x(-inc)`. `_CircularBackend.state`
   delegates only when `inc_deg != 0.0` (exact float compare, NOT near-zero).
5. RED+GREEN sourced/invariant suite:
   - `test_venus_inclined_z_nonzero` — inc=3.39467605°, lan=76.67984255°
     (SOURCED); r[2]≠0, |r[2]|≈sma·sin(inc)≈3.0e6 km at quarter-period.
   - `test_mars_inclined_z_nonzero` — inc=1.84969142°, lan=49.55953891°
     (SOURCED); r[2]≠0.
   - `test_inclined_circular_speed_preserved` — |v|==√(μ_sun/|r|). `# INVARIANT`.
   - `test_inclined_circular_angular_momentum_direction` — h=r×v ∥
     n_hat=(−sin(lan)sin(inc), cos(lan)sin(inc), cos(inc)) within 1e-12. Catches
     wrong rotation sign/order. `# INVARIANT`.
   - `test_inclined_circular_period_closes` — within 1 km of t=0 after one
     period. `# INVARIANT`.
   - `test_inclined_venus_node_in_ecliptic_plane` — at node crossing
     (longitude==lan) r[2]==0 to ~1e-3 km. SOURCED geometric consequence.
   - `test_ephemeris_model_property_still_circular` and
     `test_astropy_backend_unaffected` — smoke / interface guards.
6. Lint: `uv run ruff check src/cyclerfinder/core/constants.py
   src/cyclerfinder/core/ephemeris.py`; `uv run mypy src`.

### Risks (honest)
- i=0 short-circuit MUST be exact float compare — any tolerance silently shifts z
  for a near-zero body.
- Rotation order/sign: `R_z(-lan) @ R_x(-inc)` applied to the in-plane state;
  wrong convention flips h. The angular-momentum test is the guard.
- Do NOT add inc/lan to the astropy backend — it already reads DE440 3D states.
- frames.py / propagate.py / lambert.py need no change (all 3D-native; synodic
  frame uses ecliptic-longitude atan2(x,y) projection which is correct with z≠0).

---

## STAGE 4 — 3D Tisserand feasibility (Track B, SECOND)

### What it unlocks
Extends the Tisserand linkability predicate from coplanar to inclined orbits
(`linkable_3d`), routing `tisserand_feasible` through it when `ephem is not
None`. Completes the **VEM 3D pruning** path. For V/E/M (inc < 3.5°) the extra
pruning beyond the coplanar predicate is negligible — value is realised for
genuinely inclined spacecraft (larger `i_sc_max_deg`).

### Size / partial-ness
Medium, conservative first step. Honest limitation: the 3D routing slows the
`tisserand_feasible` tight loop ~30×/pair; the coplanar-True short-circuit is the
mitigation. **Bounded first step**, not a full 3D Tisserand.

### Files
- Modify: `src/cyclerfinder/core/constants.py` (`PlanetData` @78 — add
  `ecc: float`; populate Venus/Earth/Mars). **Shared with STAGE 2 — see conflict
  map.**
- Modify: `src/cyclerfinder/search/tisserand.py` (add `linkable_3d`; update the
  "do not silently extend to 3-D" module docstring @20–28)
- Modify: `src/cyclerfinder/search/sequence.py` (`tisserand_feasible` @265 —
  replace `del ephem` with 3D routing when `ephem is not None`)
- New: `tests/search/test_tisserand_3d.py`

### TDD task list
1. GREEN (constants): add `PlanetData.ecc` with sourced inline values — Venus
   0.00677727, Earth 0.01671123, Mars 0.09340065 (Standish & Williams). If STAGE
   2 already added inc_deg, append ecc in the same field block (single merge).
2. RED `tests/search/test_tisserand_3d.py`:
   - `test_linkable_3d_em_at_5_5_kms_true` — `linkable_3d('E','M',5.5,
     i_sc_max_deg=0.01)` True (coplanar limit agrees with 2D at Aldrin band).
     `# COMPUTED`.
   - `test_linkable_3d_em_at_0_5_kms_false` — 0.5 km/s, i_sc_max=30° → False.
     `# COMPUTED` physics.
   - `test_linkable_3d_never_raises` — vinf in (0.0, 1e-6, 1e6) → False, no
     raise. Contract.
   - `test_linkable_3d_symmetric` — `linkable_3d('E','M',v)==linkable_3d('M',
     'E',v)` for v in (3.0,5.5,8.0). `# COMPUTED`.
   - `test_linkable_3d_coplanar_agrees_with_2d_for_em` **[GATE]** — for
     i_sc_max=0.01, `linkable_3d==linkable` over linspace(0.5,12.0,40). Critical
     regression guard. `# COMPUTED`.
   - `test_linkable_3d_inclined_orbit_extends_or_preserves_linkable_range` —
     monotonicity: `linkable` True ⇒ `linkable_3d(...,i_sc_max=30)` True.
     `# COMPUTED`.
   - `test_tisserand_feasible_routes_to_3d_when_ephem_provided` —
     `tisserand_feasible(em_cell, vinf_cap=8.0, ephem=Ephemeris('circular'))`
     True for Aldrin-band cell. `# COMPUTED` routing.
   - `test_tisserand_feasible_coplanar_default_unchanged` **[REGRESSION]** —
     ephem=None path identical to pre-STAGE4 (same gate cells as test_sequence.py).
3. GREEN: `linkable_3d(body_a, body_b, vinf_kms, *, i_sc_max_deg=30.0,
   tol_cos_i=1.0e-3, a_range_au=(0.3,5.0), n_points=80, n_a_points=80) -> bool`.
   At fixed V_inf, derive cos(i_sc) from body A's Tisserand eqn at each (a,e);
   evaluate body B's eqn at the same (a,e,i_sc); 2D grid scan over (a,e) tests
   agreement within tol_cos_i. Reuse `_orbit_crosses_planet` for both bodies
   (perihelion ≤ a_p ≤ aphelion — frame-independent). Formula (sourced — Strange
   & Longuski 2002): `T_p = a_p/a + 2cos(i)√((a/a_p)(1−e²))`,
   `V_inf² = (μ_sun/a_p)(3 − T_p)`.
4. GREEN: `tisserand_feasible(cell, vinf_cap, ephem=None)` — keep coplanar path
   for ephem=None; when `ephem is not None`, **short-circuit: if
   `linkable(a,b,vinf)` True return True immediately** (coplanar True ⇒ 3D True);
   only fall through to `linkable_3d` for the False cases. Update module docstring
   to sanction `linkable_3d` as the 3D extension while `linkable` stays coplanar.
5. Lint: `uv run ruff check src/cyclerfinder/core/constants.py
   src/cyclerfinder/search/tisserand.py src/cyclerfinder/search/sequence.py`;
   `uv run mypy src`.

### Risks (honest)
- Performance: O(n_a·n_e)≈6400/call vs ~200; ~30×/pair in the tight loop.
  Coplanar-True short-circuit limits the expensive scan to rescue-only cases.
- Negligible extra pruning for V/E/M (inc < 3.5° ⇒ cos i within 0.18%, inside the
  0.01 AU tolerance). Value appears only for larger i_sc_max_deg. Honest: this is
  a conservative extension, not a search-space reducer for the current body set.
- `tol_cos_i` calibration — derived ~0.005 from 0.01 AU propagation; default 1e-3
  conservative; document as requiring calibration against the coplanar GATE test.
- ecc/inc values MUST cite Standish & Williams — never computed from astropy
  (circular for a constants module).

---

## File / conflict map across stages

| File | STAGE 3 | STAGE 1 | STAGE 2 | STAGE 4 | Conflict handling |
|---|---|---|---|---|---|
| `search/phase_match.py` | add `leg_duration_seeds`, `find_candidate_windows`, `__all__` | — | — | — | none |
| `verify/real_closure.py` | swap `_resolve_real_t_start` body | — | — | — | none |
| `search/optimize.py` | add `_resolve_t0_multi_seed`, swap call @1309 | `_ephemeris_tof_seed_and_bounds` floor, `optimise_cell_ephemeris` thread revs | — | — | **Same file, both Track-A stages.** STAGE 3 first; STAGE 1 edits different functions. Sequence the two; no parallel edits to optimize.py. |
| `search/maintain.py` | — | `_build_chain`, `_objective`, `optimise_maintenance_dv`, audit Aldrin wrapper | — | — | none |
| `core/constants.py` | — | — | add `inc_deg`, `lan_deg` to PlanetData | add `ecc` to PlanetData | **ONLY cross-track shared file.** Build STAGE 2 before STAGE 4; both add fields to the same PlanetData block + PLANETS dict. Merge the field additions in one block; keep Standish & Williams citation single. |
| `core/ephemeris.py` | — | — | add `_InclinedCircularBackend`, delegate in `_CircularBackend` | — | none |
| `search/tisserand.py` | — | — | — | add `linkable_3d`, update docstring | none |
| `search/sequence.py` | — | — | — | route `tisserand_feasible` to 3D | none |
| `tests/search/test_s1l1_real_rediscovery.py` | xfail reason → "topology confirmed" | xfail reason → "topology resolved by STAGE 1 plumbing; flip once confirmed" | — | — | **Both Track-A stages edit the xfail reason.** STAGE 1 supersedes STAGE 3's edit. Apply in order. |

New test files (no conflicts): `tests/search/test_phase_match_robust.py` (S3),
`tests/search/test_s1l1_idealised_multirev.py` (S1),
`tests/core/test_ephemeris_inclined.py` (S2), `tests/search/test_tisserand_3d.py`
(S4).

### Cross-stage conflicts found (summary)
1. **`core/constants.py` PlanetData — the one genuine cross-track conflict.**
   STAGE 2 adds `inc_deg`/`lan_deg`; STAGE 4 adds `ecc`. Both touch the same
   frozen-dataclass field list and the same PLANETS records. Resolution: build
   STAGE 2 first; STAGE 4 appends `ecc` to the existing block. Field ordering:
   all new fields have defaults and follow non-default fields (Python 3.11). One
   shared Standish & Williams citation.
2. **`search/optimize.py` — intra-track ordering.** STAGE 3 and STAGE 1 both edit
   this file but different functions. Must be sequenced (S3 → S1), not edited in
   parallel.
3. **`test_s1l1_real_rediscovery.py` xfail reason — sequential overwrite.** S3
   sets "topology confirmed"; S1 overwrites with the plumbing-resolved reason.
   Harmless if applied in order.

No other file is touched by more than one stage. Track A and Track B are
otherwise fully independent and safe to develop in parallel up to the
`constants.py` merge point.

---

## What each stage unlocks (rollup)
- **STAGE 3** → correct-basin epoch selection; precondition for any real-ephemeris
  family search. Confirms (does not fix) the S1L1 topology diagnosis.
- **STAGE 1** → multi-rev/branch plumbing end-to-end; makes the ~200 SnLm rows and
  Hollister-15 reachable on real ephemeris; provides the machinery for S1L1
  closure once topology is pinned.
- **STAGE 2** → inclined Venus/Mars 3D states; the VEM 3D foundation that lifts
  the coplanar limitation hosting real cyclers.
- **STAGE 4** → 3D Tisserand pruning (bounded first step); completes the VEM 3D
  feasibility path; conservative for the current low-inclination body set.

## Honest large/partial flags
- STAGE 1 is the largest and most error-prone (scipy args coupling, Aldrin
  bit-reproducibility) and is **likely partial** — the S1L1 idealised test stays
  xfail until the topology characterisation script names the winning
  per_leg_revs/topology.
- S1L1 closure is **not** delivered by STAGE 3 alone (topology), and only
  potentially by STAGE 1 once topology is confirmed — treat S1L1 as an open item
  spanning both Track-A stages.
- STAGE 4 delivers correct 3D machinery but **negligible practical pruning** for
  V/E/M at current inclinations; its value is forward-looking (inclined spacecraft
  orbits).
