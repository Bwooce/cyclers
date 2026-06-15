# Epoch-aware genome data model + closure framework (#297 / #289 Phase 1)

**Status:** landed 2026-06-16. Substrate only; no catalogue writeback.

## What Phase 1 delivers

A thin, dependency-free (in the sense of reusing every existing dynamical
primitive) substrate for the three epoch-locked classes admitted by schema
v4.7 (`quasi_cycler`, `precursor_mga`, `mga_tour`):

- `src/cyclerfinder/genome/epoch_aware_genome.py`
  - `EpochLockedTrajectory` (frozen dataclass) — encodes a candidate
    epoch-locked trajectory: body sequence, per-leg time-of-flight, sourced
    `|V_inf|` at each encounter, `launch_epoch_utc`, `validity_window_*`,
    `orbit_class`, `n_returns`, optional `inserts_into` (required for
    `precursor_mga`), and an optional per-encounter `periapsis_altitudes_km`
    override (so a published mission like Tito at 100 km Mars periapsis can be
    expressed without relaxing the global safe-altitude default of 300 km).
  - `EpochLockedClosure` (frozen dataclass) — result of closing a trajectory
    at a specific epoch: the worst-encounter `|V_inf|` residual, the
    worst-flyby ballistic-continuity Delta V, the per-leg Lambert solutions,
    the per-encounter `|V_inf|` (mean of in/out at intermediate bodies), an
    independent-cross-check residual, and the `converged` gate.
  - `close_epoch_locked(...)` — runs per-leg single-rev prograde Lambert
    against the real ephemeris, computes per-encounter `|V_inf|`, computes
    flyby ballistic-continuity Delta V at every intermediate body, and runs
    an independent universal-variable Kepler re-propagation cross-check on
    every leg.
  - `search_validity_window(...)` — walks `launch_epoch_utc` across a grid
    around the nominal launch and returns the converged closures.
  - `utc_to_tsec(...)` — helper matching the script's TDB-seconds-since-J2000
    convention.

- `tests/genome/test_epoch_aware_genome.py` — 12 tests, all green at 5.8 s:
  - Tito 2018 end-to-end closure within the row-documented <1.5% per-`|V_inf|`
    residual on DE440 vs DE421.
  - Tito Mars flyby ballistic-continuity Delta V below 0.05 km/s.
  - Tito independent Kepler-vs-Lambert cross-check below 0.1 km/s.
  - Class-invariant gates: strict `cycler` rejected; `precursor_mga`
    requires `inserts_into`; `n_returns` must be a positive integer;
    `leg_tofs` / `vinf_kms` shapes must match `sequence`; unknown body codes
    rejected.
  - Window search around Tito ±10 d at 1 d steps yields a non-empty,
    bounded set of converged epochs that includes the nominal launch.
  - Negative gate: a bogus `vinf_kms_at_encounters` tuple does not pass the
    closure gate.

Tito reproduction numbers (DE440):

| Quantity | Ours | Tito 2013 (DE421) | Residual |
|---|---|---|---|
| `|V_inf|` Earth depart | 6.226 km/s | 6.232 km/s | -0.10% |
| `|V_inf|` Mars (mean in/out) | 5.419 km/s | 5.417 km/s | +0.03% |
| `|V_inf|` Earth arrive | 8.900 km/s | 8.837 km/s | +0.71% |
| Mars flyby continuity | 0.0107 km/s | (ballistic) | — |
| Lambert/Kepler cross-check | ~5e-5 km/s | — | — |

The 0.71% Earth-arrival residual matches the catalogue row note: dominated by
the DE421->DE440 ephemeris difference accumulated over the 274-day return
arc, not a model error.

## What Phase 1 does NOT deliver

- No Tisserand-graph or Tisserand-Poincaré window enumerator. The closure
  framework can score a *given* `(sequence, leg_tofs_days,
  launch_epoch_utc)` triple; it does not yet propose candidates.
- No MGA-chain optimisation (no DSM placement, no leg-TOF root-finding).
- No precursor-MGA insertion-into-cycler matching (the
  `inserts_into` field is validated as non-empty but not yet referentially
  resolved against the catalogue).
- No V0-V5 gauntlet integration for the new classes. The taxonomy note
  sketches the per-class semantics (V1 within window, V2 collapses to
  "launch-to-terminal-epoch" for `mga_tour`); the gauntlet plumbing change
  comes later.
- No multi-revolution Lambert sweep. Tito and the early `mga_tour` targets
  are single-rev; this is a Phase-2 option, not a Phase-1 requirement.
- No deep-space maneuvers (DSMs) inside a leg. Phase 1 is ballistic-only.
  Petropoulos pump tours and many real precursor-MGA designs require DSMs;
  the data model intentionally has no DSM field yet — adding it is a
  separable extension that won't break existing serialisations.
- No catalogue writeback. Phase 1 is substrate + tests; admission of new
  `mga_tour` / `precursor_mga` / `quasi_cycler` rows happens later under the
  schema v4.7 + V0-V5 gauntlet rules.

## Phase 2 path of least resistance

The lowest-friction next deliverable that turns this substrate into a real
search is a **Tisserand-Poincaré graph epoch-windowed enumerator over
(Earth, Venus, Mars) `V_inf` tuples for `mga_tour` chains**. The Tisserand
graph (Strange-Russell 2007; Petropoulos-Longuski 2000) is purely geometric
— a node is an (encounter body, `V_inf`-magnitude) pair and an edge is a
ballistic leg whose flight time lands in the right Tisserand box — so the
chain-construction pass is fast (graph BFS) and the chain-validation pass
is exactly `close_epoch_locked` from Phase 1. Concrete first deliverable:

- `src/cyclerfinder/search/tisserand_mga_window.py` — `find_mga_chains(
  launch_window=(t0, t1), planet_set=("E","V","M"), max_legs=4,
  vinf_grid_kms=(...), tof_box_days_per_leg=...)`, returning a stream of
  `EpochLockedTrajectory` candidates the caller can score with
  `close_epoch_locked`.

Phase 3+: precursor MGA insertion matching (resolve `inserts_into` against
the live catalogue; verify the precursor's final `|V_inf|` matches the
cycler's insertion `|V_inf|`); quasi-cycler closure semantics (N-return
drift bound across `validity_window`); V0-V5 gauntlet wire-up; DSM
extension to the data model.

## Structural surprises in the existing stack

- The CR3BP / cyclers-only `Cycler` value type has no `launch_epoch` and no
  `validity_window`; the epoch-locked classes are a structurally distinct
  shape, so the cleanest move was a parallel dataclass family rather than
  bolting epoch fields onto `Cycler`. The `EpochLockedTrajectory`
  intentionally lives in `genome/` alongside `tulip` / `family_switch`
  rather than overloading `core/`.
- The Tito `safe_alt_km=300.0` km Mars default in
  `cyclerfinder.core.constants.PLANETS` is engineering-conservative for the
  generic flyby case; Tito 2013 published 100 km. Without a per-encounter
  override the test `flyby_continuity_max_dv_kms` would land at ~0.1 km/s
  instead of ~0.01 km/s — the trajectory is still ballistically feasible at
  the *published* periapsis, but a generic search will need a way to
  declare per-encounter altitudes when the literature value differs from
  the safe default. Hence `periapsis_altitudes_km` in
  `EpochLockedTrajectory`. Phase 2 enumerator should set this to `None`
  (use safe default) and only override when a row's `notes` records a
  published value.
- `lambert()` returns `list[LambertSolution]` with the single-rev branch at
  index 0 and optional multi-rev branches after it. The closure driver
  takes `sols[0]` so a future multi-rev option can be added without
  changing the function signature.
- The Lambert solver and the Kepler propagator share the universal-variable
  family of methods but solve different problems (BVP vs IVP) with
  different residual functions and Newton loops. Re-propagating
  `(r1, v1, tof)` and comparing to `r2` is a legitimate independent
  cross-check at the conic-solver layer; for a true integrator-class
  cross-check (per the orbit-closure discipline) Phase 4+ should add a
  DOP853 propagation in patched-conic form. Phase 1's cross-check catches
  the most common Lambert-stack failure mode (wrong branch / wrong
  prograde-sense) without the integrator-startup cost.

## References

- `docs/notes/2026-06-16-catalogue-scope-taxonomy.md` — schema v4.7 four-class
  taxonomy and V0-V5 per-class semantics.
- `docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md` — Axis 5
  build cost + first-IC scoping that motivates Phase 1.
- `docs/notes/2026-06-13-tito-maccallum-2018-free-return-reproduction.md` —
  the Tito reproduction note this Phase generalises.
- `scripts/tito_free_return_repro.py` — the reference closure script.
- `data/catalogue.yaml` row `tito-2018-mars-free-return` — the V0 sourced
  reference for the Phase-1 integration test.
