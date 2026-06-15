# Tisserand-Poincaré MGA-chain enumerator (#298 / #289 Phase 2)

**Status:** landed 2026-06-16. Search infrastructure only; no catalogue
writeback.

## What Phase 2 delivers

A pure-geometry MGA-chain enumerator that *proposes* candidates for the
Phase-1 `close_epoch_locked` driver to validate. The enumerator is the
search layer over the Tisserand-Poincaré graph (Strange & Longuski JSR
39(1):9-16, 2002; Campagnola & Russell JGCD 33(2):476-486, 2010); it owns
adjacency, BFS, TOF estimation, and ranking. It owns *nothing* dynamical
— every Lambert and every flyby check goes through Phase 1.

- `src/cyclerfinder/search/tisserand_mga_window.py`
  - `MGAChainCandidate` (frozen dataclass) — pure-geometry proposal:
    `sequence`, `vinf_tuple_kms`, `leg_tofs_days`, `launch_epoch_utc`,
    `tisserand_parameter`, `chain_score`.
  - `find_mga_chains(launch_window, planet_set, max_legs, vinf_grid_kms,
    tof_box_days_per_leg, epoch_step_days, chain_score_threshold,
    start_body_filter, a_range_au) -> Iterator[MGAChainCandidate]` — BFS
    over the (planet, V_inf-bin) graph. Hetero edges use
    `cyclerfinder.search.tisserand.linkable_3d` (3-D Tisserand admissibility,
    inclined to 30 deg by default). Same-body edges generate resonance
    bands `(M:N)` whose TOF lands inside `tof_box_days_per_leg`.
  - `validate_chain_candidate(candidate, ephemeris, ...) ->
    EpochLockedClosure | None` — wraps the proposal as an
    `EpochLockedTrajectory` and dispatches to Phase 1's
    `close_epoch_locked`. Returns `None` if any gate fails.
  - `scan_window_and_validate(...)` — convenience driver: enumerate +
    validate + return surviving pairs.

- `tests/search/test_tisserand_mga_window.py` — 10 tests, all green at
  ~4.8 s on the analytic backend (only the `@pytest.mark.slow` round-trip
  test exercises DE440):
  - Shape: candidates have the right `len(sequence) = len(vinf_tuple)`,
    `len(leg_tofs) = len(sequence) - 1`, every body is a valid PLANETS key,
    every TOF is positive.
  - V_inf conservation: every chain lives on a single V_inf bin
    (Tisserand-Poincaré invariant — pump tours that span multiple shells
    are Phase 3+).
  - Tisserand parameter equals `vinf_to_tisserand(seq[0], vinf[0])`.
  - `max_legs` bounds chain length.
  - `chain_score_threshold` caps output.
  - `start_body_filter` restricts the first body.
  - Input validation rejects bad inputs.
  - `validate_chain_candidate` returns `None` on impossibly-tight tolerance.
  - **Galileo VEEGA structural reproduction**: a (V, E, E, J) chain
    surfaces in the 1989-10-01..1989-11-15 window.
  - `scan_window_and_validate` returns equal-length candidate / closure
    lists with consistent sequences (slow, DE440).

## Galileo VEEGA reproduction verdict — STRUCTURAL only

**Sequence + launch-window: YES (4 candidates inside ±4 weeks of
1989-10-18).**

The enumerator surfaces (V, E, E, J) chains at V_inf = 10 km/s and
11 km/s, with launch epochs Oct 1 / Oct 16 / Oct 31 / Nov 15 — every grid
point is within ±4 weeks of Galileo's published 1989-10-18T16:53:40 UTC
liftoff (Diehl-Belbruno-Roberts 1986 AAS).

**Closure residual: FAR ABOVE the publishable gate.** The per-encounter
V_inf at the Phase-2 geometric proposal lands at (≈7, ≈24, ≈40, ≈10) km/s
vs Galileo's published (5.3, 8.93, 8.93, 5.6) km/s. Worst residual is
28-34 km/s; flyby continuity ΔV is 30-50 km/s. This is **NOT a bug** —
it's the *honest* calibration gap the spec asked us to document.

**Root cause: V_inf conservation is exact across one flyby, NOT across a
chain.** The Tisserand-Poincaré graph is a single-V_inf-shell pre-screen.
Galileo VEEGA is structurally a **multi-shell pump tour**: V_inf grows
from ~4 km/s at launch to ~8.9 km/s by Earth-2 (the pump segments) and
back to ~5.6 km/s at Jupiter (the Jupiter approach). The single-shell
enumerator surfaces VEEJ at V_inf = 10 km/s because that's the lowest
shell at which `linkable_3d(E, J)` is True at the default
`a_range_au=(0.3, 6.0)`; widening to `(0.3, 8.0)` admits V_inf = 8.

**This is a *structural* reproduction, NOT a closure-grade reproduction.**
The structural test asserts only that a (V, E, E, J) sequence exists with
a launch epoch in the window; the closure-grade reproduction is Phase 3's
job (DSM + per-leg TOF optimisation).

The data probe `scripts/probe_298_galileo_veega.py` writes 8 rows to
`data/scan_298_galileo_veega.jsonl` recording:

- `sequence`: `["V","E","E","J"]`
- `vinf_tuple_kms`: the bin V_inf (10 or 11 km/s)
- `leg_tofs_days`: geometric Hohmann-half-period proposals
  `[146, 365, 998]` (V→E, E→E synodic-style return, E→J)
- `closure_residual_kms`: 28.5 - 34.4 km/s (PHASE-3 GAP)
- `flyby_continuity_max_dv_kms`: 30-50 km/s (PHASE-3 GAP)
- `per_encounter_vinf_kms`: `[7, 24, 40, 10]` at the proposal

The JSONL is **DIAGNOSTIC, NOT SILVER**. None of these candidates is
catalogue-eligible; the file documents the seed-to-closure gap so Phase 3
can target the right loss function.

## Candidate-stream funnel at a Galileo-scale sweep

```
launch_window:  (1989-10-01, 1989-11-15)   span = 45 days
planet_set:     ("V","E","J")              3 bodies
vinf_grid_kms:  (8.0, 9.0, 10.0, 11.0)     4 shells
max_legs:       4
epoch_step:     15 days                    4 epochs

ENUMERATED              : 2 512   (BFS over Tisserand-admissible chains)
V_inf-coherent          : 2 512   (every chain on one shell by construction)
(V,E,E,J)               :     8   (sequence-filtered)
CLOSURE-attempted       :     8   (all 8 sent through close_epoch_locked)
CLOSURE-passing 0.5km/s :     0   (single-shell seed; multi-shell tour)
CLOSURE-passing 50km/s  :     0   (flyby ΔV > 50)
CLOSURE-passing 500km/s :     8   (diagnostic recording only)
literature-fresh        :     0   (Galileo is the archetype)
```

The funnel says the Phase-2 enumerator is *aggressive* (it emits 2 512
candidates in ~1 s) and *correctly biased* (it surfaces the canonical
archetype) but *seed-only* — closure-grade trajectories require Phase 3's
DSM + per-leg-TOF optimisation around the seed.

## What Phase 2 does NOT deliver

- **No multi-V_inf-shell pump tours.** The Tisserand-Poincaré graph is
  single-shell by construction. Pump tours (Petropoulos-Longuski 2000)
  walk across shells via Earth-Earth resonant returns; surfacing them
  requires a shell-walking BFS, which is a Phase 3 deliverable.
- **No DSMs.** Phase 2 is ballistic-only. The data model accommodates
  DSMs (the `EpochLockedTrajectory` will gain a separable DSM extension
  in Phase 3 without breaking serialisation).
- **No per-leg TOF optimisation around the seed.** The Phase-2 TOF
  estimate is *geometric* (Hohmann half-period for hetero hops; resonance
  ratio for same-body hops); the actual flight time at the published
  launch epoch is found by Phase 1's Lambert, but the seed is rarely
  close enough for the Lambert to converge to a publishable residual.
- **No catalogue writeback.** SILVER candidates would still need to pass
  the V0-V5 gauntlet (#274) before catalogue admission. The Phase-2
  output is a discovery JSONL, not a row.
- **No novelty claim.** `cyclerfinder.search.literature_check` (#272 /
  KNOWN_CORPUS at corpus rev `568d8a4`) is the mandatory gate.

## Phase 3 hand-off

The Phase-3 deliverable should resolve the Galileo VEEGA seed gap, which
fixes both the single-shell limitation AND the geometric-TOF seed
mismatch in one pass:

1. **Multi-shell BFS.** Generalise `_enumerate_sequences` to allow each
   intermediate flyby to change V_inf within a small per-flyby budget
   (the pump increment). The pump budget is bounded by the body's
   gravity-assist envelope (Strange-Longuski 2002 eq. 12); for Earth at
   V_inf = 9 km/s the maximum per-flyby V_inf change is ~2.7 km/s, which
   matches the (4 km/s at launch → 8.9 km/s at Earth-1 → 8.9 km/s at
   Earth-2) Galileo pattern.

2. **Per-leg TOF optimisation.** Around each seed (sequence + launch
   epoch + per-leg TOF estimate), run a small Nelder-Mead / SLSQP
   refinement over (`leg_tofs_days[k]`, `launch_epoch_utc`) minimising
   `closure_residual_kms + alpha * flyby_continuity_max_dv_kms`. The
   Phase-1 `close_epoch_locked` is the loss function; the search space is
   bounded by `tof_box_days_per_leg` and `validity_window`.

3. **DSM extension to `EpochLockedTrajectory`.** Add an optional
   `dsm_specs` field (one per leg, defaulting to None = no DSM); the
   closure driver picks up the DSM contribution from
   `cyclerfinder.search.dsm_leg`. Petropoulos pump tours and most real
   precursor-MGA designs need DSMs; the data model is currently
   ballistic-only.

4. **Precursor-MGA insertion matching.** Resolve `inserts_into` against
   the live catalogue: a `precursor_mga` candidate matches a cycler row
   if its terminal (`vinf_kms_at_encounters[-1]`,
   `sequence[-1]`) is within tolerance of the cycler's seed `|V_inf|` and
   anchor body. Use `cyclerfinder.search.literature_check` to gate
   "first-found" claims for the chain pair.

5. **V0-V5 gauntlet wire-up for `mga_tour`.** Once the seed → closure
   loop is grade-A, plumb the surviving candidates through the existing
   V0-V5 gauntlet (#274) — `mga_tour` collapses V2 to "launch-to-terminal
   epoch" per the schema-v4.7 taxonomy note.

## References

- `docs/notes/2026-06-16-297-289-phase1-epoch-aware-genome.md` — Phase 1
  data model + closure driver.
- `docs/notes/2026-06-16-catalogue-scope-taxonomy.md` — schema v4.7
  four-class taxonomy.
- `src/cyclerfinder/search/tisserand.py` — coplanar + 3-D Tisserand
  predicates inherited by Phase 2.
- `data/scan_298_galileo_veega.jsonl` — diagnostic probe output (8 rows).
- Strange, N. J. & Longuski, J. M., "Graphical Method for Gravity-Assist
  Trajectory Design", *J. Spacecraft and Rockets*, 39(1):9-16, 2002.
- Campagnola, S. & Russell, R. P., "Endgame Problem Part 2: Multibody
  Technique and the Tisserand-Poincaré Graph", *JGCD*, 33(2):476-486,
  2010.
- Petropoulos, A. E. & Longuski, J. M., "Shape-Based Algorithm for the
  Automated Design of Low-Thrust, Gravity Assist Trajectories",
  AAS/AIAA Astrodynamics Specialist Conf., 2000.
- Diehl, R. E., Belbruno, E. A., Roberts, R. M., "Galileo Mission Plan
  for the (V, E, E, J) Trajectory", AAS 1986 (KNOWN_CORPUS at corpus rev
  `568d8a4`).
