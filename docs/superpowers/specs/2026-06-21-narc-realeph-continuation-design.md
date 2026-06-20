# N-arc Real-Ephemeris Continuation Lane (Russell parent → DE440) — Design Spec

**Date:** 2026-06-21
**Status:** Design approved (brainstorming), pending implementation plan.
**Tracks:** #388 / #307. **Builds on:** the validated Russell ψ generic-return
generator (`2026-06-21-russell-psi-generic-return-generator-design.md`).

## Goal

Continue the validated Russell idealized generic-return parents out to the real
DE440 ephemeris via Russell's §5.4 homotopy, using the **N-arc** corrector — supplying
the in-basin multi-arc seed that prior single-ellipse continuation attempts lacked —
then attempt V0→V1 promotion on the catalogue SnLm rows. A genuine real-ephemeris
closure (emerged V∞ within the campaign gate of the sourced anchor, bend-feasible) is
a *proposed* promotion, **held** for review; no catalogue writeback.

## Why now (what changed)

Prior real-eph closure of these cyclers failed on **seed quality** — the seed was
off-basin / single-ellipse (memory `project_s1l1_realeph_closure_blocker`:
"modelling stack complete; S1L1 fails on family-selection (off-basin), not infra";
and the #388 spikes showed the DSM-Lambert lane doesn't even close on the circular
model). The new Russell generator now produces a **true, golden-validated multi-arc
circular-coplanar parent** in the correct basin — the missing input. This lane bridges
that parent into the existing real-eph continuation machinery.

## Reuse (no rebuild)

- `search/continuation.py::ramped_ephemeris(λ_e, λ_i)` — the Russell Fig-5.3 homotopy
  backend (ramps J2000-mean eccentricity then inclination; `a,Ω,ω,ν` pinned) — and its
  `LADDER` rung schedule.
- `search/correct.py::ballistic_correct(...)` — the N-arc real-eph multiple-shooting
  corrector (free vars `[t0, leg ToFs]`, one leg pinned by period; flyby V∞-continuity
  + closure residuals). Use `residual_mode="vector"` (bend-feasibility inside the
  residual → steers to the ballistic families Russell's AR/TR gate selects).
- `search/generic_return.py`, `search/cycler_assembly.py` — the Russell parent
  (`descriptor_to_phsi`, `assemble_cycler`, `Cycler`, `RussellModel`).
- Existing epoch-selection (`scripts/campaign_russell12.py`, `run_365_*.py` compute a
  `best_t0`) — reuse if it implements the phase-angle method below; else implement it
  in this lane.

## Architecture — one new module `search/narc_continuation.py` + a batch script

### Component 1 — the bridge: Russell parent → N-arc real-eph seed
`russell_parent_to_ballistic_seed(model, cycler, row) -> NarcSeed`. Converts the
idealized canonical `Cycler` to `ballistic_correct` inputs:
- `sequence` ← `row["sequence_canonical"].split("-")` (e.g. `E-E-M-M`).
- `per_leg_revs`, `per_leg_branch` ← the cycler's generic-return rev/branch per leg
  (transit legs single-rev; same-body resonant legs carry the generic return's
  `n_revs`/`branch`).
- `tof_seed_days` ← per-leg ToFs converted **canonical TU → days** (`× model.tu_days`,
  58.1324409).
- `period_sec` ← `p × REAL Earth-Mars synodic` (2.1354 yr = 1/(1/1.0 − 1/1.8808);
  the **real** synodic, NOT Russell's idealized 1.875-yr-Mars value — real-eph uses
  real periods). Express in seconds.
- `vinf_anchor_e/m_kms` ← the row's sourced Russell-table V∞ cells (for the held
  anchor-match check; emerged, never imposed).
`NarcSeed` is a frozen dataclass holding these.

### Component 2 — phase-angle epoch derivation (Russell §5.3, sourced)
`candidate_epochs(seed, *, launch_window_synodics=range(1, 22), grid=100) -> list[float]`:
1. Compute the parent's **beginning Earth-Mars relative phase angle** (from the
   idealized parent's departure geometry).
2. Over **LaunchWindow ∈ 1..21 synodic periods after J2000** (Russell's range — spans
   ~3 of the ~7-synodic true-repeat cycles), find the real-ephemeris epochs where
   Earth and Mars sit at that relative phase (propagate until the phase is achieved).
3. Refine each with Russell's **1-D grid search over 100 equally-spaced intervals**
   within the local synodic window (his density, not invented), scoring by the
   homotopy's final residual.
Returns the candidate `t0_seed_sec` list (J2000-relative), best-first.

### Component 3 — the N-arc continuation driver
`narc_continuation_correct(seed, *, ladder=LADDER, final_ephemeris=DE440, residual_mode="vector") -> NarcContinuationResult`:
For each candidate epoch (Component 2): walk the ramp — for each `(λ_e, λ_i)` rung in
the schedule call `ballistic_correct(..., ephem=ramped_ephemeris(λ_e, λ_i))`, re-seed
the next rung from the prior rung's converged `(t0, tofs)`; final step → DE440 (the
seven-cycle propagation is the implicit multi-rev initial guess via the corrector's
period closure). Keep the lowest-residual converged DE440 result across epochs and
rungs. `NarcContinuationResult` carries: converged, max_residual_kms, emerged
per-encounter V∞, the converged `(t0, tofs)`, bend-feasible flag, the winning epoch.

### Component 4 — batch + held report `scripts/narc_continuation_batch.py`
For each descriptor-bearing row → `descriptor_to_phsi` → `assemble_cycler` → bridge →
`narc_continuation_correct` on DE440. Record per row: converged / emerged real-eph V∞
(E and M) vs the sourced Russell anchor (0.5 km/s campaign gate) / bend-feasible / max
residual / winning epoch. Flag the single V0 target `mcconaghy-2006-em-k2` as a
**proposed V0→V1 (held)** iff converged AND both anchors match AND bend-feasible.
Runlog `data/runs/narc-continuation-<ts>.jsonl`. NO catalogue writeback.

## Data flow

row → Russell parent (idealized, canonical) → bridge (canonical→real units, real
synodic period) → phase-angle epoch derivation (LaunchWindow 1..21, 100-grid refine)
→ homotopy ramp (ballistic_correct vector-mode per rung) → DE440 closure → emerged V∞
vs sourced anchor → held report.

## Error handling

- A rung / epoch that fails to converge drops that epoch sample (recorded, not crashed).
- If NO epoch closes for a row on DE440, that is a **recorded honest negative** — and a
  meaningfully different one from the prior off-basin failure: it means a true Russell
  in-basin parent reached DE440 and still didn't close (a real, well-characterized
  terminal result for this lane, not a seed defect).
- Hyperbolic / over-VINF_CEILING emerged V∞ → loud flag, reject.
- `assemble_cycler`/`descriptor_to_phsi` returning None (ocampo rows) → out-of-scope,
  recorded, skipped.

## Testing (TDD)

- **bridge:** `russell_parent_to_ballistic_seed` converts a known `Cycler`: ToF in days
  (`× 58.1324409`), `period_sec` from the REAL synodic (2.1354 yr × p), v∞ anchors in
  km/s; sequence/revs/branch wired per leg.
- **epoch:** `candidate_epochs` returns epochs whose Earth-Mars relative phase matches
  the parent's to tolerance; respects LaunchWindow 1..21 and the 100-interval grid.
- **driver smoke:** `narc_continuation_correct` runs end-to-end on a ramped backend at
  `(λ_e, λ_i) = (0, 0)` (= circular) and reproduces the idealized closure (residual ~0)
  — isolates the ramp/corrector plumbing from the DE440 swap.
- **batch (empirical):** the catalogue run is the result, not a unit assertion.
- **golden honesty:** EXPECTED V∞ = the sourced Russell anchor; emerged not imposed; no
  tolerance loosening; no writeback.

## Honesty gates (orbit-closure-discipline — non-negotiable)

1. Sourced V∞ anchor is the comparison target, emerged from the converged DE440 result,
   never imposed.
2. Independent cross-check: bend-feasibility (vector residual) + the V∞ cap sanity gate.
3. No tolerance / cap loosening to force a closure.
4. No catalogue writeback; any V0→V1 is **held** for session review.
5. A clean negative (true parent reaches DE440 but doesn't close) is a success.

## Out of scope / deferred

- Out-of-plane / 3D parent (the parent is the coplanar generator; #414). The ramp's
  inclination step still walks toward real (small) inclinations from the coplanar seed,
  as Russell does — but the PARENT is coplanar.
- SNOPT elastic-mode (Russell's robustness trick for surviving intermediate
  infeasibility) — we use `scipy.least_squares` as `ballistic_correct` already does;
  the small homotopy steps substitute for elastic mode. If a row fails ONLY at an
  intermediate rung (not the DE440 endpoint), that is flagged as an elastic-mode gap,
  not silently dropped.
- Catalogue writeback (post-review only).

## References

- `docs/notes/2026-06-07-russell-2004-continuation-deepdive.md` — §5.3/§5.4 epoch
  (phase-angle + LaunchWindow 1..21 + 100-interval grid + seven-cycle propagation),
  the homotopy schedule, the 5n multiple-shooting method.
- `src/cyclerfinder/search/continuation.py` (`ramped_ephemeris`, `LADDER`),
  `search/correct.py` (`ballistic_correct`), `search/generic_return.py`,
  `search/cycler_assembly.py`.
- Memory: `project_s1l1_realeph_closure_blocker`, `project_dsm_closure_modeljump_blocker`,
  `feedback_golden_tests_sourced_only`, `feedback_orbit_closure_discipline`,
  `project_s1l1_nomenclature`.
