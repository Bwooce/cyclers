# Seed the full multiple-shooting shooter from the Russell parent — Design Spec

**Date:** 2026-06-21
**Status:** Design approved (brainstorming), pending implementation plan.
**Tracks:** #388 / #307 / #135. **Builds on:** the Russell ψ generic-return generator
(`2026-06-21-russell-psi-generic-return-generator-design.md`) and the N-arc
continuation lane (`2026-06-21-narc-realeph-continuation-design.md`).

## Goal

Feed the existing full-state multiple-shooting corrector (`nbody/shooter.py`) the
**constructed, golden-validated Russell ψ generic-return parent** as its conic seed —
the in-basin seed the #135 verdict said the shooter needs but never had (it was only
ever fed the searched #133 near-miss survey) — and run it on the real n-body model.
This is the decisive test of whether the literal Russell parent breaks the
seeding/family-selection wall that has blocked real-eph closure of the SnLm cyclers.
No new corrector — the shooter is Russell's §5.4.4 / Jones-SNOPT full-state method.

## Why (the localized blocker)

The #135 like-for-like diagnostic (recorded in `nbody/shooter.py`) found the full
shooter "closes geometrically but lands OFF-ANCHOR; verdict = seeding/basin, not
solver." The narc-continuation lane (`2026-06-21-narc-continuation-results.md`)
independently landed the same low-energy off-anchor basin, AND showed the target
model (DE440 vs Russell mean-element) is not the cause. The one lever never tried is
the seed: a CONSTRUCTED parent (the new ψ generator) rather than a searched near-miss.

## Scope

IN: a seed adapter (Russell `NarcSeed` → `shooter.ShootingSeed` via the existing
v∞-vector node machinery) + a batch that runs `shooter.shoot` on the descriptor-bearing
catalogue rows and reports held outcomes.

OUT: any change to `nbody/shooter.py`, `correct.py`, or the generator (pure reuse); the
near-miss survey path (unchanged); out-of-plane (#414); catalogue writeback.

## Architecture — one small module `search/shooter_russell_seed.py` + a batch script

### Component 1 — seed adapter
`russell_shooting_seed(narc_seed, t0_sec, ephem) -> ShootingSeed`:
1. Reconstruct the free-leg ToFs + slack leg from `narc_seed` (the `NarcSeed` from
   `narc_continuation.russell_parent_to_ballistic_seed`: `sequence`, `per_leg_revs`,
   `per_leg_branch`, `tof_seed_days`, `period_sec`). Pin the slack leg = longest ToF
   (the `correct.ballistic_correct` convention); the free ToFs are the rest.
2. `nodes = correct._vinf_nodes(sequence=…, per_leg_revs=…, per_leg_branch=…,
   t0_sec=t0_sec, free_tof_days=free_tofs, slack_leg=slack_leg,
   period_days=period_sec/DAY_S, ephem=ephem)` — the per-encounter v∞ **vector**
   mapping (`b{i}_in`/`b{i}_out`) that the shooter consumes (computed by Lambert per
   leg). This is the bridge the shooter never had from a constructed parent.
3. `return shooter.seed_from_conic(sequence=…, vinf_nodes=nodes, t0_sec=t0_sec,
   tofs_days=full_tofs, slack_leg=slack_leg, period_days=period_sec/DAY_S, ephem=ephem)`.
Raises `ValueError` (caught by the batch) if `_vinf_nodes` fails (degenerate Lambert).

### Component 2 — batch + held report `scripts/shooter_russell_batch.py`
For each descriptor-bearing catalogue row: `descriptor_to_phsi` → `assemble_cycler`
→ `russell_parent_to_ballistic_seed` → `candidate_epochs` (best K, default K=3, to
keep full n-body propagation tractable — log the cap) → for each epoch build the
shooting seed (Component 1) and run `shooter.shoot(seed, ephem=real n-body, …)` →
keep the lowest-defect `ShootResult`. Record per row: converged (defect < tol),
emerged n-body V∞ per encounter vs the sourced anchor (0.5 km/s gate, both E&M),
bend-feasible. Flag `mcconaghy-2006-em-k2` PROPOSED V0→V1 (HELD) iff converged AND
anchor-matched AND bend-feasible. Runlog `data/runs/shooter-russell-<ts>.jsonl`;
summary; "HELD — no writeback".

## Data flow

row → Russell parent (golden idealized) → NarcSeed → candidate epochs → `_vinf_nodes`
(v∞ vectors) → `seed_from_conic` (ShootingSeed) → `shooter.shoot` (full-state n-body
multiple-shooting) → emerged n-body V∞ vs sourced anchor → held report.

## Error handling

- `_vinf_nodes` / Lambert degenerate → adapter raises ValueError → batch records the
  row as un-seedable, continues.
- `descriptor_to_phsi` / `assemble_cycler` None (ocampo / unassemblable) → skip.
- `shoot()` non-convergence → recorded `ShootResult` with its defect; not a crash.
- Epoch cap (K) is logged so coverage is honest (no silent truncation).

## Testing (TDD)

- **adapter smoke:** `russell_shooting_seed` returns a `ShootingSeed` with one
  `(6,)` node state per encounter, correct epochs/tofs, for `russell-ch4-4.991gG2`.
- **shoot smoke:** `shooter.shoot` on that seed (small step / few nodes) returns a
  `ShootResult` without error (converged or not — the result, not closure, is the
  smoke assertion).
- **batch (empirical):** the catalogue run is the decisive result.
- **golden honesty:** emerged n-body V∞ compared to the SOURCED anchor; never imposed.

## Honesty gates (orbit-closure-discipline)

1. Sourced V∞ anchor is the comparison target; the shooter's n-body V∞ emerges.
2. No tolerance loosening; bend-feasibility is an independent check.
3. No catalogue writeback; any V0→V1 is held for session review.
4. A non-closure from the literal Russell parent is the decisive characterized
   negative — the strongest evidence yet that these cyclers are not landable by our
   stack from even the true parent (closing the long #388/#135/S1L1 thread), NOT a
   reason to loosen anything.

## References
- `nbody/shooter.py` (`ShootingSeed`, `seed_from_conic`, `shoot`, the #135 verdict),
  `search/correct.py` (`_vinf_nodes`), `search/narc_continuation.py`
  (`russell_parent_to_ballistic_seed`, `candidate_epochs`),
  `search/cycler_assembly.py` (`descriptor_to_phsi`, `assemble_cycler`).
- `docs/notes/2026-06-21-narc-continuation-results.md` (the seeding/family-selection
  diagnosis), `docs/notes/2026-06-06-russell12-likeforlike.md` (#135).
- Memory: `project_s1l1_realeph_closure_blocker`, `project_dsm_closure_modeljump_blocker`,
  `feedback_golden_tests_sourced_only`, `feedback_orbit_closure_discipline`,
  `feedback_never_give_up_reproducing_papers`.
