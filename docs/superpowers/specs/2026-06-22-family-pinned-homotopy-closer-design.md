# Family-pinned penalty homotopy closer — Design Spec

**Date:** 2026-06-22
**Status:** Design approved (brainstorming), pending implementation plan.
**Tracks:** #388 (SnLm V0→V1 real-ephemeris closure). **Builds on:** the STM
shooter (`nbody/shooter.py`, commits `8196f55`/`8cfe894`) and its decisive
literal-parent negative (`docs/notes/2026-06-21-shooter-stm-batch-results.md`).

## Goal

Test whether each SnLm Earth-Mars cycler row has a **ballistic family member at
the PUBLISHED V∞** by biasing the full STM multiple-shooting corrector into the
published-V∞ basin with a penalty, then **ramping the penalty to zero** so the
recorded V∞ emerges from an unpenalized ballistic solve. This is the
global / family-targeted lever the STM batch verdict pointed to: the local
corrector either fails to engage (3 rows: parent too far from any continuous
trajectory) or closes into the WRONG high-energy family (1 row, the #135
signature). Held; no catalogue writeback.

## Why this approach (and why golden-safe)

The obstruction is **family selection**, not solver performance (the STM made the
solve tractable and it still lands off-anchor) nor model fidelity (the conic N-arc
lane reached DE440 and the same off-anchor basin,
`2026-06-21-narc-continuation-results.md`). A V∞-anchor penalty biases the
corrector toward the published basin from the start — directly attacking family
selection — at low dimensional cost (one homotopy scalar `λv`), reusing the STM
shooter.

**Golden discipline (load-bearing, per `feedback_golden_tests_sourced_only` /
`feedback_orbit_closure_discipline`).** The penalty pulls toward the *sourced*
anchor but is **ramped to zero**; the recorded V∞ is read ONLY from the `λv=0`
unpenalized solve. The penalty is a basin-selector, never the recorded value.
- If the `λv=0` ballistic solve converges AND emerged V∞ is within 0.5 km/s of
  BOTH the E and M anchors AND it is bend-feasible → a ballistic family member
  exists at the published V∞ → **PROPOSED V0→V1, held, not applied**.
- If V∞ snaps off-anchor when the penalty lifts (or the `λv=0` solve fails to
  converge) → the published cycler is not ballistic in DE440 at that family → a
  **stronger characterized negative** than the current one (it now rules out the
  family-near-anchor basin explicitly, not just the literal-parent seed).

This is the right next experiment regardless of outcome: it converts "the parent
seed lands off-anchor" into "is there *any* ballistic member near the published
V∞?".

## Reuse (no rebuild)

- `nbody/shooter.py`: `defect_residual`, `_stm_jacobian`, `shoot`, `ShootingSeed`,
  `ShootResult`, `_node_vinf`, `_seed_with_states`, `_x_to_states`/`_states_to_x`.
  The penalty is an **opt-in augmentation** of the residual + STM Jacobian; the
  default (no anchors) path is byte-identical.
- `search/shooter_russell_seed.py::russell_shooting_seed` + `narc_continuation`
  (`russell_parent_to_ballistic_seed` carries `vinf_anchor_e_kms` /
  `vinf_anchor_m_kms`) — the seed + anchors, unchanged.
- `scripts/shooter_russell_batch.py` (detached/checkpointed/heartbeat) — extended
  to call the homotopy driver.
- `search/continuation.py::continuation_correct` — the **ladder + warm-start**
  pattern the homotopy driver mirrors (NOT its ephemeris ramp; this runs on DE440).

## Architecture — three units, each independently testable

### Unit 1 — Penalty-augmented residual + STM Jacobian (`nbody/shooter.py`)

A V∞-anchor penalty appended to `defect_residual`, behind an optional argument so
the existing path is unchanged.

- **API.** `defect_residual(..., vinf_anchors: Mapping[str, float] | None = None,
  vinf_weight: float = 0.0)`. `vinf_anchors` maps a body code → sourced |V∞|
  (km/s) (e.g. `{"E": 4.70, "M": 5.00}`). When `vinf_anchors` is None OR
  `vinf_weight == 0.0` the returned residual is **byte-identical** to today's.
- **Penalty rows.** For each encounter node `i` whose body has an anchor, append
  one scalar `sqrt(vinf_weight) * (|v∞_i| − anchor[body_i])`, where
  `|v∞_i| = ||state_i[3:] − v_planet(body_i, epoch_i)||`. All E-nodes use the E
  anchor, all M-nodes the M anchor (Russell cyclers carry one V∞ per body by
  Tisserand construction). Ordering: penalty rows are appended AFTER the existing
  leg-defect + hinge + wrap blocks, in node order, so the layout is deterministic.
- **STM Jacobian rows.** Extend `_stm_jacobian(..., vinf_anchors, vinf_weight)`
  with matching appended rows: `∂(sqrt(w)·|v∞_i|)/∂state_i = sqrt(w) · [0,0,0, v̂∞_i]`
  (zero on the position block, the V∞ unit vector on the velocity block of node
  `i`'s 6-state); zero on all other nodes. Degenerate `|v∞_i| < eps` → that row's
  Jacobian is zero (and the residual uses `|v∞_i|` directly), so a near-zero-V∞
  node does not produce a NaN gradient.
- **`shoot` wiring.** `shoot(..., vinf_anchors=None, vinf_weight=0.0)` threads both
  into `residual_of_x` and the STM `jac_of_x_stm`. FD path unaffected (it
  finite-differences the augmented residual automatically when anchors are set).

### Unit 2 — Homotopy driver (`search/family_pinned_shoot.py`, new)

`family_pinned_shoot(seed, *, ephem, bodies, vinf_anchors, weight_ladder,
accuracy, max_nfev, max_wall_sec, progress=None) -> FamilyPinnedResult`.

- **Ladder.** Iterate `vinf_weight` over `weight_ladder` (descending to 0, e.g.
  `(W, 0.6*W, 0.3*W, 0.1*W, 0.0)` with `W` chosen so the penalty is commensurate
  with the leg-defect scale — see Open Questions). Each rung calls
  `shoot(..., vinf_anchors=vinf_anchors, vinf_weight=w, jacobian="stm")`,
  **warm-started** from the previous rung's `corrected_states` (rebuild the seed
  via `_seed_with_states`). The final `w=0` rung is the verdict solve.
- **Result.** `FamilyPinnedResult` (frozen): the final (`w=0`) `ShootResult`, the
  per-rung `(weight, defect_norm, vinf_per_encounter_kms)` trace, and a derived
  **anchor-retention** metric = the change in best E/M anchor residual from the
  first penalized rung to the final `w=0` rung (did V∞ hold near anchor or snap
  away when the penalty lifted?).
- **Errors.** A diverging rung does not abort: record it in the trace and stop the
  ladder early, returning the best `w`-reached result flagged not-converged
  (mirror `continuation_correct`).

### Unit 3 — Batch wiring + verdict (`scripts/shooter_russell_batch.py`)

Add a `--family-pinned` mode (or a sibling script `shooter_family_pinned_batch.py`
— see Open Questions) that, per descriptor row × best-phase epoch, calls
`family_pinned_shoot` with the row's `{"E": anchor_e, "M": anchor_m}` and records
to a JSONL runlog: final-rung converged?, defect, emerged V∞ vs anchor,
anchor-match, anchor-retention, bend, plus the per-rung trace. Detached,
per-(row,epoch) append+fsync checkpoint, `--resume`, per-evaluation heartbeat —
the same resilience pattern as the STM batch. PROPOSED V0→V1 flag only under the
golden gate above; HELD, no writeback.

## Data flow

ShootingSeed + `{E,M}` anchors → ladder `λv: W → 0` [penalty-augmented STM shoot,
warm-started each rung] → final `λv=0` ShootResult + per-rung trace →
emerged V∞ vs anchor + anchor-retention → held verdict (PROPOSED or stronger
negative).

## Testing (TDD)

- **Unit 1 residual:** with `vinf_anchors=None` or `vinf_weight=0`, residual is
  byte-identical to plain `defect_residual` (array equality). With `w>0`, the
  appended rows equal `sqrt(w)*(|v∞_i|−anchor)` computed independently.
- **Unit 1 Jacobian:** augmented `_stm_jacobian` matches a finite-difference of
  the augmented residual (rel < 5e-3) on the cheap Sun-only fixture from
  `tests/nbody/test_shooter_stm_jacobian.py`, including the penalty rows.
- **Unit 2 driver:** on a constructed seed sampled from a continuous arc whose
  node V∞ is set as the anchor, the `λv→0` solve **retains** that V∞ (anchor
  residual stays small across the ladder); a smoke test that the ladder
  warm-starts and returns a `FamilyPinnedResult`. `@pytest.mark.slow`.
- **Golden honesty:** an explicit test that the recorded V∞ is read from the
  `w=0` rung (not a penalized rung), and that no test EXPECTS a self-computed V∞ —
  the only sourced number is the anchor, which is an INPUT to the penalty, never
  the asserted output.

## Honesty gates

1. Penalty ramps to zero; recorded V∞ is the `w=0` solve's emerged value, never
   imposed; no tolerance loosening; no catalogue writeback (held).
2. Anchor-match (within 0.5 of both anchors at `w=0`, converged, bend-feasible) =
   PROPOSED V0→V1, recorded not applied. Anchor snap-away or `w=0` non-convergence
   = stronger characterized negative.
3. The sourced anchor is an INPUT (basin selector); the test under test is whether
   a ballistic member survives at `w=0`.

## Out of scope / deferred

- The MBH multi-start and the hybrid (penalty-homotopy inner + MBH) escalations —
  documented as the next lever only if the penalty homotopy shows partial
  anchor-retention.
- Ephemeris (e/i) trajectory homotopy — already tried (the conic N-arc lane).
- Out-of-plane parent (#414); freeing node epochs in the shoot.

## Open questions (resolve in the plan)

1. **Penalty weight scale `W`.** The leg defects are km / km·s⁻¹ (O(1e3–1e9 on a
   bad seed); the V∞ penalty is km/s (O(1–30)). `W` must make the penalty bite at
   the top of the ladder without swamping continuity. Plan: pick `W` from a short
   calibration on row `9.353Gg2` (the one row that already drives defect to ~7e3)
   so the penalty residual is ~comparable to the converged continuity residual;
   document the chosen value. This is a tuning knob, not a sourced value.
2. **New script vs `--family-pinned` flag.** Lean to a **sibling script**
   (`shooter_family_pinned_batch.py`) to keep the STM-batch verdict script frozen
   as the published negative; decide in the plan.

## References
- `docs/notes/2026-06-21-shooter-stm-batch-results.md` (the verdict this answers),
  `docs/notes/2026-06-21-narc-continuation-results.md` (conic off-anchor basin).
- `nbody/shooter.py` (`defect_residual`, `_stm_jacobian`, `shoot`),
  `search/shooter_russell_seed.py`, `search/narc_continuation.py`,
  `search/continuation.py` (ladder/warm-start pattern).
- Memory: `project_dsm_closure_modeljump_blocker`,
  `project_s1l1_realeph_closure_blocker`, `feedback_golden_tests_sourced_only`,
  `feedback_orbit_closure_discipline`, `feedback_long_runs_acceptable`,
  `feedback_never_give_up_reproducing_papers`,
  `reference_rebound_variation_custom_force_gotcha`.
