# #423 / M7 scoping plan — per-row real-ephemeris maintenance-ΔV (horizon TCM)

**2026-06-23.** Scoping (not implementation) per the "do #424 then scope #423" decision.
M7 is the missing measurement that lets the bands→validation coupling bite at scale:
it populates `real_closure.RealClosureResult.horizon_tcm_mps` / `per_cycle_tcm_mps`
(today hardcoded `0.0`) with the real continuous maintenance ΔV over N cycles, which
then feeds the already-built consumers:
- `dv_band_acceptance.classify_dv_band` / `assign_dv_band_from_measurement` (#422) →
  acquire a computed `dv_band` (provenance `computed-v3`) for the ~215 null rows, or
  flag a mismatch vs a sourced band.
- `dv_band_acceptance.v3_class_split_verdict` (#424) → the programmatic V3
  ballistic/powered pass/fail (ballistic ≤120 m/s; powered ≤ own budget ×1.10).

## What already exists (build ON these, do not reinvent)
- **The #169 method** — "continuous-from-one-seed horizon TCM": already computes the
  per-row real-eph TCM for the handful of V3 rows (S1L1 = 62 m/s/7cyc; App-C #188 =
  163.6; #192 = 2040). Gated by `tests/nbody/test_s1l1_corrected_nbody.py` and
  `tests/nbody/test_appc_batch_nbody.py`. **M7 = generalise this from hand-run rows
  to a callable per-row stage** wired into `real_closure`.
- **The analytic-STM shooter (#388)** — block-bidiagonal multiple-shooting Jacobian
  (~8.6× faster/iter than FD). The #388 saga proved the FD-Jacobian full n-body
  shoot is compute-INFEASIBLE on multi-yr multi-rev cyclers (single shoot >400s); the
  analytic STM is what makes per-row TCM tractable. M7 MUST use it, not FD.
- **`real_closure.py`** (M6b) — already does the real-eph phase-matched propagation +
  geometric-drift closure; M7 adds the TCM-budget layer on top of its propagation.
- **`writeback.apply_v3_to_entry`** — already wires `horizon_tcm_dv_mps` into the V3
  block; it just receives 0.0 today.

## The hard part (the validation ceiling, and the #388 wall)
1. **B-plane flyby targeting.** validate.py records: "the higher-fidelity
   Mars-perturbed continuous run needs B-plane flyby targeting (naive
   patched-conic→continuous handoff diverges, an artifact not a fuel cost)." So M7's
   per-cycle TCM must target the B-plane at each flyby, NOT bill the handoff
   divergence as fuel. This is the core correctness risk — get it wrong and every
   row looks powered.
2. **Compute.** Even with the analytic STM, a multi-yr × 7-cycle shoot is seconds-to-
   minutes per row. Full-catalogue (318 rows) is a days-scale batch → must be
   detached + checkpointed (incremental runlog), per the long-runs-acceptable rule.
3. **Feasibility floor.** Some cyclers (the #388 off-anchor / high-V∞ families) may
   not shoot to convergence at all — those stay V0/null (honest), not forced.

## Proposed phasing
- **Phase 1 (golden reproduction, no scale).** Wrap #169 as a callable
  `compute_horizon_tcm(cycler, n_cycles, ephem) -> (horizon_mps, per_cycle)` and
  prove it reproduces the 3 recorded V3 numbers (S1L1 62, #188 163.6, #192 2040) as
  goldens. Wire its output into `real_closure` (replace the 0.0). Run #424's
  `v3_class_split_verdict` on the result and confirm the recorded V3 decisions fall
  out programmatically (not from hand-written evidence). **This alone closes the
  "V3 is manual" gap for the existing V3 rows.**
- **Phase 2 (banded + V2 rows).** Run M7 on the 11 banded rows + the 31 V1+ rows;
  use #422 to (a) assign computed bands to the unbanded-but-reproduced, (b) flag any
  sourced-vs-measured band mismatch for review. Reconcile.
- **Phase 3 (catalogue-wide, compute-gated).** Detached/checkpointed batch over the
  remaining ~215 null rows. Each acquires a `computed-v3` band or stays null
  (infeasible-to-shoot). This is the V0→V3 mass-promotion effort = the validation
  ceiling; expect it to take days and to leave a residue of unshootable rows.

## Acceptance / definition of done (M7)
- `real_closure.horizon_tcm_mps` is a real measurement for at least the Phase-1 rows;
  the 3 #169 goldens reproduce to tolerance.
- V3 promotion uses `v3_class_split_verdict` programmatically (the #175 manual
  convention retires into code).
- New computed bands carry `dv_band_source="computed-v3"`; sourced bands untouched;
  mismatches surfaced, never silently overwritten.
- Full `tests/data tests/search tests/verify` + the `tests/nbody` goldens green.

## Dependencies / sequencing
- #424 (done) — the verdict consumer.
- #422 (done) — the classifier/assignment consumer.
- This plan = #423. Phase 1 is the high-value, bounded slice; Phases 2-3 are the
  compute-bound ceiling. Per `project_dvband_validation_coupling_gate`, completing
  at least Phase 1 (programmatic V3) is the gate before a new discovery+validation
  campaign.

## References
- docs/notes/2026-06-22-dv-band-definitions.md; docs/notes/2026-06-08-v3-powered-classsplit.md;
  docs/notes/2026-06-08-appc-v3-batch-results.md; the #388 family-pinned results;
  src/cyclerfinder/verify/{real_closure,dv_band_acceptance}.py; project_validation_ceiling.
