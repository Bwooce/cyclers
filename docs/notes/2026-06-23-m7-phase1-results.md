# M7 Phase 1 results — per-row real-ephemeris maintenance-ΔV (horizon TCM) is built

**2026-06-23 (#423).** Phase 1 of the M7 build (plan: `2026-06-23-m7-implementation-plan.md`)
is complete: `real_closure.RealClosureResult.horizon_tcm_mps` (hardcoded `0.0` since M6b)
is now a real measurement under the opt-in `compute_tcm` flag, golden-validated on the
real S1L1 cycler.

## What was built

- **`nbody/maintenance_shoot.target_leg`** — n-body fixed-time position-targeting Newton
  on the propagator's co-integrated STM (`d r_f/d v_0 = STM[0:3,3:6]`). Solves for the
  departure velocity whose (perturbed) propagation hits the next encounter planet — the
  n-body analogue of a fixed-time Lambert / the flyby B-plane targeting as position
  targeting. Golden: Sun-only ≡ two-body `lambert()` to mm/s; Mars-perturbed sub-km.
- **`nbody/maintenance_shoot.continuous_maintenance_chain`** — walks a cycler's encounter
  sequence as ONE continuous trajectory, targets each leg, and bills each interior flyby's
  maintenance ΔV via `core.flyby.flyby_dv_for` (the part a free ballistic flyby cannot
  bridge). Honest divergence → `inf` (row stays V0, never forced).
- **`real_closure.verify_real_closure(compute_tcm=True)`** + `_compute_horizon_tcm` — builds
  the N-cycle node sequence (incl. the inter-cycle home flyby), seeds per-leg from the
  constructed cycler's rev-correct `leg.v_depart`, runs the chain, populates
  `horizon_tcm_mps`/`per_cycle_tcm_mps`, which feed `v3_class_split_verdict` (#424).

## Three design decisions (with evidence)

1. **Leg-by-leg position-targeting, not global B-plane multiple-shooting.** Lighter, reuses
   the propagator STM + `flyby_dv` + `real_closure`'s Cycler builder; does not touch
   `shooter.py`. The infra audit confirmed soundness.
2. **Per-leg perturbers = system bodies − the leg's endpoint flyby bodies.** Naive continuous
   integration of an endpoint body from its own node centre diverges (the spacecraft starts
   in the softened core) = the patched-conic handoff artifact, not fuel. Endpoint exclusion
   removes it while keeping every genuine non-endpoint third-body perturbation — so M7 works
   for moon tours / Venus-flyby chains, not just Earth-Mars (where it reduces to Sun-cruise).
   Validated: naive diverges; default converges; a non-endpoint Jupiter perturber shifts the
   TCM (real signal).
3. **Per-leg seeds from the rev-correct sourced/constructed departure velocity.** Decisive:
   seeding a multi-rev resonant cycler from a single-rev Lambert lands the Newton on a
   high-energy single-rev transfer — S1L1 read **530,000 m/s** and the powered rows diverged.
   With sourced-v∞ seeds S1L1 reads **40.1 m/s**, reproducing the corrected #169 proxy.

## Golden / validation

- S1L1 full 22-node App-C sequence, sourced-seeded, 7-cycle: **~40 m/s**, all 21 legs
  converge — matches the corrected Sun-only proxy (~40.2 m/s, post-#198) to abs 5 m/s.
  M7 position-targeting and the proxy's sourced-direction accounting measure the same
  ballistic physics on a real resonant cycler. Under the 120 m/s V3 budget → S1L1
  essentially_ballistic (unchanged).
- Wiring test: `_compute_horizon_tcm` on a real-eph E-M-E cycler returns a finite 2-cycle
  TCM (incl. the home flyby) and drives `v3_class_split_verdict` programmatically — the
  manual #175 convention is retired into code.

## Scope / limits

- **Heliocentric lane only** (μ_Sun, Sun-central propagator) — exactly the `real_closure`
  cycler lane (the ~215 V0 Earth-Mars SnLm census + heliocentric multi-planet cyclers).
  Planet-central moon-tour M7 (Saturn/Uranus central body) is a separate generalisation;
  those rows validate via their own V3/V4 lanes (`v2_moontour`/`v3_3d`/`v4_uranus`), not
  touched here.
- **Phases 2-3 remain** (not started): run M7 on the banded + V1+ rows, then a catalogue-wide
  detached/checkpointed batch (days-scale) = the V0→V3 mass-promotion / validation-ceiling
  lift. Expect a residue of unshootable (off-anchor / high-V∞) rows that stay V0 honestly.

## Gate status

Per `project_dvband_validation_coupling_gate`: the gate to resume discovery+validation
campaigns was "land #424 + ≥ M7 Phase 1". Both are now done — the gate is **cleared**.
Phases 2-3 are the ceiling-raise, not the gate.
