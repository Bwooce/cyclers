# M7 implementation plan — Mars-perturbed B-plane horizon-TCM (the #423 build)

**2026-06-23.** Decided after the infra audit (see the agent map in the session) and
the #198-epoch finding (`2026-06-23-appc-s1l1-tcm-epoch-rederivation.md`). The user
chose the **full Mars-perturbed** build (option B), not the Sun-only bounded wrap.

## Design of record: leg-by-leg n-body position-targeting (NOT global B-plane shooting)

The infra audit proposed extending `nbody/shooter.py`'s global multiple-shooting
residual with B-plane constraint rows + Jacobian blocks. We take a **lighter,
equally-faithful** formulation that reuses more and touches `shooter.py` not at all:

**Why the Sun-only proxy diverged (the #169 Mars-perturbed test):** it propagated each
leg with the *sourced* App-C v∞ direction through *continuous Mars gravity* and snapped
to the planet — under perturbation the sourced direction no longer reaches the next
encounter, so the chain blew up post-flyby. The cure is to **target the next encounter**.

**M7 measure (operational TCM).** Walk the cycler node-to-node as ONE continuous
trajectory, Mars-perturbed, no re-anchoring of velocity:

1. At node *i* the spacecraft is at the real planet position `r_i` (patched-conic flyby
   at the body; the ballistic miss is recorded as evidence, must be ≪ SOI), arriving
   with `vinf_in`.
2. **Target leg *i*→*i*+1**: solve for the departure velocity `v_dep` such that an
   n-body (Mars-perturbed) propagation from `(r_i, v_dep)` over the leg ToF arrives at
   the real planet position `r_{i+1}`. Newton on the 3×3 `∂r_f/∂v_dep` = `stm[0:3,3:6]`
   from `propagator.propagate(..., with_stm=True)`. This is an **n-body Lambert** and is
   the "B-plane targeting" expressed as position targeting (aiming the flyby = choosing
   the post-flyby v∞ = hitting the next planet).
3. `vinf_out_required = v_dep − v_planet(t_i)`. A real ballistic flyby rotates `vinf_in`
   for FREE up to `max_bend` at constant magnitude; the **maintenance ΔV** is the part
   it cannot supply:
   `dv = core.flyby.flyby_dv(vinf_in, vinf_out_required, mu_body, rp_min_body)`
   (= |Δ|v∞|| magnitude term + excess-bend term; 0.0 if ballistic-feasible). Reused
   verbatim — same decomposition as the Sun-only `continuous_chain`.
4. `horizon_tcm_mps = Σ dv over all flybys over N cycles`; `per_cycle_tcm_mps` likewise.

This reuses: `propagator` (Mars perturbers + co-integrated STM, gravity-gradient tensor
already applied to variational particles — REBOUND-variation gotcha handled), the Jones
bend kernel (`core.flyby` / `nbody.bplane`), `flyby_dv`, and `real_closure`'s
`construct_real_ephemeris_cycler` (the seed). It does **not** touch `shooter.py`.

## Tasks (TDD, incremental commits)

**Task 1 — `nbody/maintenance_shoot.py::target_leg`.** Fixed-time n-body position-
targeting Newton: `(r0,t0,t1,r_target,v_guess,perturbers,ephem,cache) -> (v_dep, miss_km,
converged)`. Uses `stm[0:3,3:6]`. Honest failure: return `converged=False` if Newton
stalls / propagation diverges (off-family legs stay unmeasured, not forced).
*Golden:* in Sun-only mode (`perturbers=()`) it must reproduce the two-body
`core.lambert.lambert` solution to tolerance (independent cross-check). + a Mars-perturbed
smoke (converges, miss ≪ SOI).

**Task 2 — `continuous_maintenance_chain`.** Walk the real-eph `Cycler` node-to-node,
Mars-perturbed, calling `target_leg` per leg, decomposing residual via `flyby_dv`,
accumulating. Returns nodes + `horizon_tcm_mps` + `per_cycle_tcm_mps` + a `diverged`
flag. *Golden:* Sun-only mode on S1L1 reproduces the corrected Sun-only proxy
(~40.2 m/s, since targeting the real next planet in two-body ≈ the sourced direction).
Mars-perturbed mode gives the real (higher) measurement — recorded, not pre-judged.

**Task 3 — wire into `real_closure.verify_real_closure`.** Opt-in (`compute_tcm=False`
default; the N-cycle Mars-perturbed shoot is minutes/row). When on, populate
`horizon_tcm_mps` / `per_cycle_tcm_mps` (replace the hardcoded 0.0 at ~644/788) and run
`verify.dv_band_acceptance.v3_class_split_verdict` (#424) on the result — the #175 manual
convention retires into code. *Test:* the verdict falls out programmatically for S1L1
(ballistic) and an Aldrin powered row.

**Task 4 — golden-pin + document.** Pin the Mars-perturbed M7 numbers for S1L1 + the two
App-C rows as computed regression anchors (labelled, not sourced goldens); doc note with
the measured per-cycle TCM and the convergence/divergence map. Full
`tests/data tests/search tests/verify` + nbody goldens green.

**Task 5 (Phase 2/3) — scale.** Banded + V1+ rows; then catalogue-wide
detached/checkpointed batch (days-scale, per the long-runs rule). Each row acquires a
`computed-v3` band or stays null (infeasible-to-shoot). This is the validation-ceiling
lift; expect a residue of unshootable (off-anchor / high-V∞) rows.

## Acceptance (M7 done)
- `real_closure.horizon_tcm_mps` is a real Mars-perturbed measurement for ≥ the Phase-1
  rows; Sun-only-mode goldens reproduce; V3 promotion uses `v3_class_split_verdict`
  programmatically. Computed bands carry `dv_band_source="computed-v3"`; sourced bands
  untouched; mismatches surfaced. No catalogue claim altered without review.

## References
- `2026-06-23-m7-scoping-plan.md`; `2026-06-23-appc-s1l1-tcm-epoch-rederivation.md`;
  `nbody/{propagator,bplane,flyby_gradients,shooter}.py`; `core/flyby.py::flyby_dv`;
  `verify/real_closure.py`; `verify/dv_band_acceptance.py` (#424).
