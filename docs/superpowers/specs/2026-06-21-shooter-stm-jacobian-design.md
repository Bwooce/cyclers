# Analytic STM Jacobian for the n-body shooter — Design Spec

**Date:** 2026-06-21
**Status:** Design approved (brainstorming), pending implementation plan.
**Tracks:** #388 / #135 / #347 (shared STM consumer). **Unblocks:** the shooter-from-
Russell-parent verdict (`2026-06-21-shooter-russell-seed-design.md`).

## Goal

Replace the n-body shooter's finite-difference Jacobian (`1+6n` full re-propagations
per Newton step) with a **state-transition-matrix (STM) Jacobian** assembled from
per-leg variational propagations (~`n-1` augmented propagations per step → ~`(1+6n)×`
fewer propagations). This makes `nbody.shooter.shoot` tractable on the multi-year,
multi-rev SnLm cyclers, so the decisive test — the full multiple-shooting corrector
seeded from the literal constructed Russell parent — can be run to a verdict. Held;
no catalogue writeback.

## Why this is the right lever (and the payoff)

The FD shooter is compute-infeasible on these cyclers (single shoot >400s no-return;
`max_nfev=2` >20 min). The dominant cost is the FD Jacobian's repeated full n-body
propagations. The STM gives the Jacobian from one augmented propagation per leg. The
full shooter (full 6-state-per-node DOF) seeded from the literal Russell parent is the
one configuration never tested and the one with the DOF to potentially hold the
sourced high-V∞ family where the conic lane (magnitude-only DOF) collapsed to the
low-energy basin. Bonus: an n-body STM is also the missing piece #347 flagged for
monodromy/eigenvalue work — a second consumer.

## Approach (settled in brainstorming)

REBOUND **co-integrated variational particles** for the STM: state and STM share one
adaptive IAS15 step path in the same simulation, so this is NOT the mixed-propagation
variable-step pitfall Pellegrini-Russell 2016 warns about — it is effectively the
paper's fixed-path case AND gives the one-propagation speedup. The FD path stays the
default and the parity oracle; a **complex-step (CSD) single-column check** is the
high-accuracy fidelity gate on the STM. (Pure CSD for every Jacobian was rejected: same
propagation count as FD → no speedup. Pure fixed-path-replay was rejected: REBOUND
co-integration already achieves no-mixing natively.)

## Reuse (no rebuild)

- `nbody/propagator.py::RestrictedNBody.propagate` — extend with `with_stm=True` via
  `rebound.Simulation.add_variation()` (REBOUND 5.0.0 native variational particles,
  confirmed available). Mirror the existing variational-STM pattern in
  `core/cr3bp.py` (`propagate(with_stm=True)`, the analytic-dynamics precedent).
- `nbody/shooter.py` — `defect_residual`, `_x_to_states`/`_states_to_x`,
  `_fd_jacobian` (kept as oracle), the `shoot` LM loop. Add the STM Jacobian as an
  opt-in alongside FD.
- `search/shooter_russell_seed.py` (the seed bridge) + `narc_continuation` /
  `cycler_assembly` (the parent) — unchanged.

## Architecture — phased (each phase independently testable)

### Phase 1 — STM propagation (`nbody/propagator.py`)
`RestrictedNBody.propagate(state0, t0, tof, bodies, *, with_stm=False)`:
- `with_stm=False` → current behaviour, byte-identical.
- `with_stm=True` → add 6 first-order variational particles (one per state component)
  to the sim, co-integrate with the state over the leg, return `(final_state,
  Phi)` where `Phi` is the 6×6 leg STM `∂final/∂state0` read from the variational
  particles' final states.
**Golden/parity:** `Phi` matches a finite-difference STM of `propagate` to a tolerance
(e.g. rel 1e-4) on a known two-body leg AND a flyby-containing leg; and a **complex-step
STM** (one column via `iℎ` perturbation propagated) matches `Phi` to ~1e-8 (the
paper's accuracy gate). Two-body STM also checked against the analytic Keplerian STM
(closed-form) as an independent oracle.

### Phase 2 — STM-assembled multiple-shooting Jacobian (`nbody/shooter.py`)
Add `jacobian: Literal["fd","stm"] = "fd"` to `shoot`. For `"stm"`: the defect
residual is full-state continuity `c_i = propagate(state_i, leg_i) − state_{i+1}`, so
the Jacobian is **block-bidiagonal** — `∂c_i/∂state_i = Phi_i` (the leg STM),
`∂c_i/∂state_{i+1} = −I` — assembled from the per-leg STMs (one variational propagation
per leg). Boundary/periodicity rows handled to match the FD residual's layout exactly.
Pass the assembled Jacobian to `least_squares(jac=…)`.
**Parity:** the STM Jacobian equals `_fd_jacobian` (the oracle) to a tolerance on a
small seed; the LM solve with `jacobian="stm"` reaches the same fixed point as `"fd"`
on a cheap case (e.g. a near-circular short leg), in far fewer propagations (assert the
propagation-count drop).

### Phase 3 — literal-parent shooter batch (detached, checkpointed, days-tolerant)
Re-run `scripts/shooter_russell_batch.py` with `jacobian="stm"` over the descriptor
rows × candidate epochs. Run rules (per `feedback_long_runs_acceptable`):
**incremental per-(row,epoch) JSONL checkpoint write** (a kill/restart loses nothing),
launched detached (`setsid`/`nohup`) so it survives agent reaping; progress reported
from the checkpoint file, NOT by blocking. No artificial epoch/iteration cap for
monitorability — run the methodologically-correct K (best phase-error epochs over
LaunchWindow 1..21) to completion even if it takes days. Record per row: converged /
defect / emerged n-body V∞ vs sourced anchor (0.5 km/s) / bend-feasible. Verdict +
held report; no writeback.

## Data flow

leg (state0, tof) → `propagate(with_stm=True)` (REBOUND variational) → (state, Phi) →
[Phase 2] block-bidiagonal Jacobian from per-leg Phi → `least_squares(jac=…)` LM solve
→ `ShootResult` → [Phase 3] literal-parent batch → emerged V∞ vs sourced anchor → held.

## Error handling

- `with_stm=True` on a degenerate/close-encounter leg where variational integration
  blows up → raise a typed error; the shooter falls back to FD for that Jacobian
  (logged), so a single bad leg doesn't abort the solve.
- Parity-gate failure (STM vs FD/CSD beyond tolerance) is a HARD test failure — the STM
  is wrong; do not proceed to Phase 3 until it passes.
- Phase 3: a non-converging (row,epoch) is checkpointed as a recorded negative, not a
  crash; the batch continues.

## Testing (TDD)

- Phase 1: variational STM vs FD-STM (rel 1e-4) and vs CSD-STM (~1e-8) on a two-body
  leg + a flyby leg; two-body vs analytic Keplerian STM. `with_stm=False` byte-identical.
- Phase 2: STM Jacobian vs `_fd_jacobian` parity; `shoot(jacobian="stm")` reaches the
  same fixed point as `"fd"` on a cheap case with a measured propagation-count drop.
- Phase 3: empirical (the batch). Smoke: one `jacobian="stm"` shoot returns a
  `ShootResult` far faster than the FD smoke (gate the speedup).
- Golden honesty: emerged V∞ compared to the SOURCED Russell anchor; never imposed.

## Honesty gates

1. The STM is the *means*; correctness is enforced by the FD oracle + CSD parity gate
   (the STM must agree before it's trusted).
2. Sourced V∞ anchor emerges from the converged shoot; never imposed; no tolerance
   loosening; no catalogue writeback (held).
3. A non-closure of the literal parent under the now-tractable full shooter is the
   decisive characterized negative; a closure + anchor-match is a PROPOSED V0→V1 (held).

## Out of scope / deferred

- 2nd-order STM (Φ²) — Phase 1 builds Φ¹ only (the Jacobian needs first order); Φ² is a
  later add if a #347 consumer needs it.
- Replacing FD everywhere / other shooter callers — FD stays default; STM is opt-in.
- Out-of-plane parent (#414); catalogue writeback.

## References

- `docs/notes/2026-06-17-digest-pellegrini-russell-2016.md` (variational vs CSD vs FD;
  variable-step pitfall; co-integration avoids mixing).
- `nbody/propagator.py` (RestrictedNBody/REBOUND), `nbody/shooter.py` (`defect_residual`,
  `_fd_jacobian`, `shoot`), `core/cr3bp.py` (`propagate(with_stm=True)` precedent),
  `search/shooter_russell_seed.py` (the seed bridge).
- `docs/notes/2026-06-21-shooter-russell-results.md` (why this lever),
  `docs/notes/2026-06-21-narc-continuation-results.md`.
- Memory: `feedback_long_runs_acceptable`, `project_dsm_closure_modeljump_blocker`,
  `project_s1l1_realeph_closure_blocker`, `feedback_golden_tests_sourced_only`,
  `feedback_never_give_up_reproducing_papers`.
