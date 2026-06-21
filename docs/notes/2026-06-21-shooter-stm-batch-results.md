# 2026-06-21 — #388 STM-Jacobian literal-parent shooter batch

**Status: RUNNING (detached, checkpointed). This note is updated as rows land.**
HELD — no catalogue writeback either way.

## What changed vs the FD-era run

The earlier full-shooter-from-Russell-parent attempt
(`2026-06-21-shooter-russell-results.md`) was **compute-infeasible**: the LM
solver's internal finite-difference Jacobian is `(6*n_nodes+1)` full-state
residuals per iteration. For the `E-E-M-M` rows (4 nodes → 24 free components)
that is 25 residual evaluations × 3 leg propagations each ≈ 75 leg-propagations
**per LM iteration**, ~1000 s serially. That forced a `K=1` / `max_nfev=2`
coverage cap that could only sample the seed basin, never run to convergence —
so the lane's verdict rested on standing two-line evidence, not on a converged
literal-parent solve.

This batch uses the **analytic block-bidiagonal STM Jacobian**
(`shooter.shoot(jacobian="stm")`, committed `8cfe894`): the multiple-shooting
residual is full-state continuity, so the Jacobian is exactly
`∂c_i/∂node_i = Φ_i` (the per-leg state-transition matrix), `∂c_i/∂node_{i+1} =
−I`, assembled from **one co-integrated variational propagation per leg** — 3
augmented leg-props per Jacobian instead of ~75 FD leg-props. The variational
particles carry the full dynamics including the rails perturbers (the
gravity-gradient-tensor fix, commit `8196f55`), so Φ is correct on the
flyby-containing legs, not just cruise. Parity vs the FD oracle: rel `5.1e-7` on
the Sun-only fixture.

Because the per-iteration cost collapses, the FD-era coverage cap is **removed**
(per `feedback_long_runs_acceptable` — monitorability is not the measure of the
best path): `K=3` best-phase-error epochs over `LaunchWindow 1..21`, `max_nfev=100`,
run to completion however long it takes. Each `(row, epoch)` record is written to
`data/runs/shooter-stm-batch.jsonl` immediately (append + fsync); `--resume`
skips completed pairs.

## STM-vs-FD speedup (measured)

_Pending — filled from `shoot_wall_sec` in the runlog vs the FD-era ~470 s
single-iteration baseline._

## Per-row results

_Pending — table populated from the runlog: per row (best epoch) converged /
defect_norm / emerged n-body V∞ vs the SOURCED E&M anchors / anchor-match /
bend-feasible / wall._

| id | V | conv | defect | emerged V∞ (km/s) | anchor E/M resid | match | bend | wall |
|----|---|------|--------|-------------------|------------------|-------|------|------|
| _pending_ | | | | | | | | |

## Verdict

_Pending._ Decision rule (held either way):
- A row that **closes + anchor-matches** (esp. `mcconaghy-2006-em-k2` or a V3
  regression row) → PROPOSED V0→V1, recorded, **not applied**.
- If none close to the sourced family → the decisive characterized negative: the
  full multiple-shooting corrector, seeded from the literal Russell parent and
  now run to convergence under the tractable STM Jacobian, does not land the
  sourced high-V∞ family — confirming the obstruction is the seed basin /
  family selection, now with a converged solve behind it rather than a
  compute-capped probe.

## References
- `docs/notes/2026-06-21-shooter-russell-results.md` (the FD-era infeasibility
  finding this run answers).
- `docs/notes/2026-06-21-narc-continuation-results.md` (conic off-anchor basin).
- `nbody/shooter.py::_stm_jacobian`, `nbody/propagator.py::_install_rails_forces`
  (the variational perturber gradient).
- Memory: `feedback_long_runs_acceptable`, `project_dsm_closure_modeljump_blocker`,
  `project_s1l1_realeph_closure_blocker`, `feedback_orbit_closure_discipline`.
