# 2026-06-21 — #388 Seed the full shooter from the Russell parent: built; shoot operationally infeasible

**Verdict: the seeding bridge is built and works; the full n-body multiple-shooting
solve is computationally infeasible to run to a verdict over these cyclers. The
underlying question ("does the literal Russell parent break the seeding wall?") rests
on the standing two-line evidence, which already establishes the
seeding/family-selection wall. HELD, no writeback.**

## What was built (committed)

- `src/cyclerfinder/search/shooter_russell_seed.py::russell_shooting_seed`
  (commit `5d62f79`): bridges the constructed, golden-validated Russell ψ
  generic-return parent → the shooter's seed format, via
  `correct._vinf_nodes` (per-leg v∞ **vectors**) → `nbody.shooter.seed_from_conic`
  → `ShootingSeed`. Adapter test passes; the node solve is non-degenerate for the
  catalogue rows. This is the in-basin seed the #135 verdict said the shooter
  needed but never had (it was only ever fed the searched #133 near-miss survey).
- `scripts/shooter_russell_batch.py` + a `@pytest.mark.slow` shoot smoke
  (commit `9ccdf09`): correct and CI-safe; pass `bodies=`, reduced `max_nfev`,
  `n_jobs=16`.

## Why no per-row verdict — operational infeasibility (measured)

`nbody.shooter.shoot` runs restricted-n-body propagation between encounter nodes
with a **finite-difference Jacobian** (≈ 1 + 6·n_nodes residual evals per LM
iteration; each residual eval propagates every multi-rev leg over the full
multi-synodic cycler in real n-body). For these `E-E-M-M` two-synodic cyclers
(~4-yr arcs, multi-rev legs):

- A single `shoot()` at default `max_nfev=200` ran **>400 s with no return** (killed).
- The full batch (K=3 epochs × 4 rows) **exceeded a 50-min timeout writing nothing**.
- A bounded probe of just two rows at **`max_nfev=2`** (two residual evals each)
  did **not complete in 20 min** — so even a single LM iteration's n-body
  propagations exceed a practical budget.

`shoot`'s `max_wall_sec` guard does not bound this effectively (it is checked
between solver iterations; one iteration already overruns). So a clean
shooter-from-Russell-parent DE440 verdict is **not obtainable in feasible time**
without a major performance effort.

## The underlying question is already answered by standing evidence

The point of this lane was to test whether the *constructed* Russell parent (vs a
searched near-miss) breaks the seeding/family-selection wall. Two independent lines
already establish that wall, and neither needs the (infeasible) full shoot:

1. **Conic continuation (this session, `2026-06-21-narc-continuation-results.md`):**
   the N-arc homotopy reaches DE440 *and* Russell's own mean-element model and lands
   the **low-energy off-anchor basin** (residuals 1.4–2.5 km/s; emerged V∞ ≈ half the
   sourced anchor) — the wrong family member — from the Russell parent.
2. **#135 like-for-like (`nbody/shooter.py` docstring):** the full shooter, seeded
   from the near-miss survey, "closes geometrically but lands OFF-ANCHOR; verdict =
   **seeding/basin, not solver deficiency**." S1L1
   ([[project_s1l1_realeph_closure_blocker]]) confirms family-selection persistence.

Both say: from a low-energy-relaxable seed, the corrector finds a *different*
(low-energy) member of the family, not Russell's published high-V∞ cycler. The
constructed parent does not cheaply change that, and the full n-body confirmation is
compute-bound.

## Honest status & next lever

`mcconaghy-2006-em-k2` stays V0; no catalogue writeback (`data/catalogue.yaml` /
`validate.py` untouched). This closes the #388 real-eph thread at: **the obstruction
is the seed basin / family selection, and the full-fidelity confirmation is
gated on shooter performance, not on more modelling.**

The single concrete next lever, if pursued: **replace the shooter's finite-difference
Jacobian with the analytic state-transition-matrix Jacobian** (Pellegrini-Russell
2016, already digested — `docs/notes/2026-06-17-digest-pellegrini-russell-2016.md`)
— the FD Jacobian is the dominant cost; an analytic/variational STM would cut the
per-iteration propagation count by ~(1+6·n_nodes)×, plausibly making the shoot
tractable. Plus `n_jobs` FD-column parallelism and vectorized propagation. That is a
shooter-performance project (its own spec), not a closure-modelling change.

## References
- `2026-06-21-narc-continuation-results.md` (conic off-anchor + model-isn't-the-cause).
- `nbody/shooter.py` (#135 verdict), `docs/notes/2026-06-06-russell12-likeforlike.md`.
- `docs/notes/2026-06-17-digest-pellegrini-russell-2016.md` (analytic STM Jacobian).
- Memory: `project_dsm_closure_modeljump_blocker`, `project_s1l1_realeph_closure_blocker`,
  `feedback_never_give_up_reproducing_papers`, `feedback_orbit_closure_discipline`.
