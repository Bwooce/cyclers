# #806 — EGGIE Gate-B Io V∞ miss: underdetermined projection wander, fixed by weight continuation

**Date**: 2026-08-09
**Task**: `#806` (data/OUTSTANDING.md) — root-cause and fix
`tests/search/test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface`'s
`Io_out` V∞ miss: `abs(8.313346200993646 - 8.38) = 0.0667` vs the `< 0.05` gate.
**Verdict**: genuine corrector-robustness defect in
`src/cyclerfinder/search/eggie_ballistic.py::refine` — **fixed in the source**, no test
tolerance changed, no golden value touched. `#584`'s "BLAS/Accelerate-vs-Linux rounding
artifact" diagnosis was directionally right about the *trigger* (backend-dependent solver
path) but wrong to classify it as irreducible noise: the machine-dependence was an
artifact of a fixable two-stage corrector design.

## The defect (measured, not assumed)

`table4_vinf_eggie()` ran a two-stage correction:

1. **Pull stage** — `least_squares` on 8 residuals (5 core ballistic + 3 Table-4-V∞
   pulls weighted `w = 0.5`) in 6 unknowns. Measured landing: `Io_out = 8.377696`
   (0.0023 km/s from the sourced 8.38), core resnorm `1.02e-3` km/s — i.e. still
   ~1e-3 *off* the ballistic manifold, because the weighted compromise splits the
   residual budget between the pulls and the core constraints (core resnorm scales as
   `~4e-3 * w^2`).
2. **Projection stage** — `least_squares` on the 5 core residuals alone: **5 equations,
   6 unknowns — underdetermined**. The on-manifold solution set is a 1-D family, and
   trf, started 1e-3 off-manifold, wandered **along the family tangent** before
   converging: measured `dx` had ToF components of ~200-1300 s and the converged point
   sat at `Io_out = 8.313346` — a **0.064 km/s tangential slide** from the stage-1
   point (37 nfev vs stage 1's 3). Where trf stops along that family is solver-path
   dependent, hence BLAS-backend dependent (Linux OpenBLAS landed inside the 0.05 gate,
   macOS Accelerate 33% outside it — `#584`'s observation, now mechanistically
   explained). With `include_seam=False` (the interior sub-tour) the projection is 4×6,
   a 2-D family — even more room to wander.

So the 0.0667 miss was **not** a model floor: the ballistic manifold demonstrably
contains a point at the Table-4 V∞ levels (the fix below reaches `Io_out = 8.377265`
with core resnorm `1.2e-11`), and the sourced values are only quoted to 2 decimals
(±0.005 quantization) anyway. The old code simply let the final solve drift away from
the point stage 1 had already found.

## The fix: geometric weight continuation

`refine()` now reduces the pull weight geometrically (`w *= 0.1` per stage, from
`target_w` down to `1e-5`) instead of solving once at `target_w` and then projecting.
Every continuation stage stays **overdetermined** (core + pull residuals), so the
iterate never has a free tangent direction to wander along, and the core resnorm
contracts as `~w^2` down to the Lambert/solver floor:

| stage w | core resnorm (km/s) | Io_out (km/s) |
|---------|--------------------:|--------------:|
| 0.5     | 1.02e-3             | 8.377696      |
| 0.05    | 1.21e-5             | 8.377270      |
| 0.005   | 1.21e-7             | 8.377265      |
| 5e-4    | 1.21e-9             | 8.377265      |
| 5e-5    | 1.19e-11 (floor)    | 8.377265      |

The final pure-projection solve (kept, so the `target_w=0` contract "core residuals
driven to zero" is verified on the same code path) now starts already on-manifold and
**does not move at all** (measured `dx = 0` exactly) — the wander mechanism is
eliminated rather than merely reduced, so the result is robust to BLAS-backend rounding
in the projection's iteration path.

## Post-fix results (all 3 constructions, all documented regimes preserved)

- **GATE B (`table4_vinf_eggie`)**: Europa 9.1247 / Ganymede 7.0714 (both, equal to
  4.6e-14) / Io 8.3773 km/s — all within 0.005 of the sourced Table-4 values (gate is
  0.05); core resnorm 1.2e-11; total ΔV 2.1e-9 m/s; altitudes still sub-surface
  (G1 −2627, G2 −2627, Io −1817, E −1559 km) — the #480 binding-constraint finding is
  unchanged, just cleaner.
- **GATE B interior (`interior_table4_eggie`)**: now hits 9.12 / 7.07 / 8.38 *exactly*
  (its seam-free interior manifold is 2-D, so both pulls are absorbable on-manifold);
  interior altitudes 1418 / 2233 / 7177 km in-window, seam defect 0.2508 km/s
  (> the 0.2 gate), Europa flyby sub-surface (−340 km) — regime unchanged.
- **GATE A (`feasible_ballistic_eggie`)**: `target_w=0` path untouched, byte-identical
  behavior.

## Invalidation sweep (`feedback_bugfix_invalidates_past_searches`)

`refine()`'s `target_w>0` path has exactly two callers, both in this module
(`table4_vinf_eggie`, `interior_table4_eggie`), both test-only constructions
documenting `#480`'s verdict — no catalogue row or search negative depended on the
wandered values, and the qualitative `#480` verdict (Table-4 V∞ reachable ballistically
but bend-infeasible in the strict 2D model) is unchanged. The sibling modules built
"same philosophy" were checked for the same defect class and are clean:
`eige_ballistic.refine` is a *single* always-overdetermined solve (soft altitude pulls
kept in the final residual, pinning its 2-DOF family — no separate underdetermined
projection); `ll2011_ballistic.refine_gipeipe` is likewise a single solve. No follow-up
tasks required.

## Verification

- `uv run pytest tests/search/test_eggie_ballistic.py` — 3/3 pass (the previously
  failing Gate-B test now passes with 0.0027 vs 0.05 margin).
- `uv run pytest tests/search -q`, `uv run ruff check .`, `uv run ruff format --check .`,
  `uv run mypy src tests` — clean (see commit).
