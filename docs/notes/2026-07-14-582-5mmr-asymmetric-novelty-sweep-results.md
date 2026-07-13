# #582: 5-MMR asymmetric novelty sweep results (empty, method-conditional)

**Date:** 2026-07-14
**Scope:** the full novelty sweep #582's build/positive-control dispatch deliberately
deferred to the coordinator. Positive control (MMR 3:2, hardest tabulated case)
passed with 0.6-4.1% error against 3-5% tolerances (see
`docs/notes/2026-07-13-582-asymmetric-3d-positive-control.md`); this note covers
the actual novelty sweep run afterward.

## What was run

All 5 tabulated interior MMRs (`3:2, 5:2, 3:1, 4:1, 5:1`), paper-scale GA budget
(population=200, generations=400, 8 workers), via
`scripts/run_582_asymmetric_3d_niching_search.py --mode ga`. Each MMR took
9-23 minutes wall time on this machine (longer than #583's Sun-Earth sweep --
this fitness function's period T~8-12 TU per orbit vs. #583's 1-year window).

Each finished population was then analyzed via `--mode analyze`: cluster into
distinct high-fitness basins (25 representatives per MMR, fitness floor 0.9,
bounds-normalized distance threshold 0.1), then run every representative
through the mandatory pipeline (refine via `correct_general_periodic_3d` ->
`classify_symmetry` -> `build_candidate_signature` -> `literature_anchors_engaged`).

## Result: 0/104 asymmetric

| MMR | Clusters | Converged | Symmetric | Asymmetric |
|---|---|---|---|---|
| 3:2 | 25 | 23 | 23 | 0 |
| 5:2 | 25 | 23 | 23 | 0 |
| 3:1 | 25 | 22 | 22 | 0 |
| 4:1 | 25 | 17 | 17 | 0 |
| 5:1 | 24 | 19 | 19 | 0 |
| **Total** | **124** | **104** | **104** | **0** |

Every converged cluster, across all 5 MMRs and 124 distinct cluster
representatives, classified as the already-known SYMMETRIC #440 family. Not
one asymmetric member was found anywhere in the sweep.

## Interpretation: very likely a search-box artifact, not a physical negative

`isolated_3d_asymmetric_fitness.py::mmr_bounds()`'s own docstring states the
box is deliberately narrow and centered on the known symmetric seed:
`x0_frac=0.15`, `ydot0_frac=0.35`, but critically `z0_abs=xdot0_abs=zdot0_abs
=0.05` -- tight ABSOLUTE bounds on exactly the state components whose
departure from zero is what breaks mirror symmetry in the first place. This
was a deliberate, documented engineering choice to keep the GA inside one
resonance's basin rather than drifting to a neighboring MMR or the exterior
1:2 family-selection trap #440 already documents.

Given that, a 0/104 asymmetric result is exactly what a narrow-box search
centered on a strongly-attracting symmetric basin would produce whether or
not a genuine asymmetric member exists nearby -- the box may simply never
give the GA room to reach it. This does NOT mean "no asymmetric members
exist at these resonances"; it means "none were found in a box this narrow,
with this budget." Per [[project_negative_results_registry]], empty is
always conditional on the method -- stamped as such in
`data/empty_regions.jsonl` (`region_id
er3bp-isolated-3d-asymmetric-mu0.001-5mmr-582-2026-07-14`).

## Natural next test (NOT run -- a scope decision, not an automatic follow-up)

Re-run with deliberately widened `z0_abs`/`xdot0_abs`/`zdot0_abs` (the
symmetry-breaking degrees of freedom), keeping `x0_frac`/`ydot0_frac`/`t_frac`
as-is to anchor the resonance and reduce MMR-drift risk. This is the fairer
test of the asymmetric-novelty hypothesis #582 was built to probe. Whether to
run it is a coordinator/user scope decision, not automatically warranted by
this result alone -- widening those bounds also widens the risk of the GA
drifting into an unrelated basin (a different MMR, the exterior trap, or a
non-family degenerate solution), so the bounds choice itself would need the
same care `mmr_bounds()`'s original calibration got.

## Reproduction

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode ga --workers 8
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode analyze
# (repeat --mmr for 5:2, 3:1, 4:1, 5:1)
```

Raw data: `data/found/582_niching_ga/ga_{p}_{q}_final.npz` (population) and
`{p}_{q}_analysis_summary.json` (per-cluster pipeline results) for each MMR.
