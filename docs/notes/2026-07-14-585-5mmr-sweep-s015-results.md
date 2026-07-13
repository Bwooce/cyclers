# #585: 5-MMR asymmetric novelty re-sweep at s=0.15 (empty, stronger evidence)

**Date:** 2026-07-14
**Scope:** the full novelty re-sweep #585's build deliberately deferred to the
coordinator. Build + per-rung positive control passed at s=0.15 (see
`docs/notes/2026-07-14-585-resonance-scaled-symmetry-breaking-bounds.md`);
this note covers the actual sweep run afterward.

## What was run

All 5 tabulated interior MMRs, paper-scale GA budget (population=200,
generations=400, 8 workers), at the Fable-reviewed, positive-control-confirmed
`symmetry_breaking_s=0.15` (a genuine 3.4-10.3x widening over #582's original
flat 0.05 box in exactly the symmetry-breaking directions, per-MMR).

## Result: 0/78 asymmetric (still)

| MMR | Clusters | Converged | Symmetric | Asymmetric | Drifted-to-neighbor |
|---|---|---|---|---|---|
| 3:2 | 25 | 2 | 2 | 0 | 0 |
| 5:2 | 25 | 22 | 22 | 0 | 0 |
| 3:1 | 25 | 20 | 20 | 0 | 0 |
| 4:1 | 25 | 20 | 20 | 0 | 0 |
| 5:1 | 25 | 14 | 14 | 0 | 0 |
| **Total** | **125** | **78** | **78** | **0** | **0** |

Every converged cluster is still the known symmetric family. Critically, the
new drift-detection check (#585's own addition) confirms **zero** converged
members drifted to a neighboring MMR — the widened box did not simply leak
into unrelated resonance territory, so the null result cannot be explained
away as "the widening was accidentally too tight/wrong-direction to matter."

Note MMR 3:2 (the hardest, most eccentric tabulated case) converged far
fewer clusters (2/25) than the others — consistent with #585's own
positive-control finding that 3:2 sits closest to genuine basin competition
at wider bounds; still, of the 2 that did converge, both are symmetric.

## Interpretation

This is a materially stronger negative than #582's original stamp
(`er3bp-isolated-3d-asymmetric-mu0.001-5mmr-582-2026-07-14`): the earlier
result was confounded by a demonstrably anisotropic, seed-centered box; this
one uses a resonance-scaled, positive-control-validated widening specifically
targeting the symmetry-breaking degrees of freedom, and the drift-detection
check rules out silent over-widening as a confound.

Two live possibilities remain (neither resolved here, per
[[project_negative_results_registry]] -- empty is always conditional on the
method):

1. **Asymmetric periodic members genuinely do not exist near these interior
   MMRs at mu=0.001** in this idealized CR3BP model -- plausible per the
   planar-CR3BP literature (asymmetric families are typically found at
   exterior resonances and higher eccentricity, not interior MMRs at this
   mass ratio).
2. **A wider box might still reach one**, but s=0.30 already failed its own
   positive control (genuine basin competition, not "box too wide" --
   see #585's build note) and is not certified. Closing this gap further
   would need either a value between 0.15 and 0.30, or a different
   fitness-landscape treatment (e.g. a larger population to hold the
   competing basins simultaneously, since deterministic crowding is
   specifically designed for that -- see #585's own note on this).

Stamped as a separate, additional method-conditional empty region in
`data/empty_regions.jsonl` (`region_id
er3bp-isolated-3d-asymmetric-mu0.001-5mmr-585-s0p15-2026-07-14`) --
does NOT overwrite #582's original stamp, which remains valid for its own
(narrower, flat-box) bounds.

## Reproduction

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode ga --symmetry-breaking-s 0.15 --workers 8
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode analyze --symmetry-breaking-s 0.15
# (repeat --mmr for 5:2, 3:1, 4:1, 5:1)
```

Raw data: `data/found/582_niching_ga/ga_{p}_{q}_s0p15_final.npz` (population)
and `{p}_{q}_s0p15_analysis_summary.json` (per-cluster pipeline results).
