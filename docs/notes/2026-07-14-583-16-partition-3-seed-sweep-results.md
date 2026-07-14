# #583: full 16-partition x 3-seed widened-domain sweep results

**Date:** 2026-07-14
**Scope:** the full multi-seed cluster-everything sweep #586's build deliberately deferred
to the coordinator. Build + small-scale validation (partition C, 2 seeds) passed (see
`docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md`, Part 5 addendum);
this note covers the actual full sweep run afterward.

## What was run

All 16 partitions from #583's redesign (14 single-family bands A-N, plus `DEEP_HILL`
and `BEYOND_HI_R`), 3 independent seeds each (48 GA runs total), paper-scale budget
(population=200, generations=400, 8 workers), harvested via #586's cluster-everything
pipeline: every seed's final population is clustered into distinct high-fitness basins,
each representative checked against the drift classifier (bounded vs divergent) and
matched against all 14 published Gurfil-Kasdin families (not just the partition's own
named target). Total wall time: ~5h51m (09:19-15:11).

## Result: mixed — machinery partially validated, large unadjudicated candidate pool

| Partition | Seeds recovering own family | Unmatched-bounded candidates |
|---|---|---|
| A | 0/3 | 15 |
| B | 0/3 | 17 |
| C | 3/3 | 21 |
| D | 3/3 | 0 |
| E | 2/3 | 2 |
| F | 0/3 | 24 |
| G | 3/3 | 11 |
| H | 3/3 | 12 |
| I | 0/3 | 48 |
| J | 0/3 | 5 |
| K | 2/3 | 29 |
| L | 0/3 | 5 |
| M | 3/3 | 6 |
| N | 0/3 | 12 |
| DEEP_HILL | n/a (no single target family) | 31 |
| BEYOND_HI_R | n/a (no single target family) | 26 |
| **Total** | **6/14 single-family partitions recovered in >=1 seed** | **264** |

**Machinery health (positive-control signal):** 6 of 14 single-family partitions
(C, D, E, G, H, M) recovered their own target family in at least one of 3 seeds —
confirming the corrected pipeline (bounds fix + theta0 wrap fix + multi-seed
harvesting) genuinely works when it works. The other 8 (A, B, F, I, J, K, L, N) did
NOT recover their own target family in any of 3 seeds — even with multi-seed
sampling, per-partition recovery remains unreliable under `gurfil_kasdin_fitness`'s
saturating boundedness objective, exactly as #586's design review anticipated
("recovery per partition is seed-conditional... not every seed should recover the
target family"). `DEEP_HILL` and `BEYOND_HI_R` have no single named target family by
design (they test genuinely uncharted/known-adjacent territory), so "0 recovered" for
them is not a failure signal.

**264 unmatched-bounded candidates** were surfaced across the whole sweep —
bounded (per the drift classifier), converged, and not matching any of the 14
published Gurfil-Kasdin families under the same criterion stage 2 used. **NONE OF
THESE ARE NOVELTY CLAIMS.** Per
[[feedback_literature_novelty_check_baseline]] and this project's standing
discipline: not-found-in-KNOWN_CORPUS is necessary-not-sufficient. Specifically NOT
yet done, and required before any claim:

1. **Deduplication.** 264 raw candidates across 48 independent GA runs almost
   certainly contain many near-identical repeats of the same underlying basin
   (the same physical orbit rediscovered by multiple seeds/partitions, especially
   given `DEEP_HILL`/`BEYOND_HI_R`'s 31/26 counts and `I`'s 48 -- these numbers are
   far too high to represent 264 genuinely distinct orbits). No dedup pass has been
   run.
2. **Live `check_literature()` search.** `literature_anchors_engaged()` only proves
   the offline `KNOWN_CORPUS` matcher structurally engages (a non-empty, correctly-
   scoped anchor pool) -- it does NOT run the actual live literature search. Every
   candidate's anchor list currently includes broad Sun-Earth/Earth-Mars cycler
   literature (Aldrin, Russell-Ocampo, Gurfil-Kasdin itself, Sun-Earth co-orbital,
   Henon family-f) because `primary="Sun", sequence=("E",)` pulls in the same wide
   anchor pool for every candidate regardless of its specific orbit -- this is
   necessary-not-sufficient exactly as designed, not a novelty verdict.
3. **Opus/Fable adjudication** of whatever survives dedup + live literature check,
   matching the pattern already established for #564/#565/#577 (Uranian/Galilean
   symmetric-closure survivor lists).

## Raw data

`data/found/583_widened_search/` -- per-partition-per-seed checkpoints, final
populations, and harvest summaries; per-partition aggregate harvest summaries
(`{partition}_aggregate_harvest_summary.json`). ~21 MB, 100 `.npz` + 67 `.json` files.

## Not yet dispatched

The dedup + live-literature-check + Opus/Fable adjudication pipeline for the 264
raw unmatched-bounded candidates -- a genuinely separate, substantial piece of work
(the candidate volume alone means dedup logic needs real design, not just a quick
pass), tracked as a new task rather than rushed through here.
