# #588 step 1: global dedup of the 264-candidate unmatched-bounded pool

**Date:** 2026-07-14
**Script:** `scripts/dedup_588_candidate_pool.py`
**Input:** `data/found/583_widened_search/{PARTITION}_aggregate_harvest_summary.json`
(16 files, `unmatched_bounded_candidates` lists, 264 raw candidates total —
see `docs/notes/2026-07-14-583-16-partition-3-seed-sweep-results.md` for
their origin).
**Output:** `data/found/583_widened_search/deduped_candidates.json`

## Problem

264 raw candidates were surfaced across 48 independent GA runs (16
partitions x 3 seeds), each partition fixing a different subset of the
7-gene genome (`[x, x', y, y', z, z', theta0]`) at 0. Per-partition counts
as high as `I`'s 48, `DEEP_HILL`'s 31, and `BEYOND_HI_R`'s 26 are far too
many to represent that many genuinely distinct physical orbits — almost
certainly the same handful of underlying basins, rediscovered repeatedly
by different seeds/partitions.

Raw `ga_genome` vectors are **not** directly comparable across partitions:
two candidates from different partitions can be the same physical orbit
sampled at different launch phases (wildly different raw components), or
different physical orbits that happen to land near each other in whichever
subspace happens to be free for that partition. Any dedup pass has to work
in a partition-independent, physical space instead.

## Method

1. **Characterize, don't compare genomes.** Every candidate's genome is
   re-integrated through `run_581_gurfil_reproduction.py::characterize()`
   (imported unmodified — the exact function both `run_581_*` and
   `run_583_*` already use for their own family-matching), producing
   `rmin_km_1yr`, `rmax_km_1yr`, and a `type` classification (DRO / DPO /
   DEO / ERO, optionally `3D `-prefixed). `characterize()` also reports
   5-year values but exposes no orbital period; only the 1-year
   rmin/rmax + type were used as the dedup feature space (the same feature
   set `match_family_in_widened_population` itself already keys its
   family-match tolerance on).
2. **Never merge across `type`.** A DRO and a DEO are different dynamics
   regardless of radial proximity.
3. **Greedy fitness-ranked clustering**, same structure as
   `cluster_representatives`/`cluster_population` elsewhere in this
   pipeline: walk candidates in descending `ga_fitness` order, accept a
   new cluster only if no existing cluster of the same `type` has both
   `rmin_km_1yr` and `rmax_km_1yr` within `DEDUP_RELATIVE_TOLERANCE`
   (default 0.10, i.e. 10% relative deviation) of the candidate.
4. **Representative selection**: highest `ga_fitness` within the cluster,
   ties broken by lower `drift.trend_fraction` (the more stably bounded
   member).

## Threshold justification (0.10 relative deviation)

* **Precedent.** This exact pipeline already uses 0.10 twice —
  `DEFAULT_ANALYZE_DISTANCE_THRESHOLD` (#582) and
  `DEFAULT_HARVEST_DISTANCE_THRESHOLD` (#583) — both documented as "coarse
  enough to collapse near-duplicate members of one niche, tight enough to
  keep genuinely distinct basins separate," applied to bounds-normalized
  genome-space distance within one seed's population. Reusing the same
  numeric philosophy for the cross-partition physical-space pass keeps the
  two dedup stages conceptually consistent.
* **Empirical chaining check.** A greedy single-link merge can silently
  chain far-apart points through a series of overlapping near-duplicates.
  Scanning tolerances 0.02–0.50 against the real 264-candidate pool:

  | tol  | clusters | worst within-cluster spread (rmin, max/min) |
  |------|----------|-----------------------------------------------|
  | 0.05 | 61       | — |
  | 0.10 | **45**   | ~1.18x |
  | 0.15 | 37       | ~1.26x |
  | 0.20 | 30       | ~1.44x (chaining starts to dominate) |
  | 0.30 | 24       | — |

  At 0.10 the worst observed within-cluster spread (~1.18x) stays close to
  what direct pairwise proximity alone would justify — chaining is not
  measurably inflating clusters. At 0.20 the worst cluster's internal
  spread (1.44x) already exceeds what a single 20%-tolerance hop should
  allow, i.e. chaining is doing real work by that point. 0.10 was kept as
  the more conservative of the two viable choices: under-merging only
  costs a handful of redundant downstream literature-check calls, while
  over-merging could silently drop a candidate that deserved its own
  check.

No period/T metric was included (not returned by `characterize()`); adding
one would need a genuine period-detection pass, out of scope for this
mechanical dedup step.

## Result

**264 raw candidates → 45 distinct clusters (5.87x dedup ratio).**

By type: DRO 15, 3D DRO 15, DEO 8, ERO 3, DPO 2, 3D ERO 2 (45 total).

Largest clusters (by `duplicate_count`):

| cluster | type | dup_count | partitions |
|---|---|---|---|
| 0 | DRO | 41 | A, B, C, DEEP_HILL, F |
| 1 | DEO | 25 | B, BEYOND_HI_R, F |
| 2 | DEO | 13 | A, B, F, I |
| 3 | DEO | 13 | A, B, C, F |
| 4 | ERO | 12 | H (single-partition) |
| 5 | 3D ERO | 12 | N (single-partition) |

9/45 clusters span more than one partition (accounting for a
disproportionate share of the raw 264 — the top 4 cross-partition
clusters alone absorb 92 of the 264 raw candidates). The most frequent
cross-partition pairings are `A<->F`, `B<->F`, `A<->B`, `A<->C`, `C<->F`
(5 shared clusters each) — i.e. partitions A, B, C, F (all "paper band,"
1-free-position-component families A-F) are frequently rediscovering the
same handful of basins from each other, as expected since they share
adjacent position-bound boxes.

`DEEP_HILL` and `BEYOND_HI_R` behave differently from each other here:
* `DEEP_HILL` is almost entirely self-contained — only 1 of its 31 raw
  candidates merges cross-partition (into cluster 0's 41-member DRO
  basin, alongside A/B/C/F); the other 30 collapse into 4 `DEEP_HILL`-only
  clusters (10/9/8/3 members). Its own basin is largely distinct from the
  paper-band partitions.
* `BEYOND_HI_R` is the opposite — 22 of its 26 raw candidates merge into
  cluster 1 (the 25-member DEO basin shared with B and F), and only 4
  remain in 2 `BEYOND_HI_R`-only clusters. Most of `BEYOND_HI_R`'s
  "novel-looking" 26 candidates are in fact repeat rediscoveries of a
  basin B and F's own seeds already found.

36/45 clusters (80%) are single-partition, confirming most of the
"implausibly high" per-partition counts really were repeat rediscoveries
within that one partition's own seeds, not genuinely new diversity —
though as the `BEYOND_HI_R` case shows, some of the inflation is
cross-partition rediscovery rather than purely within-partition.

## Explicitly not done here

Per this dispatch's scope (#588 step 1 only): no live `check_literature()`
web search was run, and no novelty claim is made about any of the 45
surviving clusters. `data/found/583_widened_search/deduped_candidates.json`
is provenance-complete (each cluster's `source_locations` lists every raw
partition/seed/cluster_rank that merged into it) and ready for the
live-literature-check + Opus/Fable adjudication steps #588 still owns.
