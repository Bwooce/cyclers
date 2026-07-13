# #582 (stage 3a of #581 followup): niching-GA population analysis pipeline

**Date:** 2026-07-13
**Scope of this note:** the analysis-step BUILD only. The coordinator is
running the actual paper-scale (population=200, generations=400) GA sweep
across all 5 interior MMRs directly (separately from this dispatch); this
work adds a `--mode analyze` step to `run_582_asymmetric_3d_niching_search.py`
that will turn each finished `ga_{p}_{q}_final.npz` into an adjudicated
per-MMR summary once those runs land. No `ga_{p}_{q}_final.npz` from the real
sweep was read or touched here -- everything below was validated against a
small, hand-constructed synthetic population built and deleted inside this
dispatch's own isolated worktree.

## What was built

* `scripts/run_582_asymmetric_3d_niching_search.py`:
  * `cluster_representatives()` -- greedy niche clustering over a finished
    GA population.
  * `analyze_ga_population()` -- clusters a population, routes each
    representative through the SAME mandatory pipeline the positive control
    already validated (`refine_ga_candidate` -> check `orbit.converged` ->
    `classify_symmetry` -> `build_candidate_signature` ->
    `literature_anchors_engaged`), and writes
    `data/found/582_niching_ga/{p}_{q}_analysis_summary.json`.
  * `--mode analyze` CLI wiring, plus `--fitness-floor`, `--distance-threshold`,
    `--max-clusters` overrides.

No new corrector, classifier, or fitness function -- this is pure
post-processing of an already-finished checkpoint, reusing every existing
#582 pipeline function unmodified.

## Design decisions

### 1. Clustering: greedy threshold, not scipy/sklearn

The niching GA (`search/niching_ga.py`, deterministic crowding) is explicitly
designed to hold multiple co-existing basins in one final population -- just
taking the single best-fitness member would silently discard that structure,
the entire point of using a niching method over plain differential evolution.

`cluster_representatives()` implements a simple greedy pass:

1. Sort the population by fitness, descending.
2. Stop once fitness drops below `fitness_floor` (population is sorted, so
   nothing further qualifies either).
3. Walk down the remaining list; accept a member as a NEW cluster
   representative only if its distance to EVERY already-accepted
   representative exceeds `distance_threshold`.

Distance is Euclidean in the bounds-normalized `[0,1]^6` space -- the SAME
metric `niching_ga.py`'s own parent-child pairing rule already uses internally
(see that module's `_dist`/`_Encoding.decode_norm`), not a new metric invented
for this step. `pyproject.toml` was checked: scikit-learn is not a project
dependency, and this project's own convention (per project memory) is to
avoid adding one for a step this simple -- a greedy threshold pass needs
nothing beyond numpy.

This is intentionally the SAME approach `run_583_widened_bounded_drift_search.py`
uses for its own family-matching step (nearest-neighbor RMS distance in
bounds-normalized space), generalized here from "distance to a known
published target" into "distance to an already-accepted representative"
since #582 has no published target to match against.

### 2. Fitness floor: 0.9 (default)

Same numeric choice #583's own `analyze_partition` makes, for the analogous
reason: `isolated_3d_asymmetric_fitness` is a bounded-reciprocal-penalty
objective (`1 / (1 + penalty)`) that peaks at 1.0 only when BOTH the
periodicity defect and the soft-Jacobi-band term are small, and falls off
fast moving away from any basin (quadratic penalty inside a reciprocal). 0.9
keeps genuinely basin-adjacent genomes (worth the corrector's cost) while
excluding population noise that never found a basin at all. It is
deliberately much looser than the corrector's own acceptance bar
(`orbit.converged`, which requires `corrector_residual < 1e-10` AND
`independent_closure_residual < 1e-6`) -- the corrector step immediately
downstream is the real accept/reject gate; this floor only controls how much
work `cluster_representatives()` does before that gate ever runs.

### 3. Distance threshold: 0.1 (default, bounds-normalized units)

Not paper-sourced (there is no published clustering recipe for this
un-published fitness landscape -- that gap is exactly what #582 exists to
probe). Chosen to be:

* Coarse enough that near-duplicate members of the SAME niche (a 200-member
  population run to 400 generations typically produces many genomes clustered
  tightly around its best point) collapse to one representative instead of
  each counting as its own "discovery".
* Tight enough that two genuinely different basins inside one MMR's box
  (which spans only +-15% to +-50% around the analytic e=0 guess per
  `mmr_bounds`, per-MMR) are not accidentally merged.

Exposed as `--distance-threshold` since it is a judgment call, not a derived
constant -- a future pass across the real 5-MMR sweep output may need to
retune it per-MMR (the box sizes differ: e.g. 3:2's `ydot0_frac=0.35` is
wider than the other 4 MMRs, per that function's own docstring).

### 4. Cluster cap: 25 (default, `--max-clusters`)

A defensive bound only -- keeps one `--mode analyze` invocation's wall time
proportional to "how many genuinely distinct basins survived", not to
population size, in the pathological case a badly-tuned `distance_threshold`
lets through far more clusters than expected. Every representative pays one
full `refine_ga_candidate()` Newton-iteration cost (bounded, same order as
the positive control's own single-candidate cost, ~seconds), so 25 keeps a
worst case bounded without needing to be tuned against the real sweep data
this dispatch never touches.

### 5. What's explicitly NOT done here

Per #582's own mandate and
[[feedback_literature_novelty_check_baseline]], `literature_anchors_engaged()`
only proves the STRUCTURAL matcher reaches a non-empty anchor pool -- it is
NOT the live literature search. `search/literature_check.py::check_literature()`
needs a real, injected `SearchFn` (live web search) and is explicitly left to
a future coordinator- or Opus/Fable-adjudication-owned pass. Every summary
JSON says this explicitly in its own `literature_check_status` field, rather
than leaving "not found" to be silently read as "novel" -- not-found is
necessary, not sufficient, per that same feedback note.

## Validation (this dispatch's own smoke test, not the real sweep)

Per this dispatch's own process requirement, the real coordinator-owned
`ga_{p}_{q}_final.npz` files were never read (they live in the main working
tree, this dispatch runs in an isolated worktree, and the instruction was
explicit not to wait on or touch them regardless). Validation instead used a
small, hand-constructed synthetic population for MMR 3:2, built and deleted
within this dispatch:

* **Cluster A** (3 genomes): small perturbations (~1e-3 to 1e-4 in each
  component) around the known #440 3:2 circular member's own corrected state
  (`x0=0.7310974, ydot0=0.4569989, T=11.9265124`, reused from the existing
  positive-control harness), fitness 0.999/0.998/0.997.
* **Cluster B** (1 genome): a genuinely different point elsewhere in the
  3:2 MMR bounds box (`x0=0.70, z0=0.02, xdot0=0.01, ydot0=0.30, zdot0=-0.01,
  T=10.0`), fitness 0.93.
* **14 junk genomes**: uniform-random within the bounds box, fitness
  uniform(0.05, 0.85) -- population noise that never found a basin.

Result (`--mode analyze --mmr 3:2`, default floor=0.9/threshold=0.1):

```
[3:2] population=18 n_above_floor(0.9)=4 clusters_selected=2
[3:2] cluster 0: CONVERGED fitness=0.999500 symmetric independent_closure=4.53e-11
      anchors=[... Antoniadou & Libert spatial resonant periodic orbits ..., ...]
[3:2] cluster 1: NOT CONVERGED fitness=0.930000 corrector_residual=1.158e+00 independent_closure=1.158e+00
[3:2] analysis: 2 cluster(s), 1 converged, 1 literature-matcher-ready; summary -> .../3_2_analysis_summary.json
```

Confirms:

* The 3 near-copies (cluster A) correctly collapsed into ONE representative
  (the highest-fitness member) rather than 3 separate "clusters" -- clustering
  is doing real work, not a no-op.
* Cluster B, a genuinely different point, was correctly kept as its OWN
  representative -- clustering does not over-merge distinct basins.
* All 14 junk genomes (fitness < 0.9) were correctly excluded before
  clustering even ran.
* Cluster A's representative refined cleanly (`orbit.converged=True`,
  independent Radau closure residual 4.5e-11) and was correctly classified
  `is_symmetric=True` -- it IS the known symmetric member, and the classifier
  says so.
* Cluster B's representative (an arbitrary in-box point, not a real seed) did
  NOT converge (`corrector_residual=1.158`) -- the pipeline correctly reports
  a non-convergent cluster rather than silently dropping or misreporting it.
* `literature_anchors_engaged()` fired non-empty for the converged, symmetric
  cluster, confirming end-to-end wiring through to the matcher-engagement
  check.

**Threshold-sensitivity spot-check:** re-running with
`--distance-threshold 0.0001` correctly split cluster A into 3 separate
clusters (`clusters_selected=4`) -- confirms the threshold parameter is
load-bearing and not silently ignored.

**Missing-checkpoint path:** running `--mode analyze --mmr 5:2` (no
`ga_5_2_final.npz` on disk in this worktree) prints a clear message and
returns cleanly instead of raising an unhandled exception.

All synthetic artifacts (`ga_3_2_final.npz`, its checkpoint/runlog, and the
resulting `3_2_analysis_summary.json`) were deleted after validation; `git
status` on `data/found/582_niching_ga/` is clean.

## Checks run before commit

* `uv run ruff check .` / `uv run ruff format --check .` -- clean.
* `uv run mypy src tests` -- 0 new errors (11 pre-existing errors in 4
  unrelated files, all optional-dependency import issues -- `matplotlib`/
  `pypdfium2` not installed locally; `scripts/` is not in mypy's checked
  paths, consistent with the rest of this project's `scripts/run_*.py`
  files).
* `uv run pytest tests/data tests/search tests/scripts -q` -- see the commit
  message / task report for the result.

## Reproduction

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode analyze
```

(once a real `ga_3_2_final.npz` exists from the coordinator's own `--mode ga`
sweep run).
