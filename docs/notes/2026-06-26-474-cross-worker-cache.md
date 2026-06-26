# #474 — Cross-worker effectiveness of the #472 per-leg caches

**Date:** 2026-06-26
**Status:** landed
**Depends on:** #472 (commit 54f8252, `functools.cache` on the deterministic per-leg cost functions)

## Problem

#472 memoized the deterministic per-leg cost functions (`vilm._vilm_dv_min_pair`,
`vilm.vilm_dv_floor`, `vilm.min_vinf_for_vilm`, `vilm.europa_endgame_dv`,
`tisserand.linkable`, `tisserand.vinf_to_tisserand`, …, `flyby.max_bend`,
`flyby.dv_*`, `correct._max_bend_deg_nominal`) with `functools.cache`, giving a
~28× single-process speedup (see `docs/notes/2026-06-26-472-memoization-benchmark.md`).

`functools.cache` is **per-process**. The discovery campaigns fan out under
joblib. The default `loky` backend spawns **fresh interpreters** (spawn +
cloudpickle), so every worker starts with a COLD cache and rebuilds it from
scratch — the #472 speedup largely evaporates the moment a campaign goes
parallel.

## Approach chosen — Option A (fork + COW + opt-in pre-warm)

Empirically verified the backend behaviour on this host (Linux, default start
method `fork`):

```
parent:  CacheInfo(currsize=6)          # warmed before fan-out
loky:             worker currsize = 0   # COLD — fresh interpreter
multiprocessing:  worker currsize = 6   # WARM — forked, COW-inherited
```

joblib's `"multiprocessing"` backend (already a `ParallelSweepConfig.backend`
`Literal`) forks the pool on Linux. Forked children share the parent's pages
copy-on-write, so **if the parent's `functools.cache` tables are populated before
the fork, every worker inherits them for free** (read-only cache pages shared,
not duplicated — ~free with the host RAM). This is strictly less new machinery
than Option B (joblib.Memory disk cache), and the cached values never leave RAM.

### Implementation

1. **`ParallelSweepConfig.prewarm: Callable[[], None] | None = None`**
   (`src/cyclerfinder/parallel/parallel_sweep.py`). `parallel_sweep` runs it
   ONCE in the parent before fanning out. `None` (default) is a no-op — every
   existing caller is byte-for-byte unchanged. The hook is the substrate's only
   change; it stays domain-agnostic.

2. **`cyclerfinder.search.cache_warm.warm_moon_leg_caches`**
   (`src/cyclerfinder/search/cache_warm.py`) drives the public composite
   `moon_prune.moon_leg_admissible` over the campaign's discrete
   `(leg, V∞, budget)` domain, which **transitively** populates every #472 cache
   — so it cannot drift out of sync when a new cached function is added
   downstream. It is pure (reads only the frozen `SATELLITES` / `PRIMARIES` /
   `PLANETS` tables) and best-effort (a missing body in one cell never aborts the
   rest).

Wire-up at a campaign call site:

```python
cfg = ParallelSweepConfig(
    backend="multiprocessing",  # forks on Linux -> COW-shared caches
    prewarm=functools.partial(
        warm_moon_leg_caches, legs=legs, vinf_grid=vinf_grid, budget_grid=budget_grid,
    ),
)
```

## Bit-identical proof

`prewarm` only reads frozen module constants and only *populates* a cache whose
values are byte-identical to the cold-built ones — it changes no result. Proven by:

- `tests/parallel/test_cross_worker_cache.py`:
  - **PARITY** — the multiprocessing+prewarm sweep, the loky+prewarm sweep, and
    the `prewarm=None` legacy sweep each return results **bit-identical** to a
    direct single-process reference over the representative + edge domain.
  - **SHARING** — under multiprocessing+prewarm each forked worker reports a
    NON-EMPTY `vilm.vilm_dv_floor.cache_info().currsize` at its first call
    (inherited via COW); the **positive control** asserts loky workers report
    `currsize == 0` (cold) with the same prewarm hook, proving the SHARING test
    measures a real effect, not a tautology.
- `scripts/bench_474_cross_worker_cache.py` carries an in-script PARITY assertion
  across all three passes + the serial reference.

## Measured parallel result

`scripts/bench_474_cross_worker_cache.py` (960-cell moon-tour skeleton, 4 workers;
runlog `out/bench_474_runlog.jsonl`):

| pass                          | wall (s) |
|-------------------------------|---------:|
| loky, no prewarm (status quo) |   ~9.70  |
| multiprocessing, no prewarm   |   ~0.51  |
| multiprocessing, prewarm      |   ~0.54  |

- **prewarm-vs-loky ≈ 18×.** The dominant status-quo cost is loky's
  fresh-interpreter spawn + re-import per worker, *on top of* the cold cache.
- **prewarm-vs-mp-cold ≈ 1.0** *for this skeleton.* The per-cell #472 cache build
  is sub-second here, so once you fork (cheap), the cold-vs-warm cache difference
  is in the noise. The COW pre-warm benefit **scales with the cache-build cost
  relative to the per-cell work** and with the number of worker re-spawns
  (joblib recycles workers): it is the safety margin that keeps the #472 speedup
  from evaporating when the cache build is expensive, while fork is what recovers
  the bulk of the parallelism.

**Practical recommendation for campaign callers:** switch CPU-bound moon-tour
sweeps to `backend="multiprocessing"` (Linux); add the `prewarm` hook whenever the
per-leg cache build is non-trivial relative to the per-cell work.

## Follow-up (not done)

- **Option B (joblib.Memory disk cache under `out/cache_472/`)** remains available
  if a campaign must survive *across re-runs* or runs on a non-fork platform
  (macOS/Windows default to spawn, where COW does not apply). Lower speedup than
  COW-RAM but robust. Noted, not implemented — Option A covers the Linux campaign
  hosts.
- Rewiring specific campaign entry points to pass the multiprocessing backend +
  prewarm hook is a mechanical follow-up; this task delivers the substrate + the
  warmer + the proof.
```
