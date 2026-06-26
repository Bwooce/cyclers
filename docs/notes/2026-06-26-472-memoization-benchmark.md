# #472 — Per-leg cost-function memoization: benchmark + scope note

**Date:** 2026-06-26
**Task:** #472 (efficiency — memoize the deterministic per-leg cost functions to
speed the moon-tour campaigns; regression safety is the primary deliverable).

## What was memoized

All caches use `functools.cache` (unbounded; ruff `UP033` prefers it over
`lru_cache(maxsize=None)`). Every memoized function is a **pure, deterministic,
side-effect-free** scalar/string function that reads only **frozen** module
constants (`SATELLITES` / `PRIMARIES` / `PLANETS` — frozen dataclasses built
once at import). No array arguments are keyed, so the cache is **exact** (no
rounding tolerance needed).

| Module | Functions cached |
|---|---|
| `search/vilm.py` | `_vc_adim`, `_vbar_vinf_adim`, `min_vinf_for_vilm`, `_v_m`, `_leverage_dv_kms`, `_vilm_dv_min_pair` (new cached core of `vilm_dv_min`), `vilm_dv_floor`, `europa_endgame_dv` |
| `search/tisserand.py` | `_a_p_km`, `vinf_to_tisserand`, `tisserand_to_vinf`, `linkable` (the expensive contour-intersection predicate — dominant cost) |
| `core/flyby.py` | `max_bend`, `dv_from_turn_deficit`, `dv_powered_flyby_periapsis` (scalar float primitives) |
| `search/correct.py` | `_max_bend_deg_nominal` (new cached `rp_factors=None` core of `_max_bend_deg`) |

## What was deliberately left UN-cached (honesty)

- **`core/lambert.py` `lambert()`, `core/kepler.py` `propagate()` / `coe_to_rv()`** —
  these take `Vec3` (ndarray) arguments and return arrays. Caching them would
  require a lossy rounding key (changing results beyond a tolerance) AND would
  share mutable array objects that callers may mutate. The task explicitly flags
  the array-key risk; the safe choice is to not cache them. The Lambert per-solve
  scalar helpers (`stumpff_*`, `_t_of_z`, ...) are called with continuously
  varying `z` *within* a single solve, so they have ~no repeated-arg hit rate —
  no caching value.
- **`search/correct.py` `_max_bend_deg(rp_factors=...)`** (the scaled path) —
  `rp_factors` is a mutable `dict` (unhashable) AND is a per-call config
  override; caching it would alias one config's result onto another. Only the
  nominal (`rp_factors is None` / body absent) path is cached, via the split-out
  `_max_bend_deg_nominal`.
- **`flyby.flyby_dv`, `bend_angle`, `is_ballistic_feasible`, `bend_decompose`** —
  ndarray (`Vec3`) arguments; same array-key hazard.

## Purity / config-version concern

The cache keys are body strings + scalars only. The body records they read
(`SATELLITES[m]`, `PLANETS[b]`) are **frozen dataclasses** in module-constant
dicts, built once at import — there is no mutable config (no per-run
flyby-altitude override reaches the cached path). The one config-override knob,
`rp_factors`, routes through the **un-cached** `_max_bend_deg` branch, so a
config change is never served stale. This is pinned by
`tests/search/test_correct_maxbend_memoization.py::test_config_change_would_not_be_stale`.

## Regression test infrastructure (the primary deliverable)

- `tests/search/test_vilm_memoization.py`
- `tests/core/test_flyby_memoization.py`
- `tests/search/test_tisserand_memoization.py`
- `tests/search/test_correct_maxbend_memoization.py`

Each provides: **parity** (`cached == .__wrapped__` over representative + edge
inputs), **key-correctness** (distinct args → distinct entries via
`cache_info().currsize`; repeated call is a VERIFIED hit via `cache_info().hits`;
`a_range` tuple keys do not alias), and **purity** (frozen-config guards; the
`rp_factors` override is honoured, not cached). The vector entrypoints are
asserted to have **no** cache wrapper.

The ultimate regression guard is the full `tests/data tests/search` ratchet
staying green/bit-identical with caching on (the cache changes no campaign or
gauntlet result).

## Benchmark (INFORMATIONAL — not a CI gate)

`scripts/bench_472_memoization.py` runs a small #468-style nested moon-tour
sweep (6 adjacent moon legs × 20 V∞ × 4 budgets × 8 repeats = 3840 leg-evals)
twice: cache-busted-per-call vs cache-warm. It writes a flushed per-unit JSONL
progress line to `out/bench_472_runlog.jsonl`
(`{item_id, sub_step, mode, elapsed_s, ts}`) — not a black box.

Result (this machine, 2026-06-26):

```
n_leg_evals      = 3840
t_uncached_s     = 101.50
t_cached_s       = 3.62
speedup_factor   = 28.05x
parity_bit_identical = true   (cached results == uncached results, exactly)
```

The ~28x is dominated by `tisserand.linkable` (the contour-intersection
brentq/grid solve), which the moon-prune gate re-evaluates on repeated
(body-pair, V∞, a_range, mu) keys. No timing assertion is in CI (machine-
dependent); the benchmark's own in-script `assert res_cached == res_uncached`
is the bit-identity check.
