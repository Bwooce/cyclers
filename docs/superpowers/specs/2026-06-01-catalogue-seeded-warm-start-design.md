# Catalogue-seeded warm-start for the idealised optimiser (#52)

**Status:** design approved (brainstorming, 2026-06-01) → implementing.
**Scope decision:** Aldrin gate only. **Coupling:** caller passes warm starts.
**ToF mapping:** absolute cumulative, clipped.

## Problem

`optimise_cell_idealized` seeds the free interior encounter epochs from
`_free_return_seed`, which places them at equispaced fractions `i·T/(N-1)`
of the period. For an `N=3` E-M-E cell that puts the single interior (Mars)
encounter at `T/2 ≈ 390 d`. The published Aldrin classic cycler has an
asymmetric `146 d` Earth→Mars transit, so the free-return seed (plus the
5-start grid and the DE pass) converges to an alternate ~10.6 km/s family
rather than Aldrin's published 9.7 km/s. The two entries
`aldrin-classic-em-k1-{outbound,inbound}` are therefore parked in
`EXPECTED_SKIPS` of `tests/test_catalogue_rediscovery.py`.

The fix is to let a caller hand the optimiser one or more *warm starts* —
interior-epoch vectors derived from a known geometry — so SLSQP polishes
from inside the correct basin.

## Design

### 1. Optimiser API (mechanism, catalogue-agnostic)

Add a keyword-only parameter to `optimise_cell_idealized`:

```python
warm_starts: Sequence[Sequence[float]] | None = None
```

Each warm start is a vector of **interior** encounter epochs in seconds,
length `N-2` (same shape `_multi_start_grid` returns — endpoints `0` and
`T` are pinned and added back by `_build_cycler_from_x`). Warm starts are
polished through the existing `_polish` path *before* the grid, then merged
into the same `records` list `_select_best` already ranks. They are clipped
to the strict interior and sorted (the same monotonicity `construct_cycler`
requires) defensively, so a slightly out-of-range caller value cannot make
the optimiser raise. A warm start of the wrong length is a caller bug and
raises `ValueError`.

**Opt-in and inert by default:** `warm_starts=None` leaves the start set,
ordering, and RNG draws bitwise-unchanged, so every existing test —
including the load-bearing M5 gate `test_2syn_em_rediscovers_5_65_kms_earth`
— is unaffected.

### 2. Pure ToF→epoch converter (the mapping the caller applies)

A new pure helper in `optimize.py` (no catalogue import):

```python
def interior_epochs_from_leg_tofs(
    leg_tofs_days: Sequence[float], target_period_sec: float
) -> tuple[float, ...]
```

Absolute cumulative mapping: encounter `j` sits at the cumulative sum of
the first `j` leg ToFs (seconds). Encounter 0 is pinned at `t=0` and the
final encounter at `t=T`, so the converter returns the interior cumulative
sums only — `cumsum(leg_tofs[:-1])` — each clipped to the optimiser's
strict interior via the existing `_clip_interior`, then sorted. For Aldrin
(`legs=[146, 634]`) this yields a single interior epoch at `146 d`, exactly
the basin the free-return seed misses.

`target_period_sec` is resolved by the caller the same way the optimiser
does, via `_target_period_sec(cell)`, so the warm start clips against the
same bounds the polish enforces.

### 3. Wiring (Aldrin gate only)

The caller for this scope is the rediscovery harness
`test_catalogue_entry_rediscovers`. It already builds the cell from the
entry and has `entry.leg_tofs_days`. It will build one warm start from those
ToFs and pass it as `warm_starts=[...]`. The two Aldrin entries are removed
from `EXPECTED_SKIPS`, so the gate now *requires* them to rediscover the
published 9.7 km/s signature within `VINF_TOL_KMS`.

Production `discover()` is **not** wired in this scope — matching cells to
catalogue entries pre-optimisation is a larger feature deferred beyond #52.

## Out of scope / unchanged

- `test_discover_em_k2_yields_known_for_2syn` stays `xfail` — #54 showed the
  5.65 km/s Russell cycler is a 4-encounter E-E-M-M topology unreachable by
  a 3-encounter cell; warm-starting cannot fix a structural mismatch.
- The `_multi_start_grid` sign-collapse hack (#53) is untouched.

## Verification

- `interior_epochs_from_leg_tofs` unit test (Aldrin maps to ~146 d interior;
  clipping pins out-of-range values).
- `optimise_cell_idealized` accepts/uses warm starts; `None` path unchanged.
- Full suite stays green; M5 gate still passes.

## Outcome / finding (2026-06-01) — premise falsified

The mechanism is implemented and unit-tested, but **it does not close the
Aldrin rediscovery gate**, and the #54 hypothesis that seeding would fix it
is empirically false:

- Building the cycler at the published interior epoch (Mars at `t = 146 d`)
  gives `V∞ ≈ 21.5 km/s`, *worse* than the free-return midpoint's 6.09 km/s.
- A 1-D sweep of the single interior epoch over `(0.02T, 0.98T)` never finds
  a closing sub-cap solution: the lowest max-`V∞` is **6.09 km/s at `T/2`**
  with closure residual **21.8 km/s**. No feasible Aldrin basin exists in the
  circular-coplanar model at *any* seed.

Root cause: Aldrin's 146-day Earth→Mars transit is only cheap with Mars's
**real eccentricity** (e ≈ 0.093). The circular-coplanar idealisation removes
exactly that, making the short leg near-hyperbolic. This is a **model**
limitation, not a search/seeding one — it needs the M6b real-ephemeris
optimiser.

### What shipped

1. Warm-start **mechanism** in `optimise_cell_idealized` + the pure
   `interior_epochs_from_leg_tofs` converter, both unit-tested. Inert by
   default (`warm_starts=None`).
2. The rediscovery harness wires warm starts opportunistically, **guarded** to
   only fully-tabulated entries (`len(leg_tofs) == N-1`) so under-tabulated
   catalogue rows (the common case — return loop not recorded) fall back to
   grid+DE instead of crashing on a malformed seed.
3. The two Aldrin entries are **kept** in `EXPECTED_SKIPS`, with the reason
   corrected from "needs seed-from-published" to the verified model-limitation
   diagnosis. The `test_aldrin_regression_anchor` / `test_2syn_…` xfail reasons
   are likewise corrected (warm-start no longer cited as their fix).

#52 (the mechanism) is done; the Aldrin gate is re-classified as an
ephemeris-mode (M6b) dependency, not a #52 deliverable.
