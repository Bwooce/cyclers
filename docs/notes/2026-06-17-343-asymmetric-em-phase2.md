# 2026-06-17 — #343 asymmetric-corrector EM scan Phase 2 (parallel + prioritizer)

Re-runs the #284 asymmetric-corrector novel-cycler scan with two substrate
upgrades that landed since Phase 1 (commit `ebb9e6c`):

* **#310 single-orbit prioritizer adapter** (`src/cyclerfinder/search/single_orbit_prioritizer.py`)
  closes the "prioritizer needs pair-shaped inputs" gap that blocked Phase 1
  from ranking surviving candidates.
* **#321 parallel-sweep substrate** (`src/cyclerfinder/parallel/parallel_sweep.py`)
  provides a joblib-backed `parallel_sweep(cells, closure, config=...)` wrapper
  with measured 5-8x speedup on embarrassingly-parallel single-cell sweeps.

Both wired in via a new top-level pickle-safe closure
(`src/cyclerfinder/search/asymmetric_novel_scan_parallel.py`) + a chunked
driver (`scripts/scan_343_asymmetric_em_phase2.py`).

## Grid coverage

| axis | Phase 1 (#284) | Phase 2 (#343) | factor |
| --- | --- | --- | --- |
| `(k1, k2)` bands | 5 | 9 | 1.80x |
| Jacobi `C` levels | 4 | 6 | 1.50x |
| `x0` seeds | 4 | 6 | 1.50x |
| `xdot0` seeds | 2 | 3 | 1.50x |
| `ydot0` signs | 1 | 1 | 1.00x |
| `half_crossings` | 2 | 2 | 1.00x |
| **cells** | **320** | **1,944** | **6.08x** |

`(k1, k2)` bands added in Phase 2: `(3, 4)`, `(4, 3)`, `(5, 2)`, `(2, 5)` —
the higher-order resonances #284 explicitly flagged as under-covered.

`C` grid widened to span both lower (`3.06`) and higher (`3.19`) than Phase 1's
`[3.08, 3.17]` window.

`x0` grid finer: 6 seeds across `[-0.92, -0.65]` vs Phase 1's 4 across
`[-0.90, -0.68]`.

`xdot0` adds the zero seed (`0.0`); both negative and positive seeds were
present in Phase 1.

`ydot0_sign` and `half_crossings` held at Phase 1 settings. A first Phase 2
trial expanded these to 2 signs × 3 hc (5,832 cells) but the corrector cost
on `(k1, k2) ∈ {(5, 2), (2, 5), (3, 4), (4, 3)}` × `hc = 8` cells under live
CPU contention ran past the harness budget. The 1,944-cell grid keeps the
`(k1, k2)` expansion — the axis Phase 1 most under-sampled — and trims the
two axes that contributed little new converged territory in Phase 1.

## Runtime

* **Wall time:** 1,483 s (~25 min) on 14 loky workers, chunk_size=4.
* **Per-cell summed cost (would-be-serial):** 18,594 s.
* **Speedup:** 12.5x. Above the 5-8x target; favoured by long-running
  hc=6 corrector calls where the parallel-IO overhead is negligible.
* **Cells/sec aggregate:** 1.31. Phase 1 was 320 cells in 4,810 s = 0.067
  cells/sec serial; Phase 2 at 1.31 cells/sec wall is 19.5x faster per cell
  on the wall clock.

## Guard chain results

| stage | count | notes |
| --- | --- | --- |
| cells attempted | 1,944 | |
| closure raised | 0 | (no per-cell exceptions) |
| converged + dedup-unique | 12 | (Phase 1: 8 after dedup) |
| topology match `(k_target == k_classified)` | 0 | (Phase 1: 0) |
| DOP853 independent closure `< 1e-6` | 12 | (Phase 1: 9) |
| literature-fresh (no known-family period match) | 12 | (Phase 1: 8) |
| ML `p_fp < 0.5` | 12 | (Phase 1: 9) |
| SILVER (topology + closure + ML) | 0 | (Phase 1: 0) |
| novelty-claimable (SILVER + lit-fresh) | 0 | (Phase 1: 0) |
| known-family collisions | 0 | (Phase 1: 1 — C21 reproduction) |

## Per-(k1, k2) band breakdown

Of the 9 requested `(k1, k2)` bands, only `(1, 1)` produced any converged
orbits at the seed grid:

| `(k1, k2)` target | cells attempted | converged (dedup) | SILVER |
| --- | --- | --- | --- |
| `(1, 1)` | 216 | 12 | 0 |
| `(1, 2)` | 216 | 0 | 0 |
| `(2, 1)` | 216 | 0 | 0 |
| `(3, 2)` | 216 | 0 | 0 |
| `(2, 3)` | 216 | 0 | 0 |
| `(3, 4)` | 216 | 0 | 0 |
| `(4, 3)` | 216 | 0 | 0 |
| `(5, 2)` | 216 | 0 | 0 |
| `(2, 5)` | 216 | 0 | 0 |

The asymmetric corrector — with seed `(x0, ±0.05, 0, 0)` — collapses to the
symmetric basin (`xdot0 = 0`, asymmetry = 0). Phase 1 observed the same
collapse. Phase 2's 6x larger cell budget and wider `(k1, k2)` band coverage
yields the same outcome: the symmetric basin attractor dominates every
converged orbit. All 12 unique converged orbits classify as `(2, 0)` or
`(3, 0)` — symmetric distant-retrograde-style orbits at progressively lower
periods as `C` rises:

| C | k_classified | period (d) | x0 |
| --- | --- | --- | --- |
| 3.06 | (2, 0) | 53.90 | -0.413 |
| 3.06 | (3, 0) | 80.84 | -0.413 |
| 3.10 | (2, 0) | 53.52 | -0.453 |
| 3.10 | (3, 0) | 80.28 | -0.453 |
| 3.13 | (2, 0) | 52.98 | -0.487 |
| 3.13 | (3, 0) | 79.48 | -0.487 |
| 3.15 | (2, 0) | 52.37 | -0.511 |
| 3.15 | (3, 0) | 78.55 | -0.511 |
| 3.17 | (2, 0) | 51.34 | -0.534 |
| 3.17 | (3, 0) | 77.02 | -0.534 |
| 3.19 | (2, 0) | 49.71 | -0.553 |
| 3.19 | (3, 0) | 74.57 | -0.553 |

## Single-orbit prioritizer attachment

The `score_single_orbit` adapter ran on every converged candidate. The
surrogate-pair strategy found `braik-ross-c32-cycler-2026` as the nearest
catalogue pair-shape neighbor for every cell (typical tuple-distance
`~3.5e-2` — same energy-level basin).

No tier produced a finite score for these candidates:

* **Tier 0 (NN reachability)** — skipped. The NN expects heliocentric SI
  state; CR3BP candidates carry only rotating-frame state, so the prefilter
  falls back to `model_available=False` and returns no `predicted_dv_kms`.
* **Tiers 1-2 / 3 / 4 / 5** — skipped. The default `FiveTierPrioritizer()`
  auto-builds only the NN prefilter; the other scorers are configuration-
  heavy and were not wired by this driver (out of scope per the task brief).

This is the honest verdict on `FiveTierPrioritizer()`'s default state for a
CR3BP single-orbit candidate. The adapter is attached and operating: it
correctly identifies the surrogate-pair neighbor and invokes the per-tier
hooks. The seam is open for future tasks to wire per-tier scorers (each is
config-heavy and worth its own gauntlet pass).

The driver records the prioritizer notes per row so the seam is auditable.

## Comparison to Phase 1

Phase 1 produced 0 novelty-claimable; Phase 2 reproduces that verdict at 6x
cell coverage and with the prioritizer adapter attached.

Phase 2 picked up 12 unique converged orbits vs Phase 1's 8 — the wider C
grid (3.06 / 3.19 endpoints not in Phase 1) reveals two new (2, 0)/(3, 0)
families at C ∈ {3.06, 3.19}, and the finer x0 grid resolves at C = 3.15
that Phase 1 missed.

Phase 1's one known-family collision (C21 at period 83.96d, C=3.08, asym=0.16)
does NOT reappear in Phase 2's converged set. Looking at Phase 1's record,
that hit was at seed `(x0=-0.82, xdot0=-0.05, sign=-1.0, hc=4)`; Phase 2's
grid passes through `x0 ∈ {-0.80, -0.85}` ≠ `-0.82` and `(C=3.06, 3.10) ≠
3.08`. The asymmetric (1, 1) → C21 corrector seems to need a basin window
narrower than the Phase 2 grid resolution — a Phase 3 zoom around the Phase
1 hit could re-discover it, but it's already in the literature (Braik-Ross),
so this is interesting only as a basin-sensitivity datapoint, not a science
gain.

The novelty-discovery conclusion is unchanged from Phase 1: the asymmetric
corrector at this seed convention has a symmetric-basin attractor for every
cell that converges, and the topology gate rejects every (2, 0) / (3, 0)
result against the requested `(k1, k2)` ≠ `(2, 0)` / `(3, 0)`. **Clean
negative on novelty at the Phase 2 coverage.**

## Phase 3 recommendation

The basin behaviour is the rate-limiting step, not the cell count. The
parallel substrate gives us 12.5x more cells/wall, but the corrector keeps
landing in the symmetric basin regardless of seed. Two distinct Phase 3
directions are available:

1. **Different seed convention.** Seed `xdot0 ≠ 0` with `|xdot0| > 0.05`
   away from the symmetric basin (e.g. `|xdot0| ∈ {0.1, 0.2, 0.3}`). Risk:
   the asymmetric basin may be very thin and the corrector may diverge.
   Trial under #287's 3D-corrector-spike approach is the structured way to
   probe this.
2. **Different corrector entirely.** The `correct_general_periodic` Newton
   loop targets the closure manifold; for genuine `(k1, k2) ≠ (m, 0)`
   asymmetric orbits, the literature uses continuation from a known
   symmetric family + an asymmetry-breaking parameter (Braik-Ross Table 1
   uses Floquet bifurcations to land asymmetric orbits). Reproducing that
   continuation pipeline is out of scope for #343 but worth scoping as a
   separate task.

**Recommendation:** no immediate Phase 3 in the asymmetric-corrector
direction. The clean negative result on novelty is conditional on the
seed convention; expanding the seed convention without expanding the
corrector basin won't change the verdict.

## Negative-results registry note

Per `feedback_negative_results_registry`, this scan's coverage should be
recorded against the asymmetric-corrector method version. The (method,
grid) tuple is:

* `method = asymmetric_corrector(23b980e) + prioritizer(310) + parallel(321)`
* `grid = 9 (k1, k2) bands x 6 C x 6 x0 x 3 xdot0 x 1 sign x 2 hc = 1944 cells`
* `coverage = (k1, k2) ∈ {(1,1), (1,2), (2,1), (3,2), (2,3), (3,4), (4,3), (5,2), (2,5)};
  C ∈ [3.06, 3.19]; x0 ∈ [-0.92, -0.65]; xdot0 ∈ {-0.05, 0.0, 0.05};
  ydot0_sign = -1.0; half_crossings ∈ {4, 6}`
* `verdict = 12 converged + 0 SILVER + 0 novelty-claimable`

A future scan with a method capability that subsumes this (different
corrector / different seed convention / different topology gate) is the
trigger for re-sweeping this grid.

## Files

* `scripts/scan_343_asymmetric_em_phase2.py` — driver (chunked parallel,
  streaming-write per batch).
* `src/cyclerfinder/search/asymmetric_novel_scan_parallel.py` — top-level
  pickle-safe per-cell closure.
* `data/scan_343_asymmetric_em_phase2.jsonl` — 12 unique converged rows
  sorted by prioritizer Tier 0 (all `None`) then by `(k1, k2)`.
* `data/scan_343_asymmetric_em_phase2.stream.jsonl` — partial per-batch
  stream (matches the canonical file post-dedup; useful for recovery if
  the driver is killed mid-scan).

## Discipline

NO catalogue writeback. NO novelty claims. The single converged k_target =
(1, 1) band's 12 orbits all classify as (2, 0) or (3, 0), so the topology
gate correctly rejects every cell. Prioritizer adapter attached and
operating; per-tier scoring requires per-tier config out of scope for #343.
