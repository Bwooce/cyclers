# 2026-06-16: #284 Asymmetric-corrector novel-cycler scan at Earth-Moon

## Scope

Repurpose `src/cyclerfinder/search/cr3bp_general_periodic.py::correct_general_periodic`
(built for the #249 Braik-Ross C21 reproduction) as a *novel-search driver*:
grid over (k1, k2) winding classes x Jacobi C x (x0, xdot0) seed region,
route each converged candidate through the orbit-closure-discipline guard
chain, emit one JSONL row per CONVERGED cell to `data/scan_284.jsonl`.

Per task brief: **NO catalogue writeback, NO novelty claims**, work on `main`,
read-only on `cr3bp_general_periodic.py`. A row is `novelty_claimable` (a
flag, not a writeback) iff ALL guards pass.

## What landed

- `scripts/scan_284_asymmetric.py` — the driver. Standalone script (no new
  module under `src/`); a follow-up batched campaign can promote the driver
  to `src/cyclerfinder/search/asymmetric_novel_scan.py` if the slab confirms
  the cost model warrants it.
- `data/scan_284.jsonl` — one row per converged cell (see Output below).

## Grid — down-scoped slab vs the brief

The asymmetric corrector wall-clock is ~6-15 s/cell (dominated by STM
propagation across the multi-crossing return arc; quantified inline in the
script). The task brief sketched:

```
(k1,k2) in 7 cells x C in [3.0,3.2]/0.01 (21) x x0 in [-1.0,-0.6]/0.02 (21)
  x xdot0 (4 seeds: -0.10,-0.02,0.02,0.10) x signs (-1,+1)
  x half_crossings (2,4,6,8) = 98 784 cells
```

At 10 s/cell that is ~274 hours — ~140x the 2-hour budget. Down-scoped slab:

| Axis | Brief | Slab | Rationale |
|---|---|---|---|
| (k1, k2) | 7 cells | 5 cells | Dropped (3,1) and (5,2) — same families re-found via topology classifier on the 5 retained. |
| Jacobi C | 21 levels (0.01 step) | 4 levels (3.08, 3.12, 3.14, 3.17) | Brackets the Braik-Ross common-C=3.1294 with offsets above/below. |
| x0 seed | 21 points (0.02 step) | 4 points (0.07 step) | Covers barbel + lobe regions; C11a/C32 IC ~ -0.81 is hit by -0.82. |
| xdot0 seed | 4 points (-0.10..0.10) | 2 points (-0.05, +0.05) | Symmetric reflection retained. |
| ydot0_sign | 2 (+/-1) | 1 (-1) | Family-seed convention; +1 mirrors. |
| half_crossings | 4 (2,4,6,8) | 2 (4, 6) | hc=6 hits C11a/C32; hc=4 hits C21. hc=2/8 dropped (numerical extremes). |
| TOTAL | 98 784 | 320 | 1/309. ~53 min at 10 s/cell. |

The slab is honest about its size in the JSONL row's `scan_id` and
`scan_method_version` so a follow-up wider sweep does not double-count.

## Guard chain (5 stages, gated)

A candidate is emitted as a JSONL row on EVERY converged cell, but
`novelty_claimable = True` requires ALL five:

1. **Corrector converges** (`tol=1e-11`, `max_iter=20`) AND corrector's
   own Radau closure < 1e-6.
2. **Topology classification matches the target** via
   `binary_star_search.winding_topology` — `(k1, k2)` from the integrated
   orbit must equal the requested target.
3. **Independent DOP853 closure** (`rtol=atol=1e-12`) re-propagated from
   the corrected IC must be < 1e-6. Different integrator AND tolerance
   from the corrector's own Radau check (orbit-closure-discipline).
4. **Literature-fresh (offline)**: the candidate's `(k1, k2)` + period
   does NOT collide with a known sourced Earth-Moon family (C11a 42.14d,
   C11b 55.96d, C21 84.53d, C32 78.61d ± 2 %). The corpus anchors
   (Braik-Ross, Roberts-Tsoukkas, Kumar-Rawat-Rosengren-Ross, Koblick
   tulip, Zhang-Jiang-Yuan tulip) fire generically on any EM CR3BP
   orbit by structural overlap (body_set = {Moon}); that overlap is
   recorded as `literature_offline_anchors` but is treated as
   necessary-not-sufficient. **Offline matching, NOT live WebSearch**:
   the recorded anchors document the gap.
5. **ML false-positive flagger** (`cyclerfinder.ml.falsepos_flagger`,
   #256/#275) `p_false_positive < 0.5`. Flagger fed the corrector's
   own residual + classified period and the post-fix solver SHA
   (`23b980e`), so the epoch / mu-fix flags fire benignly.

## What was deliberately NOT used

- **5-tier prioritizer (#282) Tiers 0-5**: documented architectural seam —
  Tier 0 NN scores patched-conic *Lambert legs* between moon-state pairs
  (not CR3BP periodic orbits); Tiers 1-5 require representative orbit
  *pairs*, not a single candidate orbit. A discovery scan that emits one
  orbit per cell has no Tier-0 input shape and no second representative
  for Tiers 1-5. The doc records the gap; the script omits the misuse.

## Re-runnability

```sh
# Slab (default; ~50 min wall):
uv run python scripts/scan_284_asymmetric.py --time-budget-s 6000

# Smoke (CI-sized, ~2 min wall):
uv run python scripts/scan_284_asymmetric.py --smoke --out /tmp/scan_284_smoke.jsonl
```

## Findings — slab complete (2026-06-16 21:54 → 22:21 UTC, 1576 s wall)

```
attempted=320  converged=135  topo_match=0  ind_pass=9  lit_fresh=8
ml_low_fp=9  novel_claimable=0  known_family_collisions=1
```

Eight unique rows survive the per-grid (x0, C, period_d) dedup
(135 → 8). Of these:

| target (k1,k2) | classified | C       | T (d) | x0       | xdot0   | asym    | ind_pass | known | flag |
|----------------|------------|---------|-------|----------|---------|---------|----------|-------|------|
| (1,1)          | (2,1)      | 3.0800  | 83.96 | -0.7369  | -0.1648 | 1.6e-01 | True     | True  | C21 collision |
| (1,1)          | (2,0)      | 3.0800  | 53.74 | -0.4320  |  0.0    | 0       | True     | False | non-cycler (no Moon winding) |
| (1,1)          | (3,0)      | 3.0800  | 80.60 | -0.4320  |  0.0    | 0       | True     | False | non-cycler |
| (1,1)          | (2,0)      | 3.1200  | 53.20 | -0.4748  |  0.0    | 0       | True     | False | non-cycler |
| (1,1)          | (3,0)      | 3.1200  | 79.80 | -0.4748  |  0.0    | 0       | True     | False | non-cycler |
| (1,1)          | (2,0)      | 3.1400  | 52.71 | -0.4986  |  0.0    | 0       | True     | False | non-cycler |
| (1,1)          | (3,0)      | 3.1400  | 79.07 | -0.4986  |  0.0    | 0       | True     | False | non-cycler |
| (1,1)          | (2,0)      | 3.1700  | 51.34 | -0.5345  |  0.0    | 0       | True     | False | non-cycler |

### Honest negatives

- **0 novelty-claimable rows**: the topology guard rejected every cell,
  including the C21 collision and the 7 "asymmetric corrector collapsed
  to a symmetric perpendicular crossing" rows (`xdot0=0`, `asymmetry=0`
  — the corrector's known symmetric-special-case fallback when the
  basin attractor wins over the perturbation).
- **No (1,1), (1,2), (3,2), (2,3) members** at any of (3.08, 3.12, 3.14,
  3.17) crossed with x0 ∈ {-0.90, -0.82, -0.75, -0.68}, xdot0 ∈
  {-0.05, +0.05}, half_crossings ∈ {4, 6}.
- The 7 symmetric (2,0) / (3,0) rows have `independent_closure_dop853 <
  1e-6` so they are real periodic orbits at the requested Jacobi, just
  **not cyclers** by the winding-topology classifier (the orbit fails
  to wind the Moon at all, classified k2=0). These are interesting in
  their own right but outside the cycler scan's scope.
- **1 inadvertent C21 reproduction**: at (target=(1,1), seed (x0,xdot0)
  = (-0.75, -0.05), C=3.08, hc=4) the corrector landed in the C21
  asymmetric family at T=83.96 d, classified (2,1), `asymmetry=0.165`.
  Caught by the known-family period filter (Braik-Ross C21 = 84.53 d at
  C=3.1294, +-2% band). NOT NOVEL.

### What this means

This is the **clean-negative outcome** the orbit-closure-discipline
memory prescribes ("a clean negative is success"). On this slab the
asymmetric corrector's quadratic-basin radius does not catch novel
(k1, k2) cycler members against the topology guard. The C21
reproduction confirms the corrector itself works (it found a known
asymmetric family from a rough seed); the absence of (1,1) / (3,2) /
etc matches at our grid resolution is **not** an emptiness claim
about the full parameter space — it is conditional on:

  * the down-scoped grid (1/309 of the brief),
  * the corrector's narrow quadratic basin (#249 reproduce-test
    documents condition number ~1e7 for C11a),
  * the topology classifier's `n_samples=4000` sampling,
  * the half_crossings (4, 6) subset (not 2, 8).

### Follow-up batched plan

The full brief grid is ~150 hours at the current per-cell cost. Three
tractable expansions:

1. **Wider half_crossings**: add hc=8 to catch high-winding members
   (e.g. (3, 2) C32 which closes at half_crossings=6 symmetric =
   12 asymmetric general-corrector). Cost: ~+50% wall.
2. **Pre-seeded by published members**: for each known cycler (C11a,
   C11b, C21, C32), perturb the IC by xdot0 ∈ {0.001, 0.01, 0.05} and
   walk the family. Bypasses the wide-grid basin radius. Cost: ~30 min.
3. **Continuation from #287 3D seeds** (sibling task #291): if the 3D
   corrector lands a broken-plane (1,1) family at C in [3.05, 3.15],
   project to z=0 and seed this scan. Coordinate with #291 once their
   Phase-1 lands.

Per the negative-results registry (#172): this slab becomes an entry
in the asymmetric-corrector capability-versioned empty-region record,
not a permanent emptiness claim. Re-running becomes mandatory once a
new capability subsumes this method (e.g. multi-shooting genome, or
the 3D continuation seeding from #291).

## Concurrent-agent boundary

This task strictly touched only:

- NEW: `scripts/scan_284_asymmetric.py`
- NEW: `data/scan_284.jsonl`
- NEW: `docs/notes/2026-06-16-284-asymmetric-novel-scan.md`

No file under `src/cyclerfinder/search/cr3bp_general_periodic.py` was
modified (READ-ONLY contract). No catalogue writeback. No edits to the
5-tier prioritizer, the literature_check corpus, or the false-positive
flagger. `#285` (Saturn/Uranus repeated-moon) and `#291`
(`cr3bp_general_periodic_3d.py`) paths untouched.

## Commit framing

> `search: asymmetric-corrector novel-cycler scan EM (#284) — N closures, K literature-fresh`

No `Co-Authored-By` line, no AI attribution, no novelty claims in the
commit message. N/K filled from the final stats summary.
