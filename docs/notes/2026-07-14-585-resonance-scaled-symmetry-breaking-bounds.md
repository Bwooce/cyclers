# #585 (#582 follow-up): resonance-scaled symmetry-breaking bounds ladder

**Date:** 2026-07-14
**Scope of this note:** the bounds-widening build + the mandatory
per-rung positive-control gate, per #585's own Fable-reviewed design
("MANDATORY positive control before any novelty claim"). No novelty sweep
was run here — that is a separate, coordinator-owned dispatch, exactly as
#582's own was.

## Origin

#582's 5-MMR asymmetric novelty sweep (commit `b0225f4`) found 0/104
asymmetric members across all 5 tabulated interior MMRs — every converged
cluster classified `symmetric` (see `empty_regions.jsonl`'s
`er3bp-isolated-3d-asymmetric-mu0.001-5mmr-582-2026-07-14` entry). A Fable
design review diagnosed this as very likely a search-box artifact:
`mmr_bounds()`'s 3 flat absolute half-widths (`z0_abs = xdot0_abs =
zdot0_abs = 0.05` for every MMR, regardless of resonance) were anisotropic
~5-8x in exactly the state components whose departure from zero breaks
mirror symmetry — `xdot0_abs=0.05` at the 3:2 MMR is only an
eccentricity-proxy budget of e≈0.03-0.044, while `ydot0_frac=0.35` already
grants the SYMMETRIC direction e≈0.2 reach.

## What was built

* `src/cyclerfinder/search/isolated_3d_asymmetric_fitness.py::mmr_bounds()`
  — replaced the 3 flat absolutes with one resonance-scaled symmetry-breaking
  fraction `s`:

  ```
  v_circ = sqrt((1 - mu) / a1)          # same quantity ydot0_guess is built from
  xdot0_abs = zdot0_abs = s * v_circ
  z0_abs = max(0.05, s * a1)
  ```

  `s` is a new keyword-only parameter with **no default** — every caller
  must pick a value explicitly (a flat 0.05 box has no single `s` that
  reproduces it identically across all 5 MMRs, since that flat-ness is
  exactly the anisotropy #585 corrects, so there is no safe "old behaviour"
  default to fall back on). `x0_frac=0.15`, `ydot0_frac=0.35`, `t_frac=0.5`
  are completely unchanged. The `max(0.05, s*a1)` floor guarantees no rung
  ever shrinks `z0_abs` below #582's own already-stamped box.

* `src/cyclerfinder/search/isolated_3d_asymmetric_fitness.py::mmr_a1_from_t0()`
  — the exact inverse of `mmr_t0()` (`a1 = (1 + 2*pi/T0)**(-2/3)`), added to
  back out a converged candidate's implied semi-major axis from its period
  for the drift-detection check below. Round-trips to machine precision on
  all 5 tabulated `a1` (verified in
  `test_mmr_a1_from_t0_inverts_mmr_t0`).

* `scripts/run_582_asymmetric_3d_niching_search.py`:
  - New required `--symmetry-breaking-s` CLI flag (float, no default),
    threaded through `run_ga` / `positive_control` / `analyze_ga_population`
    into `mmr_bounds()`. All 3 modes require it explicitly.
  - Output filenames are now tagged with the `s` value used (e.g.
    `positive_control_3_2_s0p15_summary.json`, `ga_3_2_s0p3_final.npz`,
    `3_2_s0p15_analysis_summary.json`) so runs at different rungs never
    collide, and #582's own untagged files stay untouched.
  - New `nearest_mmr_by_implied_a1()` + drift-detection wired into
    `--mode analyze`: for each converged cluster representative, the
    converged period is inverted to an implied `a1` via `mmr_a1_from_t0`
    and compared against every tabulated MMR's own `a1`
    (`MMR_SEMI_MAJOR_AXES`). If the nearest MMR by that implied `a1` is not
    this run's own target, the representative is flagged
    `drifted_to_neighboring_mmr=true` in its entry, and the summary JSON
    reports a per-run `n_drifted_to_neighboring_mmr` count. Drifted
    candidates are reported, not discarded (a drifted-but-converged orbit
    may be a genuine member of ITS OWN resonance, just not this run's
    target) — this makes over-widening detectable rather than silently
    folded into "this MMR's" cluster count.

* Tests: `tests/search/test_isolated_3d_asymmetric_fitness.py` gained
  6 new tests (`test_mmr_bounds_requires_s_explicitly`,
  `test_mmr_bounds_s_scaling_formula` (both rungs x all 5 MMRs),
  `test_mmr_bounds_z0_never_shrinks_below_582_floor`,
  `test_mmr_bounds_rung2_wider_than_rung1`,
  `test_mmr_a1_from_t0_inverts_mmr_t0`,
  `test_mmr_a1_from_t0_raises_on_nonpositive_or_nonfinite`) and all
  existing `mmr_bounds(...)` call sites were updated to pass `s` explicitly.
  New `tests/scripts/test_run_582_drift_detection.py` (3 tests) pins the
  drift-detection helper's own logic (recovers own MMR at exact T0, flags
  drift to a neighbor at a period near a different MMR's T0, `_s_tag`
  filesystem-safety).

## Numeric bounds values (MMR 3:2, `a1=0.763143`, `mu=0.001`)

| Quantity | #582 original (flat) | s=0.15 | s=0.30 |
|---|---|---|---|
| `z0_abs` | 0.05000 | 0.11447 | 0.22894 |
| `xdot0_abs = zdot0_abs` | 0.05000 | 0.17162 | 0.34324 |

`v_circ = sqrt(0.999/0.763143) = 1.14413`; `s=0.15` gives `xdot0_abs =
0.15*1.14413 = 0.17162`; `s=0.30` doubles that to `0.34324`. `s*a1 =
0.15*0.763143 = 0.11447` already exceeds the 0.05 floor at 3:2, so the
`max()` floor is inactive here (it only binds at wider MMRs with smaller
`a1`, e.g. 5:1's `a1=0.3419`: `s=0.15` gives `s*a1=0.05128`, barely above
the floor).

All 5 MMRs, both rungs (`z0_abs`, `xdot0_abs=zdot0_abs`):

| MMR | a1 | s=0.15 z0_abs | s=0.15 xdot0_abs | s=0.30 z0_abs | s=0.30 xdot0_abs |
|---|---|---|---|---|---|
| 3:2 | 0.7631 | 0.11447 | 0.17162 | 0.22894 | 0.34324 |
| 5:2 | 0.5428 | 0.08142 | 0.20350 | 0.16284 | 0.40699 |
| 3:1 | 0.4807 | 0.07211 | 0.21624 | 0.14421 | 0.43248 |
| 4:1 | 0.3968 | 0.05952 | 0.23801 | 0.11904 | 0.47601 |
| 5:1 | 0.3419 | 0.05128 | 0.25640 | 0.10257 | 0.51281 |

Every rung/MMR combination widens `xdot0_abs`/`zdot0_abs` well beyond
#582's flat 0.05 (3.4x-10.3x), consistent with Fable's ~5-8x anisotropy
finding; `z0_abs` widens too, though less dramatically at the largest-`a1`
MMR (3:2).

## Mandatory positive-control gate result (MMR 3:2, the hardest tabulated case)

Both runs used the SAME small positive-control GA budget as #582
(`_small_ga_config`: population=40, generations=60, seed=`582000+p*100+q` =
582302 for 3:2) — not the eventual paper-scale sweep budget.

### s=0.15: **PASS**

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode positive-control --symmetry-breaking-s 0.15 --workers 8
```

Wall time ~29s (8 workers). Result (`positive_control_3_2_s0p15_summary.json`):

| Quantity | Known #440 member | GA + corrector result | Tolerance | Error | Margin |
|---|---|---|---|---|---|
| x0 | 0.7310974 | 0.7333928 | 3% rel | 0.314% | 9.6x |
| ydot0 | 0.4569989 | 0.4402828 | 5% rel | 3.658% | 1.4x |
| T | 11.926512 | 11.624758 | 5% rel | 2.530% | 2.0x |
| C (Jacobi) | 3.0622660 | 3.0656329 | 0.02 abs | 0.00337 | 5.9x |

`orbit.converged=True`, `corrector_residual=3.06e-12`, independent Radau
closure `2.62e-12`. `classify_symmetry` correctly reports
`is_symmetric=True` (crossing residual 7.8e-15, machine zero) and
`degenerate_planar=True` — the widened box still lands the GA in the SAME
known symmetric basin as #582's original flat box, just via a moderately
different genome path (compare this run's converged `state0` against
#582's own note's `x0=0.7356`/`ydot0=0.4390` — both are valid points on the
1-parameter family the under-determined corrector selects from). **The
widened s=0.15 box does not break basin containment.**

### s=0.30: **FAIL** (basin-competition, not simple box-blowout)

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode positive-control --symmetry-breaking-s 0.30 --workers 8
```

Wall time ~64s (8 workers). Result (`positive_control_3_2_s0p3_summary.json`):
GA best genome converged to `T=6.713` (vs known `T=11.927`), `x0=0.7172`,
`ydot0=0.3430` — the corrector did NOT converge on this seed
(`corrector_residual=0.139`, `independent_closure=0.139`,
`degenerate_planar=False`). `matched=False`: `ydot0_err=23.5%`,
`T_err=44.1%`, `jacobi_err=2.98e-2` — all well outside tolerance.

**This is a diagnosed basin-competition effect, not a simple "box too
wide" failure.** An independent 3-seed probe (not the official run, an
ad-hoc diagnostic at the same `s=0.30` box, same small 40x60 budget) found:

| seed offset | best fitness | T | x0 | ydot0 | reads as |
|---|---|---|---|---|---|
| +0 (= official run's seed) | 0.9610 | 6.713 | 0.7172 | 0.3430 | wrong basin |
| +7919 | 0.9626 | 6.846 | 0.7066 | 0.3745 | wrong basin |
| +15838 | 0.9963 | 12.510 | 0.7168 | 0.3066 | near-known basin |

1 of 3 seeds lands near the known basin (fitness 0.9963, `T` within 5% of
known); 2 of 3 land on a DIFFERENT near-periodic basin around `T≈6.7-6.8`
with comparable-or-higher raw fitness under this small budget (0.961-0.963
vs the correct basin's 0.987 at s=0.15). This means at `s=0.30`, the widened
box makes a genuinely competing basin reachable with fitness the small
40-member population's simple argmax pick cannot reliably discriminate
against the known member's basin — exactly the failure mode Fable's own
design review flagged ("over-widening wastes GA capacity but is post-hoc
detectable, not silently misleading" — this is that detection working as
intended, just showing up at the positive-control stage instead of only in
the eventual sweep's drift-check).

**Verdict for s=0.30: the mandatory positive control did NOT pass under
this reporting scheme (`positive_control` mode's simple best-of-population
argmax).** Per #585's own gate, this is reported honestly rather than
forced to a pass. Two live possibilities, not resolved here (out of scope
per #585's own "don't run the full sweep" instruction):

1. **s=0.30 is genuinely too aggressive for this fitness landscape** at
   3:2 (the largest, most eccentric, `a1`-largest MMR) — a smaller `s`
   between 0.15 and 0.30 might be the real safe ceiling.
2. **s=0.30 is usable but the SMALL positive-control budget (40 pop x
   60 gen) is the actual limiting factor**, not the box itself: the
   niching GA's own deterministic-crowding design is meant to hold
   MULTIPLE co-existing basins in one population simultaneously (this is
   the whole reason `--mode analyze` clusters the population instead of
   taking a single argmax) — the paper-scale sweep budget (200 pop, 400
   gen, per `_small_ga_config`'s own docstring contrast) may hold BOTH the
   known-basin niche and the competing basin, and `--mode analyze`'s
   clustering (not argmax) would then still recover the known member as
   one of several cluster representatives, with the competing basin
   correctly reported as its own separate (possibly drift-flagged, since
   its very different `T` may well imply a different resonance) cluster.

**Recommendation for whoever dispatches the full sweep:** treat `s=0.15`
as the confirmed-safe rung (positive control passed cleanly, single seed,
small budget). Do NOT certify `s=0.30` from this note's evidence alone —
either (a) re-run the `s=0.30` positive control at the actual paper-scale
budget (200 pop / 400 gen) before trusting it, since the small-budget
argmax check is demonstrably basin-competition-sensitive at this rung, or
(b) run `s=0.30` only through the full `--mode ga` + `--mode analyze`
pipeline and rely on the clustering (which is designed for exactly this
multi-basin situation) rather than any single-seed argmax positive control
to confirm the known member is still reachable.

## Not done here (explicitly out of scope per #585's own instructions)

* The full 5-MMR `--mode ga` novelty sweep at either rung (a multi-hour
  job, ~10-12 runs at 9-23 min each — coordinator-owned, per #585's own
  cost estimate).
* Any `empty_regions.jsonl` entry for a rung's actual sweep result (there
  is none yet to stamp — this dispatch is the bounds-widening build + the
  fast positive-control gate only).
* `t_frac=0.5`'s own flagged gap (excludes period-doubled asymmetric
  branches, `T≈2*T0`) — noted by Fable's own review as out-of-scope for
  this dispatch, carried forward as-is.

## Reproduction

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode positive-control --symmetry-breaking-s 0.15 --workers 8
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --mode positive-control --symmetry-breaking-s 0.30 --workers 8
uv run pytest tests/search/test_isolated_3d_asymmetric_fitness.py tests/scripts/test_run_582_drift_detection.py -v
```
