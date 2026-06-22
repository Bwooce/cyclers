# #169/#170 continuous-TCM goldens re-derived after the #198 epoch fix

**Date:** 2026-06-23. Surfaced while starting M7 Phase 1 (#423), whose first step is
"reproduce the 3 recorded #169/#170 horizon-TCM goldens." They did **not** reproduce
on current `main`. Root-caused, fixed, re-pinned. No catalogue row changes.

## What broke

The Sun-only patched-conic, continuous-from-one-seed horizon-TCM proxy
(`search.s1l1_corrected.continuous_chain` / `search.appc_corrected.continuous_chain`,
sum of `dv_total_kms`) no longer produced the values recorded on 2026-06-08:

| row | recorded (pre-#198) | current (corrected) | vs 120 m/s V3 budget |
|---|---|---|---|
| S1L1 (mcconaghy-2006-em-k2 parent) | 62 m/s | **40.2** | under (more so) |
| `russell-ch4-8.049gGf2` (#188) | 163.6 m/s | **114.4** | now *under* the proxy |
| `russell-ch4-8.165Gfh-f2` (#192) | 2040.6 m/s | **2020.7** | ≫ 120 |

`tests/nbody/test_appc_batch_nbody.py` was therefore **red** (its `tcm > 120` gate
failed for #188: 114.4 < 120). It went unnoticed because both tests are
`@pytest.mark.slow` (deselected in default CI).

## Root cause (airtight bisect)

Commit **`439d279` — "core: #198 fix 63s UTC/TDB epoch offset (J2000 TDB JD +
TDB-second TimeDelta) + Horizons absolute-epoch anchor test"** (2026-06-12), a
*descendant* of the #170 golden commit `bf8d2a8`.

- Run at `439d279^` (pre-fix): #188 → **163.61 m/s**.
- Run at `439d279` (the fix): #188 → **114.37 m/s**.

The 63 s epoch correction shifts `ephem.state()`. The per-leg, *re-anchored* v∞
reproduction is robust to it (encounter v∞ still matches published to < 5 m/s — those
assertions stayed green). But the *continuous, non-re-anchored* chain accumulates the
shift across 7 cycles, moving the v∞-magnitude maintenance terms → a different total.
This is a **correctness fix** (validated against a Horizons absolute-epoch anchor), so
the pre-#198 goldens were computed under a since-fixed epoch bug. Classic
"bug-fix invalidates past searches."

Confirmed *not* the cause: `s1l1_corrected.py`, `appc_corrected.py`, `core/flyby.py`
are byte-identical to their golden commits; `constants.py` changed only by adding
dwarf planets (Earth/Mars/Sun values identical). The Sun-only IAS15 path equals the
analytic Kepler path here (RestrictedNBody Sun-only is two-body), so the propagator
choice was never the variable.

## Conclusion-level impact: none on the catalogue

- **S1L1**: 40.2 m/s, still well under 120 → essentially_ballistic / V3 verdict
  **unchanged and strengthened**.
- **#188 / #192**: stay V0 / powered. The binding V3-rejection is the **sourced**
  published App-C total Δv (`block.total_dv_kms` = 420 / 1678 m/s, Russell Table 5.5),
  which is epoch- and model-independent. The Sun-only continuous proxy was always an
  **underestimate** (Sun-only, no Mars perturbation) and is a *lower bound*, not a V3
  acceptance measure. The old `tcm > 120` gate happened to hold pre-#198 by luck; the
  epoch fix pushed #188's proxy under 120, exposing that the gate was mis-grounded.

## Fixes applied

- `tests/nbody/test_appc_batch_nbody.py`: renamed
  `..._tcm_over_v3_budget` → `..._published_dv_over_v3_budget`; the V3-rejection now
  asserts `block.total_dv_kms > 120 m/s` (sourced, binding); the Sun-only proxy is
  pinned as a labelled **computed regression anchor** (114.4 / 2020.7 m/s, `abs=1.0`),
  explicitly NOT the gate.
- `tests/nbody/test_s1l1_continuous.py`: prose 62 → 40.2 m/s with the rederivation
  note; assertions were upper-bound only, so already green.

## The M7 design lesson

The Sun-only continuous-chain proxy can make a 420 m/s **powered** cycler look almost
ballistic (114 m/s). A trustworthy per-row maintenance-ΔV (M7's `horizon_tcm_mps`)
**must** use the Mars-perturbed run with B-plane flyby targeting — confirming the
M7 scoping plan's "hard part #1." M7 Phase 1 now reproduces the **corrected** goldens
(40.2 / 114.4 / 2020.7), not the stale pre-#198 ones.

## References

- `docs/notes/2026-06-08-appc-v3-batch-results.md` (original #170 numbers, pre-#198);
- `docs/notes/2026-06-08-s1l1-continuous-v4-results.md`;
- `docs/notes/2026-06-23-m7-scoping-plan.md`;
- commit `439d279` (#198 epoch fix); `feedback_bugfix_invalidates_past_searches`.
