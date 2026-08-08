# Task #807: PC (3,3) sweep topology failure — root cause and fix

**Date:** 2026-08-09
**Task:** #807 (trigger: self-hosted CI runner now runs on this Mac, so the
`#584` "Linux CI is green → local-Mac-only artifact" closure no longer holds)
**Failing test:** `tests/search/test_504_pluto_charon_kk_sweep.py::test_504_sweep_33`
**Fix commit:** see `git log` for this note's commit (pathspec-staged with
`src/cyclerfinder/search/pluto_charon_kk_sweep.py` + test + OUTSTANDING.md)

---

## Verdict

**Not a knife-edge integer-classification wobble. Not a real (3,3) orbit.**
The (3,3) mu-continuation loses the family branch at the very FIRST mu step
and (on this Mac) the corrector then captures a stable **retrograde (7,0)
near-primary orbit** — a completely different family — which `sweep_33`
reported as `stable_found=True, topology_ok=False` because it lacked the
wrong-topology → clean-negative gate that `sweep_31` (inline, since 2026-07-01)
and `real_binary_kk_sweep._finalize_candidate` (#660) have always applied.
The fix adds that gate (`_topology_gated_result`) to the four non-control
sweeps. The (3,3) verdict is now the same **clean negative on every
platform**, with the recovered topology recorded in `note`.

## Measured evidence (direct, #774-style — not assumed)

Orbit the sweep converged to (this Mac, Accelerate BLAS):

- `x0=+0.12271887`, `ydot0=-2.19931950`, `T=4.6664` TU (~4.75 d),
  `C=3.16131`, `nu=3.3e-11`, `crosscheck_dj=2.6e-12`
- Recovered winding: `w1=-7.000000`, `w2=+0.000000` → (k1,k2)=(7,0),
  retrograde, `reaches_secondary=False`, `x_max=0.1227 < L1`
- The Table-I (3,3) anchor is `x0=-0.32248`, `T=19.504` — nothing about the
  found orbit is "approximately (3,3)". The integer winding check was doing
  exactly its job.

Mu-continuation path trace (fixed C=3.18338, 40 steps, `hc=None`):

| step | mu | topology |
|---|---|---|
| 0 | 0.012151 | (3,3) prograde — anchor OK |
| 1 | 0.014566 | (3,4) prograde — **branch already lost** |
| 2 | 0.016981 | (6,2) prograde, T jumps 19.8 → 29.4 |
| 3 | 0.019397 | **(7,0) retrograde**, x0 jumps −0.309 → +0.219, T → 4.63 |
| 4–40 | → 0.108765 | rides the (7,0) retrograde family to PC mu |

Sensitivity of the END STATE (the actual knife-edge):

| perturbation | outcome |
|---|---|
| baseline (n_steps=40) | converges, retrograde (7,0) |
| anchor x0 × (1 ± 1e-9) | **diverges** (returns None) |
| anchor x0 + 1e-7 | diverges |
| n_steps = 39 or 41 | diverges |
| n_steps = 80 | converges to **prograde (5,0)** — a third family |
| n_steps = 200 | diverges |

A 1e-9 relative seed perturbation flips converged ↔ diverged, and different
step counts land on different families. The platform-dependence is in the
**corrector's diverge-vs-capture outcome after branch loss**, not in the
winding-number computation (which returns exact integers, `w1=-7.000000`).

## Reconciliation of the two prior findings

- **2026-07-01 verdict** (`docs/notes/2026-07-01-504-pluto-charon-kk-sweep-verdict.md`):
  "(3,3): mu-continuation failed (mu-step diverges before PC mu)" — that is
  the DIVERGED side of the measured knife-edge, on the old Linux x86_64
  environment. Correct then, still the correct physics: the (3,3) family
  does not survive continuation to PC mu.
- **#584 (2026-07-14)**: found the same tree on this Mac (M3/Accelerate)
  converging to a stable, crosscheck-clean orbit that fails the topology
  check, and closed it as "local-Mac-only BLAS rounding tipping an integer
  winding check at a basin boundary", backstopped by Linux CI staying green.
  Half right: BLAS/platform rounding DOES control the outcome, but the flip
  is diverge-vs-off-branch-capture in the corrector — the classification
  itself is exact and unambiguous ((7,0) retrograde, not near-(3,3)).
  `git log` confirms neither `pluto_charon_kk_sweep.py` nor the test changed
  between 2026-07-01 (`2b2f74a6`) and #584's window — the behaviour change
  was purely the platform switch.
- Both prior observations are the SAME underlying event (branch loss at the
  first mu step); only the post-loss wreckage differs by platform. That is
  why the correct fix maps both sides to the same clean-negative verdict
  rather than tightening tolerances to chase one side of a 1e-9 knife-edge.

## Why this fix is not a fudge

- The integer topology check is untouched and still exact. No approximate
  matching, no hysteresis band — none is needed, because the recovered
  winding is not near a boundary.
- "Stable orbit of a DIFFERENT family found while sweeping family (k1,k2)"
  = "the target family was not found" = clean negative for the target.
  This is the module's own pre-existing semantics (`sweep_31` inline since
  2026-07-01; the 2026-07-01 verdict's (3,1) row is worded exactly this way)
  and `real_binary_kk_sweep._finalize_candidate`'s (#660). `sweep_33`,
  `sweep_11`, and the two grid sweeps simply never got the gate.
- `sweep_32_positive_control` is deliberately NOT gated: a topology mismatch
  there must fail the hard control loudly, never fold into a clean negative.

## Changes

- `src/cyclerfinder/search/pluto_charon_kk_sweep.py`: new
  `_topology_gated_result()` applied to `sweep_11`, `sweep_21`, `sweep_22`,
  `sweep_33` (the rare gated path recomputes the winding once to record the
  recovered topology in `note`); `sweep_33` docstring documents the branch
  loss + knife-edge measurements.
- `tests/search/test_504_pluto_charon_kk_sweep.py`: `test_504_sweep_33`
  docstring records the expected all-platform clean negative and mechanism.
- No corrector, classifier, or tolerance change anywhere.

## Bug-fix-invalidates-past-searches check

The fix changes REPORTING only (a wrong-family stable orbit is now labelled
a clean negative instead of `stable_found=True, topology_ok=False`); no
solver or classifier changed, so no past negative can become a positive.
Re-verified after the change on this machine:

- **(3,2) positive control** (catalogue row `ross-rt-pc-cycler-32-2026`,
  load-bearing): PASSES the hard assertion path — reproduces cleanly.
- (1,1), (2,1), (2,2), (3,1), (3,3): all clean negatives, unchanged verdicts
  (full `tests/search/test_504_pluto_charon_kk_sweep.py`: 6/6 pass, ~9 min).
- All other importers of the module: 40/40 pass
  (`test_660`, `test_665`, `test_656`, `test_deflated_newton`, `test_627`,
  `test_629`, `test_549`, `tests/scripts/test_633`). One transient
  `test_656` grid-seed failure during the first 8-way-parallel run was
  re-run in isolation and passes in 6 s — CPU contention with the
  self-hosted CI runner (load average ~47–50 at the time) tripping the 5 s
  SIGALRM per-call bound, per [[feedback_serialize_verification_runs]].
- The 2026-07-01 verdict document is left as the historical record; its
  (3,3) mechanism line remains accurate for the environment it was run on.
- `uv run ruff check .`, `ruff format --check .`, full `uv run mypy src tests`:
  all clean.

## Follow-up registered

- **#808**: `real_binary_kk_sweep._finalize_grid_candidate` deliberately
  (documented, #660 scope decision) does NOT gate the post-C-sweep topology
  on grid paths — the exact latent failure class that bit #807, and
  `test_549` asserts `topology_ok` on such results, so a future
  platform/rounding shift could reproduce this failure mode there. Decide
  whether to extend the gate (with the same clean-negative semantics) or
  document why grid paths are exempt.
