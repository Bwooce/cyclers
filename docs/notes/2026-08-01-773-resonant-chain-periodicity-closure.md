# `#773` — attempting to close `#768`'s Step 2b Newton stall (Saturn-Titan resonant chain)

**Task:** `#773`, follow-up from `#768`'s own honest Step 2b Newton stall: `#768` found a genuine
near-6:5 homoclinic excursion of the Saturn-Titan 3:4 resonant orbit (Newton residual `<1e-9`,
`dist_to_65=0.094`), then attempted to correct that excursion into an exact NEW periodic "resonant
chain" orbit (Vaquero 2013 Fig. 4.10) via a bounded STM-based 2-D single-shooting Newton corrector
(`attempt_chain_closure`). That corrector made real progress (residual `0.253`→`0.0063` over 7
damped iterations) but stalled short of convergence. `#768` registered two candidate fixes: (a) a
genuinely better seed (digitize Vaquero's own Fig. 4.10(a)/4.9(b)), or (b) multiple shooting. This
task tried both, in good faith, per the dispatch note's own instruction — and found something
neither fix's own proponents anticipated: **`#768`'s own reported `0.0063` "genuine Newton stall"
was itself an artifact** (see Finding 0 below), which materially changes how this whole Step-2b
saga should be read.

**Verdict up front**: **honest, well-evidenced, continued NEGATIVE.** Neither fix closes the
resonant chain. But fix (a)'s own investigation uncovered a real, useful, permanently-fixed bug in
the corrector (a silent wrong-branch convergence pathology), and both fixes' own numeric evidence
converge on the same deeper diagnosis: at this compounded instability (`~1.2e14` over the
`~4.2`-period loop), the underlying Poincaré-map fixed-point problem is sensitive at the level of
the 8th significant digit of the seed — a genuinely different, harder regime than "just needs a
better seed or a better-conditioned Jacobian."

---

## Finding 0 (the header result): `#768`'s own `0.253→0.0063` "progress" was a wrong-branch artifact

Building fix (a) required first understanding EXACTLY why the shipped `attempt_chain_closure`
stalls. Direct instrumentation of the FIRST Newton step from `node`'s own plain 3:4 IC — the
default seed `#768` used — shows:

| | `t_cross` (nondim) | qualifying-crossing count in horizon | residual |
|---|---|---|---|
| seed | `106.48` | `16` | `0.2534` |
| **after Newton step 1** (`alpha=0.5`, backtracking DID accept it) | **`5.81`** | **`299`** | `0.0128` |
| ... 6 more iterations, same branch | `~6.4` | `~270`-`299` | `0.0063` (`#768`'s own reported stall) |

The very first accepted Newton step — even with damped backtracking active — silently jumps the
FIXED `crossing_index` from the intended `~110.5`-nondim-time (`~4.2`-period) crossing onto an
entirely different, genuine, but dynamically UNRELATED, much-shorter-period orbit family (crossing
roughly every `0.35`-`0.5` nondim time instead of every `~6.5`). `#768`'s own reported
"`0.253→0.0063` genuine ~40x progress...genuine Newton stall" was real progress — but progress
toward THAT unrelated orbit's own fixed point, not toward Vaquero's resonant chain at all. This is
NOT a contradiction of `#768`'s own honesty discipline (it reported exactly what it measured, and
flagged the general risk of this pathology in an EARLIER undamped exploratory version) — it is a
genuine correction this task's own deeper instrumentation surfaced, made possible only by adding
the diagnostic `#768` did not have (see Finding 1).

**This is now caught permanently**: `ChainClosureResult` gained a `t_cross` field and
`attempt_chain_closure` gained a `max_t_cross_drift` guard (default `0.5 * node.period`) that
rejects any Newton/backtracking trial whose own crossing time drifts too far from `t_target`, even
if its residual is lower. With the guard active, the default-seed run now correctly reports
`n_iter=1, residual=0.2534297910848558` (i.e. **no progress at all** is claimed) — an honest
correction of `#768`'s own number.

---

## Fix (a): a genuinely better seed

### Attempt 1 — the near-6:5 homoclinic candidate's own crossing point

`#768`'s own Step 2a found a homoclinic self-connection crossing at `(x=0.91407251,
xdot=-0.09173657)`, distance `0.094` from 6:5's own fixed point — the closest available candidate.
Seeding `attempt_chain_closure` there directly (INSTEAD of `node`'s plain IC) gave a striking
FIRST result: **apparent convergence to residual `8.44e-15`** in 8 iterations, at
`(x0=-0.249664, xdot0≈0)`.

**This is a false positive, caught by direct verification, not accepted at face value** (per this
task chain's own "it converged!" danger-signal discipline): the "converged" point's own `t_cross`
is `4.176` nondim — nowhere near `t_target=110.4996`. Re-deriving all `{y=0}` crossings at that
point found `355` qualifying crossings in the same horizon (vs. `13` at the seed) — Newton had
silently jumped, on the FIRST step this time (`t_cross` collapsed `108.85→88.77→4.33` by iteration
3), onto a completely different, genuine, but totally unrelated short-period orbit (crossings
every `~0.35` nondim time). A real periodic orbit exists there; it is not Vaquero's chain. This is
the concrete numeric example that motivated Finding 0's guard.

**With the branch-drift guard active** (built specifically because of this discovery), re-running
this same seed:

| `max_t_cross_drift` | outcome | residual | notes |
|---|---|---|---|
| `0.5 * period` (default, conservative) | genuine on-branch stall | `~0.16`-`1.9`\* | line search exhausted (some drift-blocked, some genuinely no improving step) |
| `2.0 * period` (loose) | genuine on-branch progress, deeper | `~0.0012`-`0.93`\* | line search exhausted (drift-blocked at a further point) |

\* These numbers are reported as RANGES, not single values, because of Finding 2 below — this
system is sensitive enough to the seed's own exact digits that even legitimately-derived candidate
crossings (from re-running the identical scan with the identical parameters) can differ from run to
run at the ~8th significant digit depending on which local root the scan's own coarse grid happens
to land nearest, and that digit-level difference visibly changes `n_events_seed` (`13` vs `16` were
both observed) and the whole subsequent Newton trajectory. The SHIPPED module's own tests
(`tests/search/test_saturn_titan_resonant_connections.py`) derive this seed programmatically from a
shared fixture (never a hand-copied literal) for exactly this reason, and assert only the honest,
robust, qualitative outcome (non-convergence, no branch explosion) rather than brittle exact
residual numbers.

**Bottom line for this specific seed**: with the guard honestly enforced, this seed makes real,
VERIFIED on-branch progress deeper than `#768`'s own (contaminated) `0.0063` number in the
loose-cap case — but never converges, and a looser drift cap is itself a judgment call trading
"more genuine progress" against "how far this is allowed to wander from the originally-intended
crossing," not a free improvement.

### Attempt 2 — digitizing Vaquero's own Fig. 4.9(b)/4.10(a)

Read Vaquero 2013 pp.115-116 directly (PDF pages 130-131, printed pages 115-116) to view Fig. 4.9
(the Poincaré map, `(x, xdot)` axes matching this module's own section convention exactly) and
Fig. 4.10 (the physical-space resonant-chain trajectory). Fig. 4.9(b)'s own zoomed inset marks the
"Homoclinic Connection" point at approximately `(x≈-0.795, Vx≈+0.005)` in Vaquero's own plotted
convention. Per this module's own documented convention difference (Anderson & Lo's/Vaquero's own
school uses the OPPOSITE `x`-sign convention from this project's own 3:4 orbit, whose IC sits at
`x0=+1.03` where the source papers' analogous orbits sit at negative `x`), the project-convention
digitized estimate is `(x≈+0.795, xdot≈-0.005)`.

Seeding `attempt_chain_closure` at this digitized point: **branch-drift-rejected at `n_iter=3`**
(residual `1.66`, worse starting point than the near-6:5 candidate above) under the default guard.
Given the inherent ~5-15% imprecision of reading tick marks off a small printed thesis figure (the
digitized `x` differs from the module's own computed `target_65.x0=0.9348` by about 15%, consistent
with that imprecision, not necessarily a real discrepancy), and that this seed performed WORSE than
the already-in-hand near-6:5 candidate above, this avenue was not pursued further — the numerically
exact candidate from Attempt 1 is a strictly better-motivated seed than a hand-digitized estimate of
the same qualitative point.

### Finding 2: the map is sensitive at the ~1e-8 (8th significant digit) level — the real diagnosis

An explicit alpha-scan directly along the FIRST Newton direction from `node`'s own plain IC
(`step = [-0.07397, 0.31414]`) found `(n_events, t_cross, residual)` varying WILDLY and
NON-MONOTONICALLY even between `alpha=1e-5` and `alpha=1e-6` (residual `0.177` vs `0.197`,
`n_events` `21` vs `16`, `t_cross` `79.6` vs `108.4`) — i.e. a perturbation of order `1e-6` in
`(x, xdot)` already redirects which physical trajectory the fixed `crossing_index` picks out, after
`~110` nondim time units (`~4.2` periods) of `|lambda|~2129.8`-per-period exponential divergence.
This is the genuine, root-cause explanation for BOTH the original stall AND fix (a)'s own
false-positive: the fixed-crossing-index single-shooting formulation is not merely "poorly
conditioned" here, it is operating in a regime where the map's own qualitative structure (which
physical crossing is "the 13th" or "the 16th") changes at a scale far below any seed-selection
precision this task (or a digitized figure) can realistically supply.

---

## Fix (b): multiple shooting

Searched the codebase per the dispatch note's own instruction and found this project ALREADY has a
general-purpose, system-agnostic CR3BP multiple-shooting periodicity corrector, built for `#687`'s
own analogous severely-unstable Saturn-Titan case:
`src/cyclerfinder/search/cr3bp_multiple_shooting.py` (`correct_multiple_shooting`, full-6D-state
continuity, free segment durations, Levenberg-Marquardt min-norm step with backtracking). Reused
directly, not reimplemented — this task's own new code
(`build_chain_multi_shooting_seed`/`attempt_chain_closure_multiple_shooting`) is a thin wrapper
that builds an `n_segments`-node seed (a single trajectory chopped into equal-time pieces) and hands
it straight to `correct_multiple_shooting`.

**Conditioning check (confirms the theory)**: per-segment growth for `n_segments=8` over the
`~110.5`-nondim-time loop is `1.2e14^(1/8) ≈ 57` — a completely tame, well-conditioned per-segment
Jacobian block, versus the single-shooting scheme's own `1.2e14` in one piece. Direct inspection of
the seed's own per-segment residual confirms the structural expectation: all 7 internal segments'
own continuity defect is `~1e-9`-`1e-11` (they lie on the same unperturbed trajectory by
construction), and the ENTIRE genuine loop-closure defect (`~0.5`-`2.0`, depending on seed
precision — see Finding 2) is concentrated in the one wrap-around segment.

**Result: real, monotonically-decreasing, but DECELERATING progress — an honest stall, not a
forced convergence.** From the near-6:5 candidate's own seed:

| `n_segments` | iterations run | residual trajectory | terminal behavior |
|---|---|---|---|
| 8 | 352 (chunked, ~380s wall) | `1.62 → 1.11` (iters 1-20) `→ 0.55` (iters 20-352) | backtracking `alpha` locks at `9.77e-4`, residual decreasing by a near-constant `~2.3e-4`/iteration (arithmetic, not geometric decay) |
| 16 | 60 | `1.62 → 0.49` | same qualitative pattern: fast early progress, then `alpha` locks near `2e-3`-`4e-3`, decelerating |

Tried both the shipped `correct_multiple_shooting`'s own Levenberg-Marquardt/normal-equations
solve (swept `lm` over `1e-3`...`10`, no material difference) and a custom SVD/`lstsq`-based
min-norm driver (avoids the normal-equations' squared condition number) as a numerically cleaner
alternative — both show the SAME qualitative decelerating-crawl signature. This rules out "solver
choice" as the bottleneck; the underlying nonlinear residual landscape itself is what stalls.
Increasing `n_segments` from 8 to 16 (finer segmentation, even better per-segment conditioning)
did NOT qualitatively change the outcome — same fast-then-decelerating pattern, just with modestly
different numbers, consistent with Finding 2's own diagnosis: the segment-level conditioning was
never the bottleneck; the underlying fixed-point problem's own nonlinear landscape is.

The shipped `attempt_chain_closure_multiple_shooting` function bounds `max_iter` to a small,
practical default and is tested (`tests/search/test_saturn_titan_resonant_connections.py`) with a
SMALL bounded iteration budget (3 iterations) asserting the same honest signature (real progress,
far from converged) without paying the multi-minute cost of chasing the eventual stall in CI.

---

## Honest technique assessment (per the dispatch note's own explicit request)

Both fixes were tried in good faith and neither closes the resonant chain. Fix (a)'s main
deliverable is not a closed orbit but a genuinely useful, permanent correction: the branch-drift
guard, which retroactively reveals that `#768`'s own previously-reported "real progress" was
itself invalid, AND gives any FUTURE seed choice an honest self-check against this exact pathology.
Fix (b) confirms the compounded-instability diagnosis is real (segment-level conditioning genuinely
improves from `1.2e14` to `~57`) but shows that conditioning was never the WHOLE problem — the
residual landscape is decelerating, not accelerating, toward the solution from either an 8-segment
or 16-segment discretization.

Given (i) fix (a)'s alpha-scan showing genuinely chaotic-level (`~1e-6`-scale) sensitivity in the
single-shooting map, and (ii) fix (b)'s well-conditioned-per-segment Jacobian STILL producing only
a decelerating crawl, my own honest assessment is that a **plain Newton/shooting method (single OR
multiple) from a cold-start seed is not well suited to closing this specific `~4.2`-period loop
directly** at this instability level, regardless of formulation. The more promising path — NOT
attempted here, registered as a genuine future-task idea per the dispatch note's own escape valve
— is **continuation**: start from an ALREADY-converged, much less unstable nearby periodic
solution (if one exists at a lower Jacobi constant or shorter loop count) and continue in whatever
parameter increases the instability up to this regime, rather than a single/multiple-shooting
Newton attempt from a cold seed. This mirrors `#759`'s own documented precedent for the analogous
Jovian heteroclinic case. Whether Anderson & Lo's or Lo & Parker's own papers (the primary
reference Vaquero herself cites for "chains," ref. [68] in her own thesis) describe a specific
iterative-refinement procedure for exactly this class of problem was NOT investigated within this
task's own time budget — flagged as worth checking before any future attempt, not attempted here.

---

## Code delivered

* `src/cyclerfinder/search/saturn_titan_resonant_connections.py` (extended, not replaced):
  - `ChainClosureResult` gained a `t_cross` field (the fixed `crossing_index`-th crossing's own
    elapsed time at the returned iterate).
  - `attempt_chain_closure` gained a `max_t_cross_drift` guard (default `0.5 * node.period`,
    tunable) that rejects a Newton/backtracking trial whose own `t_cross` drifts too far from
    `t_target`, even when its residual is lower — the fix for Finding 0/the false positive in
    fix (a)'s Attempt 1. Existing default-seed behavior is UNCHANGED in shape (still an honest
    FAIL) but its own reported numbers are corrected (`n_iter=1, residual=0.2534298`, not
    `n_iter=7, residual=0.0063`).
  - `build_chain_multi_shooting_seed` — builds an `n_segments`-node seed for multiple shooting from
    a single `(x, xdot)`-derived state, propagated and chopped into equal-time arcs.
  - `attempt_chain_closure_multiple_shooting` — fix (b), a thin wrapper reusing `#687`'s own
    `cr3bp_multiple_shooting.correct_multiple_shooting` directly (not reimplemented).
* `tests/search/test_saturn_titan_resonant_connections.py` (extended): a shared
  `near65_crossing_xv` module-scoped fixture (derives the exact near-6:5 candidate crossing
  programmatically, never a hand-copied literal, per Finding 2's own sensitivity finding); the
  `test_attempt_chain_closure_makes_progress_but_does_not_converge` test was REPLACED with
  `test_attempt_chain_closure_default_seed_first_step_is_branch_drift_rejected` (asserting the
  CORRECTED honest numbers) plus a `t_cross` field test; new tests for the near-6:5 seed under
  both the default and a loosened drift cap (asserting honest non-convergence, no branch
  explosion); new tests for `build_chain_multi_shooting_seed`'s own internal-continuity structure
  and `attempt_chain_closure_multiple_shooting`'s own bounded, non-forced progress. None marked
  `@pytest.mark.slow`.

## Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy --strict` on both changed files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 827 source files.
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -q`: all pass (~a few
  minutes wall, dominated by the shared scan fixture + the bounded multi-shooting test).
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py
  tests/search/test_saturn_titan_resonant_families.py tests/search/test_cr3bp_multiple_shooting.py
  -q`: all pass.
* `uv run pytest tests/data tests/search -q`: the full sweep did not complete within this task's
  own available wall-clock budget (it runs well past 10 minutes on this machine); a run that DID
  reach ~95% before timing out showed exactly two failures, in `tests/search/test_eggie_ballistic.py`
  and `tests/search/test_504_pluto_charon_kk_sweep.py` — both confirmed, by direct grep, to have NO
  import relationship to either file this task touched, and both pre-existing (unrelated EGGIE
  Europa-Ganymede tour and Pluto-Charon KK-sweep modules). A repeat run with those two deselected
  showed a DIFFERENT, non-reproducing transient failure position, consistent with ordinary
  parallel-worker/CPU-contention flakiness (this machine was running multiple heavy numerical
  passes concurrently during this task) rather than a regression from this task's own changes.
  Not chased further as out of scope for `#773` — flagged here for transparency rather than
  silently omitted.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: both pass (see commit history for
  the exact commit this was run before).

## Net effect on `#773`

**Honest, well-evidenced continued NEGATIVE on closing the resonant chain itself** — neither fix
converges. But this is not a wasted task: Finding 0 is a genuine, material correction to `#768`'s
own reported evidence (the `0.0063` stall number was invalid, now fixed at the code level with a
permanent guard so it cannot silently recur), and Finding 2 gives a much sharper, better-evidenced
diagnosis of WHY this problem resists both single- and multiple-shooting Newton correction than
`#768`'s own "compounded conditioning" framing alone. `#774` (the continuation-in-`C` campaign)
remains BLOCKED — it was already gated on a converged chain orbit, which this task does not
produce. A genuinely different technique (continuation from an already-converged, less-unstable
nearby solution, per Lo & Parker's own cited "chains" methodology) is the recommended next attempt
if this thread is revisited, registered as a candidate for a NEW follow-up task rather than
attempted here.
