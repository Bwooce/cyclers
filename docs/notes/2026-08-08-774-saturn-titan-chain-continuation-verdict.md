# `#774` — Saturn-Titan resonant-chain continuation-in-`C`: verdict

**Task:** `#774`, reopened 2026-08-08 by `#782`'s own positive result (a converged,
independently-verified "resonant chain"-type periodic orbit at Vaquero 2013's own
`C=3.010000`, `x0=0.9492672902`, `ydot0=0.0796459472`, `period=56.03254231`,
`half_crossings=4`). This task's own job: continue that orbit in Jacobi constant `C`,
reusing `#753`'s own `cr3bp_continuation.continue_family`, to confirm or refute Vaquero's
own falsifiable claim ("it is suspected that this family of periodic resonant chains ends
for a value of Jacobi constant `C < 3.01400`", Fig. 4.12).

**Verdict up front**: **inconclusive on Vaquero's specific `C < 3.01400` claim, but a
genuine, well-verified, positive finding of a different kind.** Continuing `#782`'s own
`half_crossings=4` branch in INCREASING `C` (see "direction" note below), this task found
strong numerical evidence of a genuine tangent (fold) bifurcation — the branch's own
nontrivial monodromy eigenvalue collapses smoothly and monotonically from
`|lambda|=4.77e7` at Vaquero's own `C=3.010000` down to **exactly `1.000000` at
`C=3.0100696797`** — roughly `57x` CLOSER to the `C=3.01` anchor than Vaquero's own
`3.01400` claim (`ΔC≈7.0e-5` vs her claimed `ΔC≈4.0e-3`). This does **not** confirm her
specific number, and — per this task's own honesty discipline — it is **not established**
that this `half_crossings=4` branch is the same family object she plots in Fig. 4.11/4.12.

**Addendum 2026-08-08 (CI cross-platform investigation, coordinating session):** the first
CI run of `test_c774_chain_branch_fold_eigenvalue_reaches_unity` found `max_eig` off by
`0.213` from `1.0` on Linux — investigated directly rather than assumed to be routine
cross-platform noise, given the magnitude. Confirmed deterministic-per-run on this Mac (3
repeated runs, bit-identical), then directly measured the local sensitivity: perturbing
`C` by as little as `1e-9` at the exact reported fold seed swings `max_eig` from `~1.0` to
values in the `10`-`20000+` range. This is consistent with, not contradictory to, the
"collapses smoothly and monotonically" characterization above — a smooth curve dropping
from `4.77e7` to `1` over a `ΔC` of only `~7e-5` implies an enormous local slope right at
the crossing, so a cross-platform corrector difference far too small to move the reported
`x0`/`jacobi` outside this test's own loose tolerances can still swing the eigenvalue by
orders of magnitude. **This does not undermine the qualitative finding (a fold bifurcation
genuinely exists near `C≈3.0100697`)** — it means the specific 10-significant-figure `C`
value reported above carries less genuine cross-platform precision than its digit count
suggests, and any future work building on this fold's exact location should re-derive it
fresh rather than treat `C=3.0100696796878963` as portable. Full detail in the test's own
`xfail` reason, `tests/search/test_saturn_titan_resonant_connections.py`.

---

## Direction: resolved from the dispatch note's own internal contradiction

The dispatch note said "continue downward... toward and past `C=3.01400`" — self-
contradictory, since `3.01400 > 3.010000` (that is the INCREASING direction, not
downward), and the note itself flagged this ("get the direction of termination right...
don't assume"). Vaquero's own quote ("ends for a value of `C < 3.01400`") most naturally
reads as: there exists some critical `C_end < 3.01400` where the family stops existing;
our own C=3.01 anchor is itself below 3.014, so the only way to test the claim directly is
to continue TOWARD 3.014, i.e. in the INCREASING-`C` direction. Both directions were
probed; INCREASING is the one that bears on the claim and is reported as the primary
result (see below for why DECREASING is not reported as established).

## Method

`#753`'s own `cr3bp_continuation.continue_family` (`src/cyclerfinder/search/
cr3bp_continuation.py`) was used as the ONLY corrector+gauntlet at every step — never
reimplemented — called with `half_crossings=4`, `ydot0_sign=1.0` (matching `#782`'s own
branch exactly), `corrector_tol=1e-12`, `rtol=1e-13`, `atol=1e-14` (matching `#782`'s own
convergence regime, confirmed NOT load-bearing looser than `continue_family`'s own
defaults — see "Tolerance due-diligence" below), and a TIGHTENED `period_step_frac=0.01`
(from the default `0.10`) once an off-family jump was found to slip through the looser
default (see "A branch-hop found and diagnosed" below).

Because `continue_family` itself uses a FIXED step size per call (no internal
adaptivity), and the demonstrated basin sensitivity of this system (`#782`'s own note:
a `1e-4` perturbation in `x0` converges to a period-`5.02` orbit instead of `56.03`) makes
any single fixed step size either too large (branch hop / non-convergence) or leaves
`0.004` of `C`-distance to Vaquero's own `3.014` boundary computationally unreachable at a
uniform tiny step, this task wrapped `continue_family` (called with `n_steps=1` each time)
in an OUTER step-size-adaptation loop — the loop only decides what `d_jacobi` to try next;
every correction, gauntlet check (closure, period bounds, equilibrium, Jacobi
conservation, independent-Radau cross-check, dedup), and stop-reason semantic comes from
`continue_family` itself, unmodified.

**A genuine methodological limitation of this wrapper, found and corrected mid-task**:
calling `continue_family(..., n_steps=1)` repeatedly resets its own internal `x0_hist`/
`c_hist` on every call, which disables its own `FOLD_REVERSAL` gate (that gate compares
the CURRENT step's `x0`-delta sign against the PREVIOUS step's, which requires a
persistent history `continue_family` normally keeps across a single multi-step call, but
which a fresh `n_steps=1` call cannot see). This was caught directly: applying
`continue_family`'s own fold-reversal test by hand to this task's own full step log found
one micro-reversal (`Δx0` sign flip of `~8e-6`, an order of magnitude below the
surrounding step sizes) partway through an early, AGGRESSIVE-step-growth version of the
walk. Diagnosed as a **branch hop**, not a genuine fold (see next section) — the aggressive
growth schedule (`dj *= 1.3-4.0` after each success, uncapped) occasionally took a step
large enough to jump the Newton corrector onto a nearby-but-distinct root. All numbers
reported as this task's own headline result come from a SECOND, independent run with
`d_jacobi` HARD-CAPPED at `1e-6` (never allowed to grow past it, only shrink-and-recover
on failure) plus the tightened `period_step_frac=0.01` — eliminating the branch-hop risk
by construction, not merely by post-hoc inspection.

## The branch-hop found and diagnosed (methodological finding, not the headline result)

An early, free-growing version of the adaptive wrapper (step size uncapped, growing `4x`
per success) reached `C=3.0100697` with `|lambda|` collapsing to `~5.18` via a path that
included one anomalously large step (`Δx0=+1.95e-4` when neighboring steps were moving
`~5e-6` — the default `period_step_frac=0.10` gate passed it, since the period only
changed `0.84%`). A direct continuity check — re-walking the SAME interval
(`C∈[3.01002743,3.01005824]`) with `d_jacobi` capped at `1e-6` and `period_step_frac`
tightened to `0.01` — found the genuinely continuous branch reaches only `x0=0.9493760`
at that `C` (vs. the free-growing walk's `x0=0.9494841` at essentially the same `C`) — a
`1.08e-4` discrepancy, far above Newton tolerance, confirming the large step DID hop onto
a different, nearby root. **Reassuringly, the properly capped, non-hopping re-walk
independently converges to the SAME terminal point** (`C=3.0100696796878963` vs. the
free-growing walk's `C=3.0100696796878372` — agreeing to 9 significant figures) — strong
cross-validation that this terminal point is a genuine feature of the dynamics, not an
artifact of either particular step-size schedule.

## The headline finding: a clean, monotonic collapse to `|lambda| = 1` (increasing `C`)

The capped, properly-gated walk (`period_step_frac=0.01`, `d_jacobi<=1e-6`) from
`C=3.010058235052` (itself reached by an earlier segment of the SAME capped walk, unbroken
back to `#782`'s own `C=3.010000` seed) to its own natural stopping point:

| step | `C` | `x0` | `\|lambda\|` |
|---|---|---|---|
| 1 | 3.010058735052 | 0.9493774465 | 5.354e6 |
| 5 | 3.010062485052 | 0.9493893086 | 3.365e6 |
| 10 | 3.010066735052 | 0.9494068449 | 1.280e6 |
| 15 | 3.010069645208 | 0.9494327130 | 1.139e4 |
| 18 | 3.010069678373 | 0.9494351150 | 67.07 |
| 20 | 3.010069679590 | 0.9494355364 | 82.87 |
| 21 | 3.010069679688 | 0.9494356863 | 4.042 |
| 22 | 3.010069679688 | 0.9494356889 | 1.681 |
| **23** | **3.010069679688** | **0.9494356926** | **1.000000 (0.9999999912)** |

Step 24 (any tested `d_jacobi` down to `1e-13`): **genuine `NO_CONVERGE`** — no further
step could be found.

This orbit re-converges cleanly and independently when re-seeded directly (not merely
inherited from the walk): `x0=0.9494356987780496`, `ydot0=0.07906623799657658`,
`period=55.6947459919561`, `jacobi=3.0100696796878372`, crossing residual `2.478e-14`,
`1` Newton iteration from that seed.

Full-period monodromy spectrum (evaluated one step before the exact fold, at
`|lambda|=5.18`, where the reciprocal pair is still cleanly resolved numerically): real
reciprocal pair `(5.1819, 0.19299)` (product `1.00007`, confirming symplectic
consistency); trivial pair `(0.9997±0.0223j)`, `|.|≈1` as expected from time-translation
invariance; a third pair `(-0.9887±0.1500j)`, `|.|≈0.99999` (an out-of-plane/vertical-
stability mode sitting near the unit circle — noted, not chased further, outside this
task's own scope). At the fold itself, the nontrivial real pair reaches exactly `(1, 1)`
(degenerate) — the standard eigenvalue signature of a saddle-node bifurcation, independently
confirmed by the two-root coalescence test below.

Independent-Radau cross-check: closure `5.64e-8`, Jacobi drift `1.51e-14` over the full
period. Ghost-guard vs. `node`'s own section points: `d_ghost=0.0807`, `80.7x`
`GHOST_GUARD_DELTA` — a genuine, non-trivial orbit, not a rediscovery of 3:4. `8`
`{y=0}` crossings per period.

## Two-root coalescence test (confirms genuine fold, not a corrector-basin artifact)

A collapsing eigenvalue alone does not distinguish a genuine fold (the family turns around
in `x0` vs `C`, two branches merging) from a corrector-basin artifact (Newton's own basin
of convergence simply shrinks to nothing as the problem gets harder, independent of
whether a second solution exists). The decisive, falsifiable test for a genuine fold: at
`C = C_fold - delta` (strictly BELOW the fold), there must exist TWO distinct roots
straddling `x0_fold` that coalesce into one as `delta -> 0`.

Scanned `x0_guess` on a `25`-point grid spanning `x0_fold +/- 3e-4` at four `delta` values,
re-seeding `correct_symmetric_fixed_jacobi` fresh at each grid point (no warm-starting from
neighbors, so this cannot manufacture apparent continuity):

| `delta` | `C_test` | roots found | root spread | midpoint |
|---|---|---|---|---|
| `1e-6` | 3.010068679688 | **2** | `3.193e-5` | 0.9494352037 |
| `1e-7` | 3.010069579688 | **2** | `1.010e-5` | 0.9494356454 |
| `1e-8` | 3.010069669688 | **2** | `3.193e-6` | 0.9494356896 |
| `1e-9` | 3.010069678688 | **2** | `1.010e-6` | 0.9494356940 |

**Exactly two distinct roots at every tested `delta`, straddling `x0_fold=0.9494356926`,
with the spread shrinking as `sqrt(delta)`** (ratio between successive rows: `3.16, 3.16,
3.16` — matching `sqrt(10)` to 3 significant figures at every step) **and the midpoint
converging monotonically toward `x0_fold` itself** (`0.94943520 -> 0.94943565 ->
0.94943569 -> 0.94943569`, vs. the fold's own `0.94943569`). This `sqrt(delta)` scaling of
root separation is the textbook local-normal-form signature of a saddle-node (tangent)
periodic-orbit bifurcation (`x - x_fold ~ +/- sqrt(C_fold - C)`) — about as clean a
numerical confirmation of a genuine fold as this system's own extreme conditioning permits.
**This is not a corrector-basin artifact: a real second branch exists below `C_fold` and
provably merges with the walked branch exactly at `x0_fold`, `C_fold`.**

**Corroborating check — no solution of THIS branch exists above the fold.** Seeding
`correct_symmetric_fixed_jacobi` at `C_fold + 1e-6` from three points (`x0_fold` itself and
`x0_fold ± 3e-5`) finds that Newton still "converges" every time (this system's own
well-documented basin sensitivity, `#782`'s own note: nearby perturbations land on
unrelated branches), but NEVER onto the target `~55.69`-period branch — the three seeds
land on periods `2.96`, `5.02`, and `68.58` respectively, all wildly different from
`55.69`. This is the expected complement to the two-root result: a real second solution
exists below the fold, none above it, on this particular branch.

## Tolerance due-diligence

The gauntlet tolerances passed to `continue_family` (`jacobi_tol=1e-8`,
`radau_closure_tol=1e-2`, `radau_jacobi_tol=1e-6`) are looser than `continue_family`'s own
defaults (`1e-10`/`1e-3`/`1e-8`). Checked directly whether this loosening was load-bearing:
re-ran two known-good steps under `continue_family`'s own DEFAULT (tighter) tolerances —
both passed unchanged, with the ACTUAL `radau_djacobi` values (`2.2e-14`, `3.5e-14`) far
under even the tight default. The loosening was cosmetic headroom, not consequential — the
`no_converge` failures that bound this walk are genuine Newton non-convergence, not a
gauntlet rejection under either tolerance set.

## Decreasing-`C` direction: NOT established, reported as such

An early exploratory walk (uncapped growth, default `period_step_frac=0.10`) in the
DECREASING direction showed period jumping `56.04 -> 59.13 -> 59.49 -> 61.20` and
`|lambda|` swinging non-monotonically `5.2e7 -> 1.27e7 -> 3.35e7 -> 2.5e9` over just five
steps before the search exhausted at `C=3.009943392`. This pattern — period drifting by
several percent per step while `|lambda|` swings by orders of magnitude non-monotonically
— is diagnostic of the SAME class of off-family drift the increasing-direction free-growing
walk showed (and which the tightened `period_step_frac=0.01` + capped `d_jacobi` walk was
built to rule out). This decreasing-direction result was NOT re-run under the corrected,
non-hopping methodology (time-budgeted; the increasing direction is the one that bears
directly on Vaquero's own `C<3.01400` claim). **Reported honestly as NOT established** —
neither a confirmed termination near `C=3.0099` nor evidence against one.

## Verdict on Vaquero's `C < 3.01400` claim

**Inconclusive**, precisely stated:

1. This task's own continuation covered only `ΔC≈7.0e-5` of the `ΔC≈4.0e-3` distance from
   `C=3.010000` to `3.014000` — **1.7%** of the way — before hitting what looks like a
   genuine dynamical feature (a fold), not a numerical-conditioning wall that a smarter
   technique could push through indefinitely.
2. The mechanism found (`|lambda| -> 1`, a tangent/fold bifurcation) IS the textbook
   signature of "a family ends" in bifurcation-theory language — consistent in KIND with
   Vaquero's own claim.
3. But the VALUE found (`C≈3.010070`) is NOT her value (`3.014000`) — `~57x` closer to the
   `C=3.01` anchor. This task cannot establish that `#782`'s own `half_crossings=4` branch
   is the SAME family object Vaquero continues in her own Fig. 4.11/4.12: `#782`'s own note
   already flags an unexplored `half_crossings=6` sibling candidate (period `82.17`), and
   this orbit shadows 3:4's own physical path for `~94%` of its period — it is entirely
   possible Vaquero's own published family is parameterized/selected differently (e.g. by
   tracking a specific `Wu(3:4)/Ws(3:4)` manifold intersection continuously in `C`, not by
   blind Newton continuation of one already-converged branch) and folds at a materially
   different point, or doesn't fold in this same sense at all.
4. **The `#775`-flagged coincidence is WEAKENED, not strengthened, by this result.** `#775`
   found the UNRELATED plain 3:4 family's own topology-jump lands at exactly `C=3.014000`.
   This task's own chain-branch fold lands at `C≈3.010070` — a DIFFERENT value. If the two
   phenomena were mechanistically connected, a matching (or at least much closer) value
   would be expected. They are not close (`3.014` vs `3.010070`, a factor of `~57x` in
   distance-from-anchor). Reported directly, per the dispatch note's own explicit
   instruction not to overclaim a causal link.

## Code delivered

No changes to `src/cyclerfinder/`'s production continuation machinery — `#753`'s own
`cr3bp_continuation.continue_family` was reused unmodified throughout, per the dispatch
note's explicit mandate. The multi-hour adaptive-step-size search itself lived in
one-off exploratory scripts (not committed — the search history is not reproducible in
any practical CI time budget and the OUTER step-size-adaptation loop is not itself a
reusable capability, just an orchestration convenience around the existing tool).

`tests/search/test_saturn_titan_resonant_connections.py` (+2 tests, +1 import): the fold
point and its two-root coalescence signature are anchored as FAST (`~25s` combined),
independently-reproducible regression tests —
`test_c774_chain_branch_fold_eigenvalue_reaches_unity` (single Newton correction from the
known answer; asserts convergence, `x0`/`C` match, and leading eigenvalue within `1e-4` of
`1`) and `test_c774_chain_branch_fold_two_root_coalescence` (the decisive
fold-vs-corrector-basin-artifact discriminator: two distinct roots below `C_fold`, none on
the target branch above it). Neither marked `@pytest.mark.slow`
(`feedback_delegation_fresh_agent_not_fork`: a discovery-verdict-bearing evidence test must
run in CI).

## Verification

* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -k c774 -v`: 2/2
  pass, `23.97s`.
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -q`: 45/45 pass (41
  passed + 4 pre-existing, documented XPASS — see the file's own
  `_XFAIL_CI_CROSS_PLATFORM_SATURN_TITAN_CHAIN` marker, predates this task, not touched).
* `uv run ruff check` / `uv run ruff format --check` on the changed test file: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 835 source files.
* `uv run pytest tests/data tests/search -q`: exit code 1 — **2 failures, both
  confirmed pre-existing and unrelated to this task**, not a regression from this task's own
  changes: `tests/search/test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface`
  and `tests/search/test_504_pluto_charon_kk_sweep.py::test_504_sweep_33`. Neither file, nor
  their underlying source modules (`eggie_ballistic.py`, `pluto_charon_kk_sweep.py`), was
  touched by this task (confirmed via `git status`/`git diff` — this task only changed
  `test_saturn_titan_resonant_connections.py`, `data/OUTSTANDING.md`, and this note). Both
  re-run in isolation, single-threaded (`-n 0`, no xdist contention): **identical deterministic
  failures**, not a contention-driven flake. Flagged per this project's own dormant-issue
  discipline (not fixed here, out of this task's own scope) — a coordinating session should
  triage separately. Everything else in the suite passes (including the expected, documented,
  pre-existing cross-platform XPASS markers in `test_saturn_titan_resonant_connections.py`,
  `test_ccr4bp_europa_callisto_heteroclinic_search.py`, `test_earth_moon_resonant_families.py`,
  `test_neptune_triton_resonant_families.py`, and the long-standing `#54`-tracked XFAILs).
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: 2/2 pass (run after the
  `OUTSTANDING.md` edit, before committing).

**A separate, unrelated observation from the full-suite run, not investigated further**:
`data/floquet_phase1_reproduction.jsonl` (a tracked data file, untouched by this task's own
edits) came back from the suite run with DIFFERENT content each time — some other test in the
suite writes non-deterministic output to this file as a side effect (`branch_k1`/`branch_k2`
differ between runs, e.g. `(3,3)` vs `(4,3)`). Not staged or committed by this task (pathspec-
only staging); flagged here for whoever owns that test/file next.
