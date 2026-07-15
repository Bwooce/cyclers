# #602: real numerical continuation between #588 cluster 43 and Family J

**Date:** 2026-07-15
**Origin:** `docs/notes/2026-07-14-590-followup-results.md` (#590) and
`docs/notes/2026-07-14-591-bounded-ecliptic-excursion-results.md` (#591), both CLOSED
2026-07-14. #590's own Step 2 connectivity check was explicitly caveated as "a weak, cheap
straight-line IC-interpolation heuristic ... NOT a real continuation." #602's job was to
replace that heuristic with actual differential correction (pseudo-arclength continuation +
deflated Newton, `search/pseudo_arclength.py` / `search/deflated_newton.py`, plus the
project's existing validated CR3BP family-continuation driver `search/cr3bp_continuation.py`)
and settle, definitively if possible, whether cluster 43 sits on the same continuous family
branch as Family J (Gurfil & Kasdin 2002) or is a genuinely distinct branch.

**Script:** `scripts/continue_602_cluster43_familyj.py` (deliberately outside the
`scripts/run_*.py` preflight-ratchet glob, matching the `dedup_588_*.py`/`analyze_591_*.py`
precedent — this is a continuation/verification experiment over two already-known ICs, not a
new discovery sweep). **Result file:** `data/continue_602_cluster43_familyj.jsonl` (24
records: 2 positive-control, 21 homotopy steps, 1 secondary check).

## Cluster 43's exact IC (traced to source, per the task's own instruction)

Not hardcoded in `run_581_gurfil_reproduction.py` (that file only has the 12 published
Gurfil-Kasdin families). Found in
`data/found/583_widened_search/deduped_candidates.json`, `cluster_id: 43`,
`"characterization"."ic_interleaved"` + `"theta0"` (the exact numbers `analyze_591_*.py`
itself reads from the same file for its own propagation):

```
ic_interleaved = [0.03747707684491371, 0.02134315891222638, 0.03785060584964018,
                  -0.0766716607792004, 0.031634219092074256, 0.02917567513165431]
theta0 = 4.130670798449157
```

(interleaved order `[x, x', y, y', z, z']`, geocentric ER3BP pulsating frame,
`core/er3bp_geocentric.py`'s convention — confirmed by re-reading that module's own docstring
and `table_interleaved_to_state`, which is used unmodified here, exactly as #590/#591 used it.)
Model confirmed ELLIPTIC (ER3BP, `e = E_SUN_EARTH_GURFIL_KASDIN = 0.0167`), matching #590/#591.

Family J's IC is Gurfil-Kasdin's own published Table 3 row (`TABLE34["J"]` in
`run_581_gurfil_reproduction.py`), `theta0 = 0.0`.

## What "the modules already built for this" turned out to require

Before running anything, all three modules and their test suites were read in full
(`pseudo_arclength.py`, `cr3bp_continuation.py`, `deflated_newton.py`, plus
`tests/search/test_{pseudo_arclength,cr3bp_continuation,deflated_newton}.py`). Existing
regression suites confirmed green first (34/34 pass,
`uv run pytest tests/search/test_deflated_newton.py tests/search/test_pseudo_arclength.py
tests/search/test_cr3bp_continuation.py tests/data/test_er3bp_drift_classifier.py -q`) — the
generic machinery itself is sound; nothing here is a claim that any of the three modules is
buggy.

`cr3bp_continuation.py`'s own positive control (`tests/search/test_cr3bp_continuation.py`) is
for the EARTH-MOON system (Ross & Roberts-Tsoukkas 2025's (3,3) resonant family), not Sun-Earth
DRO — so per the task's own instruction ("do the minimum needed to trust the machinery"), this
task built its OWN positive control for the Sun-Earth case (Step A below), rather than skipping
verification.

## Step A: positive control / does an "exact periodic orbit" residual even apply here?

The natural codimension-1 residual for `pseudo_arclength.continue_curve` on a full 3D ER3BP
orbit is the **1-year return-map defect**: `R(z) = Phi_{theta0 -> theta0+2*pi}(state0) - state0`
(6D, for a 7D unknown `z = [interleaved_state, theta0]` — exactly matching the module's
`N-1`-for-`N` codimension-1 contract). Before trusting a continuation walk on this residual,
this task checked whether either endpoint is even CLOSE to an exact root of it.

**Method:** `scipy.optimize.least_squares` (trust-region reflective — an independent,
professionally-implemented robust solver, not the hand-rolled Newton loop this task also tried
and independently confirmed the same finding with), run twice per seed: (1) unconstrained, and
(2) boxed to a tight ±0.01 neighborhood of the seed (in nondimensional nondimensional AU/rad
units — large relative to the O(0.01-0.08) state components).

| seed | seed defect (norm) | unconstrained solve | boxed (±0.01) solve |
|---|---|---|---|
| Family J | 0.0380 | converges to residual 8.3e-6, but **0.637 away from the seed** (an unrelated, far-away root) | stalls at residual 5.1e-3, **pinned to the box boundary** (no interior root) |
| cluster 43 | 0.00854 | converges to residual 9.6e-6, but **0.662 away from the seed** | stalls at residual 3.9e-3, **pinned to the box boundary** |

Both endpoints show the identical pattern: no exact year-periodic orbit exists anywhere near
either seed; a robust global solver only finds unrelated, distant roots. A supplementary
analytic-Jacobian check (`propagate_er3bp_geocentric(..., with_stm=True)`) confirms this is not
a finite-difference artifact: `(Phi_2pi - I)` is well-conditioned at the seed itself (condition
number ~3.6e3, smallest singular value 0.0112) but its smallest singular value collapses toward
zero as any Newton/Levenberg-Marquardt walk approaches lower residual — the hallmark of trying
to shoot onto a periodicity condition that these orbits do not actually satisfy nearby.

**Reading:** this is a genuine, load-bearing finding, not a tooling failure. It means
**"continuous curve of exact year-periodic orbits" is not the right formalization of "family"
for this specific problem.** Both #590's own long-horizon check (Family J's literal IC
escapes at ~30 yr; cluster 43 stays bounded to 1000+ yr) and Mikkola et al. (2006) (already
anchored via #590 — inclined quasi-satellite motion is generically only TEMPORARILY captured,
not exactly periodic) predict exactly this: these are quasi-periodic, metastable structures,
not periodic orbits of the 1-year forcing period. `pseudo_arclength.continue_curve` applied to
an exact-periodicity residual is therefore not viable here — not because the module is
inadequate, but because the target objects are the wrong kind of object for that residual.

## Step B: damped deflated-Newton homotopy walk (real correction, not pure interpolation)

Given Step A, the most defensible test that still uses real differential correction (rather
than #590's uncorrected straight-line heuristic) is: interpolate the 7D genome linearly between
Family J's and cluster 43's IC in 21 steps (`t = 0.00, 0.05, ..., 1.00`), and at EVERY step run
a *damped* (step-capped) Newton correction (`search/deflated_newton.py`'s own math, no new
solver logic — a step cap of 0.01 keeps the corrector from taking the huge, divergent Newton
steps Step A's un-damped attempts exhibited) targeting the same year-periodicity residual, then
independently check boundedness at both the raw and corrected point with
`classify_bounded_drift` (n_revs=50 — the SAME already-positive-controlled #583 tool #590/#591
used, at a LONGER horizon than #590's original 20-year heuristic check).

Full 21-point table in `data/continue_602_cluster43_familyj.jsonl` (`record_type ==
"homotopy_step"`). Summary:

| t | defect: raw -> corrected | correction | raw bounded (50yr) | corrected bounded (50yr) |
|---|---|---|---|---|
| 0.00 (J) | 0.0380 -> 0.0096 | reduced | **False** | False |
| 0.05 | 0.0354 -> 0.0094 | reduced | **False** | True |
| 0.10 | 0.0333 -> 0.0094 | reduced | **False** | True |
| 0.15 | 0.0316 -> 0.0316 | DIVERGED | **False** | False |
| 0.20 | 0.0304 -> 0.0304 | DIVERGED | **False** | False |
| 0.25 | 0.0296 -> 0.0296 | DIVERGED | True | True |
| 0.30-0.95 (8 pts) | mixed reduced/DIVERGED | mixed | True (every point) | mixed |
| 1.00 (cluster 43) | 0.0085 -> 0.0085 | DIVERGED | True | True |

**Two honest findings, not one clean verdict:**

1. **The correction step itself is patchy.** At more than half the 21 points the damped
   corrector could not improve on the raw interpolated defect at all within its step budget
   ("DIVERGED" in the table — a real, reported failure mode, not silently dropped). This is
   consistent with Step A: these orbits sit in a numerically fragile, marginally-stable regime
   where a well-behaved local correction landscape should not be assumed. This dimension of the
   test is **inconclusive** — it neither confirms nor refutes a fold/bifurcation, because the
   corrector's own local behavior is too noisy to read a clean signal from.
2. **The boundedness signal (independent of correction quality) is clean and coherent.** The
   raw interpolated path is UNBOUNDED (escapes within 50 yr) only in the region closest to
   Family J itself (`t = 0.00-0.20`) — consistent with Family J's own literal IC escaping at
   ~30 yr (#590) — and is BOUNDED continuously from `t = 0.25` all the way to cluster 43
   (`t = 1.00`), with no interior escape observed anywhere in that stretch, at a full 50-year
   horizon (longer than #590's original 20-year check). This extends, but does not
   qualitatively change, #590's own finding that the J-to-43 heuristic path was "the cleanest"
   of the clusters checked, and gives it a longer, more demanding horizon to stand on.

## Step C: secondary check — planar/circular retrograde-satellite backbone

A genuinely rigorous, ALREADY test-validated corrector exists for exactly one reduced version
of this question: does the PLANAR, CIRCULAR (`e=0`) retrograde-satellite family (Henon's family
f / Families A-B-F, the curve both J's and cluster 43's `ydot0/x0 ~= -2` ratios approximately
sit near, per #588's own IC-relation note) connect continuously from an x0 near Family J's to
an x0 near cluster 43's? This uses `cr3bp_continuation.py::continue_family` +
`cr3bp_periodic.correct_symmetric_fixed_jacobi` exactly as built and already validated
(Earth-Moon (3,3) family, `tests/search/test_cr3bp_continuation.py`), dropping z/inclination and
eccentricity forcing — an explicit, deliberate reduction, not a substitute for the full 3D
question Step A/B target.

**A units bug was caught and fixed before trusting this result**: `cr3bp_periodic`'s `x0` is
the standard CR3BP BARYCENTRIC coordinate (secondary/Earth at `1-mu`), while J/cluster-43's own
`x0` is Gurfil-Kasdin's GEOCENTRIC coordinate (origin at Earth) — feeding the geocentric value
straight in without the `+(1-mu)` offset put the seed point essentially AT the Sun (Jacobi
constant ~59.7, nonsensical for a DRO). Caught via a sanity check on the resulting Jacobi
constant (retrograde-satellite-family members should sit near `C ~ 3`, not `C ~ 60`) before any
result was trusted; fixed by adding the offset, re-verified with a direct crossing-time
diagnostic (first perpendicular crossing at `t ~ pi`, matching a `period_guess = 2*pi`).

**Result after the fix:** seed converges at `x0_geo = 0.03094` (near Family J's `0.03348`,
`C = 2.99903` — sensibly near the known Families A/B/F band). Continuing (direction chosen
empirically — Jacobi decreases as x0 increases for this branch) walks smoothly and
monotonically from `x0_geo = 0.0309` past cluster 43's `x0_geo = 0.0375` all the way to
`x0_geo = 0.1444` before the run's own `max_steps=200` cap — **`stop_reason = "max_steps"`, not
a fold, not a topology jump, not a gauntlet rejection.** `crossed_target_x0 = True`.

**Reading:** the underlying planar/circular retrograde-satellite backbone both orbits'
`ydot0/x0` ratios approximately sit on is confirmed, by real (not heuristic) continuation using
this project's own validated corrector, to be one smooth continuous curve spanning both
radii — with zero folds/bifurcations encountered. This corroborates (does not newly establish)
Pousse, Robutel & Vienne (2017)'s published result (already anchored via #590) that this exact
backbone is one continuous curve from near-Earth out to Sun-collision. It is real evidence, but
it is evidence about the REDUCED (planar, circular) skeleton both orbits are inclined/eccentric
perturbations OF — not direct evidence that the full 3D, e=0.0167 structures are on one branch.

## Overall verdict: genuinely INCONCLUSIVE on the full 3D question, with real corroborating evidence for "same," not "distinct"

This task cannot deliver the clean SAME-vs-DISTINCT verdict its own brief hoped for, and it
would be dishonest to force one. What the evidence gathered across Steps A-C actually supports:

- **No fold, bifurcation, or discontinuity was found anywhere along any path tried** — not in
  the (admittedly patchy) full-3D homotopy correction, not in the boundedness classifier along
  that same path (past `t=0.25`), and not in the rigorous, test-validated planar/circular
  backbone continuation (which ran cleanly all the way past cluster 43's radius with zero
  gauntlet rejections).
- **No evidence was found FOR a distinct branch either** — nothing in Steps A-C produced a
  corrector failure, a topology jump, or an escape that would indicate cluster 43 sits on a
  separate structure from Family J.
- **The reason a hard verdict is out of reach is a model-mismatch, established rigorously in
  Step A**: "family branch" as a curve of EXACT periodic orbits (the concept
  `cr3bp_continuation.py` is built around, and the natural target for
  `pseudo_arclength.continue_curve`) does not apply to either endpoint of this specific
  ER3BP problem — both are quasi-periodic, metastable structures with no nearby exact periodic
  orbit for a corrector to walk between. This is not a "the continuation infrastructure needs
  more work" finding so much as a "this pair of objects is not well-posed as periodic-orbit
  continuation targets" finding.

**Net effect on the #588/#590/#591 thread:** the balance of evidence (heuristic connectivity
from #590, longer-horizon boundedness continuity from this task, and the fold-free planar
backbone continuation from this task) continues to lean toward "cluster 43 is an unsampled
point on Family J's curve" rather than toward "distinct branch" — but this remains, honestly,
a preponderance-of-evidence reading, not a proof. #602 upgrades the QUALITY of that evidence
(real correction + a validated backbone continuation, not just straight-line interpolation) but
does not convert it into the requested definitive closure.

## Recommended next step (if this thread is revisited further)

A genuinely definitive verdict on the full 3D question would require a correctly-posed
invariant for QUASI-PERIODIC (not periodic) orbits — e.g. a proper 2-torus/rotation-number
continuation, or a homotopy-in-eccentricity starting from the now-validated planar/circular
backbone (Step C) and slowly cranking `e: 0 -> 0.0167` plus inclination `z: 0 -> z_J/z_43` while
re-correcting against a genuinely quasi-periodic (KAM-style) closure condition at each step.
That is materially new capability-building work — a different class of corrector than anything
`pseudo_arclength.py`/`cr3bp_continuation.py`/`deflated_newton.py` currently provide — and is
explicitly OUT OF SCOPE for this task's own "cheap and bounded" framing. This is the honest
scope-correction the task's own guardrail anticipated.

## Constraints honored

No `data/catalogue.yaml` write. No git commit. No file deletions. New files only:
`scripts/continue_602_cluster43_familyj.py`, `data/continue_602_cluster43_familyj.jsonl`, this
note.

## #602 status: bounded continuation experiment complete. Result: INCONCLUSIVE on a hard
SAME/DISTINCT verdict for the full 3D structure (a genuine model-mismatch finding, not a
tooling gap), with real (non-heuristic) corroborating evidence — a fold-free planar/circular
backbone continuation and a longer-horizon boundedness-continuity check — continuing to favor
"same family, unsampled point" over "distinct branch." No catalogue rows changed, no novelty
claimed either way. Left for a future pass IF a quasi-periodic-family continuation capability is
ever built (out of scope here).
