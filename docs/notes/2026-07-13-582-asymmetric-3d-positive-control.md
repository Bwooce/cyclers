# #582 (stage 3a of #581): asymmetric/spatial-isolated 3D CR3BP niching-GA -- fitness build + positive control

**Date:** 2026-07-13
**Scope of this note:** the positive-control deliverable only, per #582's own gate
("Do NOT run the full 5-MMR asymmetric novelty sweep yourself -- that's a
longer job the coordinator will launch and own separately"). No novelty claim
is made anywhere in this note.

## What was built

* `src/cyclerfinder/search/isolated_3d_asymmetric_fitness.py` -- the new GA
  fitness function. Genome `(x0, z0, xdot0, ydot0, zdot0, T)`, `y0=0`.
  Fitness = `1 / (1 + w_period*||X(T)-X(0)||^2 + w_jacobi*((C-C_target)/sigma)^2)`
  (bounded reciprocal form, same shape as Gurfil-Kasdin 2002 Eq. 15). Death
  penalties (fitness forced to 0.0): `T <= T0/2`, primary/secondary collision
  at any point in the propagated arc (event-based, not endpoint-only), and
  non-finite defect/Jacobi/integrator failure. `mmr_bounds()` builds the GA
  search box for one interior MMR, centered on the #440 analytic e=0 guess.
* `src/cyclerfinder/search/isolated_3d_asymmetric_pipeline.py` -- the
  mandatory downstream pipeline: `refine_ga_candidate()` (routes through the
  EXISTING `cr3bp_general_periodic_3d.correct_general_periodic_3d`, no new
  corrector), `classify_symmetry()` (explicit perpendicular x-z-plane-crossing
  test -- a real function, not a comment), `build_candidate_signature()`
  (populates a `literature_check.CandidateSignature`), and
  `literature_anchors_engaged()` (proves the matcher reaches a non-empty
  anchor pool for a given signature).
* `scripts/run_582_asymmetric_3d_niching_search.py` -- driver. `--mmr p:q`
  selects one of the 5 tabulated interior MMRs; `--positive-control` runs the
  full GA -> corrector -> symmetry -> signature -> compare pipeline
  end-to-end; `--mode ga` runs (or resumes) just the checkpointed GA for a
  future, larger-budget novelty sweep. Calls `preflight_search()`
  unconditionally (task_no=582). `--workers` (default 8, this machine's
  logical core count).
* Tests: `tests/search/test_isolated_3d_asymmetric_fitness.py` (16 tests),
  `tests/search/test_isolated_3d_asymmetric_pipeline.py` (10 tests) -- 26
  new tests, all passing.

## Positive control result: PASS (MMR 3:2)

Command: `uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --positive-control --workers 8`
(wall time ~14s on this machine, 8 workers, small GA budget 40 pop x 60 gen).

The 3:2 MMR was chosen deliberately: it is the module docstring's own
"GO-gate-validated recipe" and, per the existing `test_er3bp_isolated_seeds.py`
comment, "the most eccentric member and the tightest test of the recipe" --
the hardest of the 5 tabulated cases, not the easiest.

| Quantity | Known #440 member | GA + corrector result | Tolerance | Error | Margin |
|---|---|---|---|---|---|
| x0 | 0.7310974 | 0.7356006 | 3% rel | 0.616% | 4.9x |
| ydot0 | 0.4569989 | 0.4389736 | 5% rel | 3.944% | 1.3x |
| T | 11.926512 | 11.433699 | 5% rel | 4.132% | 1.2x |
| C (Jacobi) | 3.0622660 | 3.0676539 | 0.02 abs | 0.00539 | 3.7x |

Corrector: converged=True, corrector_residual=1.27e-12, independent (Radau)
closure residual=1.23e-12 (independent-integrator cross-check per
[[feedback_orbit_closure_discipline]]). Symmetry classification:
`is_symmetric=True`, best perpendicular-crossing residual 1.3e-15 (machine
zero) -- correctly identifies that this survivor, though found via the
*unconstrained* asymmetric corrector with no symmetry constraint imposed,
sits on the known SYMMETRIC orbit (which it must, since it's the known #440
seed's own family). This is itself a positive control for
`classify_symmetry()`: it fires correctly on a case where the ground truth
(symmetric) is known. `degenerate_planar=True` on the refined orbit
independently confirms the same thing (z0, zdot0 -> ~0).
`build_candidate_signature()` + `literature_anchors_engaged()` returns a
non-empty anchor list including "Antoniadou & Libert spatial resonant
periodic orbits in the RTBP (2019)" -- confirms the matcher-engagement fix
works (see caveat below on over-inclusiveness).

Full machine-readable result:
`data/found/582_niching_ga/positive_control_3_2_summary.json`.

**Verdict: the positive control PASSES.** The new fitness function + bounds
(not the old `correct_symmetric_fixed_jacobi` recipe) correctly guides the
niching GA into the basin of the known #440 3:2 circular member, and the
existing asymmetric corrector converges it back onto that same member to
well within the stated tolerances -- on the hardest of the 5 tabulated cases.

## Supplementary data point: MMR 5:1 (near-miss, informational only)

Not the official positive control (only one MMR is required per #582's own
gate), but run for extra confidence since each run only takes ~15s. Result:
**near-miss**, `matched=False`. x0/ydot0/T errors are all comfortably inside
tolerance (1.36%, 0.89%, 2.24%), but two things pushed it over the line:

1. `orbit.converged=False` -- the corrector's own residual landed at
   1.869e-10, just above its default `tol=1e-10` (the independent Radau
   closure residual, 1.862e-10, is essentially the same number and is well
   inside `independent_tol=1e-6`). This reads as a threshold-edge miss, not
   a real divergence -- `refine_ga_candidate`'s default `tol` was not tuned
   per-MMR.
2. `jacobi_err=0.0285`, just over the 0.02 abs tolerance.

Not investigated further (out of scope for this dispatch -- one passing MMR
satisfies #582's gate, and #582 explicitly says not to run the broader
sweep). Flagged honestly rather than silently omitted or used to justify
loosening the PASSING 3:2 tolerances after the fact.

## Things worth flagging for whoever runs the full sweep

1. **`require_monotone_decrease=True` is load-bearing, not optional.**
   `correct_general_periodic_3d`'s default (`False`, blind Newton) actively
   DIVERGED on a GA-realistic seed during this build: fed the GA's best 3:2
   genome (already at fitness 0.9935, periodicity defect norm ~0.06 -- a
   genuinely good seed), the blind-Newton path landed at
   `corrector_residual=0.30`, `independent_closure=2.85`, on a completely
   different point (`x0` flipped from +0.74 to -0.85). Switching to
   `require_monotone_decrease=True` (which the corrector's OWN docstring
   recommends for "ill-conditioned long-arc closures (asymmetric
   full-6D-at-T mode)") converged the SAME seed cleanly to
   `residual=1.3e-12`. `refine_ga_candidate()` now defaults to `True`; this
   is pinned by
   `test_refine_ga_candidate_blind_newton_can_diverge_documented_reason` so a
   future change to the default is a conscious decision. Anyone hand-rolling
   a different call to `correct_general_periodic_3d` for this genome family
   should NOT use the library default.

2. **The literature-matcher label reuse is over-inclusive.** Per #582's own
   note, the only spatial-CR3BP anchor at mu=0.001
   (`genome/known_corpus_3d.py`'s Antoniadou & Libert 2019 record) is itself
   pinned under `primary="Earth"`, `body_set={"Moon"}` even though it is
   explicitly annotated mu=0.001 (NOT physical Earth-Moon). Reusing that same
   label to reach it (`build_candidate_signature`'s `primary="Earth"`,
   `sequence=("Moon",)`) works -- but it ALSO pulls in every OTHER anchor
   whose `body_set` is a superset of `{"Moon"}` under `primary="Earth"`,
   including genuinely PHYSICAL Earth-Moon cycler papers (Braik-Ross,
   Kumar-Rawat-Rosengren-Ross, Hiraiwa et al., Sanaga-Park-Howell -- all
   observed engaging in this run's output). None of those anchors' spatial
   topology (`topology_3d`) actually constrains anything here, because the
   Antoniadou-Libert anchor's own `topology_3d=None` means the (k1,k2,k_z)
   filter never activates for THIS anchor pool (see
   `literature_check._spatial_topology_matches`'s own guard: filter applies
   only when BOTH sides carry `topology_3d`). Practically, this means a
   FUTURE genuine novelty search at this mu risks `check_literature` pulling
   real-web-search hits for unrelated physical Earth-Moon papers and
   mis-scoring them against a mu=0.001 candidate's fingerprint. This did not
   block the positive control (no live search was run -- see next point) but
   is a real gap the coordinator should address (e.g. a proper "generic-mu"
   corpus label, or filtering `_candidate_anchors` on an explicit mu-band
   annotation) before trusting any `check_literature` verdict from the full
   sweep.

3. **`check_literature` itself (the live web-search verdict) was NOT
   exercised.** This positive control only proves the STRUCTURAL matcher
   ENGAGES (`literature_anchors_engaged()` returns a non-empty, correct
   anchor list) -- `check_literature()` needs an injected, live
   `SearchFn` (real web search), which is out of scope for an offline
   pytest/script positive control. The full sweep's own novelty gate must
   still run the real `check_literature()` per
   [[feedback_literature_novelty_check_baseline]] before any "novel" claim.

4. **The full-asymmetric corrector is under-determined (6 residuals / 7
   unknowns, min-norm least-squares) by design.** Different seeds in the
   same basin can converge to slightly different members of a local
   1-parameter family -- this is why the 3:2 tolerances (3-5% relative) are
   looser than a bit-exact reproduction. This is expected corrector
   behaviour (documented in `cr3bp_general_periodic_3d.py`'s own docstring),
   not a defect introduced here.

5. **GA budget used here (40 pop x 60 gen) is a small positive-control
   budget, not the eventual sweep's budget.** #581's own paper-scale
   reproduction used Gurfil-Kasdin's Table 1 constants (200 pop, 400 gen).
   The driver's `--mode ga` path defaults to those same paper-scale
   constants (`--population 200 --generations 400`) for the future sweep;
   the positive-control harness (`--positive-control`) uses its own smaller,
   separately-defined budget (`_small_ga_config`) so the control itself stays
   fast (~15s/MMR).

## Reproduction

```
uv run python scripts/run_582_asymmetric_3d_niching_search.py --mmr 3:2 --positive-control --workers 8
uv run pytest tests/search/test_isolated_3d_asymmetric_fitness.py tests/search/test_isolated_3d_asymmetric_pipeline.py -v
```
