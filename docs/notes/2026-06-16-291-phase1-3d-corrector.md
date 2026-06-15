# Track-A 3D capability — Phase 1: broken-plane corrector (#291)

**Date:** 2026-06-16
**Frame:** Phase 1 of the multi-week Track-A 3D build (#291; rank-1 frontier
from #286 scoping). Picks up the #287 spike's confirmation that the planar
restriction in our cycler search is a corrector-convention choice (not a
propagator limitation) and ships it as the production corrector that future
3D campaigns can call.

**Status: SHIPPED.** New module + tests + this doc. No catalogue writeback,
no novelty claims, no widening of `KNOWN_CORPUS` — those are Phase 2-N.

---

## What Phase 1 delivers

1. **`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`**
   — `correct_general_periodic_3d(system, state0_guess, period_guess, *, free_vars, residual_indices, is_half_period_residual, ...)`
   — Full-3D single-shooting Newton on the existing 6D CR3BP propagator + STM
     (no propagator changes; the spike confirmed the propagator was already
     3D-ready).
   — Two configurable masks: `free_vars` selects which of
     `{x0, y0, z0, ẋ0, ẏ0, ż0, T}` are Newton unknowns; `residual_indices`
     selects which of `{x, y, z, ẋ, ẏ, ż}` must match at the return time.
   — Two canonical bundles exposed:
     `FREE_VARS_SYMMETRIC_TULIP = (z0, ẏ0, T)` paired with
     `RESIDUAL_PERPENDICULAR_HALF_PERIOD = (y, ẋ, ż)` at `T/2` — the
     tulip/NRHO mode used by the #287 spike;
     `FREE_VARS_FULL_ASYMMETRIC = (all 6 IC components + T)` paired with
     `RESIDUAL_FULL_STATE_AT_T = (full 6D state)` at `T` — the new
     broken-plane asymmetric mode (6 equations in 7 unknowns; min-norm
     least-squares step).
   — Returns `Periodic3DOrbit` with both the corrector residual AND an
     **independent** Radau re-propagation closure at `rtol=atol=1e-12`. Per
     `feedback_orbit_closure_discipline`, a single-integrator residual is
     never trusted alone.
   — `degenerate_planar` flag set when the corrected IC's `|z0| < eps` AND
     `|ż0| < eps`: the planar manifold is dynamically invariant, small `z0`
     seeds correctly collapse onto it, callers searching for genuinely 3D
     orbits screen on this flag.

2. **`tests/search/test_cr3bp_3d_corrector.py` (10 tests, all green)**

   The sourced golden is the catalogue's **planar** Braik-Ross C11a state
   (`data/catalogue.yaml` row `braik-ross-c11a-cycler-2026`, traced to
   Braik-Ross 2026 Table 2, reproduced in-repo at +0.0011%). Tests assert on
   **topology + closure** for the 3D extension, not on a number our own code
   produced — per `feedback_golden_tests_sourced_only`.

   * `test_planar_seed_recovers_braik_ross_c11a` — reproduce-before-trust:
     planar seed with `z0_guess=0` recovers the catalogue state, sets
     `degenerate_planar=True`, independent closure < 1e-6.
   * `test_z0_guess_p05_finds_nontrivial_3d_orbit` — the spike's payload at
     `z0_guess=0.05` lands a 3D orbit with `|z0| > 0.05`, independent
     closure < 1e-8, period in plausible range (30-60 d), Jacobi shifted
     from the planar seed (this is a different family member, not a
     rotation of the planar one).
   * `test_small_z0_collapses_to_planar` (parametrized over
     `z0_guess ∈ {1e-4, 1e-3, 5e-3}`) — confirms the spike's observation
     that small `z0_guess` collapses to the planar manifold;
     `degenerate_planar=True` is the correct return.
   * `test_stm_z_coupling_is_exercised` — sanity check that `uzz/uxz/uyz`
     in `cr3bp_stm_eom` are integrated; a perturbation in `z0` must induce
     non-zero changes in BOTH `x` and `ẏ` at finite time.
   * `test_bogus_seed_does_not_converge_silently` — a wildly off-family
     seed fails either the corrector or the independent closure gate; the
     two-gate convergence rule catches it cleanly.
   * `test_full_asymmetric_mode_reproduces_planar_member` — the new
     7-unknown asymmetric API converges on the planar IC (it IS a valid
     periodic orbit of the full CR3BP) with `degenerate_planar=True`.
   * `test_input_validation` — malformed inputs (wrong shape, empty masks,
     out-of-range indices, non-positive period) raise `ValueError`.
   * `test_index_aliases` — `IDX_*` constants match standard 6D ordering.

3. **This doc** (`docs/notes/2026-06-16-291-phase1-3d-corrector.md`).

---

## What Phase 1 does NOT deliver

Per the #286 scoping the Track-A 3D build is ~850 LOC broken across five
buckets. Phase 1 closed bucket (a)'s prerequisite (the corrector); the
others remain Phase 2-N:

* **3D family-tracer integration** — no widening of
  `continue_nrho_family` (or a new `continue_3d_family`) for asymmetric
  broken-plane continuation in a non-`x0` parameter. The #287 spike
  walked z0 by re-calling `correct_symmetric_nrho` at each step; a
  production tracer would carry tangent prediction, fold detection, and a
  step-size controller. Phase 2 candidate.
* **Bifurcation classification UI** — no genuine 3D bifurcation
  classification (a 3D bifurcation may exchange real/imaginary
  multipliers in patterns the existing 2D bracket logic doesn't taxonomise).
* **`CandidateSignature` / `KNOWN_CORPUS` widening** — no new fingerprints
  for 3D orbits. The literature-check chain still treats z0=0 as the
  canonical IC. Anti-novelty claim discipline: until the literature check
  passes for a 3D candidate, NOTHING is called "novel" (per
  `feedback_literature_novelty_check_baseline`).
* **V-pipeline updates** — no 3D acceptance criteria for V1/V2/V3.
* **Sourced same-model 3D goldens** — none added. The 3D Antoniadou-Voyatzis
  2018 corpus (`KNOWN_CORPUS` gained that anchor at commit `568d8a4`) is the
  obvious source for the first 3D golden, but extracting a per-orbit IC from
  that paper is a separate data-acquisition task.

---

## Path of least resistance for Phase 2

**Concrete first IC:** the spike's converged 3D member at

```
x0     = -0.8116406668238195
z0     = -0.24081020834770114
ydot0  = -0.10629710963669947
T_TU   =  10.20430197041440
jacobi =  3.0270396654490
```

(planar `y0 = ẋ0 = ż0 = 0`; symmetric tulip closure).

**Recommended Phase 2 first task:** add `continue_general_3d_family` in
`nrho_continuation.py` (or a new `cr3bp_3d_continuation.py`) that walks a
continuation parameter — `z0`, `jacobi`, or `T` are the natural choices —
with secant prediction and the new `correct_general_periodic_3d` corrector
as the per-step lander. Expected outcome: the spike's 80-member family
reproduces in a single function call (`continue_general_3d_family(seed=..., parameter='z0', step=-2e-3, n_steps=80)`),
with monodromy + Floquet logging at each step so the bifurcation detector
sees 3D bifurcation brackets.

**Expected blockers:**

1. The spike's continuation in `z0` mostly stayed at `z0 ≈ -0.241` (the
   corrector tagged the same family member from many `z0_target` seeds —
   suggests the natural-parameter step in `z0` is FOLDED in this region).
   Phase 2 may need pseudo-arclength continuation (existing pattern in
   `cr3bp_jacobi_arclength.py`) or a different continuation parameter.
2. The 80-member family contains members at multiple Jacobi values
   (`C ∈ [2.78, 3.03]`), suggesting the family curve crosses Jacobi-isolines
   — `jacobi` is NOT a monotone continuation parameter. `z0` looks monotone
   in the bulk but folds near the edges. `T` may be the cleanest parameter.
3. The catalogue ratchet: Phase 2 outputs must NOT be written to
   `catalogue.yaml` until literature-check passes (the Antoniadou-Voyatzis
   2018 corpus is the gate).

---

## Honest limitations from the spike

* **Small-`z0` collapse to planar:** at `z0_guess < 0.05` the corrector
  finds the planar member (`z0 ≈ 0` returned). The 3D test asserts this
  behaviour at `z0_guess ∈ {1e-4, 1e-3, 5e-3}` so a regression would be
  caught. Mathematically correct (planar manifold is invariant), but
  callers searching for genuinely 3D orbits must use a `z0_guess` outside
  the planar basin; the spike shows `|z0_guess| ≥ 0.05` is sufficient for
  the C11a region. A more principled seed-generation strategy
  (e.g. perturbations along the planar member's center-stable eigenvector)
  is a Phase 2-N concern.

* **Family extent:** the 80-member family terminated at `z0 ≈ -0.452`
  with a corrector failure. Whether that's a true fold, a graze
  approaching either primary's r=0, or just step-size limits is unknown
  from the spike. Phase 2's family-tracer must distinguish.

* **Jacobian conditioning:** the symmetric-half-period 3-residual is
  well-conditioned (the spike converged in 13 iterations with no
  line-search); the full-asymmetric 6-residual-at-T mode is more
  ill-conditioned (caught Phase 1 in the asymmetric test — required
  `tol=1e-10` rather than `1e-11`, matching the existing `correct_periodic`
  behaviour). Phase 2 may need multi-shooting for deep low-perilune
  members or extreme broken-plane orbits.

* **Independent closure threshold:** Phase 1 uses `independent_tol = 1e-6`
  (nondim, ~384 m in EM units). This is loose enough to admit moderately
  unstable orbits whose chaos amplifies round-trip noise (e.g. the C11a
  member has Floquet `λ_max ~ 2.6e4`); a tighter threshold would
  false-reject stable members rendered noisy by integrator floor. The
  threshold is a parameter so downstream callers can tighten it for
  near-stable regimes.

---

## Structural surprises (vs. #286 LOC estimate)

The #286 scoping put Track-A 3D at ~850 LOC. Phase 1 added ~330 LOC of
module + ~290 LOC of tests + this doc, so we're tracking 620 LOC for the
corrector + tests. The remaining bucket weights (family-tracer,
bifurcation UI, `KNOWN_CORPUS` widening, V-pipeline, goldens) look
unchanged from the scoping estimate. **No change to the #286 risk
class** (LOW remains LOW — the corrector landing cleanly was the
single-biggest risk reduction, and the #287 spike already de-risked it).

One pleasant surprise: the configurable-free-vars / configurable-residuals
API design lets a single function handle BOTH the symmetric half-period
closure AND the full asymmetric closure. Phase 2's `continue_3d_family`
can be the consumer of a single function rather than branching on orbit
type. This may shave 50-100 LOC off bucket (a).

---

## Files

* New: `src/cyclerfinder/search/cr3bp_general_periodic_3d.py`
* New: `tests/search/test_cr3bp_3d_corrector.py`
* New: `docs/notes/2026-06-16-291-phase1-3d-corrector.md` (this doc)

**Unchanged** (concurrent-agent isolation): `cr3bp_general_periodic.py`
(planar asymmetric, kept backward-compat for #284); `cr3bp.py` /
`cr3bp_stm_eom` (already correct per #287's Phase A verification);
`nrho_continuation.py` (the symmetric corrector is the spike's de-risked
substrate, kept unchanged).

## References

* #286 — frontier-scoping (`docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md`).
* #287 — spike verdict (`docs/notes/2026-06-16-287-3d-aldrin-scoping-spike.md`) + data (`data/spike_287.jsonl`).
* `feedback_golden_tests_sourced_only` — test discipline.
* `feedback_orbit_closure_discipline` — independent cross-check mandatory.
* `feedback_literature_novelty_check_baseline` — novelty claims gated on
  the literature check (no 3D claim made here).
