# #437 — Fold-aware pseudo-arclength continuation in eccentricity (ER3BP)

> REQUIRED SUB-SKILL on execution: superpowers:subagent-driven-development.

**Goal:** build a pseudo-arclength continuation of ER3BP symmetric periodic
orbits in the eccentricity parameter `e` that walks THROUGH folds (turning points
in e), replacing the secant predictor that fails at folds. This is the linchpin
for #440 Phase 2 (isolated-family hunt) and the #436 re-run (whose discriminator
failed at exactly these folds).

**Why now / why it's well-specified:** two newly-digested papers confirm the
failure mode and supply same-model goldens — Peng-Bai-Xu 2017 ("we tested to
directly continue... the routine failed" at turning points; fix = pseudo-arclength
+ multiple shooting; Tables 2/3 = Sun-Mercury ME-Halo ICs + monodromy eigenvalues)
and Martínez-Cacho et al. 2025 (explicit fold at **e=0.0324**, Sun-Mars 2:3 swing
QSO; App. B 15-digit QSO ICs incl. the blank "couldn't close" entry).

**Architecture:** mirror the PROVEN pseudo-arclength implementation already in the
repo — `continue_general_3d_family` in `src/cyclerfinder/search/cr3bp_3d_family_tracer.py`
(it walks 3D CR3BP families through folds via SVD-null-tangent prediction + an
augmented residual incl. the arclength constraint). Swap: continuation parameter
= `e` (the pulsating-frame eccentricity), corrector = `correct_er3bp_periodic`
(`src/cyclerfinder/genome/er3bp_periodic.py`). The new routine lives in
`src/cyclerfinder/genome/er3bp_continuation.py` alongside the existing secant
`continue_er3bp_family_in_e`.

**Conventions:** main; `uv run` ruff+mypy; no Co-Authored-By; pathspec; imports at
top; goldens' EXPECTED side sourced (a fold eccentricity / eigenvalue from a
paper), never a value our own code computed; subagents finish through commit.

---

## Task 1: pseudo-arclength continuation-in-e core

**Files:** modify `src/cyclerfinder/genome/er3bp_continuation.py`; test
`tests/genome/test_er3bp_arclength_continuation.py`.

READ `continue_general_3d_family` (cr3bp_3d_family_tracer.py) for the arclength
pattern: at each member, form the Jacobian of the residual w.r.t. the augmented
free-vars `(IC components, e)`, get the null tangent via SVD, step along it by
`ds`, then correct back onto BOTH the periodicity residual AND the arclength
constraint `<state-prev, tangent> = ds`. Read `correct_er3bp_periodic` for the
ER3BP residual (symmetric perpendicular crossing, free_vars=(IDX_X,IDX_YDOT),
residual_indices=(IDX_Y,IDX_XDOT), is_half_period_residual).

- [ ] Implement `continue_er3bp_family_in_e_arclength(sys_base, seed_state, period_f, e_target, *, ds=..., max_steps=..., is_half_period_residual=True, tol=1e-10) -> list[ER3BPPeriodicOrbit]` — pseudo-arclength in `e`, walking forward toward `e_target` and continuing PAST a fold (where `de/ds` changes sign) rather than stalling. Return the family members with their `e`.
- [ ] **Regression golden (self-consistent):** on a SMOOTH family that does not fold below `e_target` (e.g. the Broucke Earth-Moon e=0.0549 seed used in #432), assert the arclength walk reaches `e_target` and its members match the existing secant `continue_er3bp_family_in_e` to within tol (same family, no fold → both agree). This pins correctness without a fold.
- [ ] ruff + mypy; commit `genome/#437: pseudo-arclength continuation in e for ER3BP (fold-capable core)`.

## Task 2: fold golden — reproduce a documented turning point

**Files:** test `tests/genome/test_er3bp_fold_golden.py` (slow-marked).

The sourced fold: **Martínez-Cacho et al. 2025 Sun-Mars 2:3 swing-QSO fold at
e=0.0324** (digest `2026-06-25-digest-planar-retrograde-ERTBP-2025.md`; App. B
gives the QSO seed ICs — extract the Sun-Mars 2:3 steady-QSO IC + the swing-QSO
seed). If the App. B Sun-Mars 2:3 ICs are not cleanly machine-transcribable from
the digest, fall back to the Peng Sun-Mercury ME-Halo branch (Tables 2/3: 8 M5N2
ICs + eigenvalues) and reproduce its fold/branch turn instead — note which golden
was used and why.

- [ ] Seed the documented family, run BOTH continuators toward the fold: assert the
  SECANT method stalls/dies at the fold while the ARCLENGTH method walks through
  it and reconnects past the turning point. The EXPECTED side (fold exists near
  the published e; secant cannot pass it) is sourced to the paper — assert the
  arclength member count past the fold > the secant member count, and that the
  arclength family brackets the published fold eccentricity.
- [ ] If the published fold cannot be cleanly reproduced from the available ICs,
  DO NOT fabricate one: assert instead the weaker but still-sourced property
  (arclength passes a turning point that secant cannot, on our own folding family)
  and log the gap honestly in the verdict.
- [ ] ruff + mypy; commit `genome/#437: fold golden — arclength walks through a documented ER3BP turning point`.

## Task 3: wire into #436 re-run + #440 readiness (controller)

- [ ] Expose the arclength continuator as the default for the no-CR3BP-limit test
  (the #436 `classify_no_cr3bp_limit` reverse-continuation) — its step-count
  flip in #436 was a fold artifact this fixes. Re-run a small #436 EM slice to
  confirm the classification is now step-count-STABLE (the decisive #436
  diagnostic: same orbit no longer flips cr3bp_continuous↔e_only across n_steps).
- [ ] Short verdict note: #437 capability delivered + golden-validated; #436
  discriminator now stable; #440 Phase 2 unblocked.

## Self-review
- Mirrors the proven `continue_general_3d_family` arclength pattern (not new
  numerics from scratch). ✓
- Golden EXPECTED side is a published fold eccentricity, not code-computed. ✓
- Honest fallback if the exact paper ICs aren't transcribable. ✓
- Directly unblocks #436 re-run + #440 Phase 2. ✓
