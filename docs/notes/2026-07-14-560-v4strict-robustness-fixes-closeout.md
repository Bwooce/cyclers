# #560 close-out: V4-strict Lambert branch-continuity + planet-crossing-guard robustness fixes

Date: 2026-07-14
Scope: `src/cyclerfinder/data/validation/v4_uranus_strict.py` + `tests/data/test_v4_uranus_strict.py`

## TL;DR

**The two #560 robustness fixes are DONE.** Both landed in commit `6c54bba`
("#567 (1)+(2): fix Lambert branch-continuity flip, tag (not exclude)
planet-crossing infeasibility, fix hardcoded audit-field body sampling",
2026-07-11) — #560 was folded into #567's step (1) because #567's
epoch-robustness scan could not be trusted until #560's two artifact
generators were fixed. This note closes #560 out: it records what the fixes
are, pins the two regression cases #560 named, adds the one requirement #567
did not explicitly cover (a durable pin that #312's *own canonical
single-epoch* result is UNCHANGED), and — importantly — documents ONE
deliberate deviation from #560's original wording that was made under a
later, Fable-validated correction and must NOT be "fixed back".

## Fix 1 — Lambert branch-continuity (bug 1)

**Problem (as #560 framed it).** `_cycle_v4_strict`'s rev-1 Lambert branch
selection had no continuity tracking across neighbouring epochs. The old code
picked a branch BEFORE propagating, via a bare `min()` on the
departure-velocity residual against the departing moon's own velocity. That
proxy criterion has no relationship to which branch actually flies well under
the perturbed V4-strict dynamics, so its argmin can flip between two branches
that each vary smoothly with epoch — producing a discontinuous multi-km jump
in the reported terminal miss between adjacent HOURS.

**Fix chosen.** Outcome-based selection in the new `_select_leg_transfer`:
propagate *every* rev-`n_rev` Lambert branch and select the branch with the
smallest ACTUAL terminal offset (the same quantity the pass/fail gate cares
about). Selecting on a continuous outcome makes the selected value itself
continuous — "min of continuous functions" is continuous, at worst kinked at a
true crossing, never a jump. This is a cleaner realisation of #560's fix
direction (b) (min-over-both-branches) than option (a) (previous-epoch
branch memory), and needs no cross-epoch state threaded through the driver.

**Regression pinned.** `test_567_lambert_branch_selection_is_continuous_across_adjacent_hours`
pins the 2000-09-06T03:00:00 → T04:00:00 flip on the #327 SILVER. (#559's
original ephemeral 2000-04-09 13:00/14:00 diagnostic case was not preserved in
the repo; this is a freshly-located, reproducible instance of the SAME
mechanism, confirmed by direct instrumentation before the fix landed — pre-fix
the pair jumped ~24,986 km → ~3,493 km purely from the discrete branch flip;
post-fix both epochs select the outcome-best branch and the offset varies
smoothly, hour-to-hour delta ~600 km.)

## Fix 2 — planet-crossing arc guard (bug 2)

**Problem (as #560 framed it).** A DOP853 integrator failure collapsed into a
plain validation FAIL indistinguishable from a genuine dynamical failure. At
some epochs the real root cause is a non-physical planet-crossing Lambert
"solution" (osculating periapsis inside Uranus's 25,559 km equatorial radius)
that only "fails" because the integrator chokes propagating through the
planet. There was no perijove-vs-body-radius guard at all before propagation.

**Fix chosen.** New `_leg_periapsis_km` helper computes each branch's
osculating periapsis; `_select_leg_transfer` checks it against `r_eq_km`
BEFORE propagation. A planet-crossing branch is tagged
`FAILURE_MODE_PLANET_CROSSING` (with the offending periapsis recorded on the
verdict) and is not propagated. The failure modes are now a 4-way distinction
(`converged` / `lambert_no_solution` / `planet_crossing_infeasible` /
`integrator_failure`) so a genuine geometric infeasibility is never conflated
with an unexplained solver death.

**Regression pinned.** `test_567_planet_crossing_tagged_not_silently_misclassified`
pins 2000-07-24T02:00:00 on the SILVER (leg 1 Oberon→Umbriel, both rev-1
branches have periapsis ~97 km / ~852 km, far inside `r_eq`), asserting the
cycle FAILs with `failure_mode == FAILURE_MODE_PLANET_CROSSING` and a recorded
`perijove_km` inside `r_eq`. (2000-08-17, the second epoch #560 named, is the
same mechanism on the same pair; one representative pin is sufficient.)

### DELIBERATE DEVIATION from #560's original wording — do NOT revert

#560's parking-lot allocation said: "exclude design-infeasible-geometry epochs
from the pass-rate denominator rather than counting them as dynamical FAILs."

**This direction was explicitly overridden by a later, Fable-validated #567
PIN (2026-07-11) and the code intentionally does NOT exclude them.** From
OUTSTANDING.md's #567 entry:

> the guard must NOT neutrally skip/exclude perijove-inside-Uranus epochs.
> #559's own entry documents the trigger as REAL geometry ... only the silent
> stiff-death misclassification is the artifact, not the underlying
> infeasibility. These epochs are genuine dynamical infeasibility and MUST
> count as FAIL against each candidate's PASS-band width ... a guard that
> neutrally skips them would inflate every candidate's apparent band width and
> cause the writeback validation level ... to over-claim.

So the implemented behaviour is **tag-but-count-as-FAIL**, not
exclude-from-denominator. The `#567` follow-up diagnosis
(`docs/notes/2026-07-11-567-post-fix-singleton-flip-diagnosis.md`) then showed
these planet-crossing FAILs are a genuine, fast-oscillating *synodic
feasibility boundary* (FFT-locked to each moon pair's synodic period), read as
a duty cycle — not noise to be excluded. Any future work that re-reads #560's
stale "exclude from denominator" wording in isolation and changes the code to
skip these epochs would re-introduce exactly the band-width over-claim #567
was corrected to prevent.

## #312 canonical single-epoch result: UNCHANGED (verified, not assumed)

#560 requires explicit confirmation (test or direct comparison, not
assumption) that #312's own already-valid canonical single-epoch V4-strict
result is unaffected by the fixes. Confirmed directly: the #327 SILVER
(Umbriel-Oberon-Umbriel, = #312's representative) at the canonical
2000-06-21T00:00:00 reference epoch (#338/#566 anchor) is a clean PASS under
the fixed code — all 3 cycles `converged`, drift-agreement-vs-V3 = 12,159 km,
comfortably under the 50,000 km floor, `bounded_drift_survives = True`. Now
pinned as a durable regression guard by
`test_560_silver_312_canonical_epoch_unchanged`. (The pin uses generous
headroom — `< 20,000 km` — so cross-BLAS/platform sub-km churn won't flake it
while still catching a regression that moved the result toward the floor.)

## Audit-field fix (bug 3, carried in the same #567 commit, noted for completeness)

The same commit also fixed the four audit-only e/i fields that were hardcoded
to always sample Umbriel/Oberon regardless of the candidate's sequence; pinned
by `test_567_audit_fields_track_non_umbriel_oberon_sequence`. Not part of
#560, listed here only because it shares the commit.

## Verification

`tests/data`, `tests/search`, `tests/scripts` all green; `ruff check`,
`ruff format --check`, `mypy src tests` clean. The three relevant SPICE-gated
tests (2 pre-existing #567 regressions + the new #560 canonical pin) require
the URA111 kernel and run when it is present (it is, in this environment);
they skip gracefully otherwise.
