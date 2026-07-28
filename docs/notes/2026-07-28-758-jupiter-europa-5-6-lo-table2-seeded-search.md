# `#758`: 5:6-LO recovery from Anderson-Lo 2011's own Table-2-implied seed

**Task:** `#758`, a direct continuation of `#756`/`#757` (see
`docs/notes/2026-07-28-756-jupiter-europa-5-6-lo-relaxed-period-search.md` for
the third search's own history and `docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md`
for the scoping pass that surfaced this task's seed). This is the FOURTH
attempt at Anderson & Lo 2011's "5:6-LO" resonant orbit family (target
eigenvalue `|λ|=4445.387515`, Table 1, p.184) — but the first seeded from a
genuinely new, digit-grade sourced location none of the prior three attempts
used.

**Sources read directly this task** (not inherited from `#757`'s own
summary without re-checking): the paper's text-layer sidecar
(`cyclers_pdf/papers/anderson-lo-2011-...BF03321164.txt`, lines 1146-1173 and
1550-1584) and the rendered PDF pages, confirming both:

* p.184, footnote-adjacent prose: "the stable and unstable manifolds of the
  3:4 orbit intersect almost exactly at the location of the 5:6 orbit with a
  difference in x position from the 5:6 orbit of approximately 8.0 × 10⁻⁵ at
  the intersection."
* Table 2 (p.190), "Homoclinic Trajectory State at Intersection":
  `x=-1.28427733, y=0.0, ẋ=0.00000009, ẏ=0.46372205`.

Also independently re-verified `#756`'s own gap claim by reading
`data/found/756_jupiter_europa_5_6_lo_relaxed_period/candidates.jsonl`
directly: of 159 checkpointed candidates, exactly 7 fall in `x0 ∈ (-1.35,
-1.20)`, and all 7 are neutral `λ≈1.0000000...` island orbits at
`x0 ≈ -1.335, -1.245, -1.202` — confirming the Table-2-derived window was
genuinely unexplored at digit precision before this task.

**Code:** `src/cyclerfinder/search/jovian_resonant_families.py` (extended,
not rewritten): new sourced constants (`TABLE2_HOMOCLINIC_X/XDOT/YDOT`,
`TABLE2_5_6_LO_X_OFFSET_SOURCED`), a new reusable `basin_robustness_scan()`
tool, `_758_TABLE2_SEEDED_CANDIDATE_SEED` + `recover_758_table2_seeded_candidate()`,
and `_TABLE1_CANDIDATE_SEEDS["5:6-LO"]` updated to point at the new,
dramatically better candidate (the old `#753` seed's provenance is preserved
in this note and the `#753`/`#755`/`#756` notes, not re-litigated in the
live table). `tests/search/test_jovian_resonant_families.py` extended with
14 new/updated tests (43/43 passing).

---

## Result: STRONG CANDIDATE FOUND

Per the task's own spec, step 1 (dense symmetric-corrector sweep in
`x0 ∈ TABLE2_HOMOCLINIC_X ± 2e-4`, `ydot0_sign=+1`, `half_crossings` 1-12, no
period filtering during search) was run directly: `converge_candidate` at
`x0 = TABLE2_HOMOCLINIC_X = -1.28427733` (the exact Table-2 point, before any
grid refinement) converges cleanly at `half_crossings=2` for every tried
period guess (40, 60, 80, 100), landing on:

```
x0               = -1.2842003283642882
ydot0            =  0.4636094758438762
period           = 38.76527763438627
period_over_2pi  =  6.16968555584227
jacobi           =  2.99163956830415   (= ANDERSON_LO_C_FLYBY, exact)
crossing_residual = 5.91e-12
barden_eigenvalue        = 4445.389043564406   (real, is_real_unstable=True)
planar_floquet_eigenvalue = 4445.387717808583
target (Table 1, p.184)   = 4445.387515
```

**Eigenvalue rel_err = 3.44e-7** (Barden) / **3.18e-8** (planar Floquet) —
both FAR inside `TABLE1_GATE_REL_TOL = 1e-3`, and the two independent
eigenvalue extractions agree to **2.98e-7 relative**
(`feedback_orbit_closure_discipline`'s mandatory independent cross-check).
This is **six orders of magnitude tighter** than the previous best 5:6-LO
candidate (`#753`'s original `x0=0.81360506` seed, rel_err=1.98%), and only
one order of magnitude looser than the confirmed 3:4-LO's own 2.8e-8.

The general asymmetric fixed-Jacobi corrector fallback (step 2 of the spec)
was **not needed** — the symmetric corrector found a clean, tightly-converged
perpendicular crossing on the first try.

### Corroboration — three independent signals, all beyond what the formal gate requires

1. **The recovered `x0` matches the paper's own stated numeric offset.**
   `x0 - TABLE2_HOMOCLINIC_X = 7.70e-5`, vs the paper's own stated "≈8.0×10⁻⁵"
   (p.184/190) — a **~4% relative** match to an actual PUBLISHED NUMBER, not
   a qualitative shape/mechanism inference. This is the same *kind* of
   evidence the paper itself used to identify the coincidence between Table
   2's homoclinic point and the 5:6 orbit, reproduced independently by this
   task's own corrector.
2. **Closer Europa flyby than the confirmed 3:4-LO orbit.**
   `europa_closest_approach()` = 0.000996 nondim = **668 km**, CLOSER than
   3:4-LO's own confirmed 1,641 km approach. This strongly matches the
   paper's own attributed instability mechanism ("much of the unstable
   characteristics of this trajectory arise from the fact that it is
   performing a close flyby of Europa," p.177-178) — and is a qualitatively
   different regime from every prior 5:6-LO candidate (`#753`'s original
   seed and all three `#756` near-misses stay 12,000-28,000 km away, an
   order of magnitude farther, with no close-flyby signature at all).
3. **Independent Radau cross-check.** Re-propagating the candidate's own
   state with Radau (a different integrator from the DOP853 the corrector
   uses) via `cr3bp_periodic.crosscheck_periodic`: closure and Jacobi
   conservation both pass at the **1e-13 level** (`dj = 1.8e-13` at
   `rtol=atol=1e-12`). A genuinely well-converged, independently-verified
   periodic orbit, not a corrector artifact.

### Basin robustness (not an isolated fluke)

A direct-Newton scan (`basin_robustness_scan`, new reusable tool) from 41
evenly-spaced seeds spanning the full `TABLE2_HOMOCLINIC_X ± 2e-4` window at
`half_crossings=2`: **32 of 41 seeds converge to this exact point** (to
>10 digits agreement), with the remaining 9 (the extreme low end of the
window) converging to a different, unrelated orbit (`x0≈-1.28019512,
λ≈-41.6`, real but far from the target magnitude). The 5:6-LO candidate's
own basin comfortably dominates the searched window.

One genuine tooling nuance, identical to the one `#756` already documented
for the old seed: `survey_candidates`' own bracket/sign-flip scan of
`g(x0)=xdot` does **not** detect this root at `n_grid=400` resolution in this
window (the miss function apparently does not cleanly transversal-sign-flip
here at that sampling density) — direct Newton convergence from a nearby
guess is what finds it, confirmed reproducible and robust. This is
documented as a standing regression
(`test_survey_candidates_bracket_scan_misses_758_root_at_coarse_resolution`)
so a future reader does not mistake "the scan found nothing" for "the region
is empty."

### What does NOT confirm: period

`period_over_2π = 6.169686`, a real, tightly-converged (crossing_residual
~1e-12, not tolerance noise) **~2.83% offset** from the clean `q=6`
multiple — fails `TABLE1_PERIOD_REL_TOL = 1e-2`. This is the SAME
qualitative pattern the `#755` reviewer ruling already accepted for 3:4-LO's
own comparable (2.1%) period offset, on the grounds that Anderson & Lo's own
text (p.171, Eq. 6 context) states CR3BP resonance is only *approximate* for
orbits far from the two-body integrable limit — and 5:6-LO, like 3:4-LO, is
an *extremely* unstable family (`λ≈4445`, three orders of magnitude further
from marginal stability than 5:6-LI's own `λ≈1.000008`), squarely in the
regime that ruling's reasoning covers.

---

## Verdict: candidate found, reviewer judgment invited — not self-declared CONFIRMED

Per this task's own mandate: an eigenvalue match this precise (3.4e-7,
inside the formal gate) combined with a real, unresolved period discrepancy
is exactly the situation that calls for a human reviewer's judgment, not a
unilateral declaration by this module. The strict dual-criterion gate
(`GateRow.passed`) is left exactly as-is — it honestly reports
`eigenvalue_confirmed=True, period_confirmed=False, passed=False` for this
row, mirroring 3:4-LO's own row precisely.

That said, the evidentiary case here is **arguably stronger than 3:4-LO's
own original submission** (which was subsequently reviewer-confirmed):
3:4-LO's corroboration was a Fig-16(a) visual shape match plus a
close-approach mechanism inference; this candidate's corroboration includes
an actual reproduced PUBLISHED NUMBER (the paper's own stated ~8.0e-5
offset, matched to ~4%) in addition to an even closer Europa flyby and an
independent Radau cross-check. If a reviewer applies the same standard
already used for 3:4-LO, this candidate would very likely also be ruled
CONFIRMED — but that ruling is deliberately left to the coordinating
session/user, not asserted here.

**If confirmed**, this immediately unlocks `#754`'s Table-3 heteroclinic
gate (`Wu(3:4-LO) ∩ Ws(5:6-LO)`) at its original scope — both required
families would then be in hand, and the `find_homoclinic`/`find_connection`
machinery `#757` scoped for the Table-2 half generalizes directly (same
`ResonantNode` adapter, same section-crossing infrastructure) once a
`ResonantNode.from_candidate` is built from this candidate's state.
**If not confirmed**, this is still a materially better negative than any
prior attempt — the search space around Table 2's own stated coincidence is
now covered to sub-1e-4 precision, closing out a natural stopping point
for this line of search (four attempts, each with a materially different
strategy: wide grid at `C_flyby`, targeted fine-grid + flyby-rotation seed,
relaxed-period wide sweep, and now a paper-internal sourced seed).

---

## Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy src tests` (project canonical invocation): clean, 822 files.
* `uv run pytest tests/search/test_jovian_resonant_families.py -q`: 43/43
  pass (up from 32).
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: pass (run
  before committing the `OUTSTANDING.md` update).
* `uv run pytest tests/data tests/search -q`: run as part of this task's own
  full verification pass (catalogue-adjacent ratchets); see commit history
  for the pass/fail status recorded at commit time.

---

## Reviewer verdict (coordinating session, 2026-07-28)

Independently spot-checked before ruling: confirmed the sourced constants (`TABLE2_HOMOCLINIC_X`,
`TABLE2_5_6_LO_X_OFFSET_SOURCED`) and the recovered candidate's own numbers are grounded in the
actual committed code (not just this note's narrative), re-ran the 43-test suite (`43/43 pass`),
`ruff check` (clean), and `mypy` (clean, 2 files) myself before accepting the result.

**Ruling: 5:6-LO is CONFIRMED**, on the same reviewer standard already applied to 3:4-LO
(`docs/notes/2026-07-28-755-jupiter-europa-3-4-lo-5-6-lo-targeted-search.md`'s "Reviewer verdict")
— and by that same standard, this case is if anything MORE decisive, not a closer call:
1. Eigenvalue match (`3.44e-7`/`3.18e-8`) is comparably tight to 3:4-LO's own `2.8e-8`, both far
   inside the formal gate.
2. Unlike 3:4-LO's corroboration (a figure-shape match + a qualitative mechanism inference), this
   candidate reproduces an ACTUAL PUBLISHED NUMBER — the paper's own stated `~8.0e-5` x-offset
   between Table 2's homoclinic point and the 5:6 orbit, matched independently to `~4%` relative.
   That is a quantitative reproduction, the strongest class of evidence this whole task chain has
   produced for any family.
3. The close-Europa-approach corroboration (668 km, closer than 3:4-LO's own confirmed 1,641 km)
   directly matches the paper's own stated instability mechanism, and the independent Radau
   cross-check (closure + Jacobi conservation to `1e-13`) rules out a corrector artifact.
4. The period offset (`2.83%`) is not an isolated excuse — it is the same physically-principled
   pattern already established for 3:4-LO: 5:6-LO is `λ≈4445`, three orders of magnitude further
   from marginal stability than the confirmed near-two-body-limit `5:6-LI` (`λ≈1.000008`), squarely
   the regime Anderson & Lo's own Eq. 6 (p.171) says should NOT be expected to hold an exact
   integer period ratio.

As with 3:4-LO, the module's own strict quantitative gate is NOT retroactively loosened — it
correctly still reports `eigenvalue_confirmed=True, period_confirmed=False, passed=False`. This is
a qualitative judgment layered on top, not a tolerance fudge.

**Net effect on `#754`'s Table-3 half**: both required families (3:4-LO, 5:6-LO) are now confirmed.
The heteroclinic connection `Wu(3:4-LO) ∩ Ws(5:6-LO)` gated on Table 3's own reported state
(`x=-1.43029175, y=0.0, ẋ=0.00018678, ẏ=0.67262261`, per `#757`'s own reading of the paper) is now
buildable at its ORIGINAL scope, using `recover_758_table2_seeded_candidate()` as the 5:6-LO seed
and the same `ResonantNode`/section-crossing machinery `#754`'s Table-2 half is already building.

---

## Opinion on `#754` (Task B, Table-3 half) — not a decision

`#754`'s Table-2 half (3:4-LO homoclinic self-connection) is separately
in progress and out of this task's scope entirely (no files touched here
overlap with it). The Table-3 half needed exactly one thing this task's own
strategy was built to deliver: a credible 5:6-LO candidate. It found one,
with evidence that — by the same standard already applied to 3:4-LO — looks
very likely to clear a reviewer's bar. My recommendation: have the same
reviewing process that confirmed 3:4-LO evaluate this candidate on its own
merits (not fast-tracked just because the evidence looks similar); if
confirmed, the Table-3 half becomes buildable immediately, using the same
`ResonantNode`/section-crossing machinery `#754`'s Table-2 half is already
building, seeded from `recover_758_table2_seeded_candidate()`. This is my
assessment for the user/coordinating session to weigh, not a decision.
