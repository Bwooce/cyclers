# #793 — catalogue gap execution sprint: honest results

**Date:** 2026-08-08

Dispatched as a census/hygiene sprint against four claimed near-zero-risk wins in
`data/catalogue.yaml`'s `data_gaps`. This note is the honest account: what actually closed,
what didn't, and why — per this project's own discipline that a clean negative is a legitimate
outcome, not a failure to force past.

## Item 1 — multi-rev Lambert gap closure (~415 `#54`-tagged gaps)

**Capability claim verified true.** Read `src/cyclerfinder/core/lambert.py` in full: the
`lambert(..., max_revs=N)` parameter genuinely returns `low`/`high` multi-rev branches (lines
756-782), backed by a derivative-free Illinois-method solver (`_solve_uv_branch`) distinct from
the single-rev Newton path. `tests/core/test_lambert_multirev.py` (7 cases) passes cleanly. The
module's own docstring (lines 6-7) and the `LambertGeometryError` message for 0/180-degree
transfers *did* still describe multi-rev support as "M4's future responsibility" — actively
misleading, since it has existed and been tested since M4 landed. Fixed both (commit
`efbe7239`), including correcting the 180-degree error message's own claim that multi-rev support
would fix that degeneracy (it doesn't — the 180-degree case is an in-plane geometric degeneracy,
not a revolution-count problem).

**Applying the capability to the catalogue turned out to be much narrower than the dispatch
implied.** Of the 415 `#54`-tagged gaps, 213 are `kind: derive` entries for `loop-ee`/`loop-ee-N`
(intermediate Earth-Earth loop) segments. Only 9 rows (14 arcs) have everything a per-arc Lambert
closure needs self-contained (a **sourced** `tof_days` already in the segment, or one computable
from a sourced `period.years` minus sourced adjacent-leg ToFs, plus a sourced target V-infinity) —
the rest additionally need information this simple approach can't supply (an unknown adjacent-leg
ToF, an unresolved multi-loop decomposition, etc.), which is out of `#793`'s scope by design (the
task instructions explicitly held `#54-backfill`-tagged gaps out of scope, and these deeper
dependencies are the same character of gap).

I derived all 14 candidate arcs: circular-coplanar Earth positions separated by the arc's own ToF,
fed to `lambert(r1, r2, tof, max_revs=6)`, selecting whichever branch's emergent V-infinity (at
both the departure and arrival Earth crossing) best matches the row's sourced target, computing
`(a, e)` from the winning branch's `v1` via vis-viva + angular momentum.

**Result: 1 of 14 closed cleanly.**

- `russell-ocampo-3.1.2+1`'s `loop-ee-2` (1083 d, sourced from Russell Table 3.6 cumulative
  epochs): `n_revs=1`, `branch="high"`, `a_au=1.3093`, `e=0.2371`. Emergent V-infinity 3.397 km/s
  vs. the published 3.4 km/s target — 0.09% error, decisively better than every other branch
  tried (next candidate off by an additional 0.017 km/s). Independently cross-checked against
  `lamberthub`'s `izzo2015`/`gooding1990` via `lambert_crosscheck` — machine-precision agreement
  (max component disagreement 2.3e-10 m/s), confirming the Lambert solve itself, not just the
  V-infinity match, is correct. Perihelion works out to 0.999 AU, i.e. the arc touches Earth
  almost exactly at perihelion — the physically expected shape for an Earth-departure arc.
  Written back with a full derivation citation; the row's `data_gaps` entry is narrowed (not
  removed — `ret-me` and `loop-ee-1` on the same row are still open) to record exactly this.

- **13 of 14 did not close.** The dominant reason is a genuine geometric fact, not a tooling
  gap: these `loop-ee` arcs are, by design, "V-infinity-leveraging" resonant loops with ToFs
  deliberately close to half/full/1.5 Earth-year multiples (182 d, 183 d, 365 d, 366 d, 548 d —
  and the `period.years`-derived remainder-tof cases land at 155-290 degrees, still often
  poorly-conditioned). A circular-coplanar Earth-Earth transfer at a near-half/full-year ToF sits
  at or right next to the exact 0-degree/180-degree Lambert transfer-angle singularity
  `LambertGeometryError` exists to flag — the boundary-value problem becomes numerically
  degenerate (every branch converges to the trivial `a=1 AU, e=0, V-infinity=0` non-answer) well
  before it hits the literal 1e-12-radian threshold that raises an exception. This is not a code
  bug: it is the same degeneracy the project's own `core/lambert.py` docstring already documents,
  just showing up empirically across this specific dataset because these arcs are resonant by
  construction.

  I found a candidate way past this: 14 catalogue rows carry Russell's own `free_return_arcs[]`
  arc-type descriptors (`g`/`G` generic, `f`/`F` full-rev resonant, `h`/`H` half-rev — parsed by
  `search/descriptor.py`, itself fully tested). Full-rev (`f`/`F`) arcs are closed-form solvable
  without Lambert at all: the resonance pins `a` via Kepler's third law regardless of `e`, then
  `e` follows from a 1-D vis-viva solve against the target V-infinity — no 0/180-degree issue,
  since no Lambert boundary-value problem is involved. But mapping Russell's descriptor list
  (which enumerates *candidate* arc realizations, not necessarily one-to-one with this
  catalogue's materialized `loop-ee-N` segments) onto specific segments is **not a solved problem
  even in this project's own code**: `search/cycler_assembly.py::descriptor_to_phsi`'s own
  docstring states "There is NO published crosswalk between McConaghy's per-arc descriptor... and
  Russell's p.h.s.i structure, so this is a best-effort STRUCTURAL map" with several fields
  flagged `**APPROXIMATION**`. Forcing a segment-level writeback through an already-flagged
  approximation would be exactly the kind of unforced "it closed!" this project's orbit-closure
  discipline warns against. Registered as follow-on task `#794` (not dispatched) rather than
  guessed through.

  The `ret-me` (Mars-Earth) arcs in the same rows were not attempted at all: unlike `loop-ee`
  (same body at both ends, so the transfer angle is pinned by the sourced ToF and Earth's own
  circular motion alone), `ret-me`'s transfer angle depends on the relative Mars-Earth orbital
  phase at that specific point in the cycle — information this per-arc-isolated method has no way
  to supply. A genuinely different (and larger) piece of work; also folded into `#794`.

**Honest fraction:** 1 gap-arc closed outright (with an added narrative correction on a second
row, item 3 below); ~13 explored and left open for defensible, documented reasons; ~400 more
`#54`-tagged gaps not touched at all this pass (most share the `loop-ee` near-resonance problem,
or the deeper `ret-me`/`#54-backfill` dependency chain).

## Item 2 — Russell backfill (`docs/notes/2026-07-15-596-russell-backfill-method-validated.md`)

**Already done — no execution needed.** The dispatch brief said the validated method "was never
actually run for real"; that's stale. `git log` shows it executed twice already:

- `3a439ac3` (2026-07-15, task `#596`): 161 of 197 matched Russell Table 3.4/3.9-3.11 rows
  backfilled (`trajectory.segments[out-em].a_au`/`.e` via the AR+V-infinity_Earth inversion).
- `4a46c6cd` (2026-07-16, task `#616`): the declared remainder — 3 sign-transcription renames
  (not new rows) and a confirmed genuine negative on the 36 `AR<1.0` rows (Russell's own
  footnote: these never reach Mars in the circular-coplanar model at all — a model-boundary fact,
  not a missing derivation).

Re-verified live rather than trusting the log: `out-em.a_au` is populated on 162/200
`russell-ocampo-*` rows in the current `data/catalogue.yaml` (the extra 1 over 161 is the
pre-existing independently-sourced Aldrin row). Nothing to execute; nothing was reverted.

The phrase in the dispatch brief about "sibling rows pending rev-count resolution" traces to
`scripts/backfill_russell_2004_tables.py::apply_writeback`'s own gap-narrowing note: the
remaining `ret-me`/`loop-ee-*` segments on these same rows "still need the multi-rev Lambert
rev-count resolution (task `#54`) before their own `a_au`/`e` can be safely assigned." That
boundary **is** item 1's territory, not a separate item-2 task — there is no independent item-2
work left to do beyond what item 1 already covers.

## Item 3 — corrected S1L1 topology adoption

Adopted task `#167`'s confirmed topology into `s1l1-2syn-em-cpom`, with one important scope
correction: `#167`'s actual DE440 reconstruction (`search/s1l1_corrected.py`) was performed
against a *different* catalogue row, `russell-ch4-4.991gG2`, which carries its own sourced
V-infinity anchors (4.99/5.10 km/s) — not this row's 5.65/3.05 km/s pair. Per
`[[project_s1l1_nomenclature]]`'s documented history of S1L1 nomenclature confusion, I did **not**
import `#167`'s specific numeric per-leg values (transit times, DE440 epochs, miss distances)
into `s1l1-2syn-em-cpom` — those belong to the sibling row and don't transfer.

What *did* need fixing, and was fixed: `s1l1-2syn-em-cpom`'s `ret-me` segment note, the
`trajectory.segments[ret-me]` data_gap, and the `trajectory.segments` (whole-row) data_gap all
asserted a *specific factual error* — that no direct Mars-to-Earth return leg exists in the S1L1
architecture "by construction," and that a crewed return would require a mirrored companion
cycler (`L1S1`). Task `#167` independently confirmed on real DE440 ephemeris (7 Mars encounters,
sourced Russell App-C #83 per-leg data) that this is wrong: the actual per-cycle sequence is
`E -> g(Earth-Earth free return, NO Mars) -> E flyby -> G(Earth-Mars-Earth transit) -> E`, where
the `ret-me` segment is the real Mars-to-Earth return half of a single continuous `G` arc — no
mirrored companion needed. The existing 3-segment skeleton (`out-em`/`ret-me`/`loop-ee`) was
already the right *shape*; only the claim about `ret-me`'s existence was wrong. Corrected all
three notes in place (superseding, not deleting, the prior wrong text, with a citation to `#167`
and `search/s1l1_corrected.py`), added a `notes:` paragraph explaining the correction and the
row-vs-row distinction, and left this row's own per-leg `(a,e)`/`n_revs`/`branch` derivation open
(still genuinely unresolved — task `#54`, same as item 1's territory).

No `validation_level` change, no numeric writeback on this row: the correction is entirely
structural/provenance, not a new derived value.

## Item 4 — `data/MISSING_DATA.md` refresh

Added a third dated staleness-correction block (the file already carries two, from `#595`'s
2026-07-15 pass and the 2026-07-19 pass) rather than rewriting the whole file, matching its
existing convention. Updated:

- Live counts: 383 entries, 881 `data_gaps` across 313 entries (up from the 834/291 the last
  correction recorded — organic growth from unrelated tasks between 2026-07-19 and today, not a
  regression of the Russell backfill, which was re-verified still live).
- Corrected the "multi-rev Lambert solver" framing in §5: it is NOT pending (implemented, tested,
  and its docstring now says so), but "computable by the Lambert solver" oversold how tractable
  most of the 213 `loop-ee` `kind: derive` gaps actually are — added the same near-resonance
  degeneracy finding from item 1, with a pointer to the new `#794` follow-on.

## Verification

- `uv run python -c "import yaml; yaml.safe_load(open('data/catalogue.yaml'))"` — clean.
- `uv run pytest tests/core/test_lambert.py tests/core/test_lambert_multirev.py -q` — 23 passed
  (lambert.py docstring/message fix, no behavior change).
- Full `uv run pytest tests/ -q` + `ruff check .`/`ruff format --check .` + `uv run mypy src
  tests` — see the commit message / final report for the actual run results (executed after a
  concurrent heavy-load process on this machine cleared, per this project's CPU-contention
  discipline).

## Follow-on task registered

- **`#794`** — resolve the `free_return_arcs[]`-descriptor-to-`loop-ee-N`-segment mapping
  (currently an acknowledged, unvalidated approximation even in this project's own
  `search/cycler_assembly.py::descriptor_to_phsi`) with primary-source verification per row, then
  use the closed-form resonant (`f`/`F`/`h`/`H`) vis-viva solve this note describes to close the
  remaining `loop-ee`/`loop-ee-N` `#54` gaps on the 14 rows that carry `free_return_arcs[]`.
  Separately, the `ret-me` arcs across all `#54`-tagged rows need a genuinely different
  (currently unbuilt) treatment for the unknown relative Mars-Earth orbital phase — worth scoping
  as its own piece of `#794` or a sibling task once the `loop-ee` mapping is resolved.
