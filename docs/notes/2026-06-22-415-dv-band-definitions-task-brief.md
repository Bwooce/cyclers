# Task #415 — Sourced ΔV-band definitions for the cycler taxonomy

**Status:** task brief (for the other agent). Created 2026-06-22 out of the #388
family-pinned-homotopy discussion: the catalogue distinguishes
`trajectory_regime: ballistic | powered | low-thrust` and a "near-ballistic"
net, but the **quantitative boundaries between the bands are not sourced** — they
rest on a geometric proxy (Russell's TR/AR) plus an un-anchored ΔV spread
(observed maintenance ΔV ranges 0.001–3.0 km/s, with real-ephemeris GMAT
closures from 0.175 to 7.29 km/s). This task pins the bands to the literature.

## Goal

Produce a **sourced definition of the cycler ΔV bands** by reviewing the ENTIRE
existing digested corpus, cross-referencing every paper's statements about
ballistic / near-ballistic / powered / low-thrust and any ΔV thresholds, and
flagging **alignment vs conflict** between sources. Output a synthesis note that
(a) defines each band with citations, (b) proposes quantitative thresholds where
the literature supports them, (c) recommends how to encode them in the schema /
README, and (d) flags any current catalogue rows that look mis-binned under the
sourced definition. This is **literature synthesis + a proposal**, not a schema
rewrite (the catalogue/schema area is contended — propose, don't apply; see
Concurrency).

## Why now

The #388 verdict (real-ephemeris family-pinned closer) shows the published
"ballistic" SnLm cyclers do **not** close at zero deterministic ΔV in DE440 —
but "ballistic" in those papers is a property of the *idealized circular-coplanar*
model (the geometric TR ≥ 1 turn criterion), not a verified real-ephemeris
zero-ΔV claim. Whether `mcconaghy-2006-em-k2` and the `russell-ch4-*` rows are
legitimate cyclers in the real world hinges on whether their (measured)
real-ephemeris maintenance ΔV falls inside a *sourced* "still a cycler" band.
Right now we cannot answer that because the band has no sourced ceiling.

## The four bands to define (current working understanding — refine from sources)

1. **Ballistic** — zero deterministic ΔV; closure achieved purely by gravity
   (flyby turns + coasting). Geometric criterion in the idealized model:
   required heliocentric turn ≤ max unpowered flyby turn (Russell TR ≥ 1;
   McConaghy "required turn < max possible turn"). Only stochastic nav TCMs.
2. **Near-ballistic / station-keeping** — small *bounded deterministic*
   maintenance ΔV per synodic period to hold the cycle (idealized→real gap +
   perturbations). **The band whose ceiling is unsourced.**
3. **Powered (DSM)** — deterministic deep-space / powered-flyby maneuvers to
   close a geometrically-infeasible cycle (e.g. S1L1's ~84.7° Earth turn vs
   ~72° max → powered flyby).
4. **Low-thrust / SEP** — continuous propulsion; **out of the v1 paradigm**
   (already excluded; confirm the boundary, don't re-litigate scope).

## What to review (entire corpus — start at the index)

- **Entry point:** `docs/notes/CORPUS_INDEX.md` (the master corpus list).
- **Highest-yield digests/notes for ΔV-band statements** (not exhaustive — scan
  the full index for any others):
  - `2026-06-17-digest-russell-ocampo-2003.md` (TR ≥ 1 "strictly ballistic";
    116 ballistic / 62 close / near-ballistic net; 200 km Earth flyby; "require
    no powered maneuvers").
  - `2026-06-17-digest-byrnes-longuski-aldrin-1993.md` + `2026-06-10-genova-aldrin-2015-mining.md`
    (the canonical Aldrin ballistic cycler + any station-keeping budget).
  - `2026-06-17-digest-mcconaghy-2002.md`, `2026-06-10-mcconaghy-2004-dissertation-mining.md`,
    `2026-06-10-mcconaghy-table71-reproduction.md` (ballistic vs powered "broad
    class"; required-turn vs max-turn footnote; any per-cycle ΔV tables).
  - `2026-06-07-russell-2004-member-tables-transcription.md` ("Total delta v
    over NN years (km/s)" column — the per-row ΔV the tables DO carry).
  - `2026-06-07-continuation-driver-results.md`, `2026-06-08-gmat-v4-execution-results.md`
    (real-ephemeris maintenance ΔV: 0 m/s claims, 0.175 km/s, 7.29 km/s — our
    computed real-world numbers, useful as a reality check on the bands).
  - `2026-06-13-tito-maccallum-2018-free-return-reproduction.md`,
    `2026-06-07-okutsu-tito-free-returns-mining.md` (Tito ballistic-continuity
    ΔV < 0.05 km/s — a tight near-ballistic exemplar).
  - Any low-thrust references (`digest-taheri-junkins-*`, the S3L1 SEP note in
    `catalogue.yaml`) to fix the band-4 boundary.
- Also read `data/README.md` §`trajectory_regime:` (current definitions) and the
  catalogue's `maintenance_dv_kms_per_synodic` / `delta_v_kms` field docs so the
  proposal is consistent with existing schema semantics.

## Deliverable

A note `docs/notes/2026-06-22-dv-band-definitions.md` containing:
1. **Cross-reference table** — one row per source: paper → exact statement(s) it
   makes about band membership / ΔV thresholds / station-keeping budgets, with
   page/section citations. Quote the source; do not paraphrase a number into a
   threshold it didn't state.
2. **Alignment vs conflict analysis** — where do sources agree (e.g. the
   geometric ballistic criterion), and where do they **conflict or use
   incompatible definitions** (e.g. different "near-ballistic" cutoffs, ΔV
   measured per-synodic vs total-over-mission, idealized-model ballistic vs
   real-ephemeris ballistic). Surface conflicts explicitly — do NOT smooth them
   over.
3. **Proposed band definitions** — each band with: the qualitative definition,
   the quantitative threshold(s) **where the literature supports one** (cited),
   and an explicit "UNSOURCED / project-convention" tag where it doesn't (so the
   gap is honest, per `feedback_published_rounded_values_are_display` /
   `feedback_golden_tests_sourced_only`).
4. **Encoding recommendation** — how `trajectory_regime`, the `near-ballistic`
   net, and the V-tiers (esp. V2-powered) should reference these bands; plus a
   flag of any current catalogue rows that appear mis-binned under the sourced
   definition (list them; do NOT edit the rows — that's a separate, coordinated
   catalogue task).

## Honesty gates

- Every threshold/claim traces to a cited source; conflicts surfaced, not hidden.
- No invented numbers presented as sourced. A band with no literature threshold
  is labelled project-convention, not dressed up as published.
- This is review + proposal only. No `catalogue.yaml` / `validate.py` / schema
  edits in this task (those ripple into census ratchets and the area is actively
  edited by another agent).

## Concurrency rules (binding — shared tree, my #388 batch live + active agents)

- Create ONLY the new note `docs/notes/2026-06-22-dv-band-definitions.md`
  (read-only on everything else). Do **not** touch `data/catalogue.yaml`,
  `data/catalogue.schema.json`, `validate.py`, the census ratchets, any
  `nbody/shooter*` / `search/family_pinned_shoot*` / `search/narc_continuation*` /
  `search/generic_return*` (my live batch), or `genome/qp_tori*` / `genome/er3bp*`.
- **Commit pathspec-scoped only** (`git add docs/notes/2026-06-22-dv-band-definitions.md`),
  never `git add -A` / `git add .`.
- **⛔ Do NOT run `git clean -fd`, `git stash`, or `git checkout .` / `git reset
  --hard` in the shared working tree** — these destroy other agents' untracked
  artifacts (a `git clean -fd` already wiped a running batch's runlog this
  session). Remove only your own files by name if you must.

## References
- `docs/notes/CORPUS_INDEX.md`, `data/README.md` (§trajectory_regime),
  `data/catalogue.yaml` (maintenance_dv / delta_v field docs).
- Memory: `feedback_golden_tests_sourced_only`,
  `feedback_published_rounded_values_are_display`,
  `feedback_corpus_document_policy`, `project_validation_ceiling`.
