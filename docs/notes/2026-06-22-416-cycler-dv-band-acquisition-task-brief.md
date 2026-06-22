# Task #416 — Acquire + digest cycler-ΔV-band primary references (multi-source the #415 bands)

**Status:** task brief. Created 2026-06-22 out of #415: the cycler ΔV-band
*structure* (two axes; geometric AR/TR criterion) is well multi-sourced, but the
**Axis-B real-ephemeris ΔV magnitude thresholds** (Russell 2004's < 1 / < 10 /
< 300 m/s per-7-cycle) rest on **one primary source**. This task acquires +
digests additional primaries to give Axis B **≥ 3 independent anchors** and
either corroborate or correct those cutoffs.

**Access note (gating):** this is an **acquisition** task — it needs the papers,
which are paywalled / not yet in the corpus. Run it only with paper access (the
private `Bwooce/cyclers_pdf` repo at `/home/bruce/dev/cyclers_pdf/papers/`, or
institutional / open-access copies). A 2026-06-22 web pass confirmed the leads
below exist but could not retrieve them (AIAA returned 403). If a target can't be
obtained, record it as an acquisition gap — do not fabricate its numbers.

## Acquisition targets (priority order)

1. **Genova-Aldrin, "Cycler Orbits and Solar System Pony Express," J. Spacecraft
   & Rockets (2021), DOI 10.2514/1.A35091** — reports impulsive *and* SEP
   cycler-maintenance ΔV sequences (the corpus has only the 2015 conference
   mining note; this JSR version is the citable one). HIGHEST priority — it is the
   most likely second independent Axis-B anchor.
2. **The primary behind the "~20–100 m/s/yr station-keeping" figure** — the web
   pass surfaced this magnitude for Mars cyclers via a secondary aggregator; trace
   it to a citable cycler-review / station-keeping paper (candidates: a Genova or
   McConaghy/Longuski review; the Conte-Spencer 2018 cycler survey already in the
   corpus may cite it — check `docs/notes/2026-06-17-digest-conte-spencer-2018.md`
   first, it may already contain a usable figure with a primary citation).
3. **Rauwolf-Friedlander-Nock 2002, AIAA 2002-5046** — the 1500 m/s/15-yr SEP
   maintenance figure (flagged uncaptured in #415; #384-adjacent). Lower priority
   (it's SEP/low-thrust maintenance, band 4) but it pins that boundary.
4. (Opportunistic) any **Landau-Longuski** cycler-ΔV paper or **Chen et al.**
   cycler station-keeping reference encountered while acquiring the above.

## Deliverable

1. A digest per acquired paper in `docs/notes/` (per the Corpus Document Policy:
   OCR if image-only → digest → register in `CORPUS_INDEX.md`), capturing every
   maintenance / station-keeping / DSM ΔV figure **with its basis** (per-cycle /
   per-synodic / per-year / total-over-N) and the paper's ballistic-vs-powered
   definition.
2. **Fold the sourced numbers into `docs/notes/2026-06-22-dv-band-definitions.md`**
   — update the "Multi-sourcing status" section: move corroborated thresholds from
   "provisional on acquisition" to **[sourced ×N]**, and **explicitly flag any
   source that CONFLICTS with Russell's < 1/< 10/< 300 m/s tiers** (different
   cutoff, different basis, different ballistic definition) rather than smoothing
   it. If the new sources disagree, the band note must say so and propose how to
   reconcile (or leave the cutoff explicitly contested).
3. Update `CORPUS_INDEX.md` with the new entries.

## Honesty gates

- Every ΔV figure traces to a cited paper + page; basis recorded; conflicts
  surfaced (`feedback_golden_tests_sourced_only`,
  `feedback_published_rounded_values_are_display`).
- An unobtainable target is an acquisition gap, recorded — never a fabricated
  number. Secondary aggregators (Wikipedia/Grokipedia) are leads, not citations.
- This task does NOT edit `data/catalogue.yaml` / schema / `validate.py` — it
  multi-sources the band note; a separate coordinated task encodes the bands.

## Concurrency rules (binding — shared tree)

- Create ONLY new `docs/notes/` digests + edit
  `docs/notes/2026-06-22-dv-band-definitions.md` + `docs/notes/CORPUS_INDEX.md`.
  Do not touch any source/test/catalogue files.
- **Commit pathspec-scoped only**; never `git add -A` / `git add .`.
- **⛔ No `git clean -fd` / `git stash` / `git checkout .` / `git reset --hard`
  in the shared working tree** — a `git clean -fd` already destroyed a running
  batch's runlog this session. Remove only your own files by name.

## References
- `docs/notes/2026-06-22-dv-band-definitions.md` (#415 — the bands this multi-sources),
  `docs/notes/2026-06-22-415-dv-band-definitions-task-brief.md`.
- `docs/notes/CORPUS_INDEX.md`; `feedback_corpus_document_policy`,
  `reference_cyclers_pdf_private_repo`.
