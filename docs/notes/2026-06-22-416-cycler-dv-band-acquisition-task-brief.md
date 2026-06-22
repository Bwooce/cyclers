# Task #416 — Acquire + digest cycler-ΔV-band primary references (multi-source the #415 bands)

**Status:** task brief. Created 2026-06-22 out of #415: the cycler ΔV-band
*structure* (two axes; geometric AR/TR criterion) is well multi-sourced, but the
**Axis-B real-ephemeris ΔV magnitude thresholds** (Russell 2004's < 1 / < 10 /
< 300 m/s per-7-cycle) rest on **one primary source**. This task acquires +
digests additional primaries to give Axis B **≥ 3 independent anchors** and
either corroborate or correct those cutoffs.

**CORRECTION (2026-06-22, after inventorying `/home/bruce/dev/cyclers_pdf/papers/`
— 161 PDFs):** this is **NOT primarily an acquisition task**. The original brief
wrongly assumed the targets were paywalled; a direct inventory shows the key
papers are **already on disk, just undigested for the ΔV-band question**. So #416
is mostly **digest the local holdings**, runnable now with no paper access. Only
**Rauwolf 2002 is a genuine acquisition gap** (confirmed not local).

## Targets — LOCAL (read these PDFs; no acquisition needed)

1. **`AAS-22-015-pascarella-pony-express.pdf`** — the "Pony Express" cycler paper
   (the conference version of the Genova/Pascarella work I'd mis-listed as the
   paywalled JSR 10.2514/1.A35091). **Local, no digest.** HIGHEST priority — most
   likely second independent Axis-B (impulsive + SEP maintenance ΔV) anchor.
2. **`cuevas-del-valle-2023-optimal-floquet-stationkeeping-relative-dynamics-three-body-aerospace.pdf`**
   — a dedicated **station-keeping** paper (maintenance ΔV in the 3-body problem).
   Local, no digest. Directly band-relevant (the maintenance/near-ballistic band).
3. **`chen-2002-earth-mars-cyclers.pdf`** (+ `.txt` OCR),
   **`chen-russell-ocampo-2002/2005-...free-return.pdf`** (+ `.txt`),
   **`chen-2012.pdf`** (+ `.txt`) — Chen cycler / multiple-impulse-free-return ΔV
   content. Local, OCR'd, no digest.
4. **`patel-2019-earth-mars-cycler-vehicle-conceptual-design-FIT-etd.pdf`** —
   cycler vehicle design, likely a ΔV/station-keeping budget. Local, no digest.
5. **`landau-longuski-2006-human-mars-trajectories-pt1-impulsive-JSR.pdf`** and
   **`landau-longuski-2009-comparative-assessment...`** — human-Mars cycler ΔV.
   Local; only tangentially covered by `2026-06-04-agrawal-landau-howe-mining.md`.
6. Re-mine the already-digested canonical PDFs for ΔV-*magnitude* specifically
   (the existing digests captured geometry/AR-TR, not necessarily every ΔV figure):
   `mcconaghy-landau-yam-2006-notable...`, `russell-2004-dissertation.pdf`,
   `byrnes-longuski-aldrin-1993-...`.

## Target — ACQUIRE (the one genuine gap)

7. **Rauwolf-Friedlander-Nock 2002, AIAA 2002-5046** — 1500 m/s/15-yr SEP
   maintenance figure; **not on disk**. Lower priority (band-4 / SEP boundary).
   #384-adjacent acquisition.

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
