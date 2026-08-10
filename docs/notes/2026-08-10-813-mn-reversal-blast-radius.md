# #813 — blast-radius audit of #794's f(M:N) descriptor-reversal fix: NO VERDICT FLIPS

**Date:** 2026-08-10
**Task:** `#813` (registered by `#794`, dispatched 2026-08-10) — per
`[[feedback_bugfix_invalidates_past_searches]]`, audit every past campaign that consumed the
pre-`#794` REVERSED `f(M:N)` read (M was taken as spacecraft revs and N as years; the primary
sources say M = Earth years = ToF, N = spacecraft revs) for the four affected non-1:1 rows —
`russell-ch4-5.30gGf3` (f 3:2), `-3.66gfF3` (F 3:2), `-5.30ggF3` (F 3:2), `-5.75ggF3` (F 2:1) —
and RE-RUN any plausibly-affected negative with the corrected seed.

## Headline

**No verdict flips.** Exactly ONE past campaign ever consumed the reversed M:N read for these
rows — the #125/#135 `campaign_russell12.py` lambert-genome closure campaign (2026-06-06/07) —
and it did so through its **own duplicated copy of the reversed parse**, which `#794` (fixing
only `search/descriptor.py`) did not touch. That duplicate is now fixed (delegated to the
corrected `descriptor.py` functions) and the campaign re-run A/B on the 4 rows: every outcome
class is unchanged (4× CLOSE-OFF-ANCHOR, probe 4× WALKED-AWAY). The corrected seeds are
measurably better (every seed-at-truth residual drops, one row's by 36%) but none crosses the
0.1 km/s closure floor, so no recorded negative was caused by the reversed read. The rows' V1
levels rest on the #137 free-return genome (aphelion+transit seeds — unaffected), not on this
campaign.

## 1. Consumer audit (who ever consumed the descriptor-derived M:N?)

Traced every consumer of `search/descriptor.py`'s parse AND every independent re-implementation
of the M:N read (`grep` for `parse_free_return_arcs` / `arc_tof_seed_days` /
`arc_to_leg_topology` / `resolve_seed` / `split(":")` across `src/`, `scripts/`, `tests/`):

| consumer | M:N consumed? | affected rows' negatives? |
|---|---|---|
| `search/seed_ladder.py` Rung 1 (the path #794 fixed) | yes, via `parse_free_return_arcs` | **zero production consumers** — `resolve_seed` is imported only by `tests/search/test_seed_ladder.py` and `tests/test_vem_rediscovery.py` (VEM rows, no f/F arcs). No campaign ever ran Rung-1 descriptor seeds for these rows. |
| `scripts/campaign_russell12.py` (#125/#135 lambert genome) | **yes — its OWN duplicated reversed parse** (`_full_rev_revs` read revs from M, years from N) | YES — the one real consumer; see §2 |
| `#177` self-seeding (`triage_self_seeding.py`, `validate_self_seeding_reachable.py`, `search/self_seeding.py`) | no — seeds only from g/G arcs' `tof_years` (f/F arcs carry `tof_years: null` and are filtered out) | 5.30gGf3/5.75ggF3's OFF-FAMILY-NO-CLOSE and 3.66gfF3's NO-DESCRIPTOR verdicts are pure g-arc coplanar geometry, untouched by the M:N bug |
| `search/dsm_descriptor_seed.py`, `search/multiarc_closure.py` (dsm/multiarc closure batches) | no — same `tof_years is not None` filter | unaffected |
| `search/continuation_batch.py` (V1→V3 lift, includes 5.30gGf3/5.75ggF3) | no — seeds from sourced aphelion+transit | unaffected |
| `#137` free-return genome (the source of these rows' current V1) | no — seeds `(a, e)` from sourced aphelion+transit | unaffected |
| `search/cycler_assembly.py::descriptor_to_phsi` (line 376 reads `split(":")[0]` as a rev count) | latent only — it reads `resonance` off *generic* arcs, and no catalogue generic arc carries a resonance (verified by direct scan), so `i` always falls back to 1 | dead in practice; module already flagged approximate by #794 |
| `data/gauntlet_ledger.jsonl` bronze entries for the 4 rows (2026-06-06) | no — provenance-axis (Axis C) verdicts from `sweep_gauntlet_ledger.py`, no seeds involved | unaffected |
| `data/empty_regions.jsonl` | no entries for these rows | n/a |

Other `split(":")` hits (`jovian_resonant_families.py`, `scan_558_*.py`) are different domains
(Jupiter-moon resonance labels, Uranus pair CLI args), not Russell descriptors.

## 2. The one real consumer: campaign_russell12 lambert genome — fixed + A/B re-run

`build_genome()` seeded each full-rev E-E loop with the REVERSED read: e.g. `F(2:1)` seeded
365.25 d / 2 revs instead of the correct 730.5 d / 1 rev; `f/F(3:2)` seeded 730.5 d / 3 revs
instead of 1095.8 d / 2 revs. Its recorded 2026-06-07 circular-model verdicts for the 4 rows
(`data/runs/russell12-circular-2026-06-07.jsonl`): 5.30gGf3 CLOSE-OFF-ANCHOR, 3.66gfF3
**NO-CLOSE**, 5.30ggF3 CLOSE-OFF-ANCHOR, 5.75ggF3 **NO-CLOSE**; seed-at-truth probe
(2026-06-06 like-for-like note): WALKED-AWAY ×4, with 3.66gfF3/5.75ggF3 at the 1000 km/s
"Lambert pathology at truth" sentinel — an artifact of demanding 2–3 revs inside a 1-year leg,
i.e. of the reversed read itself.

**Fix:** the duplicated parse is deleted; the full-rev branch now calls the corrected
`arc_to_leg_topology` / `arc_tof_seed_days` from `cyclerfinder.search.descriptor`. A `--rows`
CLI filter was added for targeted re-runs. Verified the fixed genomes change exactly the 4
non-1:1 rows (3:2 → 1095.8 d / 2 revs; 2:1 → 730.5 d / 1 rev) and leave every 1:1/generic row
byte-identical (mcconaghy, 8.049gGf2 checked directly), so the frozen
`tests/search/test_russell12_likeforlike_probe.py` diagnosis rows are untouched.

**A/B re-run** (today's HEAD, so the #205 Lambert fix is present in BOTH arms — the June-07
baseline is NOT directly comparable, hence the same-code control): `--model circular
--epochs 256 --workers 8 --probe-at-truth --phase-epochs 256`, 4 rows.
Runlogs: `data/runs/russell12-circular-20260810T-813-A-preM-Nfix.jsonl` (control, reversed
parse — reproduces the June-07 recorded seeds exactly) and
`...-20260810T-813-B-postM-Nfix.jsonl` (fixed parse).

| row | A outcome (reversed seed) | B outcome (fixed seed) | A truth-res km/s | B truth-res km/s | FLIP |
|---|---|---|---|---|---|
| 5.30gGf3 | CLOSE-OFF-ANCHOR (241 closed) | CLOSE-OFF-ANCHOR (204) | 10.72 | **6.89** | NO |
| 3.66gfF3 | CLOSE-OFF-ANCHOR (11) | CLOSE-OFF-ANCHOR (17) | 31.89 | **26.69** | NO |
| 5.30ggF3 | CLOSE-OFF-ANCHOR (231) | CLOSE-OFF-ANCHOR (216) | 37.55 | **26.11** | NO |
| 5.75ggF3 | CLOSE-OFF-ANCHOR (37) | CLOSE-OFF-ANCHOR (200) | 21.20 | **17.00** | NO |

(The June-07 NO-CLOSEs for 3.66gfF3/5.75ggF3 already close under today's code in BOTH arms —
that is the #205 Lambert completeness fix, an already-documented effect per
`docs/notes/2026-06-12-lambert-blast-radius-rerun.md`, not an M:N effect.)

Probe verdicts: WALKED-AWAY ×4 in both arms. The corrected seed genuinely improves the posed
geometry — every truth residual drops (and 5.75ggF3's probe corrector now converges to 0.0
where the control stalled at 9.65) — but the residual-at-truth stays 7–27 km/s above the 0.1
floor, so the standing #135 diagnosis ("the sourced geometry is not a residual-zero point of
this genome as posed") is unchanged for these rows. The reversed M:N read inflated the June
note's truth residuals for the 4 rows (e.g. 5.30ggF3's 37.55 is really 26.11, and the two
1000-sentinel "Lambert pathology" entries were pure seed artifacts), but did NOT cause any
false negative: no closed solution appears at anchor under the corrected seed either.

## 3. Why the residual stays high even with the corrected M:N (registered follow-up)

`build_genome` still carries the campaign's original structural assumption that `arcs[0]` is
the Mars free-return arc and `arcs[1:]` are the E-E loops. `#794` established from the primary
sources that the DESIGNATED (uppercase) arc is the Mars-transit leg — which is `arcs[1]` (G)
for 5.30gGf3 and `arcs[2]` (F) for the three F-rows. So even with the corrected M:N, the genome
assigns a designated arc to the E-E loop set and drops a real loop arc for all 4 rows (and the
slack-leg elimination then reconstructs a wrong loop ToF). That mis-posing is pre-existing,
affects the campaign lane only, and is registered as `#820` (re-pose with #794's semantics +
now-written-back loop segments) rather than silently widened into this audit.

## 4. Also found / registered

- **`#820`** — re-pose the russell12 lambert genome per #794 semantics (designated = uppercase
  arc, loops = the other arcs in itinerary order, seeds from #794's written-back
  `loop-ee-*` segments) and re-run; a truth-residual≈0 there would be V1-grade multi-arc
  closure evidence beyond the single-ellipse free-return genome.
- **`#821`** — `russell-ch4-5.30gGf3`'s catalogue prose note still explains f(3:2) with the
  REVERSED convention ("spacecraft 3 revs vs Earth 2 revs"); fix the prose (catalogue-ratchet
  discipline applies; deferred here to avoid colliding with `#811`'s concurrent catalogue
  writeback). The other three rows' notes don't state a rev reading; `#794`'s segments are
  correct.

## 5. Verification

- `uv run ruff check .` / `ruff format --check .` — clean.
- `uv run mypy src tests` — clean (full run).
- `uv run pytest tests/search/test_descriptor_arctype.py tests/search/test_descriptor_tof.py
  tests/search/test_descriptor_parse.py tests/search/test_descriptor_catalogue.py
  tests/search/test_russell12_likeforlike_probe.py tests/scripts -q` — see commit; the
  likeforlike probe (slow, direct consumer of the edited script) included deliberately.
- No catalogue / empty-regions edits in this task; the two new runlogs are APPENDED
  method-versioned records (negative-results-registry discipline), the 2026-06 originals are
  untouched.
