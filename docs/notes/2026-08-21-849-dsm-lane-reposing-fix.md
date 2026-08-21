# `#849`: DSM lane re-posed per #794/#820's designated-arc semantics — bug confirmed, ZERO close_row_dsm verdicts flip

**Task:** `#849`, registered 2026-08-12 (found during `#830`, not dispatched) — `dsm_descriptor_seed.
_descriptor_params` reads a row's `free_return_arcs` POSITIONALLY (`g_tofs[0]`, `g_tofs[1]`) to
identify the g (generic) vs G (designated) arc, exactly the defect class `#820` overturned for
`build_genome`. Re-pose per `#794`/`#820`'s designated-arc semantics and re-run `#388`'s canonical
`close_row_dsm` lane across every descriptor-bearing row.

## Verdict (read this first)

**The bug was real and is fixed, but it turns out NOT to be load-bearing for `close_row_dsm`'s
actual numeric output on any of the 12 `russell-ch4`/`mcconaghy` rows.** Across all 14
`free_return_arcs`-bearing catalogue rows, re-run under the corrected posing: **zero
converged/anchor_match verdicts flip**, and every previously-recorded residual/emerged-V∞ number
(`#830`'s `dsm_388_recheck.json`, and this task's own extension to the 5 rows `#830` never ran)
reproduces byte-for-byte. What DOES change: **4 rows** (both `russell-ocampo-*` rows plus
`5.30ggF3`/`5.75ggF3`) move from "spuriously fed a wrong-arc value into this lane's g/G-only shape
model, then failed downstream" to **honestly `None` (out of scope for this model)** — same ultimate
no-seed/no-closure outcome, for the structurally correct reason. **No `catalogue.yaml` writeback
warranted, and no adjudication follow-up registered** — there is nothing new to adjudicate; no row
gained a closure or an anchor match it didn't already have on record.

## 1. Why the bug existed, and why it turned out inert here

`_descriptor_params` fed `(g_tofs[0], g_tofs[1])` — the first two `free_return_arcs` entries with
a non-null `tof_years`, in raw catalogue list order — into `self_seeding.g_arc_branches`, whose
own model (`free_return_chain.free_return_chain_correct`) is documented as "Russell's
generic-return construction: two distinct Earth-to-Earth free-return arcs (arc-1 = g, arc-2 = G)".
Per `#794`'s primary-source semantics (Russell 2004 §4.8) and `#820`'s census, the designated
(uppercase) arc is `arcs[0]` for only 3 of the 12 `russell-ch4`/`mcconaghy` rows — so a positional
read swaps the g/G roles on the other 9 whenever `arcs[0]` happens to be lowercase.

Fixed by promoting `#820`'s own `_designated_arc_index` helper (uppercase `raw_descriptor` =
designated, raises unless exactly one) out of `scripts/campaign_russell12.py` into a shared
`cyclerfinder.search.descriptor.designated_arc_index`, now imported by BOTH lanes — `#820`'s
`build_genome` (thin wrapper, behaviour byte-unchanged) and this lane's `_descriptor_params`
(fixed for real). `_descriptor_params` now: (1) finds the designated arc by case; (2) requires it
be `arc_type == "generic"` (a Lambert-solvable g/G arc, not full-rev) — this g/G-only model has no
way to represent a full-rev-designated leg, so a row where the designated arc is itself `F`
(full-rev) is honestly out of scope, not force-fed a fabricated ToF; (3) finds exactly one OTHER
generic arc to play the non-designated `g` role.

**Why the residual/V∞ numbers didn't move even where the swap DID fire (`9.353Gg2`, `3.78Gg3`,
`9.94Gg3`):** `_descriptor_params`'s output feeds `g_arc_branches` for the coplanar arc SHAPE
only — `seed.arc_a_au`/`arc_e`/`transit_branch`, which ARE genuinely different under the corrected
posing (measured below) — but the corrector's actual seed vector (`x0`'s per-leg ToFs, the
`t0`/`bounds`) is built separately in `seed_dsm_chain_from_descriptor`, straight from the row's
sourced `free_return_arcs` (in raw list order, unaffected by which one is "designated") and
`invariants.transit_times_days`. The shape-fit arc is consumed only as a FALLBACK ToF for a
cross-body leg with no sourced transit time on record, and as the (audit-only, not corrector-input)
`big_g_tof_days` bound width — neither of which turned out to bind on any of these rows. So the bug
was real, the identification was wrong, but nothing downstream in `close_row_dsm` actually reads
the swapped value in a way that changes the converged outcome — a clean, verified, non-eventful
negative.

| row | designated arc position | old (g,G) yr | new (g,G) yr | arc_a_au/e OLD | arc_a_au/e NEW | per-leg seed ToF (unchanged) |
|---|---|---|---|---|---|---|
| `9.353Gg2` | `arcs[0]` | (1.7238, 2.5469) | (2.5469, 1.7238) | 1.7949 / 0.4926 | 1.7177 / 0.4772 | [629.62, 85.0, 930.26] d |
| `3.78Gg3` | `arcs[0]` | (2.9043, 3.5018) | (3.5018, 2.9043) | 1.3100 / 0.2366 | 1.3100 / 0.2366 (symmetric) | [1060.8, 171.0, 1279.03] d |
| `9.94Gg3` | `arcs[0]` | (1.7025, 4.7037) | (4.7037, 1.7025) | 1.4131 / 0.3894 | 1.4338 / 0.3934 | [621.84, 82.0, 1718.03] d |

## 2. Full 14-row re-classification (corrected posing vs `#830`'s uncorrected record)

Run: `uv run python scripts/screen_849_dsm_reposing_recheck.py` (real DE440 ephemeris,
`close_row_dsm` — `#388`'s canonical lane — under BOTH postings for direct comparison). Full
per-row record: `data/found/849_dsm_reposing_recheck/dsm_reposing_recheck.json`.

| row | sourced V∞ E/M | OLD (buggy, `#830`-matching) | NEW (corrected) | flip? |
|---|---|---|---|---|
| `mcconaghy-2006-em-k2` | 4.70/5.00 | seeded, conv=F, **match=T**, res 9.05 | identical | no |
| `russell-ocampo-4.3.1-5` | 3.1/2.5 | **seeded=T** (wrongly — fed the h-arc's 0.5yr ToF as "G"), conv=F, match=F, res 12.64 | **seeded=F** (honest: 1 generic + f + h, not a g/G pair) | reclassified, same null outcome |
| `russell-ocampo-2.5.1+0` | 7.8/9.9 | **seeded=T** (same defect: h-arc 0.5yr fed as "G"), conv=F, match=F, res 32.27 | **seeded=F** | reclassified, same null outcome |
| `russell-ch4-4.991gG2` | 4.99/5.10 | seeded, conv=F, **match=T**, res 9.18 | identical | no |
| `russell-ch4-8.049gGf2` | 8.05/10.02 | seeded=F (arc doesn't reach — NO-CLOSE either way; designated already at `arcs[1]`) | identical | no |
| `russell-ch4-9.353Gg2` | 9.35/10.52 | seeded, conv=F, match=F, res 29.47 | identical numbers (see §1) | no |
| `russell-ch4-3.64gGg3` | 3.64/4.59 | seeded=F (NO-CLOSE either way) | identical | no |
| `russell-ch4-3.78Gg3` | 3.78/4.63 | seeded, conv=F, match=F, res 20.15 | identical numbers | no |
| `russell-ch4-5.30gGf3` | 5.30/9.17 | seeded=F (NO-CLOSE either way) | identical | no |
| `russell-ch4-9.94Gg3` | 9.94/10.76 | seeded, conv=F, match=F, res 28.51 | identical numbers | no |
| `russell-ch4-3.66gfF3` | 3.66/4.66 | seeded=F (only 1 generic arc; `len(g_tofs)<2` even under the old code) | identical | no |
| `russell-ch4-5.30ggF3` | 5.30/5.44 | **seeded=F** (old `_descriptor_params` DID return a params tuple from the two loop g-arcs, but `g_arc_branches` then raised "orbit does not reach body" downstream — `#830`'s documented "off-family" case) | **seeded=F** (now fails honestly IN `_descriptor_params`: designated arc is `F`, full-rev, out of scope) | reclassified, same null outcome |
| `russell-ch4-5.75ggF3` | 5.75/9.36 | seeded=F (same downstream-raise path) | seeded=F (same honest reason) | reclassified, same null outcome |
| `russell-ch4-6.44Gg3` | 6.44/3.74 | seeded, conv=F, match=F, res 34.94 | identical (designated already `arcs[1]`) | no |

No row's `converged`/`anchor_match` triple changes. The only PASSING (`anchor_match=True`) rows
under this lane remain exactly `mcconaghy-2006-em-k2` and `russell-ch4-4.991gG2` — both already
V0/V3 respectively on other evidence, and both unaffected by the bug (designated was already
`arcs[1]` on both). `#830`'s "the high-V∞ rows collapse off-anchor, the low-energy rows keep
`anchor_match=True`" energy-selectivity finding is confirmed to survive the corrected posing
exactly as recorded, not merely re-asserted.

**Cross-check via the pre-existing driver:** `scripts/dsm_closure_batch.py` ("#404/#388 Component
4" — runs `close_row_dsm` over every row where `seed_dsm_chain_from_descriptor` returns a seed,
unmodified here) independently reproduces the same 6 seedable rows and identical numbers under
the corrected code: `data/runs/dsm-closure-20260821T205218.jsonl`.

## 3. The ocampo-row finding is worth flagging on its own

The registration text described the bug as "reads `free_return_arcs` positionally" — a swap
between two arcs of the SAME type. The two `russell-ocampo-*` rows show the bug was actually worse
for them: their `free_return_arcs` is `[g(generic), f(full-rev), h(half-rev)]` — only ONE generic
arc. The old code's `len(g_tofs) >= 2` guard checks `tof_years is not None`, and half-rev (`h`)
arcs DO carry a `tof_years` (`descriptor.py`: "`tof_years` -- Earth-Earth leg ToF in years; g/h
arcs only"). So the old code silently swept the `h`-arc's 0.5-year ToF in as if it were a second
GENERIC arc, feeding it into a model that assumes both inputs are Lambert-solvable g/G arcs. This
never mattered for any recorded verdict (both rows already fail `anchor_match` either way, and
neither is descriptor-bearing enough to have been part of any promotion claim), but it is the kind
of silent type confusion the positional read invited beyond the case-swap the registration named.

## 4. Code changes

- `src/cyclerfinder/search/descriptor.py`: added `designated_arc_index()` — the shared,
  case-based (not positional) designated-arc identification, promoted out of
  `campaign_russell12.py` (docstring/behaviour preserved verbatim).
- `scripts/campaign_russell12.py`: `_designated_arc_index` is now a thin wrapper delegating to
  the shared function (dedup only — `build_genome`'s own behaviour is byte-unchanged; `#820`'s
  pinned tests still pass).
- `src/cyclerfinder/search/dsm_descriptor_seed.py`: `_descriptor_params` re-posed per §1 — finds
  the designated arc by case, requires it be `arc_type == "generic"` (else `None`, honestly out
  of scope), finds exactly one other generic arc for the `g` role.
- `scripts/screen_849_dsm_reposing_recheck.py`: new driver, runs `close_row_dsm` under BOTH the
  verbatim pre-#849 positional posing and the corrected posing across every descriptor-bearing
  row, for direct comparison. NO catalogue writeback (compute-only, per the established
  compute/adjudicate split).
- `tests/search/test_dsm_descriptor_seed.py`: 5 new tests — designated-by-case-not-position pin
  (`9.353Gg2`), full-rev-designated rows correctly `None` (`5.30ggF3`, `5.75ggF3`,
  parametrized), non-two-generic-arc rows correctly `None` (both ocampo rows, parametrized), and
  a `close_row_dsm` regression pin locking the "numerically unchanged" finding for the 3
  swap-affected rows against `#830`'s own recorded numbers.

## 5. Follow-ups

**None registered.** No row's verdict changed in a way that bears on `validation_level` /
`our_status` / `orbit_class` — the two already-passing rows (`mcconaghy-2006-em-k2`,
`russell-ch4-4.991gG2`) were already correctly posed and are unaffected; every other row's
no-seed/no-close outcome is unchanged, only its stated REASON is now honest. Nothing
catalogue-worthy to adjudicate. (Current highest registered task number checked before
concluding this: `#859`.)

## 6. Verification

- `uv run ruff check .` / `ruff format --check .` — clean on every file this task touched
  (`src/cyclerfinder/search/descriptor.py`, `dsm_descriptor_seed.py`,
  `scripts/campaign_russell12.py`, `scripts/screen_849_dsm_reposing_recheck.py`,
  `tests/search/test_dsm_descriptor_seed.py`); two pre-existing issues in
  `resonant_atlas_stage_a.py`/`run_859_resonant_atlas_stage_a.py` belong to a concurrent `#859`
  build, not this task, and are untouched here.
- `uv run mypy src tests` — clean except the same pre-existing `#859`-lane error in
  `resonant_atlas_stage_a.py`, likewise untouched.
- `uv run pytest tests/search/test_dsm_descriptor_seed.py -q` — 13 passed (8 pre-existing + 5
  new).
- `uv run pytest tests/search -q` — sanity pass (chunked; see commit log for the exact split run
  under this session).
