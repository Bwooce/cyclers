# #587: corrected — the 6 Uranian rows needed `invariants`/`legs`, not `orbit_elements.cr3bp`

**Date:** 2026-07-14

#587's original spec (written during the earlier website-bug investigation) assumed the fix was
"populate `orbit_elements.cr3bp` (`jacobi_constant`, `period_nd`, `stability_index`) for the 6
#569 Uranian rows" because the website's placeholder text claims that block is "tabulated above."
Investigating before executing found this premise was **wrong**: per this project's own catalogue
schema (`src/cyclerfinder/data/catalog.py::CatalogueEntry.fully_defined`), the `cr3bp` identity
triple applies ONLY to `cycler_class="non-keplerian"` rows (a single rotating-frame periodic orbit).
All 6 rows are `cycler_class="multi-arc"` — a Uranus-centered tour across two distinct moon-pair
restricted-three-body systems (e.g. Uranus-Umbriel then Uranus-Oberon), which has no single Jacobi
constant/period/stability triple to tabulate at all. The schema's own completeness check for
`multi-arc` rows uses `invariants` (cycle-level descriptors: `aphelion_ratio`, `transit_times_days`,
`turn_ratio`), not `cr3bp`. This exact situation already has an established precedent in the
catalogue: the Jovian/Saturnian Tisserand-graph tours carry a note stating "Jacobi constant is not
conserved in this model."

## What was actually done

1. **Re-ran the existing `scripts/backfill_invariants.py`** (idempotent, previously-run tool that
   simply hadn't been re-executed since 51 multi-arc rows — including these 6 — were added to the
   catalogue). This inserted generic gap-flagged `invariants` blocks (all null, honest — the
   automated text-extraction found no explicit AR/TR/transit-time values in these rows' free text).
   Incidental effect: touches 45 OTHER pre-existing multi-arc rows too (unrelated to #587), all via
   the same safe, idempotent, gap-only mechanism — no fabricated data anywhere.
2. **Hand-corrected the 6 target rows' `invariants` blocks** with real sourced data:
   `transit_times_days` populated from each row's own already-cited source
   (`data/silver_327_verified.jsonl` for `umbriel-oberon`; `data/gauntlet_566_five_representatives.jsonl`
   for the other 5) — the same `tof_days` values already used to compute their `vinf_kms_at_encounters`.
   `aphelion_ratio` and `turn_ratio` are left `null` with an accurate reason (not the generic
   "Russell source" boilerplate the automated tool inserts, which doesn't even apply here).
3. **Populated `legs[]`** (previously `legs: []`) with the same sourced `tof_days` per leg (`n_revs:
   0`), fixing `legs_tof_days` from empty to real values.
4. **Added `data_gaps[]` entries** for `invariants.aphelion_ratio`, `invariants.turn_ratio`, and
   `orbit_elements.cr3bp` — the last one explicitly documenting WHY no `cr3bp` block exists (multi-arc
   tour, no single restricted-three-body system), matching this project's honest-gap convention.
5. **Fixed the website's placeholder text** (`cyclers.space/src/components/OrbitView.astro`): it
   previously showed the SAME "not renderable... CR3BP identity tabulated above" text for both genuine
   single-orbit `non-keplerian` rows (correct — e.g. the Ross Earth-Moon rows, which DO carry a real
   `cr3bp` block) and these `multi-arc` + `model_assumption=cr3bp` rows (wrong — no such tabulation
   exists or is coming). Now distinguishes the two cases (`isPatchedMultiArcTour` check on
   `cycler_class`) and shows accurate, honest text for the multi-arc case, pointing to the row's own
   `notes` for per-leg detail instead of a nonexistent tabulated identity.

## Result

None of the 6 rows are `fully_defined=True` (they still carry acknowledged, honest `data_gaps`
entries for `aphelion_ratio`/`turn_ratio`/`cr3bp` — this is CORRECT, not a bug: `fully_defined`
returns `False` whenever `data_gaps` is non-empty by design). What changed is that the gaps are now
accurately characterized (not-applicable vs. genuinely-unknown, rather than a schema mismatch), real
sourced data (`transit_times_days`, `legs[].tof_days`) replaced empty placeholders, and the website no
longer makes a false claim about tabulated data that will never exist for these rows.

## Scope note

No new search, no novelty claim, no change to any row's `validation_level` or `our_status`. Pure
catalogue-accuracy + website-honesty correction, verified against the actual catalogue schema code
rather than assumed from the original task's framing.

`#587 STATUS: CLOSED` (scope corrected from the original spec; the corrected fix is what actually
needed doing).
