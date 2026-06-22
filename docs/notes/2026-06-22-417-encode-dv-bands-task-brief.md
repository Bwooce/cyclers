# Task #417 — Encode the ΔV bands into `catalogue.yaml` + schema

**Status:** ready-to-execute brief. Created 2026-06-22 out of #415: the band note
(`2026-06-22-dv-band-definitions.md`) deliberately did NOT touch the data layer;
this task encodes its Axis-B taxonomy as a structured, validated catalogue field.

## Goal

Add a per-row **`dv_band`** field (Axis B — real-ephemeris deterministic maintenance
ΔV magnitude) to the catalogue, schema-validated, and backfill it **only where a
sourced ΔV + basis already exists** in the row. Never fabricate a band.

## Design (locked)

- **Field:** `dv_band` — string enum, **nullable** (null = unclassified, the default).
- **Allowed values** (from `2026-06-22-dv-band-definitions.md`, Axis B):
  - `strictly_ballistic`     — < 1 m/s / 7 cycles
  - `essentially_ballistic`  — < 10 m/s / 7 cycles
  - `low_maintenance`        — < 300 m/s / 7 cycles
  - `powered_dsm`            — ≥ 300 m/s / 7 cycles (impulsive/DSM)
  - `low_thrust_sep`         — SEP/low-thrust maintenance (out of the impulsive ladder)
  - `null`                   — unclassified (no sourced ΔV in the row, or basis unclear)
- **Orthogonal to `trajectory_regime`** (the existing powered/ballistic Axis-A label).
  Do NOT overwrite or merge them — `dv_band` is the *magnitude* axis, `trajectory_regime`
  is the *geometric* axis. A row can be `trajectory_regime: powered` AND
  `dv_band: powered_dsm`, or `ballistic` AND `strictly_ballistic`, etc.
- **Add a sibling `dv_band_source`** (string, nullable): the digest/paper the band
  assignment traces to (e.g. `russell-ocampo-2006`, `mcconaghy-2006`, `friedlander-1986`).
  Mandatory whenever `dv_band` is non-null (honesty gate — every band traces to a source).

## Backfill rules (honesty-gated)

1. Classify a row **only** if it already carries a sourced maintenance/deterministic
   ΔV with a known basis. Convert to the band's basis (total m/s over ~7 cycles) where
   the basis allows; if the conversion is ambiguous, leave `dv_band: null`.
   - e.g. `maintenance_dv_kms_per_synodic: 1.52` (Aldrin) → ~10.6 km/s/7-syn → `powered_dsm`.
   - A row with `delta_v_kms: null` and no maintenance figure → `dv_band: null`.
2. The known ballistic-family rows (S1L1, VISIT, n=7 cyclers, the stable EM cyclers)
   may be `essentially_ballistic`/`strictly_ballistic` **only** if a source in the row
   (or its `orbit_source`/`vinf_source` digest) states the near-zero maintenance — cite
   it in `dv_band_source`. Otherwise null.
3. Do NOT infer a band from `trajectory_regime` alone (ballistic ≠ a specific m/s tier).
4. Do NOT infer from V∞ or `v_infinity_leveraging_dv_kms` (that's establishment, not
   maintenance — see band note Conflicts section).

## Files

- `data/catalogue.schema.json` — add `dv_band` (enum + null) and `dv_band_source`
  (string + null) to the row properties. Not required (nullable/optional).
- `data/catalogue.yaml` — add the field to rows that qualify under the backfill rules.
- Possibly `src/cyclerfinder/.../spec` or README band-count surface — only if a count is
  published; otherwise skip.

## Ratchets / tests (BINDING — see `feedback_catalogue_edits_run_all_ratchets`)

- Any `catalogue.yaml` change ripples into MULTIPLE frozen-census ratchets
  (`cycler_class_census`, `validation_tier_census`, `rediscovery`, `validate`) + consumers.
- Run the **full** `uv run pytest tests/data tests/search -q` — NOT a hand-picked subset —
  before commit. Regenerate any census fixtures the new field legitimately changes.
- `uv run ruff check . && uv run ruff format --check .` pre-commit.
- Commit **through** the hooks (no `--no-verify`). Pathspec-scoped adds only.

## Honesty gates

- Every non-null `dv_band` has a `dv_band_source` tracing to a published source +
  basis. No band guessed from regime/V∞. A row with no sourced ΔV stays `null`.
- If most rows end up `null` (likely — most catalogue rows lack a per-row maintenance
  ΔV), that is the correct, honest outcome; do not pad coverage.

## References
- `2026-06-22-dv-band-definitions.md` (the bands), `2026-06-22-416-cycler-dv-universe-sweep.md`
  (evidence), `feedback_catalogue_edits_run_all_ratchets`, `feedback_golden_tests_sourced_only`.
