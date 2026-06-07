# Catalogue ingest — Rall 1970 + Russell 2004 Table 3.4 (#142)

Date: 2026-06-07
Scope: ingest the catalogue-eligible member rows transcribed in the two source
notes, in two gated commits. Source notes (verbatim transcriptions, cited by
author/title/report only — never any file path or private repo):

- `docs/notes/2026-06-07-rall-1970-appendices-transcription.md`
- `docs/notes/2026-06-07-russell-2004-member-tables-transcription.md`

All ingested rows: `validation_level: V0`, v4.2 per-segment provenance
(center / tof_days / source_ephemeris), `data_gaps[]` for every missing field,
source cited by author/title/report only.

---

## Batch 1 — Rall 1970 (MIT TE-34) free-fall periodic Earth-Mars orbits

Commit `3ebbb7a` `data/catalogue: ingest Rall 1970 (MIT TE-34) free-fall periodic Earth-Mars orbits`

**15 rows ingested** (provenance key `rall-1970-te34`, added to SOURCE_REGISTRY
and the schema `orbit_source`/`vinf_source` enums):

- Appendix E (circular-coplanar Model I.B): `rall-1970-m4-1`, `-m6-1`, `-m6-2`, `-m6-3`
- Appendix F (eccentric-inclined Model III.B): `rall-1970-m4-1a`, `-m5-1a..e`, `-m5-2a..e`

All `cycler_class: multi-arc`, `source_ephemeris: patched-conic mean-element
(Rall TE-34, 1969)` (NEVER a DE kernel).

### Conversions (units kept faithful; conversions documented per-row)
- V_inf: EMOS -> km/s at 1 EMOS = 29.77 km/s (Earth Mean Orbital Speed, per Rall
  / mining-note Table 3-1). Original EMOS value preserved in every note.
- Turn angles: degrees (App F, SOURCED where tabulated, e.g. first Earth
  encounter of each App-F orbit).
- Passing distances: local planetary radii (Earth radii at E, Mars radii at M) —
  NOT converted to km (`min_altitude_km` left null; km conversion is derived).
- App-E per-leg ToF: difference of consecutive JD-2440000 dates (the clean
  encoding), e.g. M4-1 E->M = 1194-1030 = 164 d.

### Data gaps / hazards
- **App-F absolute dates blocked**: the last-5-digits-of-JD hyphen encoding has a
  documented reconstruction ambiguity (source note §1). All 11 App-F rows keep
  per-version segment `tof_days: null` and a `data_gap` pointing at the
  resolution path (anchor App-F tokens to App-E JD-2440000). Never guessed.
- **multi-arc transit_times_days invariant**: the test
  `test_multi_arc_invariants::test_transit_times_present_for_all_multi_arc_rows`
  requires every multi-arc row to carry `invariants.transit_times_days`. The
  App-F rows' per-version ToFs are blocked by the date ambiguity, so
  `transit_times_days` is anchored to the **App-E clean JD-2440000 E->M leg of
  the same family** (sourced, not invented; source note §3):
  - M4-1a -> 164 d (M4-1 App-E)
  - M5-1a..e -> 107 d (M5-1 App-E: E1901->M2008)
  - M5-2a..e -> 160 d (M5-2 App-E: E1818->M1978)
  Per-version App-F segment `tof_days` stay null (the absolute-date ambiguity is
  not resolved by the family anchor). This was the one fix to the prior agent's
  WIP (see "Prior-WIP reconciliation").
- Mars V_inf for App-F rows varies per encounter; stored value is the
  per-encounter range midpoint with a `data_gap` (full per-encounter list in the
  source note §4).
- `max_turning_angle_deg` is COMPUTED (max_bend) and flagged not-source-attested.

### Deferred (NOT ingested) — recorded with reasons (source note §6)
- **MARGINAL (Model II only, not realistic Model III)**: M6-4, M6-5, M6-6, M6-7.
- **NOT eligible (author states the eccentric/"a" version does NOT exist /
  intersects a planet)**: M5-3, M5-5, M5-6, M6-8, M6-9, M6-10.
- **DEFER (existence verdict not printed / undetermined)**: M5-4 (no convergence
  near 180 deg), M6-11 (no existence comment).
- App-E base rows M5-1, M5-2 are not ingested as standalone rows (only their
  App-F a..e versions are eligible per the brief); their clean App-E ToFs are
  used only as the family transit-time anchors above.

---

## Batch 2 — Russell 2004 dissertation Table 3.4 circular-coplanar cyclers

Commit `5dd4742` `data/catalogue: ingest Russell 2004 Table 3.4 circular-coplanar cyclers`

**16 rows ingested** (provenance key `russell-2004-t34`, already in the registry/
schema from prior work), as `russell-ocampo-<p.h.s+/-i>` mirroring the existing
Table 3.4 row shape, `cycler_class: multi-arc`,
`source_ephemeris: circular-coplanar simple model (Russell 2004 Ch.3)`.

### De-duplication (verify by descriptor — brief requirement)
Table 3.4 (p.83) has **44** rows. Reconciled against the live catalogue:
- **27 already catalogued** as `russell-ocampo-*` (2.1.1+2, 2.3.1+1, 2.5.1+0,
  3.1.1+3, 3.1.2+1, 3.3.1+2, 3.5.1+1, 3.5.2+0, 3.7.1+1, 3.9.1+0, 4.0.3+1,
  4.1.1-5, 4.3.1-5, 4.3.1-4, 4.5.1-4, 4.5.1-3, 4.5.2-2, 4.5.3-1, 4.7.1-3,
  4.9.1-2, 4.9.2-1, 4.10.1+2, 4.11.1-2, 4.12.1+1, 4.13.1-1, 4.14.1-1, 4.14.1+0).
- **1 already catalogued under a different id**: Aldrin cycler 1.0.1.-1 ==
  `aldrin-classic-em-k1-outbound` (Table 3.4 explicitly cited on that row).
- **16 genuinely absent** -> ingested: 3.1.1+2, 3.1.3+0, 3.5.1+2, 4.1.1-6,
  4.1.1-4, 4.1.2-3, 4.1.2-2, 4.1.4-1, 4.6.1-4, 4.6.3+0, 4.7.1-2, 4.8.1+3,
  4.8.1+2, 4.9.1-3, 4.10.1-3, 4.12.1-2.

Note: the source note's prose estimated "~38 eligible-but-absent" rows; that
predates the live `russell-ocampo` coverage. The descriptor reconciliation
(27 present + Aldrin) is authoritative -> 16 absent. Brief's "never duplicating
existing members (verify by descriptor/period)" honoured.

### Conversions / derivations (flagged)
- AR / TR / E->M time (days) / Earth V_inf / Mars V_inf / per-flyby geocentric
  turn angles: verbatim Table 3.4 values.
- `aphelion_au = AR x 1.52 AU` (Mars sma; Russell §3.3) — derived.
- `max_turning_angle_deg = degrees(max_bend(mu, R+300 km, vinf))` — COMPUTED,
  not source-attested (reproduces existing T3.4 rows exactly: Earth R=6378.137,
  Mars R=3396.19, mu from core/constants).
- `period.years = k x 2.135 yr`.

### Data gaps
- Return M->E leg `tof_days`, intermediate Earth-Earth loop topology
  (n_flybys-1 loops/cycle), and per-arc / per-orbit `(a,e)`: not tabulated by
  Russell Table 3.4 -> `data_gaps[]` referencing task #54 (multi-rev Lambert).
- Near-ballistic rows (AR<1 or TR<1: 3.5.1+2, 4.1.1-6, 4.1.2-3, 4.6.1-4,
  4.8.1+3, 4.9.1-3, 4.10.1-3, 4.12.1-2) get `delta_v_kms: null` (Russell does
  not tabulate the per-cycle DV); strictly-ballistic rows (AR>=1 AND TR>=1) get
  `delta_v_kms: 0.0`.

---

## Census ratchets (moved deliberately, derived not invented)

| Census | Baseline | After batch 1 | After batch 2 |
|---|---|---|---|
| catalogue rows (`grep -c '^- id:'`) | 237 | 252 (+15 Rall) | 268 (+16 Russell) |
| cycler_class multi-arc | 203 | 218 | 234 |
| cycler_class single-ellipse | 28 | 28 | 28 |
| cycler_class non-keplerian | 6 | 6 | 6 |
| tier consistency_checked | 218 | 233 | 249 |
| tier cross_validated | 5 | 5 | 5 |
| tier unvalidated | 14 | 14 | 14 |

Tests updated in the same commit as the rows:
`tests/data/test_cycler_class_census.py` (MULTI_ARC_ALLOWLIST + counts) and
`tests/data/test_validation_tier_census.py` (EXPECTED_TIER_CENSUS).

All ingested rows are consistency_checked: each has
`orbit_source == vinf_source` at a single matching fidelity (one citation, one
fidelity) — consistent but not independently cross-validated.

---

## Prior-WIP reconciliation (what was kept vs fixed)

The prior agent (killed mid-batch-1) left the 15 Rall rows, the schema enum
addition, the SOURCE_REGISTRY entry, and the two census tests' batch-1 updates,
all uncommitted. Verification before adoption:

- **Spot-check (digit-for-digit vs the source note)**: M4-1, M4-1a, M6-1 fully
  checked; M6-2, M6-3, M5-1a, M5-2a, M5-2b sampled. All V_inf conversions
  (EMOS x 29.77), App-E ToFs (date differences), turn angles, and Mars-range
  midpoints reproduced exactly. **Kept as-is.**
- **One fix**: the 11 App-F rows had `invariants.transit_times_days: null`,
  which failed `test_multi_arc_invariants::test_transit_times_present_for_all_
  multi_arc_rows`. Populated from the App-E family E->M anchor (164/107/160 d) —
  a sourced value, not from the ambiguous App-F encoding; per-version segment
  `tof_days` kept null. This was the only change to the prior WIP.
- Concurrent-agent files (`src/cyclerfinder/nbody/shooter.py`,
  `src/cyclerfinder/search/continuation.py`, `tests/nbody/*`, `tests/search/*`)
  were deliberately NOT staged in either commit (xdist/stash-race hazard).

## Test evidence
`uv run pytest tests/data/ -q` green after each batch (only the 2 pre-existing
task-#54 XFAILs; 0 failures). Pre-commit hooks (jsonschema validation, ruff
check, ruff format, mypy) passed on both commits.
