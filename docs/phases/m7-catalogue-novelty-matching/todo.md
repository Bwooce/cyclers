# M7 — Catalogue, signature matching, novelty (todo)

Working checklist for the M7 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

The order mirrors plan.md §7: predecessor recap → `data/` subpackage scaffold → `data/catalog.py` signature half → `data/catalog.py` loader+matcher half → `data/ledger.py` → `verify/crosscheck.py` → `data/writeback.py` → `data/discover.py` → catalogue-rediscovery tagging gate → local green → CI → closeout.

## Predecessor recap

- [ ] Re-read spec §16 in full — the M7 design document (record schema, canonical signature, identity matching, attribution, discovery workflow).
- [ ] Re-read spec §13.5 (canonical signature + dedup), §13.6 (work queue + ledger), §13.8 (ledger schema + cell IDs).
- [ ] Re-read spec §12.2 (representation framework — `model_assumption` partitions the matcher pool).
- [ ] Re-read spec §14 V1 (Lambert cross-check gate) and V2 (multi-lap periodicity — M7 reads the M6a `StabilityReport`).
- [ ] Re-read spec §10 (closure-frame correctness risk — M7 inherits via M6a's `verify_long_term_stability`).
- [ ] Re-read M6a's `## Hand-off to M6b` in `phases/m6a-idealized-closure-verification/todo.md`: confirm `StabilityReport` shape (`cycler_id`, `n_laps_propagated`, `max_drift_km`, `max_drift_lap_index`, `per_lap_drift_km`, `stable`, `per_lap_dv`, `total_tcm_dv`, `frame_used`).
- [ ] Re-read M5's `## Hand-off to M6a` in `phases/m5-optimisation/todo.md`: confirm `OptimisationResult` shape (`cell`, `best_cycler`, `best_score`, `closure_residual_kms`, `optimiser_history`, `converged`, `constraints_satisfied`) and the reproduced V∞ values.
- [ ] Re-read `tests/_catalogue_loader.py` — M5-era test infrastructure; M7 reuses it without modification.
- [ ] Re-read `tests/test_catalogue_rediscovery.py` — M5-era parametrised rediscovery test; M7's tagging test runs alongside it.
- [ ] Confirm `data/seed_cyclers.yaml` row count (currently 219) — M7 binding test asserts the loader reads all rows.
- [ ] Confirm `src/cyclerfinder/search/optimize.py::optimise_cell_ephemeris` still raises `NotImplementedError("requires M6 ephemeris backend")` — M7's `discover` runner leaves `enable_v3=False` until M6b lands.
- [ ] Confirm `model/cycler.py::orbit_elements_au` exists — used by `canonical_signature` to compute `(a, e)` per leg.
- [ ] Confirm parallel M6b plan-writing agent is NOT touching `src/cyclerfinder/verify/`, `src/cyclerfinder/data/`, `tests/data/`, `tests/verify/test_crosscheck.py`, or any M7 workspace file.

## `data/` subpackage scaffold

- [ ] Create empty `src/cyclerfinder/data/__init__.py`.
- [ ] Create empty `tests/data/__init__.py`.
- [ ] Confirm `uv run mypy src` clean — empty `__init__.py` is fine.
- [ ] Confirm `uv run pytest tests/data` runs (no tests yet, exits cleanly with collected=0).

## `data/catalog.py` — signature half (lands first; spec §16.2 binding)

### Skeleton + module docstring

- [ ] Create `src/cyclerfinder/data/catalog.py` with module docstring referencing spec §16.1, §16.2, §16.3, §12.2.
- [ ] Add `CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent.parent / "data" / "seed_cyclers.yaml"` resolution (matches `tests/_catalogue_loader.py` pattern).
- [ ] Stub `CanonicalSignature` frozen dataclass with all fields per plan §3.1.3.
- [ ] Stub `canonical_signature`, `signature_distance` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src/cyclerfinder/data/catalog.py` clean on the skeleton.

### Helpers — one at a time, paired tests

- [ ] Implement `_lex_min_rotation(sequence_str)` per plan §3.1.1.
- [ ] Test `test_lex_min_rotation_em_e` — `_lex_min_rotation("E-M-E") == "E-E-M"`.
- [ ] Test `test_lex_min_rotation_singleton` — `_lex_min_rotation("E") == "E"`.
- [ ] Test `test_lex_min_rotation_already_minimal` — `_lex_min_rotation("E-M") == "E-M"`.
- [ ] Implement `_bin_vinf`, `_bin_a_au`, `_bin_e` per plan §3.1.2 (round to nearest bin).
- [ ] Test `test_bin_vinf_05_kms` — `_bin_vinf(5.674) == 5.65`; `_bin_vinf(5.676) == 5.70` (bin boundary).
- [ ] Test `test_bin_a_au_01` — `_bin_a_au(1.304) == 1.30`.
- [ ] Test `test_bin_e_01` — `_bin_e(0.257) == 0.26`.
- [ ] Implement `_canonical_json(d)` per plan §3.1.3 — `json.dumps(d, sort_keys=True, separators=(",", ":"), ensure_ascii=True)`.
- [ ] Test `test_canonical_json_deterministic` — `_canonical_json({"b":1, "a":2}) == _canonical_json({"a":2, "b":1})`.

### `canonical_signature(cycler, *, model_assumption)`

- [ ] Implement per plan §3.1.3:
  - [ ] Build `bodies` as sorted unique tuple of `cycler.bodies`.
  - [ ] Compute `sequence_canonical = _lex_min_rotation("-".join(cycler.bodies))` (handling the open-sequence closing-body convention per §5 risk #12 — dedupe consecutive repeats of the closing body).
  - [ ] Derive `sense` via `_derive_sense(cycler)` (heuristic from first-leg ToF vs period; document the limitation per §5 risk #2).
  - [ ] Read `period.k` from the cycler's owning `Cell` (passed in by caller — `find_cyclers` knows the cell).
  - [ ] Build `vinf_multiset_binned` as sorted tuple of `(body, _bin_vinf(||vinf_in||))` over `cycler.encounters`.
  - [ ] Build `leg_elements_multiset_binned` as sorted tuple of `(_bin_a_au(a_au), _bin_e(e))` over `cycler.legs` via `orbit_elements_au(enc.r, leg.v_depart)`.
  - [ ] Build the hash-input dict from the in-hash fields only (period.years EXCLUDED per §5 risk #15).
  - [ ] `hash = "sha1:" + hashlib.sha1(_canonical_json(input_dict).encode("utf-8")).hexdigest()`.
  - [ ] Return `CanonicalSignature(...)` with `model_assumption` as a tag (pool filter — not in hash).
- [ ] Test `test_signature_rotation_invariant` (**spec §16.2 BINDING GATE**) — for `aldrin_cycler`, cyclically rotate `cycler.bodies` and `cycler.encounters` and `cycler.legs`; assert all rotations produce identical `.hash`.
- [ ] Test `test_signature_binning_absorbs_noise` (**spec §16.2 BINDING GATE**) — perturb V∞ by ±0.02 km/s, `(a, e)` by ±0.005 AU / ±0.005; assert hash unchanged. Then perturb by 0.06 km/s on V∞ and assert hash *changes*.
- [ ] Test `test_canonical_signature_deterministic` — two calls produce bitwise-identical `hash`.
- [ ] Test `test_signature_hash_stable_across_python_versions` — sha1 of `"hello"` matches the hardcoded reference `"sha1:aaf4c61ddcc5e8a2dabede0f3b482cd9aea9434d"`.
- [ ] Test `test_signature_open_sequence_normalisation` — `Cycler(bodies=["E","M","E"], ...)` produces the same `sequence_canonical` as catalogue's `"E-M"`. Documents §5 risk #12.

### Module-level checks

- [ ] Confirm `uv run pytest tests/data/test_canonical_signature.py` green.
- [ ] Confirm `uv run mypy src/cyclerfinder/data/catalog.py` clean.

## `data/catalog.py` — loader + matcher half

### `CatalogueEntry` and `_extract_signature_fields`

- [ ] Implement `CatalogueEntry` frozen dataclass per plan §1.1 (full record, not just rediscovery projection).
- [ ] Implement `_extract_signature_fields(entry_dict)` per plan §1.1: build signature-input dict from YAML row's `bodies`, `sequence_canonical`, `period.k`, `vinf_kms_at_encounters`, `orbit_elements.a_au`, `orbit_elements.e`. Return `None` for family-seed / citation-only rows (missing required fields).
- [ ] Test `test_extract_signature_fields_aldrin` — Aldrin row produces a non-None dict with `bodies == ("E", "M")`.
- [ ] Test `test_extract_signature_fields_family_seed_returns_none` — `jones-2017-vem-triple-family` returns None.

### `Catalog` class + `load_catalog`

- [ ] Implement `Catalog` class per plan §1.1: `by_id`, `by_hash`, `entries`, `filter()`.
- [ ] Implement `load_catalog(path=CATALOGUE_PATH) -> Catalog` per plan §1.1: yaml.safe_load, build `CatalogueEntry`s, compute signature for each constructible entry, populate indices.
- [ ] Test `test_load_catalog_returns_all_219_entries` — `len(cat.entries) == 219`.
- [ ] Test `test_catalogue_signature_index_covers_constructible_entries` — every entry from `tests/_catalogue_loader.load_constructible_entries()` has a non-None `signature_hash` in `cat.by_id`.
- [ ] Test `test_catalogue_family_seeds_have_null_signature` — family-seed entries have `signature_hash is None`.
- [ ] Test `test_catalogue_filter_by_model_assumption_partitions_cr3bp` — `cat.filter(model_assumption="cr3bp")` includes Arenstorf, excludes circular-coplanar.
- [ ] Test `test_catalogue_filter_by_bodies` — `cat.filter(bodies=("E", "M"))` returns only Earth-Mars entries.

### `signature_distance` + `match` + `MatchResult`

- [ ] Implement `MatchResult` frozen dataclass with `outcome: Literal["known", "probable-match-NEEDS-HUMAN", "novel"]`, `entry`, `distance`.
- [ ] Implement `signature_distance(sig_a, sig_b) -> float` per plan §4.4: weighted L1 over (period_years/0.5, Σvinf/0.05, Σ(a/0.01, e/0.01)).
- [ ] Implement `TAU_NEAR: Final[float] = 0.5`.
- [ ] Implement `match(candidate, catalog) -> MatchResult` per plan §3.1.4 pseudocode.
- [ ] Test `test_match_aldrin_classic_returns_known` (**M7 GATE**) — Aldrin matches `aldrin-classic-em-k1-outbound`, priority_date 1985-10-28.
- [ ] Test `test_match_novel_synthetic_returns_novel` — synthetic V∞ multiset `[("E", 4.21), ("M", 7.83)]` returns `("novel", None)`.
- [ ] Test `test_match_partition_by_model_assumption` (**spec §12.2 BINDING GATE**) — forced V∞ collision with cr3bp entry returns `("novel", None)` for a `circular-coplanar` candidate.
- [ ] Test `test_match_result_frozen` — `result.outcome = "..."` raises `FrozenInstanceError`.

### Sensitivity (informational, not gating)

- [ ] Test `test_tau_near_sensitivity` — sweep `TAU_NEAR ∈ {0.3, 0.5, 1.0}`, count `probable-match` outcomes across constructible entries; assert count at 0.5 is between 0 and 5.

### Module-level checks

- [ ] Confirm `uv run pytest tests/data/test_catalogue_loader.py tests/data/test_match.py` green.
- [ ] Confirm `uv run mypy src/cyclerfinder/data/catalog.py` clean.

## `data/ledger.py`

### Skeleton

- [ ] Create `src/cyclerfinder/data/ledger.py` with module docstring referencing spec §13.6, §13.8.
- [ ] Stub `LedgerStatus`, `LedgerEntry`, `LedgerError`, `Ledger`, `LedgerLoader` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src/cyclerfinder/data/ledger.py` clean on skeleton.

### Implementation

- [ ] Implement `_atomic_append(path, line)` per plan §3.2.1 — single `os.open(O_APPEND | O_WRONLY)`, `os.write`, `os.close`. Assert line is below 4 KiB (PIPE_BUF on Linux).
- [ ] Implement `_parse_line(line)` strict JSONL → `LedgerEntry` parser; raise `LedgerError` on schema mismatch.
- [ ] Implement `LedgerEntry` frozen dataclass with all fields per plan §1.2.
- [ ] Implement `Ledger.__init__(path)` — open path, load existing entries into in-memory dict.
- [ ] Implement `Ledger.record(entry)` — `_atomic_append`; raise `LedgerError` on duplicate `cell_id`.
- [ ] Implement `Ledger.has(cell_id)`, `Ledger.get(cell_id)`, `Ledger.iter_pending()`, `Ledger.__len__`.
- [ ] Implement `Ledger.claim(cell_id, host)` per plan §3.2.2 with first-write-wins race resolution.

### Tests

- [ ] Test `test_ledger_round_trip` (**M7 GATE**) — record → reread → field-equal. Re-recording raises `LedgerError`.
- [ ] Test `test_ledger_persists_across_restart` — write, close, reopen, `.has(cell_id)` True.
- [ ] Test `test_ledger_atomic_append_no_partial_lines` — file ends with `\n`, no partial JSON.
- [ ] Test `test_ledger_concurrent_claim_one_wins` (smoke, `pytest.mark.flaky`) — `multiprocessing.Pool(2)`, exactly one claim returns True.
- [ ] Test `test_ledger_iter_pending_skips_solved` — mix of statuses; `iter_pending` yields only `pending`.
- [ ] Test `test_ledger_entry_frozen` — `entry.status = "..."` raises `FrozenInstanceError`.

### Module-level checks

- [ ] Confirm `uv run pytest tests/data/test_ledger.py` green.
- [ ] Confirm `uv run mypy src/cyclerfinder/data/ledger.py` clean.

## `verify/crosscheck.py`

### Skeleton + module docstring

- [ ] Create `src/cyclerfinder/verify/crosscheck.py` with module docstring referencing spec §14 V1.
- [ ] Add `from . import crosscheck` to `src/cyclerfinder/verify/__init__.py` (the only edit to an M6a-produced file).
- [ ] Stub `LambertCrosscheckResult`, `crosscheck_leg`, `crosscheck_cycler` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src/cyclerfinder/verify/crosscheck.py` clean on skeleton.

### Implementation

- [ ] Implement `LambertCrosscheckResult` frozen dataclass per plan §1.4.
- [ ] Implement `_leg_endpoints(leg, cycler, ephem)` private helper.
- [ ] Implement `crosscheck_leg(leg, cycler, ephem)` per plan §3.4: call in-house `lambert`, `lamberthub.izzo2015`, `lamberthub.gooding1990`; compute per-component velocity diff in m/s; populate `LambertCrosscheckResult`.
- [ ] Implement `crosscheck_cycler(cycler, ephem)` — apply `crosscheck_leg` across every leg; return tuple of results.

### Tests

- [ ] Test `test_v1_lambert_crosscheck_aldrin` (**M7 GATE — spec §14 V1**) — every Aldrin leg passes (`max_diff_mps < 1.0e-3`).
- [ ] Test `test_v1_lambert_crosscheck_2syn_em` — every 2-syn E-M leg passes.
- [ ] Test `test_crosscheck_leg_reports_correct_n_revs_branch` — `leg.n_revs == 1, branch == "low"` → lamberthub `M=1` low-branch.
- [ ] Test `test_crosscheck_cycler_aggregates` — aggregated `pass_` is `all(r.pass_ for r in results)`.
- [ ] Test `test_lambert_crosscheck_result_frozen`.

### Module-level checks

- [ ] Confirm `uv run pytest tests/verify/test_crosscheck.py` green.
- [ ] Confirm `uv run pytest tests/verify/test_propagate.py` (M6a) still green — no regression.
- [ ] Confirm `uv run mypy src/cyclerfinder/verify` clean.

## `data/writeback.py`

### Skeleton

- [ ] Create `src/cyclerfinder/data/writeback.py` with module docstring referencing spec §16.1, §16.4, §16.5.
- [ ] Stub `apply_v0_v1_to_entry`, `apply_v2_to_entry`, `apply_v3_to_entry`, `record_rediscovery`, `register_discovery`, `serialise_entry_yaml` with full docstrings + `NotImplementedError` bodies.
- [ ] `uv run mypy src/cyclerfinder/data/writeback.py` clean on skeleton.

### Implementation

- [ ] Implement `_highest_passing_level(gates)` helper per plan §3.3.1.
- [ ] Implement `apply_v2_to_entry(entry, report)` per plan §3.3.1 — copies M6a `StabilityReport` into `validation.gates.V2` (fields: `pass`, `max_drift_km`, `n_laps`, `per_lap_drift_km`, `frame_used`). Updates `validation.level` to highest passing.
- [ ] Implement `apply_v0_v1_to_entry(entry, v0_result, v1_result)` — copies V0 (`constraints_satisfied`) and V1 (`crosscheck_cycler` aggregate `max_diff_mps`) into the gates block.
- [ ] Implement `apply_v3_to_entry(entry, report)` — when M6b populates `report.per_lap_dv` / `total_tcm_dv`, copy into `validation.gates.V3` and `metrics.horizon_tcm_dv_mps`, `metrics.horizon_laps`, `metrics.horizon_years`. Until M6b lands, M7 tests use a synthetic report with non-zero TCM fields to exercise the path.
- [ ] Implement `record_rediscovery(entry, run_id, cell_id, date)` per plan §3.3.2 — idempotent on `(run_id, cell_id)`; never overwrites `first_published` / `priority_date`.
- [ ] Implement `register_discovery(entry_skeleton, run_id, cell_id, date, finder_version)` per plan §1.3 — produces `source="this-project"`, `our_status="candidate-novel"`, `priority_date=date`, `first_published=None`, `discovery_run={...}` entry.
- [ ] Implement `serialise_entry_yaml(entry)` per plan §3.3.3 — `yaml.safe_dump` with `sort_keys=False`, custom key order matching existing entries.

### Tests

- [ ] Test `test_v2_writeback_populates_validation_block` (**M7 GATE — M6a integration**) — synthetic `StabilityReport(stable=True, max_drift_km=12345.6, ...)`; assert `new_entry.validation["gates"]["V2"]["pass"] is True`, `max_drift_km == 12345.6`.
- [ ] Test `test_v2_writeback_promotes_validation_level` — V2 pass with V0/V1 also passing → `validation["level"] == "V2"`.
- [ ] Test `test_v2_writeback_does_not_mutate_input` — original `entry` unchanged after `apply_v2_to_entry` (frozen invariant).
- [ ] Test `test_record_rediscovery_idempotent` — recording same `(run_id, cell_id)` twice → one entry.
- [ ] Test `test_record_rediscovery_preserves_first_published` — `first_published`, `priority_date` unmodified.
- [ ] Test `test_register_discovery_sets_candidate_novel` — `our_status="candidate-novel"`, `source="this-project"`, `first_published is None`.
- [ ] Test `test_serialise_entry_yaml_round_trip` — round-trip an existing constructible entry; parse(serialise(entry)) == entry (field-equal; comments not preserved).

### Module-level checks

- [ ] Confirm `uv run pytest tests/data/test_writeback.py` green.
- [ ] Confirm `data/seed_cyclers.yaml` is NOT modified by any test (use `tmp_path` only).
- [ ] Confirm `uv run mypy src/cyclerfinder/data/writeback.py` clean.

## `data/discover.py`

### Skeleton

- [ ] Create `src/cyclerfinder/data/discover.py` with module docstring referencing spec §13.6, §13.8 + §14 auto-validation.
- [ ] Stub `discover` generator function with full docstring + `NotImplementedError` body.
- [ ] Stub `_auto_validate(result, catalog, ephem, *, enable_v3)` helper.
- [ ] `uv run mypy src/cyclerfinder/data/discover.py` clean on skeleton.

### Implementation

- [ ] Implement `_auto_validate` per plan §3.5: V0 from `OptimisationResult.constraints_satisfied`; V1 from `crosscheck_cycler`; V2 from `verify_long_term_stability` (M6a); V3 gated by `enable_v3` flag (default False — M6b's stub).
- [ ] Implement `discover` generator per plan §3.5.1 pseudocode: enumerate cells, check ledger, claim, optimise, signature, match, auto-validate, record, yield.

### Tests

- [ ] Test `test_discover_em_k2_yields_known_for_2syn` — full pipeline: `discover(("E","M"), 2, 7.0, tmp_path/"ledger.jsonl")` yields a tuple whose `match_result.outcome == "known"` and `entry.id == "s1l1-2syn-em-cpom"`.
- [ ] Test `test_discover_writes_ledger_for_every_cell_attempted` — after run, every enumerated cell has a ledger entry.
- [ ] Test `test_discover_resumes_from_existing_ledger` — pre-populate ledger with 3 cells as "solved"; second `discover` run skips those (spy on `optimise_cell_idealized` call count).
- [ ] Test `test_discover_records_signature_hash_in_ledger` — every `solved` ledger entry has non-empty `signature_hashes`.
- [ ] Test `test_discover_skips_v3_when_disabled` — `enable_v3=False` does not raise `NotImplementedError`.

### Module-level checks

- [ ] Confirm `uv run pytest tests/data/test_discover.py` green.
- [ ] Confirm `uv run mypy src/cyclerfinder/data/discover.py` clean.

## Catalogue-rediscovery tagging gate

- [ ] Create `tests/data/test_catalogue_rediscovery_tagging.py`.
- [ ] **Gate test** `test_rediscovered_2syn_em_tagged_known` (**M7 BINDING GATE — spec §8 M7**):
  - [ ] `result = find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)[0]`.
  - [ ] `sig = canonical_signature(result.best_cycler, model_assumption="circular-coplanar")`.
  - [ ] `match_result = match(sig, load_catalog())`.
  - [ ] Assert `match_result.outcome == "known"`.
  - [ ] Assert `match_result.entry.id == "s1l1-2syn-em-cpom"`.
  - [ ] Assert `match_result.entry.priority_date == "2002-08-05"`.
- [ ] Parametrised `test_catalogue_entry_tagged_known[<entry_id>]` over constructible entries (reuse `tests/_catalogue_loader.load_constructible_entries`):
  - [ ] Build cell from entry per `tests/test_catalogue_rediscovery._build_cell_from_entry` (DO NOT MODIFY that file; import the helper).
  - [ ] Run `optimise_cell_idealized` with `vinf_cap = max(target_vinfs) + 2.5`.
  - [ ] Compute signature, match against `load_catalog()`.
  - [ ] Assert `match_result.outcome == "known"` AND `match_result.entry.id == entry.id`.
  - [ ] Honour `EXPECTED_SKIPS` from `tests/test_catalogue_rediscovery`.
  - [ ] Mark `pytest.mark.slow` (the same marker M5's rediscovery test uses).

## Sanity: M0–M6a tests still pass

- [ ] `uv run pytest tests/core/` — M0/M1/M2/M3 + M6a frames green.
- [ ] `uv run pytest tests/model/` — M3/M4 green.
- [ ] `uv run pytest tests/search/` — M2/M4/M5/M6-slice green.
- [ ] `uv run pytest tests/verify/test_propagate.py` — M6a green (M7 only adds `crosscheck.py`, doesn't touch propagate).
- [ ] `uv run pytest tests/test_catalogue_rediscovery.py` — M5-era rediscovery test green (M7 does NOT modify it; `tests/_catalogue_loader.py` is also untouched).

## Local green

- [ ] `uv run pytest` → green (M0–M7 all passing).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean.
- [ ] `uv run mypy src tests` → clean under strict mode.

## CI

- [ ] Commit: `m7: catalogue loader + canonical signature matcher + novelty tagging + ledger + V1 cross-check`.
- [ ] Push; confirm GitHub Actions runs and all checks pass.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M7 status → `completed`; M8 row → `planned`.
- [ ] Update `docs/overview.md` §2 deferred-decisions table:
  - Remove "Ledger backend — SQLite vs JSONL (decided in M4/M7 when ledger is built)" from deferred.
  - Add to kept-decisions: "Ledger backend — JSONL (M7). SQLite reserved for a future scaling pass per plan §5 risk #11."
- [ ] (Optional) Update `data/README.md` to add a "Writeback" section describing the `apply_v2_to_entry` / `record_rediscovery` / `register_discovery` paths and the explicit rule that CI tests never modify `data/seed_cyclers.yaml`.
- [ ] Append a `## Hand-off to M8` section to this `todo.md` (below).

## Hand-off to M8

*(Filled in at M7 closeout. Placeholder structure below; replace bullets with measured values from the M7 test pass.)*

- **Catalogue-rediscovery tagging count over the constructible-entries set:**
  - `known`: `<TBD>` / `<TOTAL>`.
  - `probable-match-NEEDS-HUMAN`: `<TBD>` (list ids and matching candidate signatures below; operator resolves each).
  - `novel`: `<TBD>` (these are catalogue entries M5 produced a different cycler for — escalate per `tests/test_catalogue_rediscovery.EXPECTED_SKIPS` pattern).
- **`TAU_NEAR` value used:** 0.5 default; `<TBD: changed | unchanged>`. Sensitivity-test result at `{0.3, 0.5, 1.0}`: `<TBD counts per setting>`.
- **`sense` derivation ambiguities (plan §5 risk #2):** `<TBD: no | yes — entries affected: …>`. If ambiguity surfaced, M8 should extend `Cycler` or `OptimisationResult` with an explicit `sense` field; this is a known follow-up.
- **`period.years` exclusion from hash (plan §5 risk #15):** `<TBD: surfaced any matching surprises? | no>`.
- **V1 cross-check pass rate across constructible-entries set:** `<TBD>` / `<TOTAL>` legs passed at `max_diff_mps < 1.0e-3`. If 100%, M8 reporter can unconditionally show V1 status; if not, document which legs failed and why.
- **V3 status:** disabled (`enable_v3=False`) until M6b's `optimise_cell_ephemeris` body lands. M6b's plan needs to flip the flag once its gate passes.
- **Ledger size at end of `discover` run for `(("E","M"), k_synodic=2, vinf_cap=7.0, L_max=4)`:** `<TBD>` entries. Informs M8 parallel-runner sizing — typical cell-throughput per worker = ledger-entries-per-second from this measurement.
- **JSONL atomic-claim race:** `<TBD: fired in CI? | no>`. If it fired, M8 may swap to a SQLite backend per plan §6 deferred-decision.
- **Confirmation `data/seed_cyclers.yaml` was NOT modified by any M7 test:** `<TBD: confirmed | violation found>`.
- **Locked `CanonicalSignature` shape** (M8 reporter inherits):
  - `bodies: tuple[str, ...]`
  - `sequence_canonical: str`
  - `sense: str`
  - `period_k: int`
  - `period_years: float` (informational; NOT in hash)
  - `vinf_multiset_binned: tuple[tuple[str, float], ...]`
  - `leg_elements_multiset_binned: tuple[tuple[float, float], ...]`
  - `model_assumption: str` (pool filter; NOT in hash)
  - `hash: str` (sha1, prefixed `"sha1:"`)
- **Locked `MatchResult` shape:**
  - `outcome: Literal["known", "probable-match-NEEDS-HUMAN", "novel"]`
  - `entry: CatalogueEntry | None`
  - `distance: float | None`
- **Locked `LedgerEntry` shape** (M8 parallel runner inherits):
  - `cell_id: str`
  - `status: Literal["pending", "pruned", "searched", "solved", "failed"]`
  - `n_solutions: int`
  - `best_dv_kms: float | None`
  - `signature_hashes: tuple[str, ...]`
  - `validation_level: str | None` (`"V0"..."V5"`)
  - `t_done: str` (ISO-8601)
  - `host: str`
- **Locked `LambertCrosscheckResult` shape** (informs the M8 reporter's V1 render):
  - `leg_index: int`
  - `mine_v1_kms: tuple[float, float, float]`
  - `lamberthub_izzo_v1_kms: tuple[float, float, float]`
  - `lamberthub_gooding_v1_kms: tuple[float, float, float]`
  - `max_diff_mps: float`
  - `pass_: bool`
- **API surfaces M8 will consume immediately:**
  - `from cyclerfinder.data.catalog import canonical_signature, load_catalog, match, MatchResult, CanonicalSignature, CatalogueEntry, Catalog, TAU_NEAR`
  - `from cyclerfinder.data.ledger import Ledger, LedgerEntry, LedgerStatus, LedgerError, LedgerLoader`
  - `from cyclerfinder.data.writeback import apply_v0_v1_to_entry, apply_v2_to_entry, apply_v3_to_entry, record_rediscovery, register_discovery, serialise_entry_yaml`
  - `from cyclerfinder.data.discover import discover`
  - `from cyclerfinder.verify.crosscheck import crosscheck_leg, crosscheck_cycler, LambertCrosscheckResult`
- **Spec ambiguities encountered during M7 implementation:**
  - **`sense` field absence from `Cycler` / `OptimisationResult`.** M7 derives heuristically; M8 should consider extending the dataclasses (additive).
  - **`period.years` in hash vs not in hash.** Spec §16.1 example shows it in `signature_fields`; M7 excludes it from the hash subset per plan §5 risk #15 (k is exact; years is derived). If a later spec revision contradicts this, the change is one line.
  - **Pool-filter granularity.** Spec §16.3 says `filter(bodies=, k=)`; M7 adds `model_assumption=` per spec §12.2 / §16.2 final paragraph. This is consistent with the spec but more aggressive than the spec §16.3 pseudocode literally reads. Documented.
- **Recommended M8 first steps:**
  1. Read this hand-off + spec §11 (definition-of-done v1) + spec §13.7 (prioritisation) + spec §15 (dissemination).
  2. Build the CLI surface — `cyclerfinder find --bodies E,M --k 2 --vinf-cap 7` invokes `discover` and produces a per-result report.
  3. Run the VEM campaign — `discover(("V","E","M"), k_synodic=3, vinf_cap=...)` with the §13.7 prioritisation order (VEM first); M7's matcher handles the matching once VEM family-seed entries are individuated.
  4. Build the visualisations (`viz/plots.py` — spec §4) using the locked `CanonicalSignature` and `MatchResult` shapes.

M8's first task is the M8 plan doc itself; this hand-off is the input it consumes.
