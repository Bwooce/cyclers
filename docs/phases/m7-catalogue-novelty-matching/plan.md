# M7 — Catalogue loader, canonical signature matching, novelty scoring

**Spec reference:** spec.md §4 (architecture — `data/catalog.py` slot), §5 step 8 (novelty: compare closing solutions to the known-cycler catalogue and flag matches vs candidates-for-new), §6 (`find_cyclers` returns ranked cyclers — M7 consumes the result), §8 (M7 milestone: "catalogue/novelty — `catalog`; flag known vs candidate-new. *Gate: correctly tags the rediscovered E–M cyclers as known.*"), §12.2 (representation framework — `model_assumption` partitions the matching pool), §13.5 (canonical signature + persistent catalogue as the deduplication-and-novelty machinery for "running indefinitely"), §13.6 (work queue + ledger for resumability — **the ledger is M7's deliverable, deferred from M4/M5**), §13.7 (prioritisation), §13.8 (ledger schema: `cell_id, status, n_solutions, best_dv, signatures[], validation_level, t_done, host` — the M7 persistence contract), §14 (V0–V5 validation gauntlet — M7 reads `StabilityReport` from M6a as the V2 result and writes it into the catalogue `validation.gates` block), **§16 (catalogue schema, identity matching, attribution — the central design document for M7)**, §16.1 (catalogue record JSON schema + v2 additive fields), §16.2 (canonical signature definition — lexicographically-minimal rotation, V∞ multiset binned to 0.05 km/s, leg (a,e) multiset binned to 0.01 AU / 0.01), §16.3 (matcher: exact / probable-match / novel), §16.4 (attribution + literature ingest), §16.5 (discovery workflow — `source: this-project` lifecycle).

**Purpose:** stand up the **catalogue + identity** layer that turns one-shot `OptimisationResult` / `StabilityReport` outputs into a persistent, deduplicated, attribution-aware record. M0–M5 produce candidates; M6a/M6b verify them; M7 is the **identity arbiter** — it answers "is this cycler known?", "is this cycler novel?", "what attribution does it inherit?". Without M7, every finder run rediscovers the same families and emits no progress signal. With M7, every closing trajectory is reduced to a canonical signature, matched against a 219-row seeded catalogue, and either tagged `known-reproduction` (with inherited citation), flagged `probable-match-NEEDS-HUMAN`, or queued as `candidate-novel` for V4 GMAT review. M7 is also the home for the spec §13.6/§13.8 **append-only ledger** — the persistence layer that makes a long finder run resumable, non-redundant, and parallel-safe.

**Resolved nomenclature** (deferred-decision callout):

A previous conversation referenced "TCM ΔV budgeting" as M7. The canonical sources contradict that: `docs/overview.md` §4 (line 92: *"M7 — Catalogue loader, canonical signature matching, novelty scoring — planned"*) and `docs/spec.md` §8 (line 151: *"M7 — catalogue/novelty: `catalog`; flag known vs candidate-new. Gate: correctly tags the rediscovered E–M cyclers as known."*) both unambiguously assign **catalogue + novelty matching** to M7. TCM ΔV budgeting is **M6b** per `docs/overview.md` §4 (line 91: *"M6b — Ephemeris-mode TCM minimisation over 3–5 lap horizon"*). The M6a plan §1.5 (lines 72–77) lists *catalogue ingest of stability results* and *Lambert cross-check* as deferred to M7, which this plan absorbs as **sub-deliverables under the M7 catalogue/novelty umbrella**:

- **Primary deliverable** (spec §8 M7 gate): `data/catalog.py` — catalogue loader, canonical signature matcher, novelty tagging. Gate test: every catalogue entry the M5 optimiser rediscovers is tagged `known-reproduction` with attribution inherited.
- **Sub-deliverable A** (M6a hand-off, fixing the "M7 writes StabilityReport into YAML" debt): catalogue-side V2 gate writeback. The catalogue record's `validation.gates.V2.max_drift_km` field (spec §16.1) is populated by an M7 batch-validate runner that consumes `StabilityReport` from `verify_long_term_stability`. No changes to M6a's source code; M7 supplies the runner that drives it across the catalogue.
- **Sub-deliverable B** (M6a non-goal §1.5 row 5): `verify/crosscheck.py` — the spec §14 V1 single-leg Lambert-vs-lamberthub cross-check. Single-leg, independent of multi-lap propagation; M7's V0/V1 batch runner uses it. Stays out of M8 because V1 must run before V2 in the gauntlet (cheapest-first).
- **Sub-deliverable C** (spec §13.6 + §13.8, deferred from M4): the append-only **ledger** — a SQLite or JSONL store of `(cell_id, status, signature_hash, best_score, validation_level, t_done, host)` per cell. M7's `find_cyclers` is replaced by a ledger-backed `discover` runner that skips cells already in the ledger, claims pending cells atomically, and writes hits/signatures back. M5's `find_cyclers` stays as the single-process API; M7's `discover` wraps it with persistence.

The TCM-budget integration is the **catalogue read-side** of M6b: when M6b's `optimise_cell_ephemeris` populates `StabilityReport.per_lap_dv` / `total_tcm_dv`, M7's catalogue writer copies those into the entry's `metrics.horizon_tcm_dv_mps` and `validation.gates.V3.horizon_tcm_mps` fields. The *budget itself* is M6b's; M7 only persists the result. This is the clean responsibility split.

**Gate (definition of done):**

1. `tests/data/test_catalogue_rediscovery_tagging.py::test_rediscovered_2syn_em_tagged_known` (**M7 BINDING GATE — spec §8 M7 anchor**) asserts that when `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)` is run end-to-end, the best result is reduced to a canonical signature via `canonical_signature(result.best_cycler)`, matched against the loaded seed catalogue via `match(candidate, catalog)`, and returns `("known", entry)` where `entry.id == "s1l1-2syn-em-cpom"`. The match outcome must inherit `entry.first_published` as attribution; the entry's `rediscoveries[]` list gets appended (in-memory only — no YAML write during the gate test). **The catalogue entry is the fixture, not the tolerance** — no per-test V∞ widening, no signature-distance loosening.
2. `tests/data/test_canonical_signature.py::test_signature_rotation_invariant` (**spec §16.2 BINDING GATE**) asserts that for a fixed cycler, `canonical_signature(cycler)` is bitwise identical when the cycler is provided with its encounter sequence rotated by 1, 2, ..., len(encounters)-1 positions. Two `Cycler` instances whose encounter sequence is a cyclic rotation of one another MUST produce identical `signature.hash`. This is the *identity* property the entire matcher depends on.
3. `tests/data/test_canonical_signature.py::test_signature_binning_absorbs_noise` (**spec §16.2 BINDING GATE**) asserts that perturbing every encounter V∞ magnitude by ±0.02 km/s (below the 0.05 km/s bin width) and every leg `(a, e)` by ±0.005 AU / ±0.005 (below the 0.01 bin width) leaves `signature.hash` unchanged. Perturbing by 0.06 km/s on V∞ OR 0.02 AU on `a` changes the hash. This is the binding noise-tolerance property — too tight and accidental re-derivations slip past the matcher as "novel"; too loose and physically distinct cyclers collide.
4. `tests/data/test_match.py::test_match_aldrin_classic_returns_known` asserts that `match(canonical_signature(aldrin_classic_cycler), load_catalog())` returns `("known", entry)` where `entry.id == "aldrin-classic-em-k1-outbound"` and `entry.priority_date == "1985-10-28"` (the earliest published date — attribution rule per spec §16.4). The Aldrin classic is the only literature anchor in the seed catalogue with an unambiguous 1985 priority date; mistakes here would mis-attribute discovery credit.
5. `tests/data/test_match.py::test_match_novel_synthetic_returns_novel` asserts that a synthetic cycler with a fabricated V∞ signature ([("E", 4.21), ("M", 7.83)] — not within 0.05 km/s of any catalogue entry) returns `("novel", None)`. Confirms the matcher does not false-positive against the seed catalogue.
6. `tests/data/test_match.py::test_match_partition_by_model_assumption` (**spec §12.2 BINDING GATE**) asserts that a `circular-coplanar` candidate is never matched against `cr3bp` catalogue entries. The Arenstorf figure-8 (`arenstorf-em-figure8-1963`, `model_assumption: cr3bp`) lives in a non-comparable representation; if a circular-coplanar V∞ signature happened to numerically collide with its Jacobi-constant-derived signature fields, the match must still return `("novel", None)` (assuming no circular-coplanar collision). Pool partitioning is the spec §12.2 / §16.2 contract.
7. `tests/data/test_catalogue_loader.py::test_load_catalog_returns_all_219_entries` asserts the loader reads every entry in `data/seed_cyclers.yaml` (currently 219), distinguishes per-entry `source`, computes the canonical signature for every literature entry that has the requisite fields (bodies, sequence_canonical, period.k, vinf_kms_at_encounters with all non-null), and surfaces both `by_id` and `by_hash` indices.
8. `tests/verify/test_crosscheck.py::test_v1_lambert_crosscheck_aldrin` asserts that for an Aldrin-cycler leg, `crosscheck_leg(leg, ephem)` runs both the in-house `lambert` and `lamberthub.izzo` / `lamberthub.gooding` solvers and returns `max_diff_mps < 1.0e-3` per spec §14 V1. Single-leg; no multi-lap propagation. This is the V1 gate machinery the V0/V1 batch runner consumes.
9. `tests/data/test_ledger.py::test_ledger_round_trip` asserts that recording a cell → reading it back via `Ledger.has(cell_id)` returns `True`; recording two cells with the same `cell_id` raises `LedgerError` (idempotency / no-double-record); and that the JSONL persistence layer survives a process restart (`Ledger(path).has(...)` after `Ledger(path).record(...)` and a re-open).
10. `tests/data/test_writeback.py::test_v2_writeback_populates_validation_block` asserts that running the M7 batch-validate runner on a single cycler with a known `StabilityReport` produces an in-memory updated catalogue entry whose `validation.gates.V2` block has `pass: bool`, `max_drift_km: float`, and `n_laps: int` matching the report fields. The runner does NOT modify `data/seed_cyclers.yaml` during the test — writeback is to an in-memory copy or temp file, validated then discarded. The schema contract is what matters; on-disk YAML mutation belongs to operator-driven runs, not CI.
11. `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green on the M7 commit.

The 0.05 km/s V∞ bin and 0.01 AU / 0.01 (a, e) bin are the **binding signature tolerances** per spec §16.2. They are NOT test-tunable, NOT widened per-gate, and NOT changed by M7 — they are spec-fixed. See §4.3 below for why these specific values, and what they buy.

---

## 1. What this milestone delivers

Three new source modules, one new subpackage, two new files in `verify/`, one new test subpackage. M7 is **additive** for every existing source file except:

- `src/cyclerfinder/verify/__init__.py` — gains one new submodule re-export (`from . import crosscheck`).

No edits to M0–M6a source bodies. The catalogue YAML at `data/seed_cyclers.yaml` is **read-only from M7 source code** under CI; the on-disk writeback path exists and is exercised by integration tests writing to `tmp_path`, never to the canonical file.

### 1.1 `src/cyclerfinder/data/__init__.py` (NEW) and `src/cyclerfinder/data/catalog.py` (NEW)

First occupants of the `data/` subpackage (spec §4 architecture line: `data/catalog.py # known cyclers (Aldrin, McConaghy 2-synodic, Russell–Ocampo set) + novelty check`). Public surface, in dependency order:

- `CatalogueEntry` — frozen dataclass; the in-memory projection of a `data/seed_cyclers.yaml` row. Carries the full record (not just the rediscovery subset that `tests/_catalogue_loader.py` carries today) — `bodies`, `sequence_canonical`, `sense`, `period_k`, `period_years`, `vinf_kms_at_encounters`, `legs`, `orbit_elements`, `model_assumption`, `primary`, `trajectory_regime`, `first_published`, `priority_date`, `source`, `our_status`, plus optional `signature_hash` (computed on load), `validation` block (read-through), `discovery_run` block (for `source: this-project` entries).
- `CanonicalSignature` — frozen dataclass per spec §16.2. Fields: `bodies` (sorted tuple), `sequence_canonical` (lex-min cyclic rotation string), `sense` (str), `period_k` (int), `period_years` (float, NOT in hash — informational only), `vinf_multiset_binned` (sorted tuple of `(body, vinf_binned)`), `leg_elements_multiset_binned` (sorted tuple of `(a_au_binned, e_binned)`), `model_assumption` (str — pool-filter key, NOT in hash per spec §16.2 final paragraph), `hash` (str — sha1 of the canonical-JSON of the in-hash fields).
- `canonical_signature(cycler, *, model_assumption=...) -> CanonicalSignature` — computes the signature from a `Cycler` instance. The cyclic-rotation canonicalisation and the multiset binning are done here per spec §16.2 algorithm. `model_assumption` is a caller-supplied tag (the matcher pool filter); M5's `find_cyclers` defaults it to `"circular-coplanar"`; M6b ephemeris-mode results stay `"circular-coplanar"` because the signature is invariant to ephemeris realisation (only the launch window changes).
- `Catalog` — class wrapping the loaded `list[CatalogueEntry]`. Exposes `by_id: dict[str, CatalogueEntry]`, `by_hash: dict[str, CatalogueEntry]`, `entries: tuple[CatalogueEntry, ...]`, and `filter(*, bodies=None, k=None, model_assumption=None) -> tuple[CatalogueEntry, ...]` (the coarse prefilter for §16.3 matcher pool partitioning).
- `load_catalog(path=CATALOGUE_PATH) -> Catalog` — reads the YAML, builds `CatalogueEntry` instances, computes signatures for entries that have the full set of signature-required fields (skips citation-only / family-seed entries — they remain in the catalogue but have `signature_hash = None` and never participate in matching).
- `match(candidate, catalog) -> MatchResult` per spec §16.3 — exact / probable / novel. `MatchResult` is a tagged-union dataclass with fields `(outcome: Literal["known", "probable-match-NEEDS-HUMAN", "novel"], entry: CatalogueEntry | None, distance: float | None)`.
- `signature_distance(sig_a, sig_b) -> float` — weighted L1 over `(period_years, vinf_multiset_binned, leg_elements_multiset_binned)` per spec §16.3 pseudocode. Used by the probable-match stage.
- `TAU_NEAR: Final[float] = 0.5` — probable-match distance threshold. See §4.4 for derivation.

Module-internal helpers (private):

- `_lex_min_rotation(sequence_str)` — returns the lexicographically-minimal cyclic rotation of a dash-joined body sequence per spec §16.2.
- `_bin_vinf(vinf_kms)` — rounds V∞ to the 0.05 km/s bin (spec §16.2).
- `_bin_a_au(a_au)` and `_bin_e(e)` — round leg `a, e` to 0.01 AU / 0.01 (spec §16.2).
- `_canonical_json(d)` — sorted-key JSON serialiser used for the sha1 input.
- `_extract_signature_fields(entry)` — for catalogue rows, build the signature-input dict from the YAML fields (`bodies`, `sequence_canonical`, `period.k`, `vinf_kms_at_encounters`, `orbit_elements.a_au` + `orbit_elements.e`). Family-seed / citation-only rows return `None`.

### 1.2 `src/cyclerfinder/data/ledger.py` (NEW)

The spec §13.6 + §13.8 append-only ledger. Public surface:

- `LedgerStatus` — `Literal["pending", "pruned", "searched", "solved", "failed"]` per spec §13.8.
- `LedgerEntry` — frozen dataclass: `cell_id: str`, `status: LedgerStatus`, `n_solutions: int`, `best_dv_kms: float | None`, `signature_hashes: tuple[str, ...]`, `validation_level: str | None` (`"V0"`..`"V5"`), `t_done: str` (ISO-8601), `host: str`.
- `LedgerError` — module exception; raised on duplicate-`cell_id` writes (no double-record), on schema-mismatched reads, on persistence path errors.
- `Ledger` — class wrapping a JSONL file (sole backend in M7 per the deferred-decisions table in `docs/overview.md` §2; the SQLite option is reserved for a future scaling pass). Methods: `record(entry)`, `has(cell_id) -> bool`, `get(cell_id) -> LedgerEntry`, `claim(cell_id, host) -> bool` (atomic — returns `False` if another worker already claimed), `iter_pending() -> Iterator[str]`, `__len__`.
- `LedgerLoader` — read-only loader for analysis (e.g. M8 reporter consuming ledger to produce a catalogue summary).

Module-internal helpers:

- `_atomic_append(path, line)` — POSIX-atomic `O_APPEND` write of a single JSONL line. The atomicity is what makes the ledger parallel-safe; multiple workers can append concurrently without inter-line interleaving as long as each writes a single `\n`-terminated line below the kernel's atomic-write threshold (4 KiB on linux, well above our per-line budget).
- `_parse_line(line)` — strict JSONL → `LedgerEntry` parser; raises `LedgerError` on schema mismatch.

### 1.3 `src/cyclerfinder/data/writeback.py` (NEW)

The catalogue-side writeback path that consumes `StabilityReport` (M6a) and the M5/M6b optimisation results, and produces an updated `CatalogueEntry` for the matched / discovered cycler. Public surface:

- `apply_v0_v1_to_entry(entry, v0_result, v1_result) -> CatalogueEntry` — copies internal-consistency + Lambert-crosscheck fields into the entry's `validation.gates.V0`, `validation.gates.V1` blocks. Returns a new frozen entry (immutability preserved).
- `apply_v2_to_entry(entry, report) -> CatalogueEntry` — copies the M6a `StabilityReport` fields into `validation.gates.V2`: `pass: report.stable`, `max_drift_km: report.max_drift_km`, `n_laps: report.n_laps_propagated`, `per_lap_drift_km: report.per_lap_drift_km` (as a list). Returns a new entry.
- `apply_v3_to_entry(entry, report) -> CatalogueEntry` — copies the M6b TCM-budget fields (when M6b populates them) into `validation.gates.V3` and the entry's `metrics.horizon_tcm_dv_mps`, `metrics.horizon_laps`, `metrics.horizon_years` fields. **M7 reads** the M6b output; the optimisation itself is M6b's responsibility. M7 only persists.
- `record_rediscovery(entry, run_id, cell_id, date) -> CatalogueEntry` — appends a `rediscoveries[]` entry per spec §16.4 (retroactive correction path). Idempotent: a `(run_id, cell_id)` tuple already present is not re-appended.
- `register_discovery(entry_skeleton, run_id, cell_id, date, finder_version) -> CatalogueEntry` — for `source: this-project` candidates per spec §16.5. Produces a `CatalogueEntry` with `source="this-project"`, `our_status="candidate-novel"`, `priority_date=date`, `first_published=None`, `discovery_run={...}`. Catalogue writer is the consumer.
- `serialise_entry_yaml(entry) -> str` — produce the YAML representation of a `CatalogueEntry` for writing back to `data/seed_cyclers.yaml`. Round-trip stable: `serialise_entry_yaml(load_catalog().by_id["s1l1-2syn-em-cpom"])` == the on-disk YAML block, modulo comment formatting.

### 1.4 `src/cyclerfinder/verify/crosscheck.py` (NEW — sub-deliverable B)

The spec §14 V1 single-leg Lambert cross-check, deferred from M6a §1.5 row 5. Public surface:

- `LambertCrosscheckResult` — frozen dataclass: `leg_index: int`, `mine_v1_kms: tuple[float, float, float]`, `lamberthub_izzo_v1_kms: tuple[float, float, float]`, `lamberthub_gooding_v1_kms: tuple[float, float, float]`, `max_diff_mps: float`, `pass: bool` (`max_diff_mps < 1.0e-3`).
- `crosscheck_leg(leg, ephem) -> LambertCrosscheckResult` — re-solve a single `Leg` with both `lamberthub.izzo` and `lamberthub.gooding`, compare to the in-house `lambert` result, return the worst per-component velocity difference in m/s.
- `crosscheck_cycler(cycler, ephem) -> tuple[LambertCrosscheckResult, ...]` — runs `crosscheck_leg` across every leg in the cycler. The aggregated `pass` is `all(r.pass for r in results)`.

Module-internal helpers (private):

- `_leg_endpoints(leg, cycler, ephem) -> tuple[Vec3, Vec3]` — extract the leg's heliocentric endpoint positions from the cycler's encounter records.

### 1.5 `src/cyclerfinder/data/discover.py` (NEW — ledger-backed runner)

The persistent, resumable, parallel-safe replacement for the in-process `find_cyclers`. Public surface:

- `discover(bodies, k_synodic, vinf_cap, ledger_path, *, ephem=None, l_max=4, n_max=0, branch_set=("single",), n_starts=5, seed=0, use_de=True, rp_factors=None, host=...) -> Iterator[OptimisationResult]` — wraps `find_cyclers` with ledger persistence. For each cell yielded by the M4 enumerator: check ledger; if already done, skip; if pending, claim atomically; else mark pending, run `optimise_cell_idealized`, compute signature, match against catalogue, write a `LedgerEntry` with `validation_level` set per the auto-pipeline result, and yield the result.

Module-internal helpers:

- `_auto_validate(result, catalog, ephem)` — runs V0 (internal-consistency, already gated by `OptimisationResult.constraints_satisfied`), V1 (`crosscheck_cycler`), V2 (`verify_long_term_stability` from M6a), and tags the result with the highest passing level. V3 (ephemeris-mode TCM) is gated by `optimise_cell_ephemeris` availability; M7 leaves V3 disabled until M6b lands, with a feature-flag `enable_v3: bool = False` (the call would fail with `NotImplementedError` from M5's stub).

### 1.6 Test files

- `tests/data/__init__.py` — empty new test package.
- `tests/data/test_canonical_signature.py` — spec §16.2 binding tests (rotation invariance, binning), unit tests on the helpers.
- `tests/data/test_match.py` — spec §16.3 matcher tests including the partition-by-`model_assumption` invariant.
- `tests/data/test_catalogue_loader.py` — loader sanity, all-219 enumeration, signature index population.
- `tests/data/test_catalogue_rediscovery_tagging.py` — the M7 binding gate test, plus per-entry integration tests parametrised over the constructible-entries set (reuse `tests/_catalogue_loader.py` from M5 — keep that test infrastructure).
- `tests/data/test_ledger.py` — round-trip, idempotency, atomic-claim, restart-survives.
- `tests/data/test_writeback.py` — V0/V1/V2/V3 field copying, `record_rediscovery` idempotency, `serialise_entry_yaml` round-trip.
- `tests/data/test_discover.py` — ledger-backed runner end-to-end on a tiny enumeration; resume-from-checkpoint test.
- `tests/verify/test_crosscheck.py` — spec §14 V1 gate; Aldrin and 2-syn E-M test legs.

### 1.7 Explicit non-goals (M7 boundaries)

These belong to M8, the V4 GMAT stretch, or are deliberate carve-outs from spec §16. **Do not stub or partially implement them in M7** beyond what the locked interfaces above require:

| Out of M7 | Where it lands | Why deferred |
|---|---|---|
| M6b TCM budget *computation* (per-lap ΔV, total horizon TCM) | **M6b** | M7 reads the `StabilityReport.per_lap_dv` / `total_tcm_dv` fields produced by `optimise_cell_ephemeris`; the optimisation itself is M6b's body. M7's `apply_v3_to_entry` will populate the catalogue's V3 block when those fields are non-zero; until M6b lands they stay zero per M6a's lock. |
| Phase-match epoch-picking integration in `discover` | **M6 slice (already shipped) + M6b** | `phase_match.find_real_windows` exists; M7's `discover` does not call it. The launch-window list is a *derived view* per spec §12.2; the catalogue does not store calendar dates. M6b's TCM optimisation consumes the window dates internally; the catalogue records the result, not the dates. |
| GMAT V4 bridge (`verify/gmat_bridge.py`) | **Stretch** per spec §7 | The V4 gate runs externally; M7's `validation.gates.V4` block records the *result* (tool, date, pass), not the run. The external-tool integration is a separate spec §7 stretch deliverable. |
| Multi-body VEM cyclers | **M8** | The signature canonicalisation and matcher are body-agnostic in *form* — they accept `["V", "E", "M"]` cyclers — but the seed catalogue's VEM entries are family-seeds (`jones-2017-vem-triple-family`, `wittal-2022-em-cycler-family`) with `null` numeric fields. They contribute nothing to matching today. M8's VEM campaign will produce individual VEM members that M7's matcher will ingest unchanged. |
| Site / `cyclers.space` integration | **Live (deployed)** per `docs/overview.md` §4 | The public site reads `data/seed_cyclers.yaml` directly. M7's writeback path produces YAML the site already understands; no M7-side site code is needed. |
| Automatic literature ingest (web-scrape new papers) | **Out of scope v1** per spec §16.4 | Ingest is *curated* — a human commits a new YAML row with full attribution per spec §16.4. M7 supplies the *matching* machinery once a row exists; it does not source rows. |
| SQLite ledger backend | **Future scaling pass** per `docs/overview.md` §2 | The deferred-decisions table reads "Ledger backend — SQLite vs JSONL (decided in M4/M7 when ledger is built)." M7 ships the JSONL backend; SQLite is a swap-in (`Ledger` is a class — the persistence is one private helper). Decision recorded at M7 closeout. |
| Parallel-worker orchestration (subprocess pool, distributed) | **M8** | The ledger is parallel-*safe* (atomic `O_APPEND`, atomic claim); the M7 deliverable is the safety mechanism. The actual concurrent runner (worker pool, work-distribution policy) is M8's job since it depends on the VEM campaign's compute scale. |
| Cell prioritisation per spec §13.7 (VEM-first, long-period-first) | **M8** | Prioritisation matters only at scale. M7's `discover` walks cells in the M4 deepening order; M8 introduces a priority queue feeding the ledger. The §13.7 ordering is a configuration choice at the runner level, not a property of the catalogue or ledger. |

---

## 2. File tree after M7

```
cyclers/
├── … (M0..M6a layout preserved unchanged)
├── src/cyclerfinder/
│   ├── core/                            # unchanged
│   ├── search/                          # unchanged
│   ├── model/                           # unchanged
│   ├── verify/
│   │   ├── __init__.py                  # M6a — gains one re-export line for crosscheck
│   │   ├── propagate.py                 # M6a — unchanged
│   │   └── crosscheck.py                # NEW (M7 sub-deliverable B)
│   └── data/                            # NEW subpackage (M7)
│       ├── __init__.py                  # NEW
│       ├── catalog.py                   # NEW — loader, signature, matcher
│       ├── ledger.py                    # NEW — JSONL append-only persistence
│       ├── writeback.py                 # NEW — entry mutation helpers
│       └── discover.py                  # NEW — ledger-backed find_cyclers
└── tests/
    ├── verify/
    │   ├── test_propagate.py            # M6a — unchanged
    │   └── test_crosscheck.py           # NEW (M7)
    └── data/                            # NEW
        ├── __init__.py                  # NEW
        ├── test_canonical_signature.py  # NEW (M7 binding gates §16.2)
        ├── test_match.py                # NEW (M7 binding gate §16.3 + §12.2)
        ├── test_catalogue_loader.py     # NEW
        ├── test_catalogue_rediscovery_tagging.py  # NEW (M7 binding gate §8 M7)
        ├── test_ledger.py               # NEW
        ├── test_writeback.py            # NEW
        └── test_discover.py             # NEW
```

Subpackage `viz/` still remains uncreated — M8 territory. `verify/gmat_bridge.py` stays absent per spec §7 stretch deferral.

The one M6a edit is `src/cyclerfinder/verify/__init__.py`: M6a left it as an empty `__init__.py` per the M6a plan §1.2; M7 adds `from . import crosscheck` so the submodule is importable under `mypy --strict` without a path adjustment in test code. No edits to `verify/propagate.py`.

---

## 3. Module designs

This section walks the four major pieces of M7 in dependency order: (3.1) `data/catalog.py` — the signature + matcher + loader; (3.2) `data/ledger.py` — the persistence layer; (3.3) `data/writeback.py` — the catalogue-side mutation helpers; (3.4) `verify/crosscheck.py` — the V1 Lambert cross-check; (3.5) `data/discover.py` — the ledger-backed runner that composes everything. (3.6) imports/dependency graph; (3.7) API summary.

### 3.1 Canonical signature + matcher (`data/catalog.py`) — spec §16.2 + §16.3 binding

The signature is defined in spec §16.2 as **invariant to absolute epoch/phase, to loop-start choice, and to small numerical noise; nothing else.** M7's implementation has four parts: (a) the per-field canonicalisation (loop rotation, multiset sort, binning); (b) the JSON-stable sha1 hash; (c) the pool-filter partition by `model_assumption`; (d) the §16.3 exact / probable / novel state machine.

#### 3.1.1 Sequence canonicalisation — `_lex_min_rotation`

A cyclic sequence like `"E-M-E"` has three rotations: `"E-M-E"`, `"M-E-E"`, `"E-E-M"`. Spec §16.2: pick the **lexicographically-minimal** one. For `"E-M-E"` that is `"E-E-M"`. The algorithm is a single-pass Lyndon-style minimisation over the dash-joined string treated as a list of body codes:

```python
def _lex_min_rotation(sequence_str: str) -> str:
    parts = sequence_str.split("-")
    n = len(parts)
    if n <= 1:
        return sequence_str
    best = parts
    for i in range(1, n):
        rotated = parts[i:] + parts[:i]
        if rotated < best:
            best = rotated
    return "-".join(best)
```

Why `parts < best` and not Booth's O(n) algorithm: the cycler-sequence length is bounded by spec §13.2's `L_max ≈ 4..8`; O(n²) is trivially below the per-test budget and the implementation is auditable. Booth's algorithm becomes worth it at `L_max > 100`, never the cyclerfinder case.

The `sense` field (spec §16.2: "Keep a separate `sense` field so direction variants are distinguishable but groupable") is **not** folded under reversal in M7. Spec §16.2 allows optionally folding under reversal; M7 takes the conservative path and treats `outbound` and `inbound` as distinct signatures so the Aldrin classic outbound (`aldrin-classic-em-k1-outbound`) and inbound (`aldrin-classic-em-k1-inbound`) catalogue entries remain distinguishable until human review collapses them.

#### 3.1.2 V∞ and leg-element binning — `_bin_vinf`, `_bin_a_au`, `_bin_e`

Per spec §16.2: V∞ binned to 0.05 km/s, `a` to 0.01 AU, `e` to 0.01. The binning is rounding to the nearest bin centre via `round(value / bin) * bin`. The numerical-stability concern is `round` on values exactly halfway between bins: Python's `round` uses banker's rounding (round half to even), which is bitwise-stable across calls and acceptable for §16.2's purpose (the spec does not pin the rounding rule).

The V∞ multiset is **(body, vinf_binned)** tuples, **sorted by `(body, vinf_binned)` ascending** so the order is rotation-independent. The leg-element multiset is `(a_binned, e_binned)` tuples, sorted by `(a_binned, e_binned)` ascending. Both are tuples (immutable) so they hash deterministically.

#### 3.1.3 sha1 hash via canonical JSON — `_canonical_json` and `canonical_signature`

Per spec §16.2: `signature = sha1(canonical_json(signature_fields))`. The canonical JSON is the `signature_fields` dict serialised with sorted keys, no whitespace, ASCII-only — this is what `json.dumps(d, sort_keys=True, separators=(",", ":"), ensure_ascii=True)` produces. The sha1 input is the UTF-8-encoded bytes of that string; the output is the 40-char hex digest with a `"sha1:"` prefix per the spec §16.1 example.

`canonical_signature(cycler, *, model_assumption)` builds the input dict from the cycler:

1. `bodies` — sorted unique tuple of `cycler.bodies` (deduped, ordered alphabetically — body set rotation-invariant).
2. `sequence_canonical` — `_lex_min_rotation("-".join(cycler.bodies))`.
3. `sense` — the cycler's direction. M7's `Cycler` does not carry an explicit sense field today; the helper derives it from the cycler's first leg's `from_body` + the leg's `tof_days` relative to the cycler's period (short outbound first → `"outbound"`; long return first → `"inbound"`). Documented in §5 risk #2 as a precision concern; M5's `OptimisationResult` may need to surface the optimiser's chosen sense explicitly in a future refinement.
4. `period.k` — `cell.period_k` from the `OptimisationResult.cell` (M5 hand-off carries this).
5. `period.years` — informational; included in the dict so `canonical_signature` returns the years value for display, but **NOT** in the in-hash subset (the hash is over `bodies`, `sequence_canonical`, `sense`, `period.k`, `vinf_multiset_binned`, `leg_elements_multiset_binned`).
6. `vinf_multiset_binned` — for each encounter in `cycler.encounters`, `(encounter.body, _bin_vinf(||encounter.vinf_in||))`; sorted ascending.
7. `leg_elements_multiset_binned` — for each leg, `(_bin_a_au(orbit_elements_au(encounter.r, leg.v_depart)[0]), _bin_e(orbit_elements_au(...)[1]))`; sorted ascending. `orbit_elements_au` is the M3-shipped helper in `model/cycler.py`.

The `model_assumption` field is passed in by the caller; for M5 `find_cyclers` output the default is `"circular-coplanar"`. The pool-filter behaviour (spec §16.2 final paragraph, §12.2) is **enforced at match time**, not signature time — i.e. signatures across `model_assumption` are computed identically; the matcher refuses to compare them.

#### 3.1.4 Matcher state machine — `match(candidate, catalog)`

Per spec §16.3 pseudocode, three outcomes:

```python
def match(candidate: CanonicalSignature, catalog: Catalog) -> MatchResult:
    # Pool filter: only compare within same model_assumption (spec §12.2).
    pool = catalog.filter(
        model_assumption=candidate.model_assumption,
        bodies=candidate.bodies,
        k=candidate.period_k,
    )
    # 1. Exact (within binning).
    by_hash = {e.signature_hash: e for e in pool if e.signature_hash is not None}
    if candidate.hash in by_hash:
        return MatchResult(outcome="known", entry=by_hash[candidate.hash], distance=0.0)
    # 2. Probable: weighted L1 below TAU_NEAR.
    near = [
        (e, signature_distance(candidate, e.signature))
        for e in pool if e.signature is not None
    ]
    near = [(e, d) for e, d in near if d < TAU_NEAR]
    if near:
        best_entry, best_distance = min(near, key=lambda x: x[1])
        return MatchResult(
            outcome="probable-match-NEEDS-HUMAN",
            entry=best_entry,
            distance=best_distance,
        )
    # 3. Novel.
    return MatchResult(outcome="novel", entry=None, distance=None)
```

**Design notes:**

- **Pool filter is conservative.** `catalog.filter(model_assumption=...)` is mandatory — a circular-coplanar candidate is never compared to a cr3bp entry. `bodies` and `k` are coarse prefilters (spec §16.3 wording: `catalog.filter(bodies=sig.bodies, k=sig.period.k)`); they cut the search space by ~10x on a 219-row catalogue but never produce a false-negative.
- **Distance is monotone with binning.** A candidate whose V∞ differs from an entry's by 0.04 km/s falls into the same bin and produces hash-equality → `"known"`. A difference of 0.06 km/s falls into a different bin; the matcher then computes signature distance and reports `"probable-match"` if within `TAU_NEAR`, else `"novel"`. The bin-vs-distance interplay is what §4.4 calibrates.
- **The matcher does not write anywhere.** The `MatchResult` is consumed by the caller (the `discover` runner) which decides whether to append a `rediscoveries[]` entry, register a `candidate-novel`, etc. M7's discover is the writer; `match` is pure.

### 3.2 Append-only ledger (`data/ledger.py`) — spec §13.6 + §13.8 binding

Per spec §13.6 the ledger makes the finder run "resumable, non-redundant, and parallel". The atomic-claim invariant (spec §13.8: "Workers pull cells from the queue; the ledger coordinates") is what makes parallel safe. M7 implements JSONL backend.

#### 3.2.1 JSONL format

One JSON object per line. Schema:

```jsonc
{
  "cell_id": "EM|E-M-E|k2|r00|bss",
  "status": "solved",
  "n_solutions": 1,
  "best_dv_kms": 0.012,
  "signature_hashes": ["sha1:5e2f…"],
  "validation_level": "V2",
  "t_done": "2026-06-01T12:34:56+10:00",
  "host": "ci-worker-3"
}
```

Lines are written via `_atomic_append` — a single `os.open(path, O_APPEND | O_WRONLY)`, `os.write(fd, line.encode("utf-8"))`, `os.close(fd)`. POSIX guarantees `O_APPEND` writes below `PIPE_BUF` (4 KiB on Linux) are atomic; our per-line budget is ~300 bytes, well below the threshold.

#### 3.2.2 Atomic claim — `Ledger.claim(cell_id, host)`

Claim is the parallel-safety primitive. Two workers must never both run the same cell. The implementation:

1. Worker reads the ledger end-to-end into an in-memory `dict[str, LedgerEntry]` (incremental read; cached after first call).
2. Worker checks `cell_id` in the dict. If present with `status in ("solved", "pruned", "failed")` → return `False` (already done).
3. If present with `status == "pending"` and the `host` matches the calling host → return `True` (resuming own work).
4. If present with `status == "pending"` and different `host` → return `False` (another worker has it; the calling worker moves on).
5. If absent → atomically append a `LedgerEntry(cell_id, status="pending", host=host, t_done=now)` and return `True`.

**This is not a perfect distributed-lock primitive.** Two workers calling `claim(cell_id)` simultaneously *could* both append `pending` entries (two `O_APPEND` writes interleave at the line level but the in-memory dict each worker built is stale). The mitigation: after appending, the worker re-reads the ledger and checks whether *its* `pending` entry is the first one — if not, it backs off (yields the cell to the earlier claimer). The first-write-wins rule is what makes the protocol correct. Documented in §5 risk #5.

For M7's single-process gate test (`test_ledger_round_trip`) this complexity is unexercised; the test just asserts the round-trip works. The race-resolution path is exercised by `test_ledger_concurrent_claim_one_wins` (a smoke-only test using `multiprocessing.Pool(2)` writing to a temp ledger).

#### 3.2.3 Re-open safety

`Ledger(path).has(cell_id)` after a process restart must return `True` for previously-recorded cells. The implementation reads the JSONL file on `__init__` and rebuilds the in-memory index. Tested by `test_ledger_persists_across_restart`.

### 3.3 Writeback (`data/writeback.py`) — spec §16.1 v2 catalogue-record-mutation

The writeback helpers are pure functions: `entry → entry'` with the gate result merged in. The new entry is a fresh frozen `CatalogueEntry` instance (immutability preserved per the M3 / M6a pattern).

#### 3.3.1 `apply_v2_to_entry(entry, report)` — the M6a hand-off integration

```python
def apply_v2_to_entry(entry: CatalogueEntry, report: StabilityReport) -> CatalogueEntry:
    new_validation = dict(entry.validation)
    gates = dict(new_validation.get("gates", {}))
    gates["V2"] = {
        "pass": report.stable,
        "max_drift_km": report.max_drift_km,
        "n_laps": report.n_laps_propagated,
        "per_lap_drift_km": list(report.per_lap_drift_km),
        "frame_used": report.frame_used,
    }
    new_validation["gates"] = gates
    # Update the top-level validation.level to the highest passing gate.
    new_validation["level"] = _highest_passing_level(gates)
    return dataclasses.replace(entry, validation=new_validation)
```

Spec §16.1 carries the `validation.gates.V2: {max_drift_km}` field. M7 adds `pass`, `n_laps`, `per_lap_drift_km`, and `frame_used` for diagnostic completeness — these are not in the spec §16.1 example but are additive (consumers ignore unknown fields per the v2 additive-fields convention).

The `_highest_passing_level` helper inspects the gates dict and returns the largest `"V0"..."V5"` whose `pass` is `True` (`V0` is implicit from `report` validity — any `report` instance has `V0=True` by construction; `V1` is from `crosscheck_cycler`; `V2` is from `report.stable`; `V3..V5` per future gate runs). For M7's M6a-only writeback, `level == "V2"` when `report.stable` is `True`.

#### 3.3.2 `record_rediscovery(entry, run_id, cell_id, date)` — the spec §16.4 retroactive-correction path

When the matcher returns `("known", entry)`, the caller (`discover`) calls `record_rediscovery(entry, ...)` to append to `entry.discovery.rediscoveries[]`. The list is keyed by `(run_id, cell_id)` for deduplication:

```python
def record_rediscovery(entry, run_id, cell_id, date):
    rediscoveries = list(entry.discovery.get("rediscoveries", []))
    key = (run_id, cell_id)
    if any((r["run_id"], r["cell_id"]) == key for r in rediscoveries):
        return entry  # idempotent — already recorded.
    rediscoveries.append({"run_id": run_id, "cell_id": cell_id, "date": date})
    new_discovery = {**entry.discovery, "rediscoveries": rediscoveries}
    return dataclasses.replace(entry, discovery=new_discovery)
```

Spec §16.4 / §16.5: attribution (`first_published`, `priority_date`) is **never overwritten** by a rediscovery — the rediscovering run only appends to the audit trail. The implementation enforces this by `dataclasses.replace` which substitutes one field at a time; `first_published` and `priority_date` are not touched.

#### 3.3.3 `serialise_entry_yaml(entry)` round-trip

The catalogue YAML is the *single source of truth* per `data/README.md`. M7 must produce YAML compatible with both the existing seed catalogue and the `cyclers.space` static-site consumer. The implementation uses `yaml.safe_dump(entry_as_dict, sort_keys=False, allow_unicode=True, default_flow_style=False)` with a custom key order matching the existing entries (`id`, `name`, `source`, `trajectory_regime`, ..., `notes`, `source_quotes`).

The round-trip test `test_serialise_entry_yaml_round_trip` reads an entry, serialises, and reads back; the resulting `CatalogueEntry` must compare field-equal to the original. Comments in the original YAML are not preserved (PyYAML does not round-trip comments); this is acceptable because M7 writeback is for *new* entries (M7 discoveries) and *gate updates* — both are programmatic operations on programmatic fields, never on human-written comment metadata.

### 3.4 V1 Lambert cross-check (`verify/crosscheck.py`) — spec §14 V1 binding

Per spec §14 V1: "every leg re-solved with **lamberthub izzo + gooding**, agreement < 1e-3 m/s; full trajectory re-propagated with the **Kepler** propagator (not the Lambert that built it), planet positions met < tol."

M7 implements the Lambert-cross-check half; the Kepler-re-propagation half is folded into the V2 path (M6a's `propagate_lap` does exactly this — multi-lap propagation uses `core.kepler.propagate`, not Lambert). So V1 in M7 is the lamberthub agreement check.

```python
def crosscheck_leg(leg: Leg, cycler: Cycler, ephem: Ephemeris) -> LambertCrosscheckResult:
    r1, r2 = _leg_endpoints(leg, cycler, ephem)
    tof_sec = leg.t_arrive - leg.t_depart
    # In-house solution.
    mine = lambert(r1, r2, tof_sec, MU_SUN, prograde=True, max_revs=leg.n_revs)[0]
    # External cross-checks.
    izzo_v1, _ = lamberthub.izzo2015(MU_SUN, r1, r2, tof_sec, M=leg.n_revs, prograde=True)
    gooding_v1, _ = lamberthub.gooding1990(MU_SUN, r1, r2, tof_sec, M=leg.n_revs, prograde=True)
    diff_izzo = np.linalg.norm(mine.v1 - izzo_v1)
    diff_gooding = np.linalg.norm(mine.v1 - gooding_v1)
    max_diff_mps = max(diff_izzo, diff_gooding) * 1000.0
    return LambertCrosscheckResult(
        leg_index=...,
        mine_v1_kms=tuple(mine.v1),
        lamberthub_izzo_v1_kms=tuple(izzo_v1),
        lamberthub_gooding_v1_kms=tuple(gooding_v1),
        max_diff_mps=max_diff_mps,
        pass_=max_diff_mps < 1.0e-3,
    )
```

`lamberthub` is already an existing dep (the M1 gate test uses it). No new dependency.

**Design notes:**

- **Single-leg.** `crosscheck_leg` operates on one `Leg` at a time; `crosscheck_cycler` is the aggregator. This keeps each leg independently reportable in the catalogue's `validation.gates.V1` block.
- **`max_revs` not `n_revs`.** The in-house `lambert` returns a list of `LambertSolution` indexed by `(n_revs, branch)`; the cross-check picks the same `n_revs` and branch the cycler's leg actually uses (`leg.n_revs`, `leg.branch`) and compares against lamberthub's same `M=n_revs` solution. Mismatched branches would falsely fail.
- **m/s, not km/s.** Spec §14 V1 tolerance is `< 1e-3 m/s`. M7 reports in m/s for direct comparison.

### 3.5 Ledger-backed runner (`data/discover.py`) — composition layer

`discover` wraps `find_cyclers` with per-cell ledger persistence + auto-validation + matching. The runner is the M7 integration point that exercises everything else.

#### 3.5.1 Pseudocode

```python
def discover(
    bodies, k_synodic, vinf_cap, ledger_path,
    *, ephem=None, l_max=4, n_max=0, branch_set=("single",),
    n_starts=5, seed=0, use_de=True, rp_factors=None,
    host=None, enable_v3=False, finder_version="0.7.0",
):
    ephem = ephem or Ephemeris("circular")
    ledger = Ledger(ledger_path)
    catalog = load_catalog()
    host = host or socket.gethostname()
    for cell in feasible_cells(bodies, l_max, k_synodic, n_max, vinf_cap, ephem, branch_set):
        if ledger.has(cell.id):
            continue
        if not ledger.claim(cell.id, host):
            continue
        try:
            result = optimise_cell_idealized(
                cell, ephem,
                vinf_cap=vinf_cap, n_starts=n_starts, seed=seed,
                use_de=use_de, rp_factors=rp_factors,
            )
        except (ValueError, RuntimeError) as exc:
            ledger.record(LedgerEntry(cell.id, "failed", 0, None, (), None, _now(), host))
            continue
        if not result.constraints_satisfied:
            ledger.record(LedgerEntry(cell.id, "searched", 0, None, (), "V0", _now(), host))
            continue
        signature = canonical_signature(result.best_cycler, model_assumption="circular-coplanar")
        match_result = match(signature, catalog)
        level = _auto_validate(result, catalog, ephem, enable_v3=enable_v3)
        ledger.record(LedgerEntry(
            cell.id, "solved", 1, result.closure_residual_kms,
            (signature.hash,), level, _now(), host,
        ))
        yield result, match_result, level
```

**Design notes:**

- **Yields `(result, match_result, level)`** — the caller (M8 reporter, future CLI) sees the optimiser output, the match outcome, and the validation level.
- **`_auto_validate`** runs V0 (`constraints_satisfied`), V1 (`crosscheck_cycler` — if `max_diff_mps < 1e-3`), V2 (`verify_long_term_stability` — if `report.stable`), and (optionally) V3. Returns the highest passing level as a string. The wall-clock budget per cell is ~5s (V0/V1 ~1s, V2 ~3s — per M6a hand-off note runtime).
- **No YAML writeback by default.** The runner appends to the ledger only. YAML writeback for a discovered candidate (`source: this-project`) is a separate explicit step the operator runs (`scripts/promote-discovery.py` — out of scope for M7's CI tests; mentioned here for completeness).
- **`enable_v3=False` default.** Until M6b lands, `optimise_cell_ephemeris` raises `NotImplementedError`; calling V3 from the auto-pipeline would crash. M7 leaves the feature flag off; M6b flips it.

### 3.6 Imports / dependency graph after M7

```
constants.py                (M0)
ephemeris.py                (M1 + M6 slice)
lambert.py                  (M1)
kepler.py                   (M1)
flyby.py                    (M2)
frames.py                   (M3 + M6a)
tisserand.py                (M2)
resonance.py                (M2)
model/cycler.py             (M3)
search/construct.py         (M3)
search/sequence.py          (M4)
model/score.py              (M4)
search/optimize.py          (M5)
search/phase_match.py       (M6 slice)
verify/propagate.py         (M6a)
verify/crosscheck.py        (M7) ← lambert, model/cycler, ephemeris      [new]
                                  ← lamberthub (external)
data/catalog.py             (M7) ← model/cycler, search/optimize.Cell,   [new]
                                  pyyaml, hashlib, json
data/ledger.py              (M7) ← os, json, dataclasses, typing         [new]
data/writeback.py           (M7) ← data/catalog (CatalogueEntry),         [new]
                                  verify/propagate (StabilityReport),
                                  verify/crosscheck (LambertCrosscheckResult)
data/discover.py            (M7) ← data/catalog, data/ledger,             [new]
                                  data/writeback, search/optimize,
                                  search/sequence, verify/propagate,
                                  verify/crosscheck, core/ephemeris,
                                  socket, datetime
```

No cycles. The two new external dependencies are `pyyaml` (already in deps — `tests/_catalogue_loader.py` uses it) and `hashlib`+`json` (stdlib). `lamberthub` is already a dep.

### 3.7 API summary

| Symbol | Purpose | Where defined |
|---|---|---|
| `CatalogueEntry` | Frozen dataclass; in-memory projection of a catalogue YAML row | `data/catalog.py` (NEW) |
| `CanonicalSignature` | Frozen dataclass; spec §16.2 identity object | `data/catalog.py` (NEW) |
| `canonical_signature(cycler, *, model_assumption) -> CanonicalSignature` | Compute signature from `Cycler` | `data/catalog.py` (NEW) |
| `signature_distance(sig_a, sig_b) -> float` | Weighted L1 distance per §16.3 | `data/catalog.py` (NEW) |
| `Catalog` | Loaded-catalogue container with `by_id`, `by_hash`, `filter()` | `data/catalog.py` (NEW) |
| `load_catalog(path) -> Catalog` | Load + index `data/seed_cyclers.yaml` | `data/catalog.py` (NEW) |
| `match(candidate, catalog) -> MatchResult` | Spec §16.3 exact/probable/novel | `data/catalog.py` (NEW) |
| `MatchResult` | Tagged result of `match()` | `data/catalog.py` (NEW) |
| `TAU_NEAR: Final[float] = 0.5` | Probable-match distance threshold | `data/catalog.py` (NEW) |
| `LedgerStatus`, `LedgerEntry`, `LedgerError`, `Ledger` | JSONL append-only ledger | `data/ledger.py` (NEW) |
| `apply_v0_v1_to_entry`, `apply_v2_to_entry`, `apply_v3_to_entry` | Writeback helpers | `data/writeback.py` (NEW) |
| `record_rediscovery`, `register_discovery`, `serialise_entry_yaml` | Discovery-workflow helpers | `data/writeback.py` (NEW) |
| `crosscheck_leg`, `crosscheck_cycler`, `LambertCrosscheckResult` | Spec §14 V1 gate machinery | `verify/crosscheck.py` (NEW) |
| `discover(...) -> Iterator[(OptimisationResult, MatchResult, level)]` | Ledger-backed runner | `data/discover.py` (NEW) |

---

## 4. Tests + gates

Tests live under `tests/data/` and `tests/verify/test_crosscheck.py`. Tolerances are spec-fixed (binning per §16.2) or module-level constants (`TAU_NEAR`); none are per-test.

### 4.1 Gate tests (spec §8 M7 + §16.2 + §16.3 + §12.2 + §14 V1 binding)

| Test | Assertion | Tolerance |
|---|---|---|
| `test_rediscovered_2syn_em_tagged_known` (**M7 BINDING GATE — spec §8 M7**) | `result = find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)[0]`; `sig = canonical_signature(result.best_cycler, model_assumption="circular-coplanar")`; `match_result = match(sig, load_catalog())`. Assert `match_result.outcome == "known"`, `match_result.entry.id == "s1l1-2syn-em-cpom"`, `match_result.entry.priority_date == "2002-08-05"`. | bool/string exact |
| `test_signature_rotation_invariant` (**spec §16.2 BINDING GATE**) | For a fixed cycler `c`, build `c_rot_i` for `i ∈ 1..len(c.encounters)-1` by cyclically rotating `c.bodies`, `c.encounters`, `c.legs`. Assert `canonical_signature(c_rot_i).hash == canonical_signature(c).hash` for all `i`. | bitwise exact |
| `test_signature_binning_absorbs_noise` (**spec §16.2 BINDING GATE**) | Perturb every V∞ by ±0.02 km/s (random within ±0.02), every leg `(a, e)` by ±0.005 AU / ±0.005. Assert `canonical_signature(perturbed).hash == canonical_signature(original).hash`. Then perturb by 0.06 km/s on V∞ and assert the hash *changes*. | bitwise exact (both directions) |
| `test_match_aldrin_classic_returns_known` | Build `aldrin_cycler = build_aldrin_seed(Ephemeris("circular"))`; `match(canonical_signature(aldrin_cycler), load_catalog())` returns `("known", entry)` with `entry.id == "aldrin-classic-em-k1-outbound"` and `entry.priority_date == "1985-10-28"`. | bool/string exact |
| `test_match_novel_synthetic_returns_novel` | Build a synthetic `Cycler` with V∞ multiset `[("E", 4.21), ("M", 7.83)]` (not in any catalogue entry within 0.05 km/s). Assert `match(...).outcome == "novel"`. | bool exact |
| `test_match_partition_by_model_assumption` (**spec §12.2 BINDING GATE**) | Synthetic candidate with `model_assumption="circular-coplanar"`. Force its V∞ signature to numerically equal the Arenstorf cr3bp entry's (which would be a near-impossible numerical coincidence but the test fabricates it). Assert `match(...).outcome != "known"` — the pool filter excludes the cr3bp entry. | bool exact |
| `test_v1_lambert_crosscheck_aldrin` (**spec §14 V1 gate**) | For each leg of `build_aldrin_seed(Ephemeris("circular"))`, `crosscheck_leg(leg, cycler, ephem).max_diff_mps < 1.0e-3`. | `< 1.0e-3` m/s |

### 4.2 Catalogue loader tests

| Test | Assertion |
|---|---|
| `test_load_catalog_returns_all_219_entries` | `cat = load_catalog()`; `len(cat.entries) == 219` (or current row count; pulled live from the YAML so it tracks). |
| `test_catalogue_signature_index_covers_constructible_entries` | For every entry in `tests/_catalogue_loader.load_constructible_entries()`, `cat.by_id[entry.id].signature_hash is not None`. |
| `test_catalogue_family_seeds_have_null_signature` | Family-seed entries (`jones-2017-vem-triple-family`, `wittal-2022-em-cycler-family`, etc.) have `signature_hash is None` — they cannot be matched against, by design. |
| `test_catalogue_filter_by_model_assumption_partitions_cr3bp` | `cat.filter(model_assumption="cr3bp")` includes `arenstorf-em-figure8-1963` and `genova-aldrin-2015-em-3petal-cycler`; excludes every circular-coplanar entry. |
| `test_catalogue_filter_by_bodies` | `cat.filter(bodies=("E", "M"))` returns only entries whose `bodies == ["E", "M"]` (set-equal). |

### 4.3 Catalogue-rediscovery parametrised tagging tests

Reuse the M5-era `tests/_catalogue_loader.py` infrastructure. For each constructible entry in the seed catalogue:

1. Build the M4 cell from the entry.
2. Run `optimise_cell_idealized` with `vinf_cap = max(target_vinfs) + 2.5`.
3. Compute the signature of the optimised cycler.
4. Match against `load_catalog()`.
5. Assert: the matched entry's `id` equals the parametrised entry's `id` (with `EXPECTED_SKIPS` honoured).

This is the spec §8 M7 gate at scale — every literature entry M5 can reconstruct must be tagged `known`, never `novel`. Failures are diagnostic for both M5 (wrong basin) and M7 (binning too tight, sense ambiguity).

| Test | Assertion |
|---|---|
| `test_catalogue_entry_tagged_known[<entry_id>]` (parametrised) | For each constructible non-skipped entry, the rediscovered signature `match(sig, catalog).outcome == "known"` AND `match(sig, catalog).entry.id == <entry_id>`. Marked `pytest.mark.slow`. |

### 4.4 The 0.5 `TAU_NEAR` and 0.05/0.01/0.01 bin widths — derivation

Spec §16.2 fixes the bin widths (0.05 km/s for V∞, 0.01 AU for `a`, 0.01 for `e`). These are not M7's to set; they are **spec-binding**. The derivation cited in the spec:

1. **0.05 km/s V∞ bin:** the literature reports V∞ to ≤ 0.05 km/s precision (e.g. McConaghy 2006 reports 5.65 ± 0.05). A finer bin would cause every literature entry to "miss" by 1 bin against a M5 rediscovery whose convergence noise is ~0.01–0.05 km/s. The bin is exactly the literature precision.
2. **0.01 AU `a` bin / 0.01 `e` bin:** literature reports `a, e` to 2 decimal places (Rogers 2012 Table 1: a=1.30 / e=0.257; bin = 0.01 → exact). Real-ephemeris `a, e` vary by ~0.005 AU / ~0.005 across launch epochs (eccentricity-driven); 0.01 is loose enough to absorb this.
3. **`TAU_NEAR = 0.5`:** the weighted-L1 signature distance is computed as `|Δperiod_yr / 0.5| + Σ |Δvinf_kms / 0.05| + Σ |Δa_AU / 0.01| + Σ |Δe / 0.01|`. A "probable match" is one whose total weighted distance is below 0.5 — i.e. one bin's worth of mismatch across the entire signature, accumulated. This corresponds to a candidate that misses every bin by 25–50% of the bin width, which is the regime where a literature entry might have a different rounding convention but be the same cycler. Below 0.5 → human review; above 0.5 → genuinely different cycler.

The 0.5 value is the new M7 constant. §5 risk #4 covers the sensitivity test (what if it's 0.3? 1.0?); a sensitivity-sweep test (`test_tau_near_sensitivity`) is part of §4.5 — informational, not gating.

### 4.5 Sensitivity / informational tests

| Test | Assertion |
|---|---|
| `test_tau_near_sensitivity` | Sweep `TAU_NEAR ∈ {0.3, 0.5, 1.0}`; count `("probable-match-NEEDS-HUMAN")` outcomes across the constructible-entries set; assert the count at `TAU_NEAR=0.5` is between 0 and 5 (the seed catalogue has no near-duplicates by design). Informational — diagnoses tuning. |
| `test_canonical_json_deterministic` | `_canonical_json({"b": 1, "a": 2}) == _canonical_json({"a": 2, "b": 1})` — sorted-keys invariant. |
| `test_signature_hash_stable_across_python_versions` | Hash of a fixed input string `"hello"` equals a hardcoded sha1 reference — guards against an accidental switch from sha1 to another hash. |

### 4.6 Ledger tests

| Test | Assertion |
|---|---|
| `test_ledger_round_trip` (**M7 GATE**) | Record entry → reread → field-equal. Re-recording the same `cell_id` raises `LedgerError`. |
| `test_ledger_persists_across_restart` | Write via `Ledger(path)`; close; reopen `Ledger(path)`; `.has(cell_id)` returns `True`. |
| `test_ledger_atomic_append_no_partial_lines` | Open ledger, write a 200-byte line; assert the file ends with `\n` and has no partial JSON. |
| `test_ledger_concurrent_claim_one_wins` (smoke) | `multiprocessing.Pool(2)` both call `claim(cell_id)`; exactly one returns `True`. Marked `pytest.mark.flaky` because race timing is OS-dependent. |
| `test_ledger_iter_pending_skips_solved` | Mix of `solved` and `pending` entries; `iter_pending` yields only the `pending` cell ids. |

### 4.7 Writeback tests

| Test | Assertion |
|---|---|
| `test_v2_writeback_populates_validation_block` (**M7 GATE — M6a integration**) | Build entry; build `StabilityReport(stable=True, max_drift_km=12345.6, ...)`; `new_entry = apply_v2_to_entry(entry, report)`; `new_entry.validation["gates"]["V2"]["pass"] == True`, `["max_drift_km"] == 12345.6`. Entry is a new instance — `entry` itself is unmodified (frozen). |
| `test_v2_writeback_promotes_validation_level` | If V0 and V1 pass and V2 passes, `new_entry.validation["level"] == "V2"`. |
| `test_record_rediscovery_idempotent` | Recording the same `(run_id, cell_id)` twice produces one entry, not two. |
| `test_record_rediscovery_preserves_first_published` | `first_published` and `priority_date` are unmodified after `record_rediscovery`. |
| `test_register_discovery_sets_candidate_novel` | `register_discovery(skeleton, run_id, ...).our_status == "candidate-novel"`, `.source == "this-project"`, `.first_published is None`. |
| `test_serialise_entry_yaml_round_trip` | Round-trip an existing constructible entry: `parse(serialise(entry)) == entry`. Field-equality only; comments are not preserved. |

### 4.8 Crosscheck tests

| Test | Assertion |
|---|---|
| `test_v1_lambert_crosscheck_aldrin` (**M7 GATE — spec §14 V1**) | Aldrin classic; every leg passes (`max_diff_mps < 1.0e-3`). |
| `test_v1_lambert_crosscheck_2syn_em` | 2-syn E-M cycler; every leg passes. |
| `test_crosscheck_leg_reports_correct_n_revs_branch` | If `leg.n_revs == 1, leg.branch == "low"`, the lamberthub call uses `M=1` and selects the low-energy branch. |
| `test_crosscheck_cycler_aggregates` | `crosscheck_cycler(aldrin).pass_ == all(r.pass_ for r in result)`. |

### 4.9 Discover-runner integration tests

| Test | Assertion |
|---|---|
| `test_discover_em_k2_yields_known_for_2syn` | `list(discover(("E","M"), 2, 7.0, tmp_path/"ledger.jsonl"))`; the result whose `cycler.max_vinf()` is near 5.65 km/s has `match_result.outcome == "known"` and `match_result.entry.id == "s1l1-2syn-em-cpom"`. |
| `test_discover_writes_ledger_for_every_cell_attempted` | After the run, every cell yielded by the M4 enumerator has a `LedgerEntry` in the ledger with `status` in `{"solved", "pruned", "failed", "searched"}`. |
| `test_discover_resumes_from_existing_ledger` | First run writes 3 ledger entries, simulated as "solved"; second run picks up the same cells, skips them (no `optimise_cell_idealized` calls), yields nothing new for already-done cells. Run with a spy on `optimise_cell_idealized` to assert the call count is zero on resume. |
| `test_discover_records_signature_hash_in_ledger` | After running, every `solved` ledger entry has `signature_hashes != ()`. |
| `test_discover_skips_v3_when_disabled` | `enable_v3=False` (default); no `NotImplementedError` is raised even though `optimise_cell_ephemeris` is the M5 stub. |

### 4.10 Determinism, frozenness, and type-strict tests

| Test | Assertion |
|---|---|
| `test_catalogue_entry_frozen` | `entry.id = "..."` raises `FrozenInstanceError`. |
| `test_canonical_signature_frozen` | Same. |
| `test_match_result_frozen` | Same. |
| `test_ledger_entry_frozen` | Same. |
| `test_lambert_crosscheck_result_frozen` | Same. |
| `test_canonical_signature_deterministic` | Two calls produce bitwise-identical `hash`. |
| `test_load_catalog_deterministic` | Two calls produce field-equal `Catalog` instances. |

### 4.11 Tolerance summary

| Layer | Quantity | Tolerance | Source |
|---|---|---|---|
| V∞ bin | km/s | 0.05 | spec §16.2 |
| `a` bin | AU | 0.01 | spec §16.2 |
| `e` bin | dimensionless | 0.01 | spec §16.2 |
| Period years bin | years | 0.5 (informational only; not in hash) | M7 §4.4 |
| `TAU_NEAR` (probable-match) | weighted L1 | 0.5 | M7 §4.4 |
| V1 Lambert agreement | m/s | < 1.0e-3 | spec §14 V1 |
| V2 drift (M6a-defined) | km | < 50,000 | spec §14 V2 (M6a) |
| Catalogue rediscovery tagging | bool | exact ("known") | spec §8 M7 |

---

## 5. Risks

| # | Risk | Likelihood | Impact | Mitigation in M7 |
|---|---|---|---|---|
| 1 | **Signature binning collides two physically distinct cyclers.** Two cyclers with similar V∞ multisets and `(a, e)` could land in the same bin and produce identical signatures, falsely returning `("known", wrong_entry)`. | medium | **high — wrong attribution** | The §16.2 bin sizes are the literature precision; collisions in this band correspond to cyclers that *should* be considered the same identity per spec §16.2 (which is what the canonicalisation is for). The defence is the seed catalogue's coverage: if two distinct families collide, both will appear in the catalogue (or one will be a `probable-match` after the other lands), and human review per §16.4 is the conflict resolver. §5 risk #4 covers the dial. |
| 2 | **`sense` field derivation from `Cycler` is ambiguous.** The M3 / M5 `Cycler` does not carry an explicit `sense` field; M7 derives it from the first leg's `tof_days` vs the cycler's period (short outbound first → `"outbound"`). For cyclers where the leg ToFs are nearly equal (e.g. some 2-syn families), the derivation may flip and produce a different signature. | medium | medium | The derivation is a single helper `_derive_sense(cycler)`; if it produces ambiguous results on a real run, the fix is to extend `Cycler` or `OptimisationResult` with an explicit `sense: Literal["outbound", "inbound"]` field (additive change). M7 ships the helper as-is and documents the limitation; if the M7 gate fails due to sense ambiguity, escalate via the spec §16.2 "optionally also fold under reversal" path — set sense to `"either"` and accept the wider matching. |
| 3 | **Rotation canonicalisation is wrong for long sequences.** `_lex_min_rotation` is O(n²) over the sequence parts; M4's `L_max=8` is fine; M8's VEM at `L_max=10..12` would slow but never fail. | low | low | The O(n²) is unconditional but bounded; for `L=12` it's 144 string comparisons per signature, microseconds. If profiling ever shows this is a hot path, swap to Booth's O(n) — but only if measurements demand it. Documented as M8-future-work in the helper's docstring. |
| 4 | **`TAU_NEAR = 0.5` is wrong.** Too tight: every literature near-rediscovery becomes `("novel", None)` and the gate fails. Too loose: every distinct cycler becomes `("probable-match-NEEDS-HUMAN", wrong_entry)` and human review is swamped. | medium | high — gate / human-time impact | `test_tau_near_sensitivity` (informational) sweeps over `{0.3, 0.5, 1.0}` and reports the probable-match count on the constructible-entries set. The seed catalogue is small enough (219 entries) that the operator can audit the probable-matches; if the count at 0.5 exceeds 5, retune to 0.3 in a follow-up commit. The dial is exposed as a module constant. |
| 5 | **Ledger atomic-claim race window allows duplicate work.** Two workers calling `claim(cell_id)` simultaneously could both append `pending` entries. | low (in M7) | low (correctness preserved; only redundant work) | The first-write-wins protocol resolves the race after-the-fact: both workers detect the duplicate `pending` on re-read; the later-timestamped one backs off. Redundant computation is bounded by the worker count (at most N workers redo a cell in the worst case). Test `test_ledger_concurrent_claim_one_wins` is marked `flaky` because the race timing is OS-dependent; the more important `test_ledger_round_trip` is the binding gate. M8's parallel runner may revisit with file-lock-based claim if redundancy becomes operationally expensive. |
| 6 | **Catalogue YAML write corrupts the canonical file.** `serialise_entry_yaml` followed by an on-disk write that loses a comment or reorders a field would silently break the seed catalogue. | low | **high — data loss** | M7's CI tests NEVER write to `data/seed_cyclers.yaml`. Writeback paths use `tmp_path` only. Operator-driven runs (e.g. `scripts/promote-discovery.py`) write to a staged file the operator manually diffs against the canonical, then commits. The diff is the human gate. Documented in `data/README.md` (M7 will append a writeback section to it as part of closeout). |
| 7 | **`lamberthub` API drift breaks V1 cross-check.** Spec §14 V1 binds against lamberthub's izzo + gooding; an upstream API change (function rename, return-tuple change) breaks the test. | low | medium | `lamberthub` is already pinned in `uv.lock`; the M7 commit doesn't bump it. If a future dependency-upgrade PR breaks the test, the failure is loud (the V1 gate fires) and the fix is localised to `crosscheck.py`. The cross-check call is a single function per solver; the adapter surface is small. |
| 8 | **`model_assumption` pool filter excludes catalogue entries silently.** If a future M5-produced cycler tags itself with `model_assumption` other than `"circular-coplanar"`, the matcher silently sees a smaller pool and may return `("novel", None)` for genuine-known cyclers. | low (in M7) | medium | M5's `find_cyclers` produces `OptimisationResult` instances tagged `"circular-coplanar"` by `canonical_signature`'s default. M6b's ephemeris-mode optimiser will need to set the right tag — likely also `"circular-coplanar"` because the signature is *invariant* to ephemeris realisation per spec §12.2 (the cycler's identity is its idealised geometry; the launch window is a derived view). Documented in `canonical_signature`'s docstring; M6b plan will clarify. |
| 9 | **Frozen-dataclass round-trip with `dataclasses.replace` doesn't update nested dicts.** `apply_v2_to_entry`'s `validation` field is a nested dict; `replace(entry, validation=new_dict)` works because the field is the top-level reference, but mutating the dict in-place would not. | low | low | The implementation deliberately builds a new dict via `dict(...)` and `dict(new_validation.get("gates", {}))` (shallow copies) before mutating. Tests assert immutability of the input `entry` after the call. |
| 10 | **Catalogue grows beyond what an in-memory `Catalog` can index efficiently.** At 1000+ entries the `by_hash` dict is fine but `filter()` does linear scans. | low (current scale) | low | At 219 entries linear-scan `filter` is microseconds. At 10,000 it's still milliseconds. If matching becomes hot, switch `filter()` to a multi-key index — a one-class change with no API impact. Documented as M8-future-work. |
| 11 | **Sub-deliverable C (ledger) is the most novel piece and adds the most risk.** Spec §13.6 / §13.8 wording is concise; implementation details (file locking, restart semantics) are M7's call. | medium | medium | The ledger's M7 surface is *minimal*: append, has, claim, get, iter_pending. Each method has a focused unit test. The atomic-claim protocol is the only non-trivial logic; tested by `test_ledger_concurrent_claim_one_wins`. If the M8 parallel runner needs more, the ledger interface is a class — swap-in of a SQLite implementation is a one-class change. |
| 12 | **The M5 `find_cyclers` returns a cycler whose `cell` is in a different rotation from the catalogue's `sequence_canonical`.** E.g. M5 returns `("E", "M", "E")` for a `sequence_canonical: "E-M"` entry (one extra Earth from the open-sequence convention). | medium | medium — signature mismatch | `canonical_signature` rotates and dedupes; an `("E","M","E")` open sequence in `cycler.bodies` becomes `_lex_min_rotation("E-M")` (after dedup of the closing body) which matches the catalogue's `"E-M"`. The implementation must handle the open-sequence closing-body convention correctly. Tested in `test_signature_open_sequence_normalisation` — explicitly assert that `canonical_signature(Cycler(bodies=["E","M","E"], ...))` produces the same sequence-canonical string as the catalogue's `"E-M"`. |
| 13 | **`tests/_catalogue_loader.py` is shared with M5; modifying it would break M5's gate.** | low | low | M7 does NOT modify `tests/_catalogue_loader.py`. M7's new tests live under `tests/data/`. If M7's tests need a broader projection of catalogue rows than `CatalogueEntry` (M5's minimal-projection dataclass), M7 imports from `cyclerfinder.data.catalog` and reads catalogue rows directly. The two projections coexist. |
| 14 | **Parallel agent on M6b is editing `verify/` concurrently.** | low | low (per parent-agent coordination) | M7's edit to `verify/__init__.py` is a single `from . import crosscheck` line; M6b's plan is separately written and does not touch `verify/__init__.py`. The conflict surface is one line; the parent agent's commit ordering resolves it (M6a → M7 → M6b, or M6a → M6b → M7, both fine). M7's plan flags this. |
| 15 | **Signature includes `period_years` but the spec §16.2 hash subset isn't explicit.** Spec §16.1 example has `period: {pair, k, years}` in the `signature_fields` block. §16.2 says the hash is over the canonical-JSON of those fields. So `years` IS in the hash. But two literature entries describing the same cycler with slightly different `years` (4.27 vs 4.28) would mismatch. | medium | medium | M7 **excludes `period.years` from the hash subset** but keeps it in the `CanonicalSignature` dataclass for display. Rationale: `period.k` is exact (an integer); `period.years` is `k * synodic_period_years` which itself depends on which value of the synodic period the entry uses (2.135 from spec §9, vs 2.13 from Hollister). Including `years` in the hash would silently invalidate matches against literature entries with slightly different synodic-period conventions. The spec §16.1 example shows `years: 4.27` but the canonical identity is `k` — `years` is derived. Documented explicitly in `canonical_signature`'s docstring as a deliberate spec-interpretation. **Risk-resolution: M7 ships with `period.years` excluded from the hash; if a follow-up review of spec §16.2 contradicts this, the change is a one-line edit (add `years` to the canonical-JSON input).** |

---

## 6. Dependency additions

**None.** M7 uses only:

- `numpy` — already in deps.
- `scipy` — already in deps.
- `astropy` — already in deps (M6 slice).
- `pyyaml` — already in deps (used by `tests/_catalogue_loader.py`).
- `lamberthub` — already in deps (used by M1 gate).
- stdlib: `hashlib`, `json`, `dataclasses`, `pathlib`, `os`, `socket`, `datetime`, `typing`, `multiprocessing`.
- in-house M0–M6a modules.

No edits to `pyproject.toml`. No `uv.lock` regeneration.

---

## 7. Order of work

The `todo.md` mirrors this with checkboxes.

1. **Re-read predecessor docs.** Re-read spec §16 (full — the M7 design document), §13.5–§13.8 (ledger and signature), §12.2 (representation framework), §14 V1 / V2 (the gates M7 wires together), §10 (closure-correctness risk — M7 inherits the binding-frame requirement via M6a). Re-read M6a's hand-off note in `docs/phases/m6a-idealized-closure-verification/todo.md` (the `StabilityReport` shape). Re-read M5's `## Hand-off to M6a` (the `OptimisationResult` shape, the M5-reproduced V∞ values). Re-read `tests/_catalogue_loader.py` and `tests/test_catalogue_rediscovery.py` — M7 extends, does NOT replace this infrastructure.
2. **Create the `data/` subpackage scaffolding.** `src/cyclerfinder/data/__init__.py` empty; `tests/data/__init__.py` empty. `mypy --strict` clean.
3. **Implement `data/catalog.py` — the signature half first.**
   - Module docstring referencing spec §16.2, §16.3.
   - `CanonicalSignature` frozen dataclass.
   - `_lex_min_rotation`, `_bin_vinf`, `_bin_a_au`, `_bin_e`, `_canonical_json` helpers.
   - `canonical_signature(cycler, *, model_assumption)` public function.
   - Land tests `test_signature_rotation_invariant`, `test_signature_binning_absorbs_noise`, `test_canonical_signature_deterministic`, `test_canonical_json_deterministic`, `test_signature_hash_stable_across_python_versions`, `test_signature_open_sequence_normalisation`.
4. **Implement `data/catalog.py` — the loader + matcher half.**
   - `CatalogueEntry` dataclass.
   - `_extract_signature_fields(entry)` helper.
   - `Catalog` class with `by_id`, `by_hash`, `entries`, `filter()`.
   - `load_catalog(path=CATALOGUE_PATH) -> Catalog`.
   - `signature_distance`, `match`, `MatchResult`, `TAU_NEAR`.
   - Land tests in `tests/data/test_catalogue_loader.py` and `tests/data/test_match.py`.
   - **M7 BINDING GATE** `test_match_aldrin_classic_returns_known` and `test_match_partition_by_model_assumption` land here.
5. **Implement `data/ledger.py`.**
   - `LedgerStatus`, `LedgerEntry`, `LedgerError`, `Ledger`, `LedgerLoader`.
   - `_atomic_append`, `_parse_line` helpers.
   - Land tests in `tests/data/test_ledger.py` including `test_ledger_round_trip` (gate), `test_ledger_persists_across_restart`, the flaky `test_ledger_concurrent_claim_one_wins`.
6. **Implement `verify/crosscheck.py`.**
   - `LambertCrosscheckResult`, `crosscheck_leg`, `crosscheck_cycler`.
   - Add `from . import crosscheck` to `verify/__init__.py`.
   - Land tests in `tests/verify/test_crosscheck.py` including the **M7 GATE** `test_v1_lambert_crosscheck_aldrin` (`< 1e-3` m/s per spec §14 V1).
7. **Implement `data/writeback.py`.**
   - `apply_v0_v1_to_entry`, `apply_v2_to_entry`, `apply_v3_to_entry`.
   - `record_rediscovery`, `register_discovery`.
   - `serialise_entry_yaml`.
   - `_highest_passing_level` helper.
   - Land tests in `tests/data/test_writeback.py` including **M7 GATE** `test_v2_writeback_populates_validation_block`, `test_serialise_entry_yaml_round_trip`.
8. **Implement `data/discover.py`.**
   - `discover` generator function.
   - `_auto_validate` helper.
   - Land tests in `tests/data/test_discover.py`.
9. **Land the catalogue-rediscovery tagging gate.** `tests/data/test_catalogue_rediscovery_tagging.py`:
   - **M7 BINDING GATE** `test_rediscovered_2syn_em_tagged_known` (spec §8 M7 anchor).
   - Parametrised `test_catalogue_entry_tagged_known[<entry_id>]` over the constructible-entries set (reuse `tests/_catalogue_loader.load_constructible_entries`); marked `pytest.mark.slow`.
10. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
11. **Commit** as `m7: catalogue loader + canonical signature matcher + novelty tagging + ledger + V1 cross-check`. Push; confirm CI green.
12. **Update `docs/overview.md`.** §4 milestone table: M7 → completed; M8 → planned. §2 deferred-decisions table: cross out "Ledger backend — SQLite vs JSONL (decided in M4/M7 when ledger is built)" and move to kept-decisions: "Ledger backend: JSONL (M7). SQLite reserved for a future scaling pass."
13. **Hand-off note** appended to `todo.md` under `## Hand-off to M8` covering: the actual rediscovery rate across the 219-row catalogue (how many `known`, how many `probable-match`, how many `novel`), any spec §16.2 sense / period.years ambiguities encountered, ledger size at run end (informs M8 parallel-runner sizing), `TAU_NEAR` sensitivity result, and the locked `CanonicalSignature` / `MatchResult` / `LedgerEntry` shapes.

The order is "signature → loader/matcher → ledger → crosscheck → writeback → discover → catalogue gate" deliberately: each step depends only on what came before; the spec §16.2 identity gates land first (signature correctness is load-bearing); the spec §8 M7 gate (rediscovery tagging) lands last as the composition test.

---

## 8. Exit checklist (the gate, restated)

Before declaring M7 done:

- [ ] `uv run pytest tests/data/test_canonical_signature.py` green — spec §16.2 rotation invariance and binning pass at bitwise exact.
- [ ] `uv run pytest tests/data/test_match.py` green — including the spec §12.2 pool-filter test and the §16.3 exact/probable/novel state machine.
- [ ] `uv run pytest tests/data/test_catalogue_loader.py` green — all 219 entries loaded; signature index populated for constructible entries; family-seed entries have `signature_hash is None`.
- [ ] `uv run pytest tests/data/test_catalogue_rediscovery_tagging.py::test_rediscovered_2syn_em_tagged_known` green — **the M7 binding gate** (spec §8 M7) passes with `match.entry.id == "s1l1-2syn-em-cpom"`.
- [ ] `uv run pytest tests/verify/test_crosscheck.py::test_v1_lambert_crosscheck_aldrin` green — spec §14 V1 gate passes at `max_diff_mps < 1.0e-3`.
- [ ] `uv run pytest tests/data/test_ledger.py::test_ledger_round_trip` green; restart-survives test passes.
- [ ] `uv run pytest tests/data/test_writeback.py::test_v2_writeback_populates_validation_block` green — M6a `StabilityReport` integration works.
- [ ] `uv run pytest tests/data/test_discover.py::test_discover_em_k2_yields_known_for_2syn` green — full pipeline integration.
- [ ] `uv run pytest` green overall (no regression of M0–M6a tests).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true` — including the new dataclasses, the `Catalog.filter` keyword args, the JSONL `LedgerEntry` round-trip, `Literal["known", "probable-match-NEEDS-HUMAN", "novel"]` on `MatchResult.outcome`.
- [ ] CI green on the M7 commit.
- [ ] `docs/overview.md` updated: M7 status = `completed`; M8 row marked `planned`; the JSONL-vs-SQLite deferred decision moved to kept-decisions.
- [ ] `## Hand-off to M8` section appended to `phases/m7-catalogue-novelty-matching/todo.md` covering:
  - The catalogue-rediscovery tagging count (how many entries produced `("known", ...)`, how many `("probable-match", ...)`, how many `("novel", ...)`). Informs M8 reporter's catalogue summary.
  - Any `("probable-match-NEEDS-HUMAN", entry)` outcomes — list the candidate signatures and the matching entry; the operator must resolve each by either tightening the entry's source-data or accepting the match as known.
  - The actual `TAU_NEAR` value used (0.5 default; document if changed) and the sensitivity-test result.
  - The locked `CanonicalSignature` shape (so M8 knows what to display in catalogue tables) and `MatchResult` shape (so M8 reporter can render the three outcomes).
  - The locked `LedgerEntry` shape — M8's parallel-runner inherits this.
  - Whether `sense` derivation per risk #2 fired ambiguously on any catalogue entry; documents whether M8 needs to extend `Cycler` / `OptimisationResult` with an explicit `sense` field.
  - Whether `period.years` exclusion from the hash subset per risk #15 surfaced any matching surprises (it should not — all spec §9 anchors use 2.135 yr synodic).
  - The V1 cross-check pass rate across the constructible-entries set; if every entry passes, V1 can be unconditionally invoked in the M8 reporter; if some fail with `lamberthub` API surprises, the V1 gate becomes opt-in.
  - V3 stays disabled until M6b lands; M8 will turn it on after M6b's `optimise_cell_ephemeris` body is in place.
  - Whether the JSONL ledger's concurrent-claim race ever fired in CI (the test is `flaky`); if so, M8 may swap to a SQLite backend.
  - Confirmation that `data/seed_cyclers.yaml` was NOT modified by any M7 test (the CI gate writes only to `tmp_path`).

(Writing the M8 plan doc is the first task of M8, not an M7 exit criterion.)
