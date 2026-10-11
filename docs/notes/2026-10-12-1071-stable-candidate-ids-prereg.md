# #1071: stable candidate ids (`cid`) in the two-working-body gauntlet JSONs, pre-registration (2026-10-12)

Written before any code. Source of the task: #1024's first run was keyed on the gauntlet `key` field, which is
not unique (papercut review, twobody-ext's finding).

## 1. The problem, measured

`key` is the cycle template of the best member of a merged group; two different merged groups (different
V_inf, turns or dates) can share a template. In the committed gauntlet JSONs (checked 2026-10-12, all files
under `data/` named `*gauntlet*.json` with a `candidates` list):

| File | Candidates | Duplicate keys | Duplicate (key, worst ratio) |
|---|---|---|---|
| `data/998_pluto_smallmoons/ps_gauntlet.json` | 183 | 2 | 0 |
| `.../pn_gauntlet.json` | 149 | 4 | 0 |
| `.../pk_gauntlet.json` | 62 | 6 | 0 |
| `.../ph_gauntlet.json` | 63 | 5 | 0 |
| `data/942_cell_ev_gauntlet.json` | 31 | 4 | 0 |
| `data/942_cell_vm2n_gauntlet.json` | 23 | 5 | 0 |
| `data/942_cell_vmn_gauntlet.json` | 14 | 1 | 0 |
| `data/943_cell_ge_gauntlet.json` | 22 | 2 | 0 |
| `data/973_ev/k4_gauntlet.json` | 12 | 1 | 0 |
| `data/973_gc/k4_gauntlet.json` | 51 | 4 | 0 |
| `data/1065_barycentric/phb_gauntlet.json` | 52 | 5 | 0 |

(The lead's counts "ps 183 with 181 distinct keys, pk 62 with 56" agree: 2 and 6 duplicates.) The other files
have no duplicate keys. In every file the pair (key, ideal worst ratio) is unique, so it is a safe fallback.

## 2. The change

- `src/cyclerfinder/data/candidate_ids.py` (no literature import): `candidate_cid(cell, index, x_days)` =
  `<cell>-<4-digit zero-based index>-<8 hex of sha1 of the dates rounded to 6 decimals>`;
  `record_ratio(rec)` (the record's `worst_ratio`, else the largest `ratio` in its `flyby_table`);
  `record_id(rec)` = `cid` when present, else `"<key>@<ratio to 9 decimals>"`; `index_records(cands)` (dict by
  `record_id`, raises on a duplicate id); `index_by_key(cands)` (dict by `key`, raises on an ambiguous key
  instead of silently keeping the last); `same_record(a, b)` (cid equality when both have one, else key and
  ratio to 1e-9).
- Writers: `scripts/gauntlet_942.py` and `scripts/gauntlet_998.py` add `cid` to each new candidate record
  (index = position in the written list). No existing field, order of records or verdict changes.
- Readers: those that build a dict keyed by `c["key"]` on gauntlet records (`analyse_943_gc2_gancal_relation`,
  `run_968_control`, `run_968_control2`) use `index_by_key`, which only differs from the old dict when a key is
  ambiguous. `run_1024_barycentric.py` matches members with `same_record` and stores `cid` in its new records;
  `realeph_998.py` copies the record id into its output records.
- Not regenerated: no committed gauntlet or derived JSON is rewritten. A committed file has no `cid`; it is
  addressed by its INDEX in `candidates` (the address #1024's `index` field already uses) and, for lookups across
  files, by `record_id` (key plus ideal worst ratio).
- Not touched (not the two named writers; flagged to the lead): the #973 gauntlet writer in
  `scripts/run_973_enumerate.py` and `scripts/run_1065_barycentric_enumerate.py` write candidate records with the
  same non-unique `key`. One line each would add `cid`.

## 3. Expected effect (pre-registered)

- Every new record from the two writers carries a `cid`; two candidates with the same key get distinct cids.
- Re-running `gauntlet_998.py` for a cell on the same shards would give the same fields as the committed file plus
  `cid`; byte-identical otherwise (checked in the test on a small synthetic case, not by regenerating).
- The three dict-building readers are verified not to have collapsed a duplicate in any past use: their keys
  (`GC2_KEY`, `GANCAL5_KEY`, the 968 `KEY`s) occur once in `943_cell_gc_gauntlet.json` and once in
  `943_cell_ge_gauntlet.json`'s relevant entry (checked 2026-10-12: gc has no duplicate key; the ge duplicates are
  other keys). So no past result changes.
- Regression test: two records with the same key get different cids; a committed gauntlet JSON (ps) loads and
  `index_records` matches its records one-to-one (183 distinct fallback ids) and the lookup of every record by its
  own `record_id` returns itself.
