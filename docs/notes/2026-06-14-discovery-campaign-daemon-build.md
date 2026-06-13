# Discovery-campaign daemon build (#253, Track C) — 2026-06-14

Generalizes the corpus-search daemon (`scripts/search_campaign_daemon.py`) into a
**discovery** daemon: a plain, resumable, quota-proof process that sweeps a
pluggable search target, dedups against the catalogue + the negative registry, and
routes every outcome to the unchanged SILVER→gauntlet artefacts — **never** the
catalogue. The first wired target is the #254 repeated-moon multi-rev genome over
Jupiter.

This is a **built + tested + ready-to-launch** deliverable. The open-ended Jovian
hunt is the coordinator's background run (the launch command is at the bottom); it
was **NOT** run here (anti-hang discipline — only a tiny bounded smoke ran).

## Files

- `src/cyclerfinder/search/discovery_campaign.py` — the engine (testable).
- `tests/search/test_discovery_campaign.py` — 13 tests (all green).
- `scripts/discovery_campaign_daemon.py` — the thin driver (the coordinator runs this).
- `scripts/_smoke_discovery_campaign.py` — one-shot bounded smoke (untracked, like the `_probe*`/`scratch_*` convention).

## Architecture

```
enumerate ─▶ shard (index % n_workers) ─▶ skip done (checkpoint)
          ─▶ skip known (catalogue dedup) ─▶ close ─▶ skip known (closed-sig dedup)
          ─▶ residual < gate? ─┬─ yes ─▶ SILVER → data/review_queue.jsonl
                               └─ no  ─▶ (no row)
          ─▶ (sweep yields 0 SILVER) ─▶ method-versioned EMPTY → data/empty_regions.jsonl
```

- **Region-level capability gate** (before the loop): `should_sweep` consults the
  negative registry; if a prior empty-region over this region carries a method that
  *subsumes* ours, the whole sweep is skipped (learns nothing). A new/incomparable
  method re-sweeps (the #163-reopens-#137 lesson).
- **Resumability**: each evaluated candidate's `index` is appended to a gitignored
  `out/discovery_campaign/checkpoint_<target>_w<wid>.txt`; a restart loads the
  done-set and skips them. Enumeration order is deterministic (sorted moon set ×
  fixed resonance/`n_rev` grids), so the index↔candidate map is stable across runs.
- **Multi-worker**: `index % n_workers == worker_id` shards the deterministic
  stream; N plain processes use N cores with no shared-file contention (each has its
  own checkpoint shard; `data/*.jsonl` are append-only JSONL so concurrent
  one-line appends interleave cleanly).
- **No catalogue writeback, ever**: SILVER → `review_queue.jsonl` (a human
  promotes); negatives → `empty_regions.jsonl`. `catalogue.yaml` is read-only.

## The search-target interface (pluggable)

`SearchTarget` is a `runtime_checkable` `Protocol` — the engine knows nothing about
a target's genome, so the multi-arc harness or any future genome can be hosted by
implementing it:

- `target_id: str` / `primary: str` — identity + catalogue `primary` bucket.
- `method_capability() -> MethodCapability` — the capability tag set stamped on any
  empty-region record (so "empty" is always method-conditional).
- `enumerate_candidates() -> Iterator[Candidate]` — deterministic, resumable order.
- `close(candidate) -> ClosureResult` — one canonical residual (km/s).

The concrete `RepeatedMoonTarget` wraps the #254 genome:

- **Enumerate**: every length-k product over the sorted moon set whose consecutive
  bodies differ and which uses ≥2 distinct moons, × the per-leg `n_rev` grid.
- **Close** (generalized from `cge_scaffold.reproduce_member` to arbitrary
  sequences, body-agnostic registry-derived mean motions): sweep an initial-phase
  grid × a ToF-resonance grid; place circular-coplanar moons at the cumulative
  flyby epochs; solve the planet-frame multi-rev Lambert leg at the requested
  `n_rev`; the canonical residual is the **worst per-flyby V∞-magnitude continuity
  defect** over the cycle (the same quantity the #254 corrector minimises). The
  phasing that minimises it is returned. A leg with no Lambert solution at the
  requested `n_rev` marks the phasing infeasible; if no phasing is feasible the
  candidate is non-converged (routed as a non-hit, not a discovery).
- **Dedup signature** (`moon_cycler_signature_hash`): rotation-canonical sequence +
  V∞ multiset binned to 0.05 km/s, sorted-key ASCII JSON → sha1 — the same
  spec-§16.2 recipe, keyed on the fields a Jovicentric catalogue row stores
  (`primary` / `sequence_canonical` / `vinf_kms_at_encounters`), so a closed
  candidate dedups against the catalogue and against itself.

## Bounded smoke result (real numbers, TEMP paths)

`scripts/_smoke_discovery_campaign.py`, Jupiter, seq_len=3, `n_rev`∈{0,1,2},
12 phase samples × 4 ToF scales, gate 0.05 km/s, **first 10 candidates only**:

- **324** candidates enumerable in this bounded grid; **10** evaluated.
- Per-candidate close (first 10): **5 closed, 5 failed-close** (no Lambert
  solution at the requested `n_rev` — correctly handled). Example closed rows:
  - `idx=0 Callisto-Europa-Callisto n_rev=[0,0] resid=2.4457 km/s vinf=[4.52,6.68,3.89]`
  - `idx=1 Callisto-Europa-Callisto n_rev=[0,1] resid=1.4303 km/s vinf=[5.94,4.54,5.18]`
  - `idx=9 Callisto-Europa-Ganymede n_rev=[0,0] resid=0.0506 km/s vinf=[4.52,6.68,5.23]`
- **Routing**: 0 SILVER, **1 EMPTY** record written (to TEMP). The best of the
  first 10 closed to **0.0506 km/s** — a near-miss just above the 0.05 gate, so the
  bounded smoke routes the no-hit sweep to a method-versioned empty-region record
  (tags `{ballistic, coplanar, leveraging, multi-arc, patched-conic}`, centre
  Jupiter). The SILVER routing path is exercised by the unit test
  `test_silver_hit_routes_to_review_queue`.
- **Real `data/` untouched**: the smoke wrote only to a `tempfile.TemporaryDirectory`.

Interpretation: the 0.0506 km/s near-miss is honest — the smoke's coarse 12×4 grid
does not fully minimise the residual; it proves the pipeline end-to-end (enumerate →
close → dedup → route, both branches reachable), it is **not** a search result. The
real run uses a finer phase grid and is the coordinator's job.

## The exact command the coordinator runs (real Jovian hunt)

Launch one worker per core (4-core example); plain background processes, resumable,
quota-proof:

```sh
for w in 0 1 2 3; do
  uv run python scripts/discovery_campaign_daemon.py \
    --worker-id "$w" --n-workers 4 \
    --primary Jupiter --seq-lengths 3,4 --max-rev 3 --phase-samples 24 \
    --gate-residual-kms 0.05 \
    >> out/discovery_campaign/w"$w".log 2>&1 &
done
```

- SILVER survivors append to `data/review_queue.jsonl` (the gauntlet governs; a
  human promotes — the daemon never writes the catalogue).
- A no-hit sweep appends one method-versioned record to `data/empty_regions.jsonl`.
- Checkpoints live in `out/discovery_campaign/` (gitignored); kill and restart any
  worker to resume mid-stream.

Tighten `--gate-residual-kms` and widen `--seq-lengths` / `--phase-samples` to
deepen the hunt; the capability-subsumption gate ensures a re-sweep with a stronger
method is not blocked by a weaker method's prior empty record.

## Status

- ruff check + ruff format: clean. mypy `src tests`: clean (349 files).
- `tests/search/test_discovery_campaign.py`: 13 passed.
- The open-ended search was **NOT** run (only the 10-candidate smoke).
