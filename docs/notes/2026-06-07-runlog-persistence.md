# Runlog persistence — committing the (seed → basin) labels (#152)

**Date:** 2026-06-07
**Origin:** `docs/notes/2026-06-07-ml-surrogate-investigation.md`, Verdict B
precondition 1 ("First (cheap, do anyway): make the campaign **persist** its
per-seed outcome labels"). The campaign already *computes* every row's
`CLOSE-AND-MATCH` / `CLOSE-OFF-ANCHOR` / `NO-CLOSE` label with achieved-vs-sourced
V∞ per body, then prints and discards it. This task persists those labels.
Standalone value: a versioned, grep-able audit trail of every campaign closure
and a regression baseline — independent of whether a Build-B classifier is ever
trained.

## What landed

- **`src/cyclerfinder/data/runlog.py`** — a small JSON-lines writer/reader,
  mirroring `data/ledger.py`. `RunRecord` (frozen dataclass) is the schema;
  `RunLog.append` / `.extend` / `.read` are the I/O; `default_runlog_path`
  composes `<dir>/<tag>-<timestamp>.jsonl`.
- **`scripts/campaign_russell12.py`** — additive `--runlog-dir` flag (default ON,
  writes `data/runs/`) plus `--runlog-timestamp`. Writes one record per row for
  BOTH genome paths (lambert and free-return). Scientific behaviour and printed
  output are unchanged; the only new stdout line is `wrote runlog …`.
- **`tests/data/test_runlog.py`** — round-trip, schema validation (reject each
  missing required field, malformed JSON, non-object lines), append-across-reopen,
  `extend` ordering, blank-line skipping, and a `@pytest.mark.slow` smoke that the
  campaign actually writes ≥1 valid record.

## Schema (one JSON-lines record per row)

Required (a line missing any is rejected on read):
`row_id`, `genome` (`lambert`/`free-return`), `outcome`, `model`
(`circular`/`astropy`/…), `code_version` (`git rev-parse --short HEAD`).

Optional/additive (default empty so leaner callers and older lines round-trip):
`achieved_vinf_kms`, `sourced_vinf_kms`, `sourced_anchors`, `seed`,
`residual_kms`, `solver_audit` (checks / close counts / mode / feasibility),
`t_written`.

## Decisions

- **`code_version` is captured by the CALLER** (`scripts/...` shells out to git)
  and passed in; the `runlog` module never shells out. Keeps the data layer pure
  and testable, and avoids a subprocess dependency.
- **`data/runs/` IS COMMITTED** (not gitignored). The records are small, they are
  evidence artefacts, and the repo already commits `data/gauntlet_ledger.jsonl`.
  `.gitignore` only excludes `out/` (the throwaway run sink, spec §7) — there is
  no existing convention saying campaign labels are ignored. Committing is what
  lets the (seed → basin) corpus accumulate across runs, which is the whole point
  of the Build-B precondition.
- **Plain buffered append, not the ledger's `PIPE_BUF` atomic path.** Run records
  carry V∞ lists / per-check strings / a solver-audit object and can exceed 4 KiB,
  and the campaign writer is single-process (one file per invocation), so the
  atomic single-line guarantee is unnecessary.

## First committed runlog (backfill)

Re-ran ONLY the circular like-for-like campaign (the fast #106 gate):

```
uv run python scripts/campaign_russell12.py --model circular --runlog-timestamp 2026-06-07
```

→ `data/runs/russell12-circular-2026-06-07.jsonl`, **12 records**
(9 `CLOSE-OFF-ANCHOR`, 3 `NO-CLOSE` — matching the printed summary). Nothing
heavy was re-run.
