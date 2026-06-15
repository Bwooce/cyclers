#!/usr/bin/env python3
"""Schema v4.7 (task #294): one-time migration — annotate every existing row
with the orbit-class taxonomy fields under the catalogue-scope expansion.

The scope expansion (2026-06-15) admits four orbit classes: ``cycler``,
``quasi_cycler``, ``precursor_mga``, ``mga_tour``. Every pre-v4.7 row is a
strict cycler (the original cyclers-only scope), so this migration adds the
trivial backward-compatible defaults to every row that lacks them:

* ``orbit_class: cycler``
* ``epoch_locked: false``
* ``n_returns: infinite``

This is a structural annotation: no row's physics, source, or provenance
changes — only the three new explicit fields appear. The Tito 2018 row (the
single ``mga_tour`` admission) is added by a *separate* commit (admit Tito
2018) so this migration's diff stays purely "every existing row labelled
cycler".

The migration preserves comments and ordering by editing the YAML
*line-by-line* (the same idiom as ``scripts/backfill_validation_level.py``
and ``scripts/backfill_provenance_tags.py``) rather than reflowing through
PyYAML — PyYAML would strip the catalogue's inline ``# ...`` comments and
collapse the block style. The 280 ``cycler_class:`` lines are the stable
anchor: the three new fields are inserted *immediately after* each row's
``cycler_class:`` line.

Idempotent: a row that already carries ``orbit_class`` is skipped untouched
on a re-run; the script reports ``rows_already_tagged`` vs ``rows_inserted``
so a second invocation is a no-op.

After the line-insert step the script reparses the result with
``yaml.safe_load`` and validates against the updated JSON Schema
(``data/catalogue.schema.json``); a schema failure aborts the write. The
write itself is atomic (``cyclerfinder.data.catalog.atomic_write_text``) so
a process death mid-write cannot truncate the catalogue.

Usage::

    uv run python scripts/migrate_catalogue_scope_2026-06-15.py            # apply
    uv run python scripts/migrate_catalogue_scope_2026-06-15.py --dry-run  # report only
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import jsonschema  # type: ignore[import-untyped]
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import atomic_write_text

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
SCHEMA_PATH = REPO_ROOT / "data" / "catalogue.schema.json"

_INSERT_LINES = (
    "  orbit_class: cycler   "
    "# schema v4.7 (2026-06-15, task #294); scope expansion default — every "
    "pre-v4.7 row is a strict cycler (scripts/migrate_catalogue_scope_2026-06-15.py); "
    "see data/README.md\n"
    "  epoch_locked: false   "
    "# schema v4.7 (task #294); strict cyclers are NOT epoch-locked\n"
    "  n_returns: infinite   "
    "# schema v4.7 (task #294); strict cyclers complete an unbounded number of returns\n"
)


def _validate_against_schema(rows: list[dict]) -> None:
    schema = json.loads(SCHEMA_PATH.read_text())
    jsonschema.validate(rows, schema)


def main() -> None:
    dry_run = "--dry-run" in sys.argv

    raw_text = CATALOGUE_PATH.read_text(encoding="utf-8")
    rows = yaml.safe_load(raw_text)

    already_tagged: list[str] = []
    to_tag: list[str] = []
    for row in rows:
        rid = row["id"]
        if "orbit_class" in row:
            already_tagged.append(rid)
            continue
        to_tag.append(rid)

    print(f"rows total: {len(rows)}")
    print(f"rows already tagged (skipped): {len(already_tagged)}")
    print(f"rows to tag this run: {len(to_tag)}")

    if dry_run:
        print("\n[dry-run] no file written.")
        return

    if not to_tag:
        print("\nNothing to insert (every row already tagged with orbit_class).")
        return

    # Line-by-line insertion: after each row's `  cycler_class:` line, append
    # the three new field lines. `cycler_class` is universally present in the
    # current catalogue (280/280 rows) so it is a stable anchor.
    lines = raw_text.splitlines(keepends=True)
    out: list[str] = []
    current_id: str | None = None
    inserted = 0
    to_tag_set = set(to_tag)
    for line in lines:
        if line.startswith("- id:"):
            slug = line.split(":", 1)[1].strip()
            if "#" in slug:
                slug = slug[: slug.index("#")].strip()
            current_id = slug.strip().strip('"')
        out.append(line)
        if (
            line.startswith("  cycler_class:")
            and current_id is not None
            and current_id in to_tag_set
        ):
            out.append(_INSERT_LINES)
            inserted += 1

    if inserted != len(to_tag):
        print(
            f"ERROR: expected to insert {len(to_tag)} blocks, inserted {inserted}. "
            "Aborting (no file written).",
            file=sys.stderr,
        )
        sys.exit(1)

    new_text = "".join(out)

    # Reparse + schema-validate the proposed file BEFORE writing.
    try:
        new_rows = yaml.safe_load(new_text)
    except yaml.YAMLError as exc:
        print(f"ERROR: migrated catalogue does not parse as YAML: {exc!r}", file=sys.stderr)
        sys.exit(1)

    if not isinstance(new_rows, list) or len(new_rows) != len(rows):
        print(
            "ERROR: migrated catalogue row count mismatch "
            f"(before={len(rows)}, after={len(new_rows) if isinstance(new_rows, list) else 'NA'})",
            file=sys.stderr,
        )
        sys.exit(1)

    missing_after = [r["id"] for r in new_rows if "orbit_class" not in r]
    if missing_after:
        print(
            f"ERROR: {len(missing_after)} rows still lack orbit_class after migration "
            f"(first 5: {missing_after[:5]}). Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        _validate_against_schema(new_rows)
    except jsonschema.ValidationError as exc:
        print(f"ERROR: migrated catalogue fails JSON Schema: {exc.message}", file=sys.stderr)
        sys.exit(1)

    atomic_write_text(CATALOGUE_PATH, new_text)
    print(f"\nInserted orbit_class/epoch_locked/n_returns into {inserted} rows -> {CATALOGUE_PATH}")


if __name__ == "__main__":
    main()
