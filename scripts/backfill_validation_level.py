#!/usr/bin/env python3
"""Schema v4.5 (spec §16.7.12 / §14): back-fill the ``validation_level`` tag.

Applies spec §14's V0-V5 gauntlet definitions MECHANICALLY to ``data/catalogue.yaml``
from RECORDED in-repo test evidence — never aspirationally. The level a row earns
is the highest §14 gate for which a green, teeth-bearing test exists; when the
recorded evidence does not mechanically justify a higher gate, the row stays V0
(the internal-consistency floor). NO new physics value is introduced and NO
external source is consulted; this is a derived-metadata stamp keyed to the
evidence registry below.

Evidence applied (verified green against the live suite, 2026-06-06)
-------------------------------------------------------------------
* **V1** — ``aldrin-classic-em-k1-outbound``: the real-DE440 Aldrin cycler clears
  spec §14 V1 — every leg re-solved with ``lamberthub`` izzo2015 + gooding1990
  agrees to < ``V1_TOLERANCE_MPS``, AND the Kepler forward re-propagation residual
  passes (the §14 "re-propagated with the Kepler propagator, planet positions met
  < tol" half). Demonstrated with teeth by
  ``tests/verify/test_agreement_lamberthub.py::test_report_includes_lamberthub_path``
  and ``...::test_real_eph_paths_a_and_c_pass_b_flags_model_mismatch``.

* **V0** — the rows that have actually been exercised by the V0/real-closure
  gauntlet machinery (the M6b regression set + the two CONSTRUCTIBLE rows). Each
  carries internal-consistency evidence (constructed + propagated), but no
  recorded evidence reaches a higher §14 gate:
    - the Aldrin INBOUND row is powered and no test builds/cross-checks it on real
      ephemeris (its outbound twin is the binding gate) ⇒ V0;
    - the other regression rows are EXPECTED_SKIPS (incomplete leg data / wrong
      topology) and do not pass real-closure ⇒ V0;
    - the two CONSTRUCTIBLE rows are the Aldrin pair (already covered).

Everything NOT in the registry is left UNTAGGED — an absent ``validation_level``
is the explicit V0 floor for downstream views. Tagging only evidence-backed rows
keeps the stamp auditable and avoids asserting a gauntlet pass for the ~230 rows
the auto-pipeline has not run.

Idempotent: a row that already carries ``validation_level`` is skipped.

Usage::

    uv run python scripts/backfill_validation_level.py            # apply
    uv run python scripts/backfill_validation_level.py --dry-run  # report only
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"

# id -> validation_level, applied mechanically from recorded test evidence.
# V1 for the Aldrin outbound (real-DE440 lamberthub + Kepler reprop); V0 for the
# rows the gauntlet machinery has exercised at the internal-consistency floor.
_LEVEL_BY_ID: dict[str, str] = {
    "aldrin-classic-em-k1-outbound": "V1",
    "aldrin-classic-em-k1-inbound": "V1",  # #125 Part 1: real-DE440 inbound clears §14 V1
    "mcconaghy-2006-em-k2": "V0",
    "russell-ocampo-2.1.1+2-case2": "V0",
    "russell-ocampo-2.5.1+0": "V0",
}


def _render_line(level: str) -> str:
    return (
        f"  validation_level: {level}"
        "   # schema v4.5 (2026-06-06, spec §16.7.12 / §14); MECHANICAL gauntlet "
        "level from recorded test evidence (scripts/backfill_validation_level.py); "
        "see data/README.md\n"
    )


def main() -> None:
    dry_run = "--dry-run" in sys.argv
    raw_text = CATALOGUE_PATH.read_text(encoding="utf-8")
    rows = yaml.safe_load(raw_text)

    to_tag: dict[str, str] = {}
    already: set[str] = set()
    for row in rows:
        rid = row["id"]
        if rid not in _LEVEL_BY_ID:
            continue
        if "validation_level" in row:
            already.add(rid)
            continue
        to_tag[rid] = _LEVEL_BY_ID[rid]

    print(f"rows total: {len(rows)}")
    print(f"evidence-registry rows: {len(_LEVEL_BY_ID)}")
    print(f"already tagged (skipped): {sorted(already)}")
    print(f"rows to tag this run: {len(to_tag)}")
    for rid, lvl in sorted(to_tag.items()):
        print(f"  {rid}: {lvl}")

    if dry_run:
        print("\n[dry-run] no file written.")
        return
    if not to_tag:
        print("\nNothing to insert (all evidence-registry rows already tagged).")
        return

    # Line-by-line insertion after each row's cycler_class: line (matches the
    # provenance backfill's insertion site; preserves all comments/formatting).
    lines = raw_text.splitlines(keepends=True)
    out: list[str] = []
    current_id: str | None = None
    inserted = 0
    for line in lines:
        stripped = line.lstrip()
        if line.startswith("- id:"):
            slug = stripped.split(":", 1)[1].strip()
            if "#" in slug:
                slug = slug[: slug.index("#")].strip()
            current_id = slug.strip().strip('"')
        out.append(line)
        if line.startswith("  cycler_class:") and current_id is not None and current_id in to_tag:
            out.append(_render_line(to_tag[current_id]))
            inserted += 1

    if inserted != len(to_tag):
        print(
            f"ERROR: expected to insert {len(to_tag)} lines, inserted {inserted}. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    CATALOGUE_PATH.write_text("".join(out), encoding="utf-8")
    print(f"\nInserted validation_level into {inserted} rows -> {CATALOGUE_PATH}")


if __name__ == "__main__":
    main()
