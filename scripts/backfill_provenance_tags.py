#!/usr/bin/env python3
"""Phase 0 Task 3: back-fill source/fidelity provenance tags onto catalogue rows.

This attaches the four optional, top-level provenance tags consumed by
``cyclerfinder.data.validate.validate_provenance_tags`` —

    orbit_source    vinf_source     orbit_fidelity     vinf_fidelity

plus an optional declared ``validation_tier`` — to every ``data/catalogue.yaml``
row whose *existing metadata* can support them.  These are provenance metadata
DERIVED MECHANICALLY from fields already on the row; NO external source is
consulted and NO new physics value is introduced.  A row whose metadata cannot
support a tag is left untagged (the validator treats an absent tag as the
"unknown" marker → ``UNVALIDATED``); a source key is NEVER invented when the
row's own citations do not name a paper present in ``SOURCE_REGISTRY``.

Derivation rules (mechanical, idempotent)
-----------------------------------------
Fidelity (``orbit_fidelity`` / ``vinf_fidelity``) — from ``model_assumption``:

    circular-coplanar  -> "circular-coplanar"
    analytic-ephemeris -> "analytic-ephemeris"
    cr3bp              -> (no fidelity tag; "cr3bp" is not a Fidelity tier —
                          these are non-keplerian rows, out of the source/V_inf
                          cross-validation regime)

Both fields take the SAME fidelity because the catalogue states one
``model_assumption`` per row; the per-side split exists so a future row that
mixes fidelities can be represented, but no current row does.

Source (``vinf_source`` / ``orbit_source``) — by priority:

  1. An explicit Russell table reference in the relevant field's sourced text
     (the per-encounter ``note`` / ``source_quotes['vinf...']`` for V_inf; the
     ``orbit_elements.note`` / ``source_quotes['a_au'|'e'|'*perihelion*'|
     '*aphelion*'|'orbit_elements*']`` for the orbit).  Table number -> key:
         Table 3.4              -> russell-2004-t34
         Table 3.9 / 3.10 / 3.11-> russell-2004-t39_311
         Table 4.9 .. 4.13      -> russell-2004-t49_413
  2. A non-Russell author token in that same text mapped to a registry key:
         "Rogers"     -> rogers-2012-t1
         "Hollister"  -> hollister-1970-t3
         "McConaghy"  -> mcconaghy-2002 / mcconaghy-2006 (disambiguated by the
                         row's first_published.doi: 10.2514/1.15215 -> 2006,
                         else 2002)
         "Friedlander"-> friedlander-1986
         "spec.md"    -> spec-9
  3. Fallback to the row's ``first_published.doi`` when it maps to a registry
     key (DOI table below).
  4. Otherwise: no tag (the paper is not in SOURCE_REGISTRY — e.g. Jones 2017,
     Sanchez Net 2022, Arenstorf, Genova, Wittal, Hernandez, Russell-Strange).

When a field's own sourced text names a table/author DIFFERENT from
first_published (e.g. a row whose first_published is McConaghy 2002 but whose
V_inf note cites "Russell 2004 Table 4.9"), the per-field text wins — that is
the whole point of per-field provenance, and it is what makes a genuine
two-source cross-validation visible.

``validation_tier`` is NOT written by this script.  The tier is a *computed*
classification (``classify_validation``); writing it as data would duplicate
derived state and risk drift.  The live-row tier census (Task 4 ratchet) reads
the tags this script writes and classifies on the fly.

Insertion
---------
Tags are written as top-level row keys immediately after the ``cycler_class:``
line, matching where ``validate_provenance_tags`` reads them.  Line-by-line
text insertion (same approach as backfill_invariants.py) preserves all
comments and formatting.  Idempotent: a row that already carries any of the
four tags is skipped.

Usage::

    uv run python scripts/backfill_provenance_tags.py            # apply
    uv run python scripts/backfill_provenance_tags.py --dry-run  # report only
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"

# model_assumption -> Fidelity tier (cr3bp has no fidelity tier).
_MA_TO_FIDELITY: dict[str, str] = {
    "circular-coplanar": "circular-coplanar",
    "analytic-ephemeris": "analytic-ephemeris",
}

# first_published.doi -> registry key (fallback only).
_DOI_TO_SOURCE: dict[str, str] = {
    "10.2514/6.2012-4746": "rogers-2012-t1",
    "10.2514/3.30134": "hollister-1970-t3",
    "10.2514/6.2002-4420": "mcconaghy-2002",
    "10.2514/1.15215": "mcconaghy-2006",
    "10.2514/6.2002-4423": "mcconaghy-2002",
    "10.2514/1.5571": "russell-2004-t39_311",  # JSR venue for the Russell-Ocampo catalogue
}

_RUSSELL_TABLE_RE = re.compile(r"Table\s+(\d+)\.(\d+)")


def _table_to_source(text: str) -> str | None:
    """Return the Russell registry key for an explicit ``Table X.Y`` ref in
    *text*, or ``None``. Only Russell tables are mapped this way; the chapter
    digit selects the family (3.4 / 3.9-3.11 / 4.9-4.13)."""
    keys: set[str] = set()
    for chap_s, sub_s in _RUSSELL_TABLE_RE.findall(text):
        chap, sub = int(chap_s), int(sub_s)
        if chap == 3 and sub == 4:
            keys.add("russell-2004-t34")
        elif chap == 3 and sub in (9, 10, 11):
            keys.add("russell-2004-t39_311")
        elif chap == 4 and 9 <= sub <= 13:
            keys.add("russell-2004-t49_413")
    # Only return when the reference is unambiguous (a single family). Mixed
    # references (e.g. a "Table 4.9 ... Table 5.5" cross-note) fall through to
    # the author/DOI rules rather than guessing.
    return next(iter(keys)) if len(keys) == 1 else None


def _author_to_source(text: str, doi: str | None) -> str | None:
    """Return a registry key for a non-Russell author token in *text*, or None.

    McConaghy is disambiguated by *doi* (2006 JSR DOI -> mcconaghy-2006).
    """
    if re.search(r"Rogers", text):
        return "rogers-2012-t1"
    if re.search(r"Hollister", text):
        return "hollister-1970-t3"
    if re.search(r"Friedlander", text):
        return "friedlander-1986"
    if re.search(r"McConaghy", text):
        return "mcconaghy-2006" if doi == "10.2514/1.15215" else "mcconaghy-2002"
    if re.search(r"spec\.md", text):
        return "spec-9"
    return None


def _doi_to_source(doi: str | None) -> str | None:
    return _DOI_TO_SOURCE.get(doi) if doi else None


# Manually resolved (id, field) source overrides for the handful of rows whose
# field text names BOTH a Russell decimal table AND another author, where the
# table reference is a *secondary* "value also derivable from ..." aside rather
# than the primary citation. Each is resolved by reading the prose; the default
# table-first rule is correct for every other conflict (verified 2026-06-05).
_SOURCE_OVERRIDES: dict[tuple[str, str], str] = {
    # orbit elements are "per Rogers/Hughes/Longuski/Aldrin 2012 Table 1";
    # the "Russell 2004 Table 3.4 ... same value derivable" clause is a
    # corroboration aside, not the primary orbit-element source.
    ("aldrin-classic-em-k1-outbound", "orbit"): "rogers-2012-t1",
    # the row's own prose states "Russell 2004 Table 3.4 does not pair 'Case 1'
    # to a Russell-numbered row"; the steady-state V_inf source is McConaghy
    # 2002 (AIAA 2002-4420), matching first_published.doi.
    ("mcconaghy-2005-em-case1", "vinf"): "mcconaghy-2002",
}


def _resolve_source(field_text: str, doi: str | None, override_key: tuple[str, str]) -> str | None:
    """Apply the per-(id,field) override, then the three-tier priority
    (table ref -> author token -> DOI)."""
    if override_key in _SOURCE_OVERRIDES:
        return _SOURCE_OVERRIDES[override_key]
    return _table_to_source(field_text) or _author_to_source(field_text, doi) or _doi_to_source(doi)


def _vinf_text(row: dict) -> str:  # type: ignore[type-arg]
    parts: list[str] = []
    for enc in row.get("vinf_kms_at_encounters") or []:
        if isinstance(enc, dict) and isinstance(enc.get("note"), str):
            parts.append(enc["note"])
    sq = row.get("source_quotes") or {}
    for k, v in sq.items():
        if isinstance(v, str) and "vinf" in k:
            parts.append(v)
    return " || ".join(parts)


def _orbit_text(row: dict) -> str:  # type: ignore[type-arg]
    parts: list[str] = []
    oe = row.get("orbit_elements") or {}
    if isinstance(oe.get("note"), str):
        parts.append(oe["note"])
    sq = row.get("source_quotes") or {}
    for k, v in sq.items():
        if not isinstance(v, str):
            continue
        if (
            k == "a_au"
            or k == "e"
            or "perihelion" in k
            or "aphelion" in k
            or k.startswith("orbit_elements")
        ):
            parts.append(v)
    return " || ".join(parts)


def _tags_for_row(row: dict) -> dict[str, str]:  # type: ignore[type-arg]
    """Return the provenance tags derivable from *row*'s existing metadata.

    Keys present only when a value could be derived (untagged = unknown).
    """
    tags: dict[str, str] = {}
    ma = str(row.get("model_assumption") or "")
    fidelity = _MA_TO_FIDELITY.get(ma)
    doi = (row.get("first_published") or {}).get("doi")

    rid = str(row.get("id"))
    vinf_src = (
        _resolve_source(_vinf_text(row), doi, (rid, "vinf"))
        if row.get("vinf_kms_at_encounters")
        else None
    )
    orbit_src = _resolve_source(_orbit_text(row), doi, (rid, "orbit"))

    if vinf_src is not None:
        tags["vinf_source"] = vinf_src
        if fidelity is not None:
            tags["vinf_fidelity"] = fidelity
    if orbit_src is not None:
        tags["orbit_source"] = orbit_src
        if fidelity is not None:
            tags["orbit_fidelity"] = fidelity
    return tags


# Field emission order (deterministic YAML output).
_TAG_ORDER = ("orbit_source", "orbit_fidelity", "vinf_source", "vinf_fidelity")


def _render_tag_block(tags: dict[str, str]) -> str:
    lines: list[str] = []
    for i, key in enumerate(_TAG_ORDER):
        if key not in tags:
            continue
        comment = ""
        if i == 0 or (key == "vinf_source" and "orbit_source" not in tags):
            comment = (
                "   # schema v4.4 (2026-06-05); provenance tag back-filled from this row's "
                "existing source metadata (scripts/backfill_provenance_tags.py); see data/README.md"
            )
        lines.append(f"  {key}: {tags[key]}{comment}")
    return "\n".join(lines) + "\n" if lines else ""


def main() -> None:
    dry_run = "--dry-run" in sys.argv
    raw_text = CATALOGUE_PATH.read_text(encoding="utf-8")
    rows = yaml.safe_load(raw_text)

    tag_data: dict[str, dict[str, str]] = {}
    already: set[str] = set()
    for row in rows:
        rid = row["id"]
        if any(k in row for k in _TAG_ORDER):
            already.add(rid)
            continue
        tags = _tags_for_row(row)
        if tags:
            tag_data[rid] = tags

    # Coverage report.
    from collections import Counter

    src_counts: Counter[str] = Counter()
    for tags in tag_data.values():
        for k in ("orbit_source", "vinf_source"):
            if k in tags:
                src_counts[tags[k]] += 1
    n_orbit = sum(1 for t in tag_data.values() if "orbit_source" in t)
    n_vinf = sum(1 for t in tag_data.values() if "vinf_source" in t)
    print(f"rows total: {len(rows)}")
    print(f"already tagged (skipped): {len(already)}")
    print(f"rows to tag this run: {len(tag_data)}")
    print(f"  with orbit_source: {n_orbit}")
    print(f"  with vinf_source:  {n_vinf}")
    untagged = [r["id"] for r in rows if r["id"] not in tag_data and r["id"] not in already]
    print(f"rows left fully untagged (unknown marker): {len(untagged)}")
    print("source-key usage (orbit+vinf occurrences):")
    for key, c in src_counts.most_common():
        print(f"  {key}: {c}")
    if untagged:
        print("untagged ids:")
        for rid in untagged:
            print(f"  {rid}")

    if dry_run:
        print("\n[dry-run] no file written.")
        return
    if not tag_data:
        print("\nNothing to insert (all rows already tagged or untaggable).")
        return

    # Line-by-line insertion after each row's cycler_class: line.
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
        if line.startswith("  cycler_class:") and current_id is not None and current_id in tag_data:
            out.append(_render_tag_block(tag_data[current_id]))
            inserted += 1

    if inserted != len(tag_data):
        print(
            f"ERROR: expected to insert {len(tag_data)} tag blocks, inserted {inserted}. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    CATALOGUE_PATH.write_text("".join(out), encoding="utf-8")
    print(f"\nInserted provenance tags into {inserted} rows -> {CATALOGUE_PATH}")


if __name__ == "__main__":
    main()
