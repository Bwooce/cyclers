"""Task #228: validation tests for the assumed-errata ledger (data/errata.yaml).

Checks:

- every entry validates against ``data/errata.schema.json`` (the same
  jsonschema pattern as ``test_jsonschema.py`` uses for the catalogue);
- entry ids are unique;
- every ``affected_catalogue_rows`` id resolves to an existing row in
  ``data/catalogue.yaml`` (referential integrity);
- every ``evidence_refs`` entry that points into the repo tree (``docs/``,
  ``tests/``, ``src/``) names a file that actually exists — pinning-test node
  ids included (path taken before ``::``);
- crafted invalid entries (unknown field, bad confidence enum) fail validation
  (``additionalProperties: false`` and the enums are live).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import jsonschema
import pytest
import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
ERRATA_PATH = REPO_ROOT / "data" / "errata.yaml"
SCHEMA_PATH = REPO_ROOT / "data" / "errata.schema.json"
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"

# evidence_refs entries pointing into the repo tree start with one of these.
_REPO_PATH_RE = re.compile(r"^(docs|tests|src)/\S+")


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def _load_errata() -> list[dict[str, Any]]:
    entries = yaml.safe_load(ERRATA_PATH.read_text())
    assert isinstance(entries, list)
    return entries


def _minimal_valid_entry() -> dict[str, Any]:
    """A minimal entry passing the schema, for crafted-failure tests."""
    return {
        "id": "example-erratum",
        "paper": "A. Author, 'Title', Venue, 2026.",
        "doi_or_url": "https://doi.org/10.0000/example",
        "location": "Table 1 (p. 1)",
        "printed_value": "1.0",
        "derived_value": "2.0",
        "reasoning": "The printed value appears inconsistent with X; we derive Y by Z.",
        "evidence_refs": ["docs/notes/example.md"],
        "confidence": "unresolved-discrepancy",
        "affected_catalogue_rows": [],
        "vor_status": "final-form",
        "status": "open",
    }


def test_errata_matches_jsonschema() -> None:
    """Live data/errata.yaml validates against data/errata.schema.json."""
    jsonschema.validate(_load_errata(), _load_schema())


def test_errata_ids_unique() -> None:
    ids = [entry["id"] for entry in _load_errata()]
    assert len(ids) == len(set(ids)), f"duplicate errata ids: {sorted(ids)}"


def test_affected_catalogue_rows_exist() -> None:
    """Every affected_catalogue_rows id resolves to a live catalogue row."""
    catalogue_ids = {row["id"] for row in yaml.safe_load(CATALOGUE_PATH.read_text())}
    for entry in _load_errata():
        for row_id in entry["affected_catalogue_rows"]:
            assert row_id in catalogue_ids, (
                f"errata entry {entry['id']!r} references unknown catalogue row {row_id!r}"
            )


def test_evidence_ref_files_exist() -> None:
    """Repo-tree evidence refs (incl. pytest node ids) name files that exist."""
    for entry in _load_errata():
        for ref in entry["evidence_refs"]:
            match = _REPO_PATH_RE.match(ref)
            if match is None:
                continue  # e.g. "commit f4138d8 (...)" — not a repo path
            # Strip a pytest node-id suffix; the file part is before "::".
            file_part = match.group(0).split("::")[0]
            assert (REPO_ROOT / file_part).is_file(), (
                f"errata entry {entry['id']!r} evidence ref {ref!r}: "
                f"file {file_part!r} does not exist"
            )


def test_unknown_field_fails_schema() -> None:
    """additionalProperties: false is live."""
    bad = _minimal_valid_entry()
    bad["surprise"] = "field"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([bad], _load_schema())


def test_bad_confidence_enum_fails_schema() -> None:
    bad = _minimal_valid_entry()
    bad["confidence"] = "pretty-sure"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([bad], _load_schema())


def test_minimal_entry_passes_schema() -> None:
    jsonschema.validate([_minimal_valid_entry()], _load_schema())
