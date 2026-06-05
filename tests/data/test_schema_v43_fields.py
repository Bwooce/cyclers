"""Schema v4.3 field tests (spec §16.7.10).

Covers the row-supersession link pair:
- superseded_by (CCSDS NEXT_MESSAGE_ID analogue)
- supersedes    (CCSDS PREVIOUS_MESSAGE_ID analogue)

Both validation layers are exercised: the JSON Schema structural gate (array of
non-empty strings) and the Python semantic gate (validate_schema_invariants),
which adds the cross-row referential-integrity check JSON Schema cannot express:
every link target must resolve to an existing row id and must not self-reference.
The live-catalogue backfill (vem-emeeve-3syn -> attested Jones 2017 members) is
asserted to keep its links intact and resolvable.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_schema_invariants

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def _two_rows(
    *,
    superseded_by: Any = None,
    supersedes: Any = None,
) -> list[dict[str, Any]]:
    """Build a pair of minimal rows; the first may carry link fields to the second."""
    head: dict[str, Any] = {"id": "head-row", "bodies": ["E", "M"]}
    tail: dict[str, Any] = {"id": "tail-row", "bodies": ["E", "M"]}
    if superseded_by is not None:
        head["superseded_by"] = superseded_by
    if supersedes is not None:
        head["supersedes"] = supersedes
    return [head, tail]


# ---------------------------------------------------------------------------
# schema version is bumped to 4.3
# ---------------------------------------------------------------------------


def test_schema_version_is_4_3() -> None:
    assert _load_schema()["version"] == "4.3"


# ---------------------------------------------------------------------------
# (a) valid resolvable link passes both layers (both directions)
# ---------------------------------------------------------------------------


def test_valid_superseded_by_passes_both_layers() -> None:
    rows = _two_rows(superseded_by=["tail-row"])
    jsonschema.validate(rows, _load_schema())  # must not raise
    assert validate_schema_invariants(rows) == []


def test_valid_supersedes_passes_both_layers() -> None:
    rows = _two_rows(supersedes=["tail-row"])
    jsonschema.validate(rows, _load_schema())  # must not raise
    assert validate_schema_invariants(rows) == []


# ---------------------------------------------------------------------------
# (b) dangling target fails the semantic gate (referential integrity)
# ---------------------------------------------------------------------------


def test_dangling_superseded_by_fails_invariants() -> None:
    rows = _two_rows(superseded_by=["does-not-exist"])
    errs = validate_schema_invariants(rows)
    assert any("does not resolve" in m for m in errs), f"Expected dangling error, got: {errs}"


def test_dangling_supersedes_fails_invariants() -> None:
    rows = _two_rows(supersedes=["ghost-row"])
    errs = validate_schema_invariants(rows)
    assert any("does not resolve" in m for m in errs), f"Expected dangling error, got: {errs}"


# ---------------------------------------------------------------------------
# (c) self-reference fails
# ---------------------------------------------------------------------------


def test_self_reference_fails_invariants() -> None:
    rows = _two_rows(superseded_by=["head-row"])
    errs = validate_schema_invariants(rows)
    assert any("own id" in m for m in errs), f"Expected self-reference error, got: {errs}"


# ---------------------------------------------------------------------------
# (d) non-string / empty entries fail both layers
# ---------------------------------------------------------------------------


def test_empty_string_target_fails_invariants() -> None:
    rows = _two_rows(superseded_by=[""])
    errs = validate_schema_invariants(rows)
    assert any("non-empty row-id" in m for m in errs), f"Expected empty-id error, got: {errs}"


def test_empty_string_target_fails_schema() -> None:
    rows = _two_rows(superseded_by=[""])
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(rows, _load_schema())


def test_non_list_link_fails_schema() -> None:
    rows = _two_rows(superseded_by="tail-row")  # bare string, not a list
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(rows, _load_schema())


# ---------------------------------------------------------------------------
# (e) absent links are clean
# ---------------------------------------------------------------------------


def test_absent_links_are_clean() -> None:
    rows = _two_rows()
    jsonschema.validate(rows, _load_schema())  # must not raise
    assert validate_schema_invariants(rows) == []


# ---------------------------------------------------------------------------
# (f) live-catalogue backfill: vem-emeeve-3syn supersession resolves
# ---------------------------------------------------------------------------


def test_live_catalogue_passes_v43_invariants() -> None:
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    jsonschema.validate(rows, _load_schema())
    assert validate_schema_invariants(rows) == []


def test_vem_emeeve_3syn_supersession_backfill() -> None:
    """The premise-invalidated row points at its two attested Jones 2017 members,
    and both targets exist in the live catalogue (referential integrity holds)."""
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    by_id = {r["id"]: r for r in rows}
    row = by_id["vem-emeeve-3syn"]
    links = row.get("superseded_by")
    assert links == ["jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound"]
    for target in links:
        assert target in by_id, f"supersession target {target!r} missing from catalogue"
