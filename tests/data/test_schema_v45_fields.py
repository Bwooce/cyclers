"""Schema v4.5 field tests (spec §16.7.12).

Covers the additive top-level ``validation_level`` enum (V0-V5, the spec §14
gauntlet level) back-filled mechanically by
``scripts/backfill_validation_level.py``.

Both gates are exercised: the JSON Schema structural gate (enum membership) and
the Python semantic gate (``validate_validation_level``), which adds the rule
JSON Schema cannot express — a row may declare a level above V0 only when the
recorded mechanical evidence justifies it (the over-claim guard).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_validation_level

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def _load_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_4_5() -> None:
    assert _load_schema()["version"] == "4.5"


def test_valid_level_passes_schema() -> None:
    schema = _load_schema()
    jsonschema.validate([{"id": "ok", "validation_level": "V1"}], schema)  # must not raise


def test_unknown_level_fails_schema() -> None:
    schema = _load_schema()
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([{"id": "bad", "validation_level": "V9"}], schema)


def test_live_catalogue_validates() -> None:
    jsonschema.validate(_load_rows(), _load_schema())


def test_semantic_gate_rejects_unjustified_level_above_v0() -> None:
    """A row claiming V1+ without recorded mechanical evidence is refused."""
    rows = [{"id": "some-random-row", "validation_level": "V1"}]
    errors = validate_validation_level(rows)
    assert errors, "expected an over-claim violation"
    assert any("validation_level" in e for e in errors)


def test_semantic_gate_accepts_v0_anywhere() -> None:
    """V0 (internal-consistency floor) needs no evidence pointer."""
    rows = [{"id": "anything", "validation_level": "V0"}]
    assert validate_validation_level(rows) == []


def test_semantic_gate_accepts_known_justified_row() -> None:
    """The Aldrin outbound row's V1 is in the evidence registry."""
    rows = [{"id": "aldrin-classic-em-k1-outbound", "validation_level": "V1"}]
    assert validate_validation_level(rows) == []


def test_live_aldrin_outbound_is_v1_everything_else_v0() -> None:
    """Mechanical §14 application: exactly one row (Aldrin outbound) is V1;
    every other tagged row is V0; no row claims V2+."""
    rows = _load_rows()
    byid = {r["id"]: r.get("validation_level") for r in rows}
    assert byid.get("aldrin-classic-em-k1-outbound") == "V1"
    above_v0 = {rid: lvl for rid, lvl in byid.items() if lvl not in (None, "V0")}
    assert above_v0 == {"aldrin-classic-em-k1-outbound": "V1"}, above_v0
    # No row may carry V2 or higher (no mechanical evidence exists today).
    assert not any(lvl in ("V2", "V3", "V4", "V5") for lvl in byid.values())


def test_live_catalogue_validation_level_semantic_clean() -> None:
    """The whole live catalogue passes the validation_level semantic gate."""
    assert validate_validation_level(_load_rows()) == []
