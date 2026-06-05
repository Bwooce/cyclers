"""Task 1.4: JSON Schema structural validation tests.

Tests that:
- The live catalogue validates against data/catalogue.schema.json
- The schema version is 4.1 (schema v4.1 adds free_return_arcs per spec §16.7.7)
- A crafted invalid row (multi-arc with a_au) fails validation
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]
import pytest
import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_4_2() -> None:
    """The schema carries version == '4.2' (v4.2 adds center, tof_days_bounds, source_ephemeris)."""
    schema = _load_schema()
    assert schema["version"] == "4.2"


def test_catalogue_matches_jsonschema() -> None:
    """Live data/catalogue.yaml validates against catalogue.schema.json."""
    schema = _load_schema()
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    # Raises jsonschema.ValidationError on failure
    jsonschema.validate(rows, schema)


def test_multi_arc_with_a_au_fails_schema() -> None:
    """A multi-arc row with orbit_elements.a_au violates the conditional gate."""
    schema = _load_schema()
    # A single bad row as the whole array
    bad_rows = [
        {
            "id": "bad-multi-arc",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "cycler_class": "multi-arc",
            "orbit_elements": {"a_au": 1.6, "e": 0.3},
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


def test_single_ellipse_row_passes_schema() -> None:
    """A minimal single-ellipse row with a_au passes the schema."""
    schema = _load_schema()
    ok_rows = [
        {
            "id": "ok-single-ellipse",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_elements": {"a_au": 1.5, "e": 0.3},
        }
    ]
    # Should not raise
    jsonschema.validate(ok_rows, schema)


def test_legacy_row_without_cycler_class_passes_schema() -> None:
    """A v3 row with no cycler_class key still passes the schema (permissive)."""
    schema = _load_schema()
    ok_rows = [
        {
            "id": "legacy-row",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_elements": {"a_au": 1.3, "e": 0.257},
        }
    ]
    jsonschema.validate(ok_rows, schema)
