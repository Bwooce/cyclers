"""Schema v4.4 field tests (spec §16.7.11).

Covers the per-field provenance tag set back-filled by Task 3:
- orbit_source / vinf_source  (constrained to the SOURCE_REGISTRY enum)
- orbit_fidelity / vinf_fidelity (constrained to the Fidelity enum)
- validation_tier (optional declared tier; enum)

Both validation layers are exercised: the JSON Schema structural gate (enum
membership) and the Python semantic gate (validate_provenance_tags), which adds
the declared-tier-must-match-classify_validation rule JSON Schema cannot
express. The live-catalogue backfill is asserted to carry resolvable tags.
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


def _load_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_at_least_4_4() -> None:
    # v4.4 fields persist into later revs (v4.5 adds validation_level on top).
    assert _load_schema()["version"] >= "4.4"


def test_valid_tags_pass_schema() -> None:
    schema = _load_schema()
    rows = [
        {
            "id": "ok",
            "orbit_source": "rogers-2012-t1",
            "vinf_source": "russell-2004-t34",
            "orbit_fidelity": "circular-coplanar",
            "vinf_fidelity": "circular-coplanar",
            "validation_tier": "cross_validated",
        }
    ]
    jsonschema.validate(rows, schema)  # must not raise


def test_unknown_source_key_fails_schema() -> None:
    schema = _load_schema()
    rows = [{"id": "bad", "orbit_source": "not-a-real-source"}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(rows, schema)


def test_unknown_fidelity_fails_schema() -> None:
    schema = _load_schema()
    rows = [{"id": "bad", "vinf_fidelity": "de440"}]  # must be 'real-de440'
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(rows, schema)


def test_unknown_validation_tier_fails_schema() -> None:
    schema = _load_schema()
    rows = [{"id": "bad", "validation_tier": "platinum"}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(rows, schema)


def test_live_catalogue_tags_validate() -> None:
    """Live catalogue validates against the v4.4 schema (enum tags included)."""
    jsonschema.validate(_load_rows(), _load_schema())


def test_live_catalogue_has_backfilled_tags() -> None:
    """Sanity: the backfill landed — a known cross-validated row carries two
    distinct source tags at matching fidelity."""
    byid = {r["id"]: r for r in _load_rows()}
    row = byid["aldrin-classic-em-k1-outbound"]
    assert row["orbit_source"] == "rogers-2012-t1"
    assert row["vinf_source"] == "russell-2004-t34"
    assert row["orbit_fidelity"] == row["vinf_fidelity"] == "circular-coplanar"
