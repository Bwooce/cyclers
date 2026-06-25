"""Schema v4.10 bcr4bp_provenance field tests (task #305).

v4.10 adds ONE additive, optional, nullable nested row property:

* ``bcr4bp_provenance`` ({mu_sun, sun_commensurate_n, sun_phase_drift}) — the
  BCR4BP (Sun-Earth-Moon bicircular) periodic-orbit provenance for a future
  known-reproduction (Andreu / Rosales-Jorba lineage) row admitted via the
  BCR4BP V0-V5 gauntlet (src/cyclerfinder/data/validation/v{0,1,2,3}_bcr4bp.py).

These tests pin the schema-level contract (the nested shape + field types +
constraints), that the field is fully optional/backward-compatible, and the
CRUCIAL census invariant: NO existing live row carries it, so it changes no
census count (the additive-not-mutative guarantee the design demands).
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


def test_schema_version_is_at_least_5_0() -> None:
    # v5.0 carries the bcr4bp_provenance block (task #305). The 4.9 -> 5.0 bump
    # (not 4.10) keeps float(version) monotone for the >= ratchets.
    assert float(_load_schema()["version"]) >= 5.0


def test_bcr4bp_provenance_property_present_and_nullable() -> None:
    props = _load_schema()["items"]["properties"]
    assert "bcr4bp_provenance" in props
    prov = props["bcr4bp_provenance"]
    assert prov["type"] == ["object", "null"]
    sub = prov["properties"]
    assert set(sub) == {"mu_sun", "sun_commensurate_n", "sun_phase_drift"}
    assert sub["mu_sun"]["type"] == ["number", "null"]
    assert sub["sun_commensurate_n"]["type"] == ["integer", "null"]
    assert sub["sun_commensurate_n"]["minimum"] == 1
    assert sub["sun_phase_drift"]["type"] == ["number", "null"]
    assert sub["sun_phase_drift"]["minimum"] == 0


def test_valid_bcr4bp_provenance_passes_schema() -> None:
    schema = _load_schema()
    row = [
        {
            "id": "andreu-pol1-bcr4bp-known-reproduction",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "bcr4bp_provenance": {
                "mu_sun": 328900.5423094043,
                "sun_commensurate_n": 1,
                "sun_phase_drift": 0.0,
            },
        }
    ]
    jsonschema.validate(row, schema)  # must not raise


def test_null_bcr4bp_provenance_passes_schema() -> None:
    schema = _load_schema()
    jsonschema.validate([{"id": "x", "bcr4bp_provenance": None}], schema)


def test_legacy_row_without_field_still_passes() -> None:
    """Backward compat: a row with no bcr4bp_provenance is admitted."""
    schema = _load_schema()
    jsonschema.validate([{"id": "legacy", "bodies": ["E", "M"]}], schema)


def test_non_integer_commensurate_n_rejected() -> None:
    schema = _load_schema()
    bad = [{"id": "bad", "bcr4bp_provenance": {"sun_commensurate_n": 1.5}}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_zero_commensurate_n_rejected() -> None:
    schema = _load_schema()
    bad = [{"id": "bad", "bcr4bp_provenance": {"sun_commensurate_n": 0}}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_negative_phase_drift_rejected() -> None:
    schema = _load_schema()
    bad = [{"id": "bad", "bcr4bp_provenance": {"sun_phase_drift": -0.1}}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_live_catalogue_validates_against_v410() -> None:
    jsonschema.validate(_load_rows(), _load_schema())


def test_no_live_row_carries_bcr4bp_provenance() -> None:
    """The census invariant: the additive field is on ZERO existing rows.

    A passing BCR4BP family member is a known-reproduction flagged for human
    review, NEVER self-admitted — so until a human admits one, no live row
    carries the field, and the field changes no census count.
    """
    rows = _load_rows()
    carriers = [r["id"] for r in rows if "bcr4bp_provenance" in r]
    assert carriers == [], f"unexpected bcr4bp_provenance carriers: {carriers}"
