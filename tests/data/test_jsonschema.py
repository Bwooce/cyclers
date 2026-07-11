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

import jsonschema
import pytest
import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_current() -> None:
    """The schema carries version == '5.1' (v5.1 additively names three
    synodic-feasibility sub-fields of validity_window — synodic_duty_cycle_pct,
    synodic_boundary_period_days, synodic_period_days — formalizing the free-form
    extras task #569's six Uranian moon-pair quasi_cycler rows already carried
    under validity_window's additionalProperties:true; no existing row becomes
    invalid and no census count changes. v5.0 added the additive optional
    nullable nested bcr4bp_provenance block {mu_sun, sun_commensurate_n,
    sun_phase_drift} for a future Andreu/Rosales-Jorba known-reproduction row
    admitted via the BCR4BP V0-V5 gauntlet, task #305 — absent on every existing
    row so no census count changes. Bumped 4.9 -> 5.0 (not 4.10) because the
    version ratchets compare float(version) and float('4.10')==4.1 < 4.9 would
    regress the >= 4.x gates; the change itself is purely additive. v4.9 added
    the resonant_po orbit_class enum value for stable resonant/libration POs with
    no transport utility, task #453; v4.8 added the Axis-B dv_band enum + its
    mandatory dv_band_source companion for the real-ephemeris maintenance-ΔV band
    taxonomy, task #417; v4.7 added the four-class orbit_class taxonomy for the
    catalogue-scope expansion, task #294)."""
    schema = _load_schema()
    assert schema["version"] == "5.1"


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
