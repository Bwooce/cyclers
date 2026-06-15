"""Schema v4.7 orbit-class taxonomy tests (task #294).

v4.7 admits three additive optional top-level row fields under the
catalogue-scope expansion (cycler / quasi_cycler / precursor_mga / mga_tour):

* ``orbit_class`` (enum, default ``cycler``)
* ``epoch_locked`` (bool, default ``false``)
* ``n_returns`` (integer or the literal string ``"infinite"``, default
  ``"infinite"``)

…plus three optional epoch-locked fields:

* ``validity_window`` ({start, end} ISO-8601 timestamps)
* ``launch_epoch`` (ISO-8601 timestamp)
* ``inserts_into`` (catalogue row id; required when
  ``orbit_class='precursor_mga'`` — enforced by the Python semantic gate,
  not the JSON Schema's cross-row referential integrity)

These tests pin the schema-level contract: the enum + bool + oneOf shapes,
their defaults, and that the live (migrated) catalogue is uniformly tagged
``cycler`` / ``false`` / ``"infinite"`` (the one-time migration in
``scripts/migrate_catalogue_scope_2026-06-15.py``). A separate commit
admits Tito 2018 as the first ``mga_tour`` row; that one row is the only
non-cycler in the catalogue when this test runs.
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


def _load_catalogue() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_at_least_4_7() -> None:
    assert float(_load_schema()["version"]) >= 4.7


def test_orbit_class_enum_present() -> None:
    schema = _load_schema()
    props = schema["items"]["properties"]
    assert "orbit_class" in props
    enum = props["orbit_class"]["enum"]
    assert set(enum) == {"cycler", "quasi_cycler", "precursor_mga", "mga_tour"}
    assert props["orbit_class"]["default"] == "cycler"


def test_epoch_locked_bool_present() -> None:
    schema = _load_schema()
    props = schema["items"]["properties"]
    assert props["epoch_locked"]["type"] == "boolean"
    assert props["epoch_locked"]["default"] is False


def test_n_returns_oneof_int_or_infinite() -> None:
    schema = _load_schema()
    props = schema["items"]["properties"]
    one_of = props["n_returns"]["oneOf"]
    # Exactly two branches: positive integer OR the string literal "infinite".
    int_branch = next(b for b in one_of if b.get("type") == "integer")
    str_branch = next(b for b in one_of if b.get("type") == "string")
    assert int_branch["minimum"] == 1
    assert str_branch["enum"] == ["infinite"]
    assert props["n_returns"]["default"] == "infinite"


def test_optional_epoch_locked_fields_present() -> None:
    schema = _load_schema()
    props = schema["items"]["properties"]
    for key in ("validity_window", "launch_epoch", "inserts_into"):
        assert key in props, f"{key} missing from schema properties"


# --- Live catalogue invariants ---


def test_live_catalogue_every_row_carries_orbit_class() -> None:
    """Every row in the migrated catalogue declares all three new fields."""
    rows = _load_catalogue()
    missing = [r["id"] for r in rows if "orbit_class" not in r]
    assert not missing, f"{len(missing)} rows missing orbit_class: first 5 {missing[:5]}"
    missing_locked = [r["id"] for r in rows if "epoch_locked" not in r]
    assert not missing_locked, (
        f"{len(missing_locked)} rows missing epoch_locked: first 5 {missing_locked[:5]}"
    )
    missing_returns = [r["id"] for r in rows if "n_returns" not in r]
    assert not missing_returns, (
        f"{len(missing_returns)} rows missing n_returns: first 5 {missing_returns[:5]}"
    )


def test_cycler_rows_are_not_epoch_locked() -> None:
    """orbit_class=cycler MUST have epoch_locked=false and n_returns='infinite'."""
    rows = _load_catalogue()
    bad = [
        r["id"]
        for r in rows
        if r.get("orbit_class") == "cycler"
        and (r.get("epoch_locked") is not False or r.get("n_returns") != "infinite")
    ]
    assert not bad, f"cycler rows must be epoch_locked=false / n_returns=infinite: {bad}"


def test_non_cycler_rows_are_epoch_locked_with_finite_returns() -> None:
    """quasi_cycler / precursor_mga / mga_tour rows MUST be epoch_locked with a finite n_returns."""
    rows = _load_catalogue()
    bad = []
    for r in rows:
        cls = r.get("orbit_class")
        if cls in ("quasi_cycler", "precursor_mga", "mga_tour"):
            if r.get("epoch_locked") is not True:
                bad.append((r["id"], cls, "epoch_locked must be true"))
            n = r.get("n_returns")
            if not isinstance(n, int) or n < 1:
                bad.append((r["id"], cls, f"n_returns must be a positive integer, got {n!r}"))
    assert not bad, f"epoch-locked classes failed invariants: {bad}"


def test_catalogue_validates_against_v47_schema() -> None:
    """The migrated catalogue still validates end-to-end against the v4.7 schema."""
    schema = _load_schema()
    rows = _load_catalogue()
    jsonschema.validate(rows, schema)


# --- Fault injection: bad enum values must be REJECTED by JSON Schema ---


@pytest.mark.parametrize(
    "bad_class",
    [
        "Cycler",  # wrong case
        "tour",  # not the full mga_tour
        "free-return",  # made up
        "",  # empty
    ],
)
def test_bad_orbit_class_rejected(bad_class: str) -> None:
    schema = _load_schema()
    bad_rows = [
        {
            "id": "bad-orbit-class",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_class": bad_class,
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


def test_n_returns_zero_rejected() -> None:
    schema = _load_schema()
    bad_rows = [
        {
            "id": "bad-n-returns-zero",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_class": "mga_tour",
            "epoch_locked": True,
            "n_returns": 0,  # must be >= 1
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


def test_n_returns_string_other_than_infinite_rejected() -> None:
    schema = _load_schema()
    bad_rows = [
        {
            "id": "bad-n-returns-str",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_class": "cycler",
            "n_returns": "forever",  # only "infinite" is admitted
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


def test_legacy_row_without_orbit_class_still_passes() -> None:
    """Schema must remain backward-compatible: a v4.6-shape row (no orbit_class) is admitted.

    Important for staged rollouts and for ad-hoc rows in fault-injection tests.
    """
    schema = _load_schema()
    legacy_row = [
        {
            "id": "legacy-no-orbit-class",
            "bodies": ["E", "M"],
            "sequence_canonical": "E-M",
            "orbit_elements": {"a_au": 1.5, "e": 0.3},
        }
    ]
    jsonschema.validate(legacy_row, schema)
