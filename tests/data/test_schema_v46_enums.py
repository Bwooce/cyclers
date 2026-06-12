"""Schema v4.6 controlled-vocabulary enum tests (task #194).

v4.6 hardens four previously-freeform discriminator string fields with
controlled-vocabulary enums:
  - model_assumption       (top-level)
  - trajectory_regime      (top-level)
  - sense                  (top-level)
  - data_gaps[].kind       (nested)

Discovered during #184: model_assumption='bicircular' and
data_gaps[].kind='not-applicable' were accepted with zero validation, so a typo
(cr3pb, unkown) would pass silently and corrupt the catalogue.

Fault-injection pattern: each enum must REJECT a fabricated known-bad value AND
ACCEPT every value currently present in the live catalogue.
"""

from __future__ import annotations

import collections
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


def test_schema_version_is_4_6() -> None:
    assert _load_schema()["version"] == "4.6"


# --- Fault injection: known-bad values must be REJECTED ---


@pytest.mark.parametrize(
    "field,bad_value",
    [
        ("model_assumption", "cr3pb"),  # transposed cr3bp
        ("model_assumption", "bicirular"),  # dropped 'c'
        ("trajectory_regime", "balistic"),  # dropped 'l'
        ("sense", "outbund"),  # dropped 'o'
    ],
)
def test_top_level_enum_rejects_typo(field: str, bad_value: str) -> None:
    schema = _load_schema()
    bad_rows = [{"id": "fault-injection", field: bad_value}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


def test_data_gaps_kind_rejects_typo() -> None:
    schema = _load_schema()
    bad_rows = [
        {
            "id": "fault-injection",
            "data_gaps": [{"kind": "unkown"}],  # the canonical #184 typo
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad_rows, schema)


# --- Acceptance: every live value must PASS ---


def test_enums_accept_every_live_value() -> None:
    """The new enums must be the exact superset of the live catalogue: every
    distinct value currently present validates (a single bad enum would reject
    a real row)."""
    schema = _load_schema()
    rows = _load_rows()

    ma = sorted({r.get("model_assumption") for r in rows}, key=str)
    tr = sorted({r.get("trajectory_regime") for r in rows}, key=str)
    se = sorted({r.get("sense") for r in rows}, key=str)
    kinds: set[Any] = set()
    for r in rows:
        gaps = r.get("data_gaps")
        if isinstance(gaps, list):
            for g in gaps:
                if isinstance(g, dict):
                    kinds.add(g.get("kind"))

    # Each live value, exercised through the schema in isolation, must validate.
    for v in ma:
        jsonschema.validate([{"id": "live", "model_assumption": v}], schema)
    for v in tr:
        jsonschema.validate([{"id": "live", "trajectory_regime": v}], schema)
    for v in se:
        jsonschema.validate([{"id": "live", "sense": v}], schema)
    for v in sorted(kinds, key=str):
        jsonschema.validate([{"id": "live", "data_gaps": [{"kind": v}]}], schema)


def test_full_live_catalogue_still_validates() -> None:
    """The whole live catalogue validates end-to-end against v4.6 (the binding gate:
    if any enum is too narrow, this fails)."""
    schema = _load_schema()
    rows = _load_rows()
    jsonschema.validate(rows, schema)


def test_live_value_sets_are_subset_of_enums() -> None:
    """Belt-and-braces: assert the live distinct sets are subsets of the declared
    enums, so the enums are documented to track reality."""
    schema = _load_schema()
    rows = _load_rows()
    props = schema["items"]["properties"]

    ma_enum = set(props["model_assumption"]["enum"])
    tr_enum = set(props["trajectory_regime"]["enum"])
    se_enum = set(props["sense"]["enum"])
    kind_enum = set(props["data_gaps"]["items"]["properties"]["kind"]["enum"])

    live_ma = collections.Counter(r.get("model_assumption") for r in rows)
    live_tr = collections.Counter(r.get("trajectory_regime") for r in rows)
    live_se = collections.Counter(r.get("sense") for r in rows)
    live_kind: collections.Counter[Any] = collections.Counter()
    for r in rows:
        gaps = r.get("data_gaps")
        if isinstance(gaps, list):
            for g in gaps:
                if isinstance(g, dict):
                    live_kind[g.get("kind")] += 1

    assert set(live_ma) <= ma_enum, set(live_ma) - ma_enum
    assert set(live_tr) <= tr_enum, set(live_tr) - tr_enum
    assert set(live_se) <= se_enum, set(live_se) - se_enum
    assert set(live_kind) <= kind_enum, set(live_kind) - kind_enum
