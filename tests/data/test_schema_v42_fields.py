"""Schema v4.2 field tests (spec §16.7.9).

Covers the three additive optional fields:
- trajectory.segments[].center (free string, absent => "Sun")
- trajectory.segments[].tof_days_bounds (published [min, max] days range)
- source_ephemeris (ephemeris model the source's numbers were computed against)

Both validation layers are exercised: the JSON Schema structural gate and the
Python semantic gate (validate_schema_invariants). The non-containment rule
(tof_days outside its bounds still passes) is regression-critical and tested
explicitly.
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


def _row_with_v42(
    *,
    center: str = "Earth",
    tof_days_bounds: Any = (161, 172),
    tof_days: Any = 146,
    source_ephemeris: Any = "DE430",
) -> dict[str, Any]:
    """Build a single-ellipse row carrying all three v4.2 fields."""
    seg: dict[str, Any] = {"id": "out", "from": "E", "to": "M", "tof_days": tof_days}
    if center is not None:
        seg["center"] = center
    if tof_days_bounds is not None:
        seg["tof_days_bounds"] = list(tof_days_bounds)
    row: dict[str, Any] = {
        "id": "v42-row",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-M",
        "orbit_elements": {"a_au": 1.5, "e": 0.3},
        "trajectory": {"segments": [seg]},
    }
    if source_ephemeris is not None:
        row["source_ephemeris"] = source_ephemeris
    return row


# ---------------------------------------------------------------------------
# (a) valid row with all three v4.2 fields passes both layers
# ---------------------------------------------------------------------------


def test_valid_v42_row_passes_both_layers() -> None:
    """A row with center, tof_days_bounds, source_ephemeris passes schema + invariants."""
    row = _row_with_v42()
    jsonschema.validate([row], _load_schema())  # must not raise
    assert validate_schema_invariants([row]) == []


# ---------------------------------------------------------------------------
# (b) tof_days_bounds with min > max fails
# ---------------------------------------------------------------------------


def test_tof_days_bounds_min_gt_max_fails_invariants() -> None:
    """tof_days_bounds with min > max triggers a semantic error."""
    row = _row_with_v42(tof_days_bounds=(172, 161))
    errs = validate_schema_invariants([row])
    assert any("min <= max" in m for m in errs), f"Expected min<=max error, got: {errs}"


# ---------------------------------------------------------------------------
# (c) bounds with wrong length fails
# ---------------------------------------------------------------------------


def test_tof_days_bounds_wrong_length_fails_invariants() -> None:
    """tof_days_bounds with three items triggers a semantic error."""
    row = _row_with_v42(tof_days_bounds=(100, 150, 200))
    errs = validate_schema_invariants([row])
    assert any("exactly 2" in m for m in errs), f"Expected length error, got: {errs}"


def test_tof_days_bounds_wrong_length_fails_schema() -> None:
    """tof_days_bounds with three items violates the JSON Schema maxItems."""
    row = _row_with_v42(tof_days_bounds=(100, 150, 200))
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([row], _load_schema())


# ---------------------------------------------------------------------------
# (d) bounds with non-positive value fails
# ---------------------------------------------------------------------------


def test_tof_days_bounds_nonpositive_fails_invariants() -> None:
    """tof_days_bounds with a zero/negative value triggers a semantic error."""
    row = _row_with_v42(tof_days_bounds=(0, 172))
    errs = validate_schema_invariants([row])
    assert any("> 0" in m for m in errs), f"Expected positivity error, got: {errs}"


def test_tof_days_bounds_nonpositive_fails_schema() -> None:
    """tof_days_bounds with a non-positive value violates exclusiveMinimum."""
    row = _row_with_v42(tof_days_bounds=(-5, 172))
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([row], _load_schema())


# ---------------------------------------------------------------------------
# (e) empty-string source_ephemeris fails
# ---------------------------------------------------------------------------


def test_empty_source_ephemeris_fails_invariants() -> None:
    """An empty-string source_ephemeris triggers a semantic error."""
    row = _row_with_v42(source_ephemeris="")
    errs = validate_schema_invariants([row])
    assert any("source_ephemeris" in m for m in errs), f"Expected se error, got: {errs}"


def test_whitespace_source_ephemeris_fails_invariants() -> None:
    """A whitespace-only source_ephemeris triggers a semantic error."""
    row = _row_with_v42(source_ephemeris="   ")
    errs = validate_schema_invariants([row])
    assert any("source_ephemeris" in m for m in errs), f"Expected se error, got: {errs}"


def test_absent_source_ephemeris_is_clean() -> None:
    """A row with no source_ephemeris key is clean."""
    row = _row_with_v42(source_ephemeris=None)
    assert validate_schema_invariants([row]) == []


# ---------------------------------------------------------------------------
# (f) REGRESSION-CRITICAL: tof_days OUTSIDE its bounds still PASSES
# ---------------------------------------------------------------------------


def test_tof_days_outside_bounds_still_passes() -> None:
    """Non-containment: tof_days=146 outside bounds [161,172] is valid (Aldrin/Rogers case).

    Different model framings of the same leg are both sourced; containment must
    NOT be enforced. This is the regression-critical guarantee of v4.2.
    """
    row = _row_with_v42(tof_days=146, tof_days_bounds=(161, 172))
    jsonschema.validate([row], _load_schema())  # must not raise
    assert validate_schema_invariants([row]) == []


# ---------------------------------------------------------------------------
# center field
# ---------------------------------------------------------------------------


def test_segment_center_is_free_string() -> None:
    """A planet-centric segment center passes both layers (no enum)."""
    row = _row_with_v42(center="Jupiter")
    jsonschema.validate([row], _load_schema())  # must not raise
    assert validate_schema_invariants([row]) == []


# ---------------------------------------------------------------------------
# (g) the real catalogue still validates cleanly under v4.2
# ---------------------------------------------------------------------------


def test_live_catalogue_passes_v42_invariants() -> None:
    """Live data/catalogue.yaml passes both v4.2 layers (no backfill yet, must stay clean)."""
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    jsonschema.validate(rows, _load_schema())
    assert validate_schema_invariants(rows) == []
