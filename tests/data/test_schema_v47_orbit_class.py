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

import jsonschema
import pytest
import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def _load_catalogue() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def _is_epoch_free_cr3bp_corridor_subclass(row: dict[str, Any]) -> bool:
    """Schema v5.2 (task #684): the epoch-free CR3BP KAM-corridor ``quasi_cycler``
    subclass -- ``orbit_class='quasi_cycler'`` AND ``model_assumption='cr3bp'`` AND
    ``validity_window`` null/absent. A KAM torus around a linearly-stable CR3BP periodic
    orbit (#682's cycler-corridor census) winds quasi-periodically FOREVER with no real-
    ephemeris calendar window, so this narrow, conditional carve-out is required instead
    to be ``epoch_locked=False`` / ``n_returns='infinite'`` -- like ``cycler``/
    ``resonant_po`` -- rather than fabricating a window and a finite return count.
    ``precursor_mga``/``mga_tour`` are NEVER eligible for this exception (only
    ``quasi_cycler``), and a ``quasi_cycler`` row that DOES carry a real
    ``validity_window`` (e.g. every #569 Uranian moon-pair row) is NEVER eligible either,
    even if it happens to also be ``model_assumption='cr3bp'`` -- both conditions must
    hold simultaneously.
    """
    return (
        row.get("orbit_class") == "quasi_cycler"
        and row.get("model_assumption") == "cr3bp"
        and row.get("validity_window") is None
    )


def test_schema_version_is_at_least_4_7() -> None:
    assert float(_load_schema()["version"]) >= 4.7


def test_orbit_class_enum_present() -> None:
    schema = _load_schema()
    props = schema["items"]["properties"]
    assert "orbit_class" in props
    enum = props["orbit_class"]["enum"]
    assert set(enum) == {
        "cycler",
        "quasi_cycler",
        "precursor_mga",
        "mga_tour",
        "resonant_po",
        "torus_homoclinic",
    }
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
    """orbit_class in {cycler, resonant_po} MUST have epoch_locked=false and n_returns='infinite'.

    resonant_po (schema v4.9, #453) is a stable resonant/libration periodic orbit: it shares
    the cycler reachability/periodicity invariants (any epoch, infinite returns) but has no
    transport utility. Both classes are pinned to the not-epoch-locked side here.
    """
    rows = _load_catalogue()
    not_epoch_locked = {"cycler", "resonant_po"}
    bad = [
        r["id"]
        for r in rows
        if r.get("orbit_class") in not_epoch_locked
        and (r.get("epoch_locked") is not False or r.get("n_returns") != "infinite")
    ]
    assert not bad, (
        f"cycler/resonant_po rows must be epoch_locked=false / n_returns=infinite: {bad}"
    )


def test_non_cycler_rows_are_epoch_locked_with_finite_returns() -> None:
    """quasi_cycler / precursor_mga / mga_tour / torus_homoclinic rows MUST be
    epoch_locked with a finite n_returns.

    Exception (schema v5.2, task #684): a row matching
    :func:`_is_epoch_free_cr3bp_corridor_subclass` (the epoch-free CR3BP KAM-corridor
    ``quasi_cycler`` subclass) is instead required to be ``epoch_locked=False`` /
    ``n_returns='infinite'`` -- see that helper's docstring. Every other quasi_cycler row,
    and every precursor_mga/mga_tour row regardless of model_assumption, is unaffected by
    this carve-out and still hits the original hard requirement below.

    ``torus_homoclinic`` (schema v5.3, task #707/#708) joins this group: a CCR4BP
    torus-homoclinic connection is real-epoch-anchored (epoch_locked=true) and is a
    ONE-SHOT departure-and-return opportunity per recurring synodic window
    (n_returns=1) -- it is NEVER eligible for the epoch-free CR3BP carve-out above
    (that carve-out is orbit_class='quasi_cycler'-scoped only).
    """
    rows = _load_catalogue()
    bad = []
    for r in rows:
        cls = r.get("orbit_class")
        if cls not in ("quasi_cycler", "precursor_mga", "mga_tour", "torus_homoclinic"):
            continue
        if _is_epoch_free_cr3bp_corridor_subclass(r):
            if r.get("epoch_locked") is not False or r.get("n_returns") != "infinite":
                bad.append(
                    (
                        r["id"],
                        cls,
                        "epoch-free CR3BP quasi_cycler subclass must be "
                        "epoch_locked=false / n_returns='infinite'",
                    )
                )
            continue
        if r.get("epoch_locked") is not True:
            bad.append((r["id"], cls, "epoch_locked must be true"))
        n = r.get("n_returns")
        if not isinstance(n, int) or n < 1:
            bad.append((r["id"], cls, f"n_returns must be a positive integer, got {n!r}"))
    assert not bad, f"epoch-locked classes failed invariants: {bad}"


# --- Schema v5.2 (#684): the epoch-free CR3BP KAM-corridor quasi_cycler carve-out ---


def _base_quasi_cycler_row(**overrides: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "id": "synthetic-quasi-cycler",
        "bodies": ["E", "Moon"],
        "sequence_canonical": "E-Moon-E",
        "orbit_class": "quasi_cycler",
        "model_assumption": "cr3bp",
        "epoch_locked": False,
        "n_returns": "infinite",
        "validity_window": None,
    }
    row.update(overrides)
    return row


def test_epoch_free_cr3bp_quasi_cycler_subclass_permitted() -> None:
    """The new subclass (cr3bp + no validity_window) IS recognised, and its natural,
    honest field values (epoch_locked=false / n_returns='infinite') are what the
    subclass predicate is designed to admit -- confirms the amendment actually accepts
    #682's corridor rows rather than merely failing to reject them."""
    row = _base_quasi_cycler_row(epoch_locked=False, n_returns="infinite", validity_window=None)
    assert _is_epoch_free_cr3bp_corridor_subclass(row)
    # A row satisfying the predicate with these exact field values is precisely what
    # test_non_cycler_rows_are_epoch_locked_with_finite_returns treats as compliant
    # (it only flags a predicate-match row when epoch_locked/n_returns deviate).
    assert row["epoch_locked"] is False
    assert row["n_returns"] == "infinite"


def test_real_ephemeris_cr3bp_quasi_cycler_still_epoch_locked() -> None:
    """REGRESSION GUARD: a quasi_cycler row that carries a real validity_window must still
    be forced epoch_locked=true even if model_assumption='cr3bp' -- the carve-out requires
    BOTH conditions (cr3bp AND no validity_window), not cr3bp alone. Mirrors what would
    happen if a future epoch-robustness scan (like #567/#568) were run on a corridor row."""
    row = _base_quasi_cycler_row(
        epoch_locked=False,
        n_returns="infinite",
        validity_window={"start": "2000-01-01T00:00:00Z", "end": "2050-01-01T00:00:00Z"},
    )
    assert not _is_epoch_free_cr3bp_corridor_subclass(row), (
        "a quasi_cycler row with a real validity_window must NOT be treated as the "
        "epoch-free subclass, even when model_assumption='cr3bp'"
    )


def test_non_cr3bp_quasi_cycler_still_epoch_locked() -> None:
    """REGRESSION GUARD: a quasi_cycler row with NO validity_window but a non-cr3bp
    model_assumption (e.g. a future analytic-ephemeris family-seed row) must NOT slip
    through the carve-out -- both conditions are required, not either."""
    row = _base_quasi_cycler_row(model_assumption="analytic-ephemeris", validity_window=None)
    assert not _is_epoch_free_cr3bp_corridor_subclass(row), (
        "a non-cr3bp quasi_cycler row must NOT be treated as the epoch-free subclass"
    )


def test_precursor_mga_and_mga_tour_never_eligible_for_epoch_free_carveout() -> None:
    """REGRESSION GUARD: precursor_mga/mga_tour rows are NEVER exempt from
    epoch_locked=true + finite n_returns, even when model_assumption='cr3bp' and no
    validity_window is present -- the carve-out is scoped to quasi_cycler ONLY."""
    for cls in ("precursor_mga", "mga_tour"):
        row = _base_quasi_cycler_row(
            orbit_class=cls, model_assumption="cr3bp", validity_window=None
        )
        assert not _is_epoch_free_cr3bp_corridor_subclass(row), (
            f"{cls} must never be treated as the epoch-free CR3BP quasi_cycler subclass"
        )


def test_live_epoch_free_quasi_cycler_rows_are_cr3bp_and_windowless() -> None:
    """Sanity check on the live catalogue: any row currently using the v5.2 carve-out
    (epoch_locked=false under orbit_class=quasi_cycler) really is cr3bp + windowless, not
    an accidental mistagging of a real-ephemeris row."""
    rows = _load_catalogue()
    bad = [
        r["id"]
        for r in rows
        if r.get("orbit_class") == "quasi_cycler"
        and r.get("epoch_locked") is False
        and not _is_epoch_free_cr3bp_corridor_subclass(r)
    ]
    assert not bad, (
        f"quasi_cycler rows with epoch_locked=false must satisfy the v5.2 carve-out "
        f"condition (model_assumption='cr3bp' AND validity_window null): {bad}"
    )


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
