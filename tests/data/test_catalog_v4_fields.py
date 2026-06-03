"""Task 1.1 + 1.2: v4 field parsing and fully_defined dispatch tests.

Spec §16.7 (schema v4, 2026-06-03): cycler_class, orbit_elements frame/center,
invariants, cr3bp, period_basis — all additive, optional, default to v3 behaviour.
"""

from __future__ import annotations

import textwrap
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CatalogueEntry, _entry_from_yaml, load_catalog

# ---------------------------------------------------------------------------
# Task 1.1 — v4 fields parse correctly and default properly
# ---------------------------------------------------------------------------


def test_cycler_class_defaults_to_single_ellipse() -> None:
    """A row with no cycler_class key reads as single-ellipse."""
    cat = load_catalog()
    entries = cat.by_id
    # s1l1-2syn-em-cpom has no cycler_class in current catalogue
    e = entries["s1l1-2syn-em-cpom"]
    assert e.cycler_class == "single-ellipse"


def test_v4_fields_present_and_defaulted() -> None:
    """All v4 fields are present on entries and default correctly."""
    cat = load_catalog()
    e = cat.by_id["s1l1-2syn-em-cpom"]
    assert e.orbit_elements_reference_frame == "heliocentric-inertial"
    assert e.orbit_elements_center == "Sun"
    assert e.invariants is None  # only set for multi-arc
    assert e.cr3bp is None  # only set for non-keplerian
    assert e.period_basis is None  # only set for n-body


def _make_entry_from_yaml_str(yaml_str: str) -> CatalogueEntry:
    """Parse a minimal YAML row string into a CatalogueEntry."""
    row: dict[str, Any] = yaml.safe_load(textwrap.dedent(yaml_str))
    return _entry_from_yaml(row)


def test_cycler_class_explicit_value_parsed() -> None:
    """An explicit cycler_class is read correctly."""
    e = _make_entry_from_yaml_str(
        """
        id: test-multi-arc
        bodies: ["E", "M"]
        sequence_canonical: "E-M"
        cycler_class: multi-arc
        """
    )
    assert e.cycler_class == "multi-arc"


def test_orbit_elements_frame_parsed() -> None:
    """orbit_elements.reference_frame is read from the nested key."""
    e = _make_entry_from_yaml_str(
        """
        id: test-nk
        bodies: ["E", "Moon"]
        sequence_canonical: "E-Moon"
        orbit_elements:
          reference_frame: planetcentric-inertial
          center: Earth
        """
    )
    assert e.orbit_elements_reference_frame == "planetcentric-inertial"
    assert e.orbit_elements_center == "Earth"


def test_invariants_parsed_when_present() -> None:
    """invariants dict is loaded when present."""
    e = _make_entry_from_yaml_str(
        """
        id: test-multi-arc-inv
        bodies: ["E", "M"]
        sequence_canonical: "E-M"
        cycler_class: multi-arc
        invariants:
          aphelion_ratio: 1.44
        """
    )
    assert e.invariants == {"aphelion_ratio": 1.44}


def test_cr3bp_parsed_from_orbit_elements() -> None:
    """cr3bp is read from orbit_elements.cr3bp when present."""
    e = _make_entry_from_yaml_str(
        """
        id: test-nk-cr3bp
        bodies: ["E", "Moon"]
        sequence_canonical: "E-Moon"
        cycler_class: non-keplerian
        orbit_elements:
          cr3bp:
            jacobi_constant: 3.0
            period_nd: 6.2
            stability_index: 1.0
        """
    )
    assert e.cr3bp is not None
    assert e.cr3bp["jacobi_constant"] == pytest.approx(3.0)
    assert e.cr3bp["period_nd"] == pytest.approx(6.2)
    assert e.cr3bp["stability_index"] == pytest.approx(1.0)


def test_period_basis_parsed_when_present() -> None:
    """period.basis list becomes a tuple on the entry."""
    e = _make_entry_from_yaml_str(
        """
        id: test-vem
        bodies: ["V", "E", "M"]
        sequence_canonical: "E-M-V"
        period:
          pair: "E-M"
          k: 3
          basis:
            - pair: "E-M"
              k: 3
            - pair: "E-V"
              k: 4
        """
    )
    assert e.period_basis is not None
    assert len(e.period_basis) == 2
    assert e.period_basis[0] == {"pair": "E-M", "k": 3}


def test_period_basis_none_when_absent() -> None:
    """period.basis defaults to None when not present in the row."""
    e = _make_entry_from_yaml_str(
        """
        id: test-simple
        bodies: ["E", "M"]
        sequence_canonical: "E-M"
        period:
          pair: "E-M"
          k: 2
        """
    )
    assert e.period_basis is None


# ---------------------------------------------------------------------------
# Task 1.2 — fully_defined dispatches by cycler_class
# ---------------------------------------------------------------------------


def _minimal_single_ellipse_row(
    *,
    a_au: float | None = 1.5,
    e: float | None = 0.3,
    data_gaps: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Minimal dict that can pass fully_defined for single-ellipse."""
    row: dict[str, Any] = {
        "id": "test-se",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-M",
        "cycler_class": "single-ellipse",
        "orbit_elements": {"a_au": a_au, "e": e},
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 5.0},
            {"body": "M", "vinf_kms": 3.0},
        ],
        "legs": [{"tof_days": 200, "n_revs": 0}],
    }
    if data_gaps is not None:
        row["data_gaps"] = data_gaps
    return row


def test_fully_defined_single_ellipse_with_ae() -> None:
    """single-ellipse with a/e + vinf + legs is fully_defined."""
    e = _entry_from_yaml(_minimal_single_ellipse_row())
    assert e.fully_defined is True


def test_fully_defined_single_ellipse_missing_a_is_false() -> None:
    """single-ellipse missing a_au is NOT fully_defined."""
    e = _entry_from_yaml(_minimal_single_ellipse_row(a_au=None))
    assert e.fully_defined is False


def test_fully_defined_single_ellipse_missing_e_is_false() -> None:
    """single-ellipse missing e is NOT fully_defined."""
    e = _entry_from_yaml(_minimal_single_ellipse_row(e=None))
    assert e.fully_defined is False


def test_fully_defined_multi_arc_uses_invariants_not_ae() -> None:
    """multi-arc with null a/e but populated invariants + vinf + legs is fully_defined."""
    row: dict[str, Any] = {
        "id": "test-ma",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-M",
        "cycler_class": "multi-arc",
        "orbit_elements": {"a_au": None, "e": None},
        "invariants": {"aphelion_ratio": 1.44},
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 5.0},
            {"body": "M", "vinf_kms": 3.0},
        ],
        "legs": [{"tof_days": 200, "n_revs": 0}],
    }
    e = _entry_from_yaml(row)
    assert e.fully_defined is True


def test_fully_defined_multi_arc_without_invariants_is_false() -> None:
    """multi-arc with no invariants (None) is NOT fully_defined."""
    row: dict[str, Any] = {
        "id": "test-ma-no-inv",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-M",
        "cycler_class": "multi-arc",
        "orbit_elements": {"a_au": None, "e": None},
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 5.0},
            {"body": "M", "vinf_kms": 3.0},
        ],
        "legs": [{"tof_days": 200, "n_revs": 0}],
    }
    e = _entry_from_yaml(row)
    assert e.fully_defined is False


def test_fully_defined_non_keplerian_uses_cr3bp() -> None:
    """non-keplerian with complete cr3bp identity tuple is fully_defined."""
    row: dict[str, Any] = {
        "id": "test-nk-full",
        "bodies": ["E", "Moon"],
        "sequence_canonical": "E-Moon",
        "cycler_class": "non-keplerian",
        "orbit_elements": {
            "a_au": None,
            "e": None,
            "cr3bp": {
                "jacobi_constant": 3.0,
                "period_nd": 6.2,
                "stability_index": 1.0,
            },
        },
        "primary": "Earth",
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 1.0},
        ],
        "legs": [{"tof_days": 10, "n_revs": 0}],
    }
    e = _entry_from_yaml(row)
    assert e.fully_defined is True


def test_fully_defined_non_keplerian_missing_cr3bp_is_false() -> None:
    """non-keplerian without cr3bp dict is NOT fully_defined."""
    row: dict[str, Any] = {
        "id": "test-nk-no-cr3bp",
        "bodies": ["E", "Moon"],
        "sequence_canonical": "E-Moon",
        "cycler_class": "non-keplerian",
        "orbit_elements": {"a_au": None, "e": None},
        "primary": "Earth",
    }
    e = _entry_from_yaml(row)
    assert e.fully_defined is False


def test_fully_defined_non_keplerian_partial_cr3bp_is_false() -> None:
    """non-keplerian with incomplete cr3bp (missing stability_index) is NOT fully_defined."""
    row: dict[str, Any] = {
        "id": "test-nk-partial",
        "bodies": ["E", "Moon"],
        "sequence_canonical": "E-Moon",
        "cycler_class": "non-keplerian",
        "orbit_elements": {
            "cr3bp": {
                "jacobi_constant": 3.0,
                "period_nd": 6.2,
                # stability_index missing
            }
        },
        "primary": "Earth",
    }
    e = _entry_from_yaml(row)
    assert e.fully_defined is False


def test_fully_defined_with_data_gaps_is_false() -> None:
    """Any row with data_gaps is NOT fully_defined regardless of class."""
    e = _entry_from_yaml(
        _minimal_single_ellipse_row(
            data_gaps=[{"path": "orbit_elements.a_au", "kind": "missing", "note": "test"}]
        )
    )
    assert e.fully_defined is False
