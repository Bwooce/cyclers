"""Task 7: combined validate_catalogue() CI gate (closes pending #73's
in-repo semantic half).

validate_catalogue() runs the schema-shape, physical-consistency, and
provenance-tag layers together and is the single in-Python gate the suite
ratchets the live catalogue against. The two-layer loader validation
(JSON Schema 4.3 + validate_schema_invariants) already shipped; this binds
the new physical + provenance layers into the same entry point.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import (
    validate_catalogue,
    validate_provenance_tags,
)

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"


def _live_rows() -> list[dict[str, Any]]:
    return cast("list[dict[str, Any]]", yaml.safe_load(CATALOGUE_PATH.read_text()))


# ---------------------------------------------------------------------------
# Live ratchet: the real catalogue passes the combined gate
# ---------------------------------------------------------------------------


def test_live_catalogue_passes_combined_gate() -> None:
    errs = validate_catalogue(_live_rows())
    assert errs == [], "validate_catalogue violations:\n" + "\n".join(errs)


# ---------------------------------------------------------------------------
# Combined gate surfaces both shape AND physics violations
# ---------------------------------------------------------------------------


def test_combined_gate_catches_shape_violation() -> None:
    """A multi-arc row with top-level a_au (shape layer) is caught."""
    bad = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"a_au": 1.6}}
    errs = validate_catalogue([bad])
    assert any("a_au" in m for m in errs), errs


def test_combined_gate_catches_physics_violation() -> None:
    """A V∞ unit error (physics layer) is caught by the same entry point."""
    bad = {"id": "x", "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": 9700.0}]}
    errs = validate_catalogue([bad])
    assert any("unit error" in m for m in errs), errs


# ---------------------------------------------------------------------------
# Provenance-tag layer (forward-compatible; no-op on current catalogue)
# ---------------------------------------------------------------------------


def test_provenance_tags_noop_on_untagged_catalogue() -> None:
    """Current catalogue carries no per-field provenance tags -> no errors."""
    assert validate_provenance_tags(_live_rows()) == []


def test_bad_source_key_flagged() -> None:
    bad = {"id": "x", "orbit_source": "not-a-real-source"}
    errs = validate_provenance_tags([bad])
    assert any("registry key" in m for m in errs), errs


def test_bad_fidelity_flagged() -> None:
    bad = {"id": "x", "orbit_fidelity": "de440"}  # must be 'real-de440'
    errs = validate_provenance_tags([bad])
    assert any("fidelity tier" in m for m in errs), errs


def test_overclaimed_cross_validated_tier_flagged() -> None:
    """A row claiming cross_validated while sharing one source is rejected."""
    bad = {
        "id": "x",
        "orbit_source": "russell-2004-t34",
        "vinf_source": "russell-2004-t34",  # SAME source -> only consistency
        "orbit_fidelity": "circular-coplanar",
        "vinf_fidelity": "circular-coplanar",
        "validation_tier": "cross_validated",
    }
    errs = validate_provenance_tags([bad])
    assert any("only support" in m for m in errs), errs


def test_correctly_claimed_cross_validated_tier_clean() -> None:
    """Two different sources, same fidelity, declaring cross_validated: clean."""
    ok = {
        "id": "x",
        "orbit_source": "rogers-2012-t1",
        "vinf_source": "russell-2004-t34",
        "orbit_fidelity": "circular-coplanar",
        "vinf_fidelity": "circular-coplanar",
        "validation_tier": "cross_validated",
    }
    assert validate_provenance_tags([ok]) == []


def test_cross_fidelity_overclaim_flagged() -> None:
    """Two sources but different fidelity cannot be cross_validated (S1L1 class)."""
    bad = {
        "id": "x",
        "orbit_source": "rogers-2012-t1",
        "vinf_source": "russell-2004-t34",
        "orbit_fidelity": "circular-coplanar",
        "vinf_fidelity": "real-de440",
        "validation_tier": "cross_validated",
    }
    errs = validate_provenance_tags([bad])
    assert any("only support" in m for m in errs), errs


def test_unknown_tier_value_flagged() -> None:
    bad = {"id": "x", "validation_tier": "platinum"}
    errs = validate_provenance_tags([bad])
    assert any("not a recognised Tier" in m for m in errs), errs
