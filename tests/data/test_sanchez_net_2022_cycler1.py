"""Data-integrity (V0) gate for the Sanchez Net 2022 EEM Cycler 1 catalogue row.

Asserts that the sourced V-infinity multiset in the catalogue row matches the
values read directly from Sanchez Net 2022, Fig. 2a (p. 862):

  Earth flybys (events 2,4,5,7,8): 4.283, 3.605, 3.709, 5.225, 5.234 km/s
  Mars flybys (events 3,6,9): 6.207, 5.210, 7.276 km/s

Source: Sanchez Net, Pellegrini, Parker, Vander Hook, Woollands (2022).
"Cycler Orbits and Solar System Pony Express." Journal of Spacecraft and
Rockets, Vol. 59, No. 3, pp. 861-870. DOI: 10.2514/1.A35091. Fig. 2a.

This is a traceability test: it guards against the catalogue row drifting
from its sourced numbers. The expected values below are the ONLY legitimate
source for this test — they were read directly from Fig. 2a and must not be
replaced with values computed by this codebase (golden-test discipline).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

_REPO = Path(__file__).resolve().parent.parent.parent
_CATALOGUE = _REPO / "data" / "catalogue.yaml"

_ROW_ID = "sanchez-net-2022-eem-cycler1"

# ---------------------------------------------------------------------------
# SOURCED EXPECTED values — Sanchez Net 2022 Fig. 2a, p. 862.
# Earth-flyby V∞ values (events 2, 4, 5, 7, 8 in Fig. 2a):
_EXPECTED_EARTH_VINF_KMS: list[float] = [4.283, 3.605, 3.709, 5.225, 5.234]
# Mars-flyby V∞ values (events 3, 6, 9 in Fig. 2a):
_EXPECTED_MARS_VINF_KMS: list[float] = [6.207, 5.210, 7.276]
# ---------------------------------------------------------------------------


def _load_row() -> dict[str, Any]:
    rows: list[dict[str, Any]] = yaml.safe_load(_CATALOGUE.read_text())
    matches = [r for r in rows if r["id"] == _ROW_ID]
    assert matches, f"catalogue row {_ROW_ID!r} not found"
    return matches[0]


def test_row_exists_and_is_multi_arc() -> None:
    """The row is present and tagged multi-arc (EEM patched-conic)."""
    row = _load_row()
    assert row["cycler_class"] == "multi-arc", (
        f"expected cycler_class=multi-arc, got {row['cycler_class']!r}"
    )


def test_row_model_assumption_is_analytic_ephemeris() -> None:
    """Model assumption must be analytic-ephemeris (real-date patched-conic)."""
    row = _load_row()
    assert row["model_assumption"] == "analytic-ephemeris", (
        f"expected analytic-ephemeris, got {row['model_assumption']!r}"
    )


def test_row_orbit_elements_a_e_are_null() -> None:
    """Multi-arc schema invariant: a_au and e must be null."""
    row = _load_row()
    oe = row.get("orbit_elements") or {}
    assert oe.get("a_au") is None, "a_au must be null for multi-arc entry"
    assert oe.get("e") is None, "e must be null for multi-arc entry"


def test_earth_vinf_multiset_matches_fig2a() -> None:
    """The catalogue row's Earth V∞ values match Fig. 2a events 2,4,5,7,8.

    SOURCED EXPECTED: Sanchez Net 2022, JSR 59(3):861-870, Fig. 2a, p. 862.
    Events 2 (4.283), 4 (3.605), 5 (3.709), 7 (5.225), 8 (5.234) km/s.
    """
    row = _load_row()
    earth_vinfs = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"] if v["body"] == "E")
    expected_sorted = sorted(_EXPECTED_EARTH_VINF_KMS)
    assert len(earth_vinfs) == len(expected_sorted), (
        f"expected {len(expected_sorted)} Earth V∞ entries, got {len(earth_vinfs)}"
    )
    for got, exp in zip(earth_vinfs, expected_sorted, strict=True):
        assert got == pytest.approx(exp, abs=0.001), (
            f"Earth V∞ mismatch: catalogue has {got} km/s, Fig. 2a has {exp} km/s"
        )


def test_mars_vinf_multiset_matches_fig2a() -> None:
    """The catalogue row's Mars V∞ values match Fig. 2a events 3,6,9.

    SOURCED EXPECTED: Sanchez Net 2022, JSR 59(3):861-870, Fig. 2a, p. 862.
    Events 3 (6.207), 6 (5.210), 9 (7.276) km/s.
    """
    row = _load_row()
    mars_vinfs = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"] if v["body"] == "M")
    expected_sorted = sorted(_EXPECTED_MARS_VINF_KMS)
    assert len(mars_vinfs) == len(expected_sorted), (
        f"expected {len(expected_sorted)} Mars V∞ entries, got {len(mars_vinfs)}"
    )
    for got, exp in zip(mars_vinfs, expected_sorted, strict=True):
        assert got == pytest.approx(exp, abs=0.001), (
            f"Mars V∞ mismatch: catalogue has {got} km/s, Fig. 2a has {exp} km/s"
        )


def test_near_ballistic_delta_v() -> None:
    """Max per-flyby Δv ≤ 0.005 km/s (5 m/s) per Sanchez Net Fig. 2a event 4."""
    row = _load_row()
    assert row.get("delta_v_kms") == pytest.approx(0.005, abs=0.0001), (
        "delta_v_kms should be 0.005 km/s (5 m/s) — sourced from Fig. 2a event 4"
    )


def test_period_k2_years_4_28() -> None:
    """Period k=2 (2-synodic), years≈4.28 — sourced from Fig. 2a dates."""
    row = _load_row()
    period = row.get("period") or {}
    assert period.get("k") == 2, f"expected period.k=2, got {period.get('k')!r}"
    assert period.get("years") == pytest.approx(4.28, abs=0.05), (
        f"expected period.years ≈ 4.28, got {period.get('years')!r}"
    )


def test_sequence_canonical_eem() -> None:
    """sequence_canonical = 'E-E-M' (EEM structure per Sanchez Net Fig. 2a caption)."""
    row = _load_row()
    assert row.get("sequence_canonical") == "E-E-M", (
        f"expected sequence_canonical='E-E-M', got {row.get('sequence_canonical')!r}"
    )
