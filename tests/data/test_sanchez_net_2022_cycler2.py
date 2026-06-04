"""Data-integrity (V0) gate for the Sanchez Net 2022 EM Cycler 2 catalogue row.

Asserts that the sourced V-infinity multiset in the catalogue row matches the
values read directly from Sanchez Net 2022, Fig. 2b (p. 862):

  Earth flybys (events 3, 4, 6): 3.090, 5.285, 5.721 km/s
  Mars flybys (events 2, 5): 6.466, 6.871 km/s

Source: Sanchez Net, Pellegrini, Parker, Vander Hook, Woollands (2022).
"Cycler Orbits and Solar System Pony Express." Journal of Spacecraft and
Rockets, Vol. 59, No. 3, pp. 861-870. DOI: 10.2514/1.A35091. Fig. 2b.

This is a traceability test: it guards against the catalogue row drifting
from its sourced numbers. The expected values below are the ONLY legitimate
source for this test — they were read directly from Fig. 2b and must not be
replaced with values computed by this codebase (golden-test discipline).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

_REPO = Path(__file__).resolve().parent.parent.parent
_CATALOGUE = _REPO / "data" / "catalogue.yaml"

_ROW_ID = "sanchez-net-2022-em-cycler2"

# ---------------------------------------------------------------------------
# SOURCED EXPECTED values — Sanchez Net 2022 Fig. 2b, p. 862.
# Earth-flyby V∞ values (events 3, 4, 6 in Fig. 2b):
_EXPECTED_EARTH_VINF_KMS: list[float] = [3.090, 5.285, 5.721]
# Mars-flyby V∞ values (events 2, 5 in Fig. 2b):
_EXPECTED_MARS_VINF_KMS: list[float] = [6.466, 6.871]
# ---------------------------------------------------------------------------


def _load_row() -> dict[str, Any]:
    rows: list[dict[str, Any]] = yaml.safe_load(_CATALOGUE.read_text())
    matches = [r for r in rows if r["id"] == _ROW_ID]
    assert matches, f"catalogue row {_ROW_ID!r} not found"
    return matches[0]


def test_row_exists_and_is_multi_arc() -> None:
    """The row is present and tagged multi-arc (EM patched-conic)."""
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


def test_earth_vinf_multiset_matches_fig2b() -> None:
    """The catalogue row's Earth V∞ values match Fig. 2b events 3, 4, 6.

    SOURCED EXPECTED: Sanchez Net 2022, JSR 59(3):861-870, Fig. 2b, p. 862.
    Events 3 (3.090), 4 (5.285), 6 (5.721) km/s.
    """
    row = _load_row()
    earth_vinfs = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"] if v["body"] == "E")
    expected_sorted = sorted(_EXPECTED_EARTH_VINF_KMS)
    assert len(earth_vinfs) == len(expected_sorted), (
        f"expected {len(expected_sorted)} Earth V∞ entries, got {len(earth_vinfs)}"
    )
    for got, exp in zip(earth_vinfs, expected_sorted, strict=True):
        assert got == pytest.approx(exp, abs=0.001), (
            f"Earth V∞ mismatch: catalogue has {got} km/s, Fig. 2b has {exp} km/s"
        )


def test_mars_vinf_multiset_matches_fig2b() -> None:
    """The catalogue row's Mars V∞ values match Fig. 2b events 2, 5.

    SOURCED EXPECTED: Sanchez Net 2022, JSR 59(3):861-870, Fig. 2b, p. 862.
    Events 2 (6.466), 5 (6.871) km/s.
    """
    row = _load_row()
    mars_vinfs = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"] if v["body"] == "M")
    expected_sorted = sorted(_EXPECTED_MARS_VINF_KMS)
    assert len(mars_vinfs) == len(expected_sorted), (
        f"expected {len(expected_sorted)} Mars V∞ entries, got {len(mars_vinfs)}"
    )
    for got, exp in zip(mars_vinfs, expected_sorted, strict=True):
        assert got == pytest.approx(exp, abs=0.001), (
            f"Mars V∞ mismatch: catalogue has {got} km/s, Fig. 2b has {exp} km/s"
        )


def test_near_ballistic_delta_v() -> None:
    """Max per-flyby Δv ≤ 0.007 km/s (7 m/s) per Sanchez Net Fig. 2b event 3."""
    row = _load_row()
    assert row.get("delta_v_kms") == pytest.approx(0.007, abs=0.0001), (
        "delta_v_kms should be 0.007 km/s (7 m/s) — sourced from Fig. 2b event 3"
    )


def test_period_k4_years_7_87() -> None:
    """Period k=4 (nearest synodic integer), years≈7.87 — derived from Fig. 2b dates.

    Caption source: Sanchez Net 2022, Fig. 2b: '8 years between Mars 1 and 2'.
    Mars-1 (12/19/2037) to Mars-2 (11/02/2045) = 2875 days = 7.87 yr.
    """
    row = _load_row()
    period = row.get("period") or {}
    assert period.get("k") == 4, f"expected period.k=4, got {period.get('k')!r}"
    assert period.get("years") == pytest.approx(7.87, abs=0.05), (
        f"expected period.years ≈ 7.87, got {period.get('years')!r}"
    )


def test_sequence_canonical_eeeem() -> None:
    """sequence_canonical = 'E-E-E-M' (lex-min cyclic rotation of 4-node E-M-E-E cycle).

    Sanchez Net labels this 'EM' at macro scale (one Mars visit per cycle).
    The actual encounter loop has two intermediate Earth flybys per cycle:
    E(launch) -> M -> E -> E -> M -> E -> E, giving per-cycle node list
    [E, M, E, E]. Lex-min cyclic rotation of {E-M-E-E, M-E-E-E, E-E-E-M,
    E-E-M-E} is 'E-E-E-M'.
    """
    row = _load_row()
    assert row.get("sequence_canonical") == "E-E-E-M", (
        f"expected sequence_canonical='E-E-E-M', got {row.get('sequence_canonical')!r}"
    )
