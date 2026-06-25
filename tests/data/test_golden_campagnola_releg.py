"""#449 Task 3: sourced Campagnola-Russell endgame-releg golden.

The golden EXPECTED side traces ONLY to Campagnola & Russell 2010 JGCD 33(2)
printed tables / worked scalars (Part-1 Table 1 no-GA, Table 2 with-GA, the
Europa-endgame 154/147 m/s), transcribed in
``docs/notes/2026-06-05-endgame-tisserand-mining.md``. Never a value our own
releg solver computed (feedback_golden_tests_sourced_only). This test locks the
sourced values + their provenance in place; Tasks 2/5 consume them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

_GOLDEN = Path(__file__).resolve().parents[2] / "data" / "golden" / "campagnola_endgame_releg.yaml"


def _load() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(_GOLDEN.read_text())
    return data


def test_golden_exists_and_has_provenance() -> None:
    data = _load()
    assert "Campagnola" in data["source"] and "Russell" in data["source"]
    assert data["digest"].endswith("2026-06-05-endgame-tisserand-mining.md")


def test_table1_no_ga_has_at_least_two_pairs_with_sources() -> None:
    """Table 1 (no-GA) ΔV_min targets for >= 2 moon pairs, each sourced."""
    data = _load()
    rows = data["table1_no_ga"]
    assert len(rows) >= 2
    # The Ganymede-Europa 1.71 km/s anchor (the releg golden's headline floor).
    by_transfer = {r["transfer"]: r for r in rows}
    assert by_transfer["Ganymede-Europa"]["dv_min_kms"] == 1.71
    for r in rows:
        assert r["dv_min_kms"] < r["dv_max_kms"]  # ΔV_min is the floor
        assert "source" in r and "Table 1" in r["source"]


def test_table2_with_ga_present_and_cheaper_than_no_ga() -> None:
    """Table 2 (with-GA) ΔV_min < the same-ends Table 1 (no-GA) ΔV_min."""
    data = _load()
    no_ga = {r["transfer"]: r["dv_min_kms"] for r in data["table1_no_ga"]}
    rows = data["table2_with_ga"]
    assert len(rows) >= 1
    for r in rows:
        assert "source" in r and "Table 2" in r["source"]
        assert r["via"]  # routed through an intermediate moon
    # Callisto-G-Europa (1.61) < Callisto-Europa no-GA (1.94): a GA reduces ΔV.
    cg_europa = next(r for r in rows if r["transfer"] == "Callisto-G-Europa")
    assert cg_europa["dv_min_kms"] < no_ga["Callisto-Europa"]


def test_europa_endgame_scalars_sourced() -> None:
    """Europa endgame 154 m/s (3-VILM) / 147 m/s (CR3BP) + 46-day phasing."""
    data = _load()
    e = data["europa_endgame"]
    assert e["dv_discrete_3vilm_ms"] == 154.0
    assert e["dv_cr3bp_long_transfer_ms"] == 147.0
    assert e["duration_days"] == 46.0
    assert e["vinf_high_kms"] == 1.8 and e["vinf_low_kms"] == 0.77
    assert "A6" in e["source"]


def test_disjoint_contour_pairs_marked_unbridgeable() -> None:
    """The structural-emptiness golden: Uranian pairs are bridgeable=False.

    Each disjoint-contour pair traces to the structural negative
    uranus-neptune-regular-moon-endgame-vilm-2026-06-23. Ariel-Umbriel (the
    tightest adjacent Uranian pair) MUST appear and be unbridgeable.
    """
    data = _load()
    pairs = data["disjoint_contour_pairs"]
    assert len(pairs) >= 1
    by_pair = {(p["moon_a"], p["moon_b"]): p for p in pairs}
    assert ("Ariel", "Umbriel") in by_pair
    for p in pairs:
        assert p["bridgeable"] is False
        assert "uranus-neptune-regular-moon-endgame-vilm" in p["source"]
