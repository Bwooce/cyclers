"""Data-integrity (V0) gate for the two Aldrin V-infinity-leveraging
establishment entries, sourced from Rogers et al. 2012 (AIAA 2012-4746)
Table 4 (analytic-ephemeris STOUR results, p.7).

These are establishment (leveraging) trajectories, not closed ballistic
cyclers, so this is a *traceability / drift guard*, not a reproduction:
it pins the Rogers Table-4 numbers the catalogue carries and asserts the
Mars-side V-infinity stays null (Table 4 tabulates only Earth launch /
flyby V-infinity for these Earth-to-Earth leveraging arcs — a distinct
Mars V-infinity is not defined for them). If any pinned value changes,
this test fails so the change is re-verified against the PDF in
docs/refs/rogers-2012-vinf-leveraging-cyclers-AIAA-2012-4746.pdf.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

_CATALOGUE = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

# Rogers et al. 2012, AIAA 2012-4746, Table 4, p.7 (corroborated in the body
# text: the 4:3(2)- Aldrin cycler has "a maneuver delta-V of 0.568 km/s" and a
# launch V_inf of 2.608 km/s with a ~2.90 km/s leveraging increase -> flyby
# V_inf ~= 5.51 km/s).
_EXPECTED: dict[str, dict[str, float]] = {
    "aldrin-4-3-2-establishment": {
        "vinf_e": 5.509,
        "leveraging_dv": 0.568,
        "perihelion_au": 0.983,
        "aphelion_au": 2.229,
    },
    "aldrin-3-2-1-establishment": {
        "vinf_e": 6.546,
        "leveraging_dv": 0.530,
        "perihelion_au": 0.964,
        "aphelion_au": 2.229,
    },
}


def _entries() -> dict[str, dict[str, Any]]:
    return {r["id"]: r for r in yaml.safe_load(_CATALOGUE.read_text())}


@pytest.mark.parametrize("cid", sorted(_EXPECTED))
def test_establishment_entry_matches_rogers_table4(cid: str) -> None:
    exp = _EXPECTED[cid]
    entry = _entries()[cid]

    vinf = {v["body"]: v["vinf_kms"] for v in entry["vinf_kms_at_encounters"]}
    assert vinf["E"] == pytest.approx(exp["vinf_e"], abs=0.001), f"{cid} Earth V_inf"
    # Mars V_inf is not defined for an Earth-side leveraging establishment arc.
    assert vinf["M"] is None, f"{cid} Mars V_inf should stay null (N/A)"

    assert entry["v_infinity_leveraging_dv_kms"] == pytest.approx(
        exp["leveraging_dv"], abs=0.001
    ), f"{cid} leveraging dV"

    oe = entry["orbit_elements"]
    assert oe["perihelion_au"] == pytest.approx(exp["perihelion_au"], abs=0.001)
    assert oe["aphelion_au"] == pytest.approx(exp["aphelion_au"], abs=0.001)
