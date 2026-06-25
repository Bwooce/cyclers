"""#465 Task 4: sourced in-band multi-rev-leveraging golden.

The multi-rev leveraging releg (#465) closes a flyby moon-tour cycle by walking
each leg's arrival V_inf down to a common flyby target with a CHAIN of resonant
hops — paying ONLY the leveraging (begingame + endgame) quadratures, NOT the
escape/capture insertions an orbiter mission pays (a cycler flies by every moon,
design §1.3). This golden locks the sourced DECOMPOSITION:

    full Table-1 ΔV_min  ==  leverage-only  +  escape  +  capture

so the cycler's leverage-only per-transfer cost (the design-draft §6 numbers) is
the published Campagnola-Russell Table-1 ΔV_min MINUS the sourced escape/capture
insertions. The EXPECTED side is the published Part-1 Table 1 ΔV_min
(``data/golden/campagnola_endgame_releg.yaml``, transcribed from the paper),
NEVER a number our own code computed in isolation
(``feedback_golden_tests_sourced_only``) — the test asserts that our
``vilm``-derived leverage-only + escape/capture RECONSTRUCTS the published total.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.vilm import (
    _escape_dv,
    _hohmann_vinf,
    _leverage_dv_kms,
    min_vinf_for_vilm,
)

_GOLDEN = Path(__file__).resolve().parents[2] / "data" / "golden" / "campagnola_endgame_releg.yaml"


def _load() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(_GOLDEN.read_text())
    return data


def _order(moon_a: str, moon_b: str) -> tuple[str, str]:
    """Outer, inner by about-primary SMA."""
    if SATELLITES[moon_a].sma_km >= SATELLITES[moon_b].sma_km:
        return moon_a, moon_b
    return moon_b, moon_a


def _leverage_only_kms(moon_a: str, moon_b: str) -> float:
    """The flyby-cycler leverage-only ΔV (begingame at outer + endgame at inner).

    The ``_leverage_dv_kms`` halves of ``vilm.vilm_dv_min`` — the cycler pays
    these but NOT the escape/capture (it flies by, design §1.3).
    """
    outer, inner = _order(moon_a, moon_b)
    vinf_o, vinf_i = _hohmann_vinf(outer, inner)
    return _leverage_dv_kms(outer, vinf_o, exterior=True) + _leverage_dv_kms(
        inner, vinf_i, exterior=False
    )


def _escape_capture_kms(moon_a: str, moon_b: str) -> float:
    """The escape + capture insertion cost an ORBITER pays (a cycler does not)."""
    outer, inner = _order(moon_a, moon_b)
    return _escape_dv(outer, min_vinf_for_vilm(outer)) + _escape_dv(
        inner, min_vinf_for_vilm(inner, exterior=False)
    )


# Design-draft §6 leverage-only per-transfer costs (km/s). These are NOT the
# sourced EXPECTED side — they are what the test PROVES reconstructs the published
# Table-1 total; the sourced anchor is the YAML ``dv_min_kms`` below.
_DESIGN_LEVERAGE_ONLY_KMS = {
    "Ganymede-Europa": 0.305,
    "Europa-Io": 0.409,
    "Rhea-Dione": 0.196,
}


def test_golden_decomposition_reconstructs_published_table1() -> None:
    """leverage-only + escape + capture == the published Table-1 ΔV_min (sourced).

    For every Table-1 pair in the golden YAML, the ``vilm``-derived leverage-only
    cost plus the sourced escape/capture insertions must reconstruct the PUBLISHED
    ΔV_min within the documented transcription tolerance (the YAML carries the
    paper's 2-3 sig-fig rounded values). This is the non-circular anchor: the
    EXPECTED side is the printed table, and the cycler's leverage-only number is
    the published total minus the (also-decomposed) escape/capture.
    """
    data = _load()
    rows = {r["transfer"]: r for r in data["table1_no_ga"]}
    # The pairs whose moons are all in the registry (the golden has a couple of
    # Saturnian inner moons we exercise for the Saturnian positive control).
    checked = 0
    for transfer in ("Ganymede-Europa", "Europa-Io", "Rhea-Dione", "Callisto-Ganymede"):
        if transfer not in rows:
            continue
        row = rows[transfer]
        a, b = row["moon_a"], row["moon_b"]
        full = _leverage_only_kms(a, b) + _escape_capture_kms(a, b)
        # Reconstructs the published ΔV_min to the paper's printed precision
        # (2 sig figs ⇒ ~0.01 km/s rounding band).
        assert abs(full - row["dv_min_kms"]) <= 0.02, (
            f"{transfer}: reconstructed {full:.3f} vs published {row['dv_min_kms']}"
        )
        checked += 1
    assert checked >= 3


def test_leverage_only_costs_match_design_section6() -> None:
    """The cycler leverage-only costs reproduce the design-draft §6 table.

    The §6 numbers are the begingame+endgame quadratures (the cycler's actual
    per-transfer cost, no escape/capture); each must match to the documented tol.
    """
    for transfer, expected in _DESIGN_LEVERAGE_ONLY_KMS.items():
        a, b = transfer.split("-")
        assert _leverage_only_kms(a, b) == _approx(expected)


def test_galilean_cycle_leverage_only_in_band() -> None:
    """The Io-Europa-Ganymede-Io leverage-only cycle cost is IN-BAND (sourced).

    The cycle sums Europa-Io + Ganymede-Europa + Ganymede-Io leverage-only costs;
    the total is comfortably under the 3.5 km/s/cycle powered ceiling — the in-band
    closure the single-leg relegs (13.18 km/s) missed. Each per-transfer term is
    the sourced §6 quadrature.
    """
    cycle = (
        _leverage_only_kms("Io", "Europa")
        + _leverage_only_kms("Europa", "Ganymede")
        + _leverage_only_kms("Ganymede", "Io")
    )
    assert cycle < 3.5
    # Sanity: the headline ~1.5 km/s/cycle (design §6), not a degenerate near-zero.
    assert 1.0 < cycle < 2.0


def _approx(value: float, tol: float = 1.0e-3) -> Any:
    import pytest

    return pytest.approx(value, abs=tol)
