"""M-3D Phase 4: solve_at_fidelity resolves circular-inclined (plan §4).

HONEST SCOPE: no sourced inclined-closure anchor exists. Asserted: the rung
resolves (no FidelityRungUnavailableError) and its scalar ToF/V_inf equal the
coplanar rung's (inc-only lift is frame geometry; consistency invariant).
"""

from __future__ import annotations

import pytest

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.sequence import Cell
from cyclerfinder.verify.fidelity import solve_at_fidelity

_ALDRIN_ID = "aldrin-classic-em-k1-outbound"


def _aldrin_sourced_ae() -> tuple[float, float]:
    """Read the SOURCED Aldrin (a, e) from the catalogue row (no magic numbers)."""
    entry = load_catalog().by_id[_ALDRIN_ID]
    oe = entry.raw["orbit_elements"]
    return float(oe["a_au"]), float(oe["e"])


def _aldrin_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )


def test_inclined_rung_resolves_and_matches_coplanar_scalars() -> None:
    cell = _aldrin_cell()
    a_au, e = _aldrin_sourced_ae()
    coplanar = solve_at_fidelity(cell, "circular-coplanar", a_au=a_au, e=e)
    inclined = solve_at_fidelity(cell, "circular-inclined", a_au=a_au, e=e)
    assert inclined.fidelity == "circular-inclined"
    assert inclined.converged is True
    assert inclined.outbound_tof_days == pytest.approx(coplanar.outbound_tof_days)
    assert inclined.vinf_kms == pytest.approx(coplanar.vinf_kms)
