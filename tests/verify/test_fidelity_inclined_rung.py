"""M-3D Phase 4: solve_at_fidelity resolves circular-inclined (plan §4).

HONEST SCOPE: no sourced inclined-closure anchor exists. Asserted: the rung
resolves (no FidelityRungUnavailableError) and its scalar ToF/V_inf equal the
coplanar rung's (inc-only lift is frame geometry; consistency invariant).
"""

from __future__ import annotations

import pytest

from cyclerfinder.search.sequence import Cell
from cyclerfinder.verify.fidelity import solve_at_fidelity

# Concrete Earth-Mars cycler (a, e) used ONLY as solver inputs to exercise the
# inclined-rung consistency invariant below (inclined ToF/V_inf must equal the
# coplanar rung — inc-only lift is pure frame geometry, so the invariant holds for any
# plausible a, e). These are the SAIC Aldrin (a, e) that #368 retired from the
# catalogue as non-authoritative; the test no longer reads them from the row and they
# carry NO sourced claim here — they are just a concrete cell to resolve.
_FIXTURE_A_AU = 1.60
_FIXTURE_E = 0.393


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
    a_au, e = _FIXTURE_A_AU, _FIXTURE_E
    coplanar = solve_at_fidelity(cell, "circular-coplanar", a_au=a_au, e=e)
    inclined = solve_at_fidelity(cell, "circular-inclined", a_au=a_au, e=e)
    assert inclined.fidelity == "circular-inclined"
    assert inclined.converged is True
    assert inclined.outbound_tof_days == pytest.approx(coplanar.outbound_tof_days)
    assert inclined.vinf_kms == pytest.approx(coplanar.vinf_kms)
