"""Fast unit tests for the Axis-B persistence classifier (Forge phase 1).

These are synthetic-input tests of :func:`fidelity_persistence` — no physics,
no ephemeris, no catalogue. They prove the classifier's *teeth*: an undocumented
shift is flagged ``SHIFTS_UNDOCUMENTED``, while a tolerated value PERSISTS and a
source-documented shift is recognised. Provenance: all inputs are synthetic
contract fixtures (``# COMPUTED``-equivalent — they assert the classifier's
logic, not any sourced physical magnitude).
"""

from __future__ import annotations

import pytest

from cyclerfinder.search.sequence import Cell
from cyclerfinder.verify.fidelity import (
    FidelityRungUnavailableError,
    PersistenceClass,
    fidelity_persistence,
    solve_at_fidelity,
)


def _cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )


def test_persists_within_tolerance() -> None:
    rep = fidelity_persistence("vinf", 6.50, 6.55, abs_tol=0.1)
    assert rep.classification is PersistenceClass.PERSISTS
    assert rep.within_tol is True
    assert rep.delta == pytest.approx(0.05)


def test_shift_in_expected_direction_is_documented() -> None:
    # Documented to INCREASE with fidelity; it does, beyond tol.
    rep = fidelity_persistence("tof", 146.0, 165.0, abs_tol=2.0, expected_direction=+1)
    assert rep.classification is PersistenceClass.SHIFTS_DOCUMENTED
    assert rep.within_tol is False
    assert rep.delta == pytest.approx(19.0)


def test_shift_against_expected_direction_is_undocumented() -> None:
    # Documented to INCREASE, but it DECREASED — the red flag (teeth).
    rep = fidelity_persistence("tof", 146.0, 120.0, abs_tol=2.0, expected_direction=+1)
    assert rep.classification is PersistenceClass.SHIFTS_UNDOCUMENTED


def test_shift_into_band_is_documented() -> None:
    rep = fidelity_persistence("tof", 146.0, 165.0, abs_tol=2.0, documented_band=(161.0, 172.0))
    assert rep.classification is PersistenceClass.SHIFTS_DOCUMENTED


def test_shift_toward_band_but_not_inside_is_documented() -> None:
    # 151 is below the band [161,172] but closer to it than 146 was.
    rep = fidelity_persistence("tof", 146.0, 151.0, abs_tol=2.0, documented_band=(161.0, 172.0))
    assert rep.classification is PersistenceClass.SHIFTS_DOCUMENTED


def test_shift_away_from_band_is_undocumented() -> None:
    # 130 is further below the band than 146 was — moves away, undocumented.
    rep = fidelity_persistence("tof", 146.0, 130.0, abs_tol=2.0, documented_band=(161.0, 172.0))
    assert rep.classification is PersistenceClass.SHIFTS_UNDOCUMENTED


def test_no_documentation_means_any_shift_is_undocumented() -> None:
    rep = fidelity_persistence("tof", 146.0, 200.0, abs_tol=2.0)
    assert rep.classification is PersistenceClass.SHIFTS_UNDOCUMENTED


def test_overshooting_past_band_far_side_is_undocumented() -> None:
    # Band [161,172]; high=300 is far past it AND further than low=146 → away.
    rep = fidelity_persistence("tof", 146.0, 300.0, abs_tol=2.0, documented_band=(161.0, 172.0))
    assert rep.classification is PersistenceClass.SHIFTS_UNDOCUMENTED


# ---------------------------------------------------------------------------
# solve_at_fidelity dispatch contract (fast paths only)
# ---------------------------------------------------------------------------


def test_unknown_fidelity_raises() -> None:
    with pytest.raises(ValueError, match="unknown fidelity"):
        solve_at_fidelity(_cell(), "nonsense")  # type: ignore[arg-type]


def test_analytic_ephemeris_rung_is_unavailable() -> None:
    with pytest.raises(FidelityRungUnavailableError):
        solve_at_fidelity(_cell(), "analytic-ephemeris")


def test_coplanar_requires_a_and_e() -> None:
    with pytest.raises(ValueError, match="requires sourced"):
        solve_at_fidelity(_cell(), "circular-coplanar")


def test_coplanar_solves_from_sourced_elements() -> None:
    # Closed-form construction; fast. a/e are SOURCED Aldrin elements (Rogers
    # 2012 Table 1, mirrored from data/catalogue.yaml). We assert only the
    # solution SHAPE and convergence flag here — no golden magnitude (the
    # sourced-value gate lives in test_fidelity_gate.py).
    sol = solve_at_fidelity(_cell(), "circular-coplanar", a_au=1.60, e=0.393)
    assert sol.fidelity == "circular-coplanar"
    assert sol.converged is True
    assert sol.outbound_tof_days > 0.0
    assert set(sol.vinf_kms) == {"E", "M"}
