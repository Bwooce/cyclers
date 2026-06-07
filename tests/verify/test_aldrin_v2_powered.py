"""§14 V2-powered evidence — the Aldrin E-M cycler under the class-split amendment.

This module records the RECORDED, MECHANICAL outcome of the 2026-06-07 §14 V2
class-split amendment applied to the classic Aldrin Earth-Mars cycler pair. It is
the positive counterpart to ``test_aldrin_v2_v3_campaign.py`` (the #134 finding
that the *original* single V2 gate is structurally unsatisfiable for a powered,
per-cycle-retargeted cycler — drift ~4.14e8 km / 3 laps, ≈2072x tolerance).

The headline result:

* **``aldrin-classic-em-k1-outbound`` PASSES V2-powered.** Over 3 consecutive
  in-family cycles, (a) every planned encounter is achieved with the per-cycle
  maintenance applied (Mars-flyby V∞ continuity ≤1e-6 km/s, ≤ the 0.5 km/s
  encounter tolerance; strictly-positive in-family maintenance ΔV ≈2.76-2.91
  km/s/cycle) and (b) the intra-cycle drift vs the planned trajectory is bounded
  (per-leg Kepler forward-reprop residual ≤0.002 km, ≤ the 1 km bound). This is
  the recorded ``validation_level: V2`` evidence for the outbound row.

* **``aldrin-classic-em-k1-inbound`` stays V1 (honest exclusion).** Its
  real-window optimiser lands on a *ballistic ΔV≈0 neighbour* rather than the
  in-family powered solve (the recorded off-family resolver issue, #134). A
  ΔV≈0 "maintenance" is not "the documented per-cycle maintenance applied", so
  clause (a) is not met and the inbound is NOT promoted.

Discipline: these tests are ``slow`` (real-DE440 BVP solves, minutes per case).
The maintenance-ΔV magnitude is OUR computed value (McConaghy 2002 defers it) and
is only sanity-bounded, never matched against a sourced number. The gated
quantities are code-path self-consistency checks (V∞ continuity, Kepler-reprop
residual), not rediscovered published values.
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.bvp import solve_powered_periodic_cycler
from cyclerfinder.verify.v2_powered import (
    ENCOUNTER_VINF_TOL_KMS,
    INTRA_CYCLE_DRIFT_TOL_KM,
    verify_aldrin_v2_powered,
)
from tests.data._catalogue_loader_m6b import load_m6b_entries

ALDRIN_PRIORITY = datetime(1985, 10, 28, tzinfo=UTC)
ALDRIN_PERIOD_YEARS = 2.135  # catalogue Aldrin cycler repeat (= 1 E-M synodic)


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


@pytest.mark.slow
def test_aldrin_outbound_passes_v2_powered(astropy_ephem: Ephemeris) -> None:
    """V2-powered PASS: 3 consecutive in-family cycles each achieve their
    encounters with maintenance applied AND keep intra-cycle drift bounded."""
    result = verify_aldrin_v2_powered(
        astropy_ephem,
        priority_date=ALDRIN_PRIORITY,
        period_years=ALDRIN_PERIOD_YEARS,
        n_cycles=3,
    )
    assert result.n_cycles == 3
    assert len(result.per_cycle) == 3
    for i, c in enumerate(result.per_cycle):
        # In-family on the SOURCED Aldrin anchors (a≈1.59 AU, e≈0.393).
        assert c.converged, f"cycle {i} did not converge"
        assert c.a_au == pytest.approx(1.59, abs=0.05), (i, c.a_au)
        assert c.e == pytest.approx(0.393, abs=0.03), (i, c.e)
        # Clause (a): encounter achieved with maintenance applied.
        assert c.encounter_vinf_continuity_kms <= ENCOUNTER_VINF_TOL_KMS, (i, c)
        assert 0.0 < c.maintenance_dv_kms < 3.5, (i, c.maintenance_dv_kms)
        assert c.encounter_ok, (i, c)
        # Clause (b): bounded intra-cycle drift vs the planned trajectory.
        assert c.intra_cycle_drift_km <= INTRA_CYCLE_DRIFT_TOL_KM, (i, c)
        assert c.drift_ok, (i, c)
    assert result.v2_powered_passed, result.detail


@pytest.mark.slow
def test_aldrin_inbound_not_promoted_off_family_zero_dv(astropy_ephem: Ephemeris) -> None:
    """V2-powered EXCLUSION (no promotion): the inbound twin's real-window solve
    lands on a ballistic ΔV≈0 neighbour, not the in-family powered solve, so the
    'maintenance applied' half of clause (a) is not met. The inbound stays V1.
    """
    entry = next(e for e in load_m6b_entries() if e["id"] == "aldrin-classic-em-k1-inbound")
    assert entry["trajectory_regime"] == "powered"
    sol = solve_powered_periodic_cycler(
        entry,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
    )
    # The recorded off-family issue: the optimiser resolves to ΔV≈0 (a ballistic
    # neighbour) rather than the in-family powered Aldrin maintenance solve.
    assert sol.total_maintenance_dv_kms == pytest.approx(0.0, abs=1e-3), (
        f"inbound expected ballistic ΔV≈0 off-family neighbour; got "
        f"{sol.total_maintenance_dv_kms} km/s — if this is now in-family, "
        f"re-evaluate the inbound for V2-powered promotion"
    )
