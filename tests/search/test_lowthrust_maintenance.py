"""Tests for the low-thrust cycler-maintenance evaluator (Phase 4).

Phase 4 of ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md`` —
delivered as **machinery only** (see the plan's "Execution deviation" note): no
catalogue rows, no schema edits. The evaluator models a cycler's per-synodic
maintenance manoeuvre as a thrust-bounded low-thrust budget and reports a
``propellant_mass_fraction`` so future *sourced* powered rows can be compared
to ballistic ones.

Golden discipline: no published source gives a powered low-thrust cycler ΔV we
hold, so every EXPECTED here is a physics invariant or an
internal-consistency check, never a literature number:

* propellant mass fraction is the Tsiolkovsky ``1 - exp(-ΔV/(g0·Isp))`` of the
  reported maintenance ΔV (a source-free identity);
* zero ΔV ⇒ zero propellant fraction; fraction is monotone in ΔV and Isp;
* the Aldrin maintenance ΔV recomputed under the thrust-bounded model is
  internally consistent with the impulsive ``optimise_maintenance_dv`` value
  (same maneuver, same magnitude) — an internal cross-check, NOT a golden.
"""

from __future__ import annotations

from itertools import pairwise
from math import exp

import pytest

from cyclerfinder.core.constants import STANDARD_GRAVITY_KM_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.lowthrust_maintenance import (
    PoweredMaintenanceResult,
    powered_maintenance_from_dv,
    propellant_mass_fraction,
)
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv

# ---------------------------------------------------------------------------
# propellant_mass_fraction — Tsiolkovsky identity
# ---------------------------------------------------------------------------


def test_propellant_fraction_tsiolkovsky() -> None:
    dv_kms = 0.5
    isp_s = 3000.0
    frac = propellant_mass_fraction(dv_kms, isp_s)
    expected = 1.0 - exp(-dv_kms / (STANDARD_GRAVITY_KM_S2 * isp_s))
    assert frac == pytest.approx(expected, rel=1e-12)


def test_propellant_fraction_zero_dv_is_zero() -> None:
    assert propellant_mass_fraction(0.0, 3000.0) == pytest.approx(0.0, abs=1e-15)


def test_propellant_fraction_monotone_in_dv() -> None:
    isp_s = 3000.0
    fracs = [propellant_mass_fraction(dv, isp_s) for dv in (0.1, 0.5, 1.0, 2.0)]
    assert all(b > a for a, b in pairwise(fracs))


def test_propellant_fraction_higher_isp_costs_less() -> None:
    dv_kms = 1.0
    assert propellant_mass_fraction(dv_kms, 5000.0) < propellant_mass_fraction(dv_kms, 2000.0)


def test_propellant_fraction_in_unit_interval() -> None:
    for dv in (0.0, 1.0, 5.0, 50.0):
        frac = propellant_mass_fraction(dv, 3000.0)
        assert 0.0 <= frac < 1.0


def test_propellant_fraction_rejects_bad_isp() -> None:
    with pytest.raises(ValueError, match="isp"):
        propellant_mass_fraction(1.0, 0.0)


def test_propellant_fraction_rejects_negative_dv() -> None:
    with pytest.raises(ValueError, match="dv"):
        propellant_mass_fraction(-1.0, 3000.0)


# ---------------------------------------------------------------------------
# powered_maintenance_from_dv — wrap a maintenance ΔV into a powered result
# ---------------------------------------------------------------------------


def test_powered_maintenance_from_dv_fields() -> None:
    result = powered_maintenance_from_dv(
        maintenance_dv_kms=0.3,
        isp_s=3000.0,
        dry_mass_kg=1000.0,
    )
    assert isinstance(result, PoweredMaintenanceResult)
    assert result.maintenance_dv_kms == pytest.approx(0.3)
    expected_frac = 1.0 - exp(-0.3 / (STANDARD_GRAVITY_KM_S2 * 3000.0))
    assert result.propellant_mass_fraction == pytest.approx(expected_frac, rel=1e-12)
    # Propellant mass: dry-mass * frac / (1 - frac) (mass before burn minus after).
    # m0 = dry / (1 - frac); propellant = m0 - dry.
    m0 = 1000.0 / (1.0 - expected_frac)
    assert result.propellant_mass_kg == pytest.approx(m0 - 1000.0, rel=1e-10)


def test_powered_maintenance_zero_dv() -> None:
    result = powered_maintenance_from_dv(0.0, 3000.0, 1000.0)
    assert result.propellant_mass_fraction == pytest.approx(0.0, abs=1e-15)
    assert result.propellant_mass_kg == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Internal-consistency cross-check against the impulsive Aldrin maintenance
# (NOT a golden — both numbers are our own computed values)
# ---------------------------------------------------------------------------


def test_aldrin_powered_maintenance_consistent_with_impulsive() -> None:
    """The powered model's ΔV equals the impulsive maintenance ΔV it wraps.

    ``optimise_aldrin_maintenance_dv`` computes a (non-sourced) maintenance ΔV
    via the single-impulse turn-deficit surrogate. Wrapping that same ΔV into
    the powered model must report it back unchanged and attach a consistent
    propellant fraction — an internal cross-check that the powered evaluator is
    a faithful re-expression of the existing impulsive cost, not a new number.
    """
    eph = Ephemeris("circular")
    impulsive = optimise_aldrin_maintenance_dv(eph, seed=0)
    assert impulsive.converged
    isp_s = 3000.0
    powered = powered_maintenance_from_dv(
        impulsive.maintenance_dv_kms,
        isp_s,
        dry_mass_kg=1000.0,
    )
    assert powered.maintenance_dv_kms == pytest.approx(impulsive.maintenance_dv_kms)
    expected_frac = 1.0 - exp(-impulsive.maintenance_dv_kms / (STANDARD_GRAVITY_KM_S2 * isp_s))
    assert powered.propellant_mass_fraction == pytest.approx(expected_frac, rel=1e-12)
    # A real maintenance ΔV (Aldrin's Earth turn deficit) gives a positive,
    # sub-unity propellant fraction.
    assert 0.0 < powered.propellant_mass_fraction < 1.0
