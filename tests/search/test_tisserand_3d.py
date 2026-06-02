"""Tests for the 3-D Tisserand linkability predicate (:func:`linkable_3d`)
and the ephemeris-routed :func:`tisserand_feasible` path (STAGE 4).

The 3-D predicate extends the coplanar :func:`linkable` to inclined
spacecraft orbits: at a fixed V∞, body A's Tisserand equation fixes
``cos(i_sc)`` at each ``(a, e)``; body B's equation is tested at the same
``(a, e, i_sc)``. A 2-D grid scan over ``(a, e)`` asks whether both
equations agree.

GOLDEN-TEST DISCIPLINE: every EXPECTED value below is either a published
sourced anchor (Aldrin band V∞ ≈ 5.5 km/s, E-M; spec §9 / Rogers 2012) or
explicitly labelled ``# COMPUTED`` / ``# INVARIANT``. No fabricated numbers.

Formula (sourced — Strange & Longuski 2002 JSR 39(1):9-16):
    T_p = a_p/a + 2 cos(i) sqrt((a/a_p)(1 - e^2))
    V_inf^2 = (mu_sun / a_p) (3 - T_p)
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.sequence import Cell, tisserand_feasible
from cyclerfinder.search.tisserand import linkable, linkable_3d

# Aldrin-band E-M anchor: V∞ ≈ 5.5 km/s is the linkable band (spec §9,
# Rogers et al. 2012 Table 1). Used as a sourced INPUT.
ALDRIN_BAND_VINF_KMS: float = 5.5
HIGH_VINF_CAP_KMS: float = 8.0


# ---------------------------------------------------------------------------
# linkable_3d coplanar-limit and physics
# ---------------------------------------------------------------------------


def test_linkable_3d_em_at_5_5_kms_true() -> None:
    """# COMPUTED consistency: coplanar-limit 3-D agrees with 2-D at Aldrin band."""
    assert linkable_3d("E", "M", ALDRIN_BAND_VINF_KMS, i_sc_max_deg=0.01) is True


def test_linkable_3d_em_at_0_5_kms_false() -> None:
    """# COMPUTED physics: 0.5 km/s is far below the E-M linkable threshold;
    inclination (even 30 deg) cannot rescue it."""
    assert linkable_3d("E", "M", 0.5, i_sc_max_deg=30.0) is False


@pytest.mark.parametrize("vinf", [0.0, 1.0e-6, 1.0e6])
def test_linkable_3d_never_raises(vinf: float) -> None:
    """Contract: pathological V∞ returns False without raising, mirroring
    the :func:`linkable` contract."""
    assert linkable_3d("E", "M", vinf, i_sc_max_deg=30.0) is False


@pytest.mark.parametrize("vinf", [3.0, 5.5, 8.0])
def test_linkable_3d_symmetric(vinf: float) -> None:
    """# COMPUTED consistency: the predicate is symmetric in its body args."""
    assert linkable_3d("E", "M", vinf) == linkable_3d("M", "E", vinf)


def test_linkable_3d_coplanar_agrees_with_2d_for_em() -> None:
    """# COMPUTED consistency [GATE] — the critical regression guard.

    For ``i_sc_max_deg`` ~ 0 the 3-D predicate must reproduce the coplanar
    :func:`linkable` over the full V∞ scan band.
    """
    for vinf in np.linspace(0.5, 12.0, 40):
        v = float(vinf)
        assert linkable_3d("E", "M", v, i_sc_max_deg=0.01) == linkable("E", "M", v), (
            f"3-D coplanar limit disagrees with 2-D at vinf={v:.3f}"
        )


def test_linkable_3d_inclined_orbit_extends_or_preserves_linkable_range() -> None:
    """# COMPUTED monotonicity: a coplanar-True pair stays True when
    inclination is allowed (inclination opens options, never closes coplanar
    ones)."""
    for vinf in np.linspace(0.5, 12.0, 40):
        v = float(vinf)
        if linkable("E", "M", v):
            assert linkable_3d("E", "M", v, i_sc_max_deg=30.0) is True, (
                f"inclined predicate closed a coplanar-linkable option at vinf={v:.3f}"
            )


# ---------------------------------------------------------------------------
# tisserand_feasible routing
# ---------------------------------------------------------------------------


def _em_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )


def test_tisserand_feasible_routes_to_3d_when_ephem_provided() -> None:
    """# COMPUTED routing: with an ephemeris the 3-D path is active and does
    not regress the Aldrin-band feasibility result."""
    assert (
        tisserand_feasible(_em_cell(), vinf_cap=HIGH_VINF_CAP_KMS, ephem=Ephemeris("circular"))
        is True
    )


def test_tisserand_feasible_coplanar_default_unchanged() -> None:
    """# REGRESSION: ephem=None path is identical to pre-STAGE4 behaviour.

    Uses the same gate cells/caps as the existing test_sequence.py suite.
    """
    cell = _em_cell()
    assert tisserand_feasible(cell, vinf_cap=2.0) is False  # below Aldrin band
    assert tisserand_feasible(cell, vinf_cap=HIGH_VINF_CAP_KMS) is True  # covers band
