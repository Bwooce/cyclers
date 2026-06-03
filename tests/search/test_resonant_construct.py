"""Tests for :func:`cyclerfinder.search.resonant_construct.construct_resonant_cycler`.

Golden cross-check (non-circular)
---------------------------------
The construction takes S1L1's *sourced* orbit elements (a=1.30 AU, e=0.257;
Rogers 2012 Table 1 / spec §9) and computes the per-body V_inf. We assert the
COMPUTED V_inf against *independently sourced* COPLANAR V_inf anchors:

  * Russell & Ocampo 2004 tabulate 4.99 / 5.10 km/s (Earth / Mars) for the
    McConaghy "Notable Two-Synodic" cycler.
  * McConaghy 2006's abstract gives 4.7 / 5.0 km/s for the same family.

Orbit elements and V_inf are independently published, and they agree only in
the circular-coplanar model used here, so this is a legitimate golden test.

NOTE: the spec §9 values 5.65 / 3.05 km/s are a DIFFERENT, higher-fidelity
(real-ephemeris) figure and are NOT reproduced by this coplanar construction —
in particular the Mars 3.05 km/s requires an *eccentric* Mars orbit, which is
Task 3's real-ephemeris domain. We deliberately do not test against 5.65/3.05.
"""

from __future__ import annotations

import pytest

from cyclerfinder.search.resonant_construct import construct_resonant_cycler

# Sourced S1L1 heliocentric orbit (Rogers 2012 Table 1 / spec §9).
_S1L1_A_AU = 1.30
_S1L1_E = 0.257


def test_s1l1_earth_vinf_matches_russell_coplanar() -> None:
    """Computed Earth V_inf matches the Russell 2004 sourced coplanar anchor."""
    res = construct_resonant_cycler(a_au=_S1L1_A_AU, e=_S1L1_E)
    # Russell & Ocampo 2004 coplanar Earth V_inf = 4.99 km/s.
    assert abs(res.vinf_kms["E"] - 4.99) < 0.3


def test_s1l1_mars_vinf_matches_russell_coplanar() -> None:
    """Computed Mars V_inf matches the Russell 2004 sourced coplanar anchor."""
    res = construct_resonant_cycler(a_au=_S1L1_A_AU, e=_S1L1_E)
    # Russell & Ocampo 2004 coplanar Mars V_inf = 5.10 km/s.
    assert abs(res.vinf_kms["M"] - 5.10) < 0.3


def test_crossings_and_legs_are_self_consistent() -> None:
    """Crossing anomalies are outbound (in [0, pi]) and leg ToFs are positive.

    EXPECTED COMPUTED — structural contract, not a sourced value.
    """
    res = construct_resonant_cycler(a_au=_S1L1_A_AU, e=_S1L1_E)
    for body in ("E", "M"):
        assert 0.0 <= res.crossing_true_anom_rad[body] <= 3.1416
    for tof in res.leg_tofs_days.values():
        assert tof > 0.0
    # Mars (farther out) crosses at a larger true anomaly than Earth.
    assert res.crossing_true_anom_rad["M"] > res.crossing_true_anom_rad["E"]


def test_unreachable_body_raises() -> None:
    """A near-circular orbit cannot reach Mars -> ValueError."""
    with pytest.raises(ValueError, match="does not reach"):
        # a=1.0 AU, e=0.05 -> apoapsis 1.05 AU < Mars 1.524 AU.
        construct_resonant_cycler(a_au=1.0, e=0.05, bodies=("E", "M"))
