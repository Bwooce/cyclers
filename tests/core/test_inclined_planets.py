"""M-3D Phase 2: opt-in inclined planet table (plan §2). SOURCED anchor: the
inc/Ω values are Standish & Williams Table 1 (the EXPECTED side IS the source).
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import inclined_planets


def test_inclined_planets_carries_sourced_inc_lan() -> None:
    inc = inclined_planets()
    assert inc["V"].inc_deg == pytest.approx(3.39467605)
    assert inc["V"].lan_deg == pytest.approx(76.67984255)
    assert inc["M"].inc_deg == pytest.approx(1.84969142)
    assert inc["M"].lan_deg == pytest.approx(49.55953891)


def test_inclined_planets_does_not_mutate_live_PLANETS() -> None:  # noqa: N802
    _ = inclined_planets()
    assert PLANETS["V"].inc_deg == 0.0
    assert PLANETS["M"].inc_deg == 0.0  # live coplanar default UNTOUCHED


def test_inclined_circular_state_sits_on_node_with_z_excursion() -> None:
    """At t=0 each inclined body is on its ascending node (z==0); off-node it has
    a z-component of order a*sin(inc). SOURCED: the inclination is the anchor."""
    from cyclerfinder.core.ephemeris import Ephemeris

    ephem = Ephemeris.inclined_circular()
    r0, _ = ephem.state("V", 0.0)
    assert abs(float(r0[2])) < 1.0  # on the node at t=0 (km)
    # Quarter period later the z-excursion ~ a * sin(i). Sign/scale only.
    r_q, _ = ephem.state("V", 0.0 + 56.0 * 86400.0)  # ~quarter of Venus year
    assert abs(float(r_q[2])) > 1.0e6  # genuine out-of-ecliptic motion
